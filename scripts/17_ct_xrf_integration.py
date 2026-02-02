#!/usr/bin/env python3
"""
CT-XRF Data Integration for Pebas Formation Cores
==================================================
Reads DICOM CT data and integrates with XRF measurements for multi-proxy analysis.
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path
import pydicom
from pydicom.pixel_data_handlers.util import apply_voi_lut
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# CONFIGURATION
# ==============================================================================

BASE_PATH = Path("/Users/isaac/Documents/GitHub/PhD_Projects/Pebas-XRF")
CT_PATH = BASE_PATH / "Pebas-CT" / "ST1"
OUTPUT_PATH = BASE_PATH / "output"
FIG_PATH = OUTPUT_PATH / "figures" / "ct_integration"

FIG_PATH.mkdir(parents=True, exist_ok=True)

# ==============================================================================
# 1. READ ALL DICOM FILES
# ==============================================================================

def read_dicom_series(series_path):
    """Read all DICOM files in a series and return sorted by slice location."""
    dicom_files = []

    for f in sorted(series_path.iterdir()):
        if f.name.startswith("IM"):
            try:
                dcm = pydicom.dcmread(str(f))
                dicom_files.append({
                    'file': f.name,
                    'dcm': dcm,
                    'slice_location': float(dcm.SliceLocation),
                    'instance_number': int(dcm.InstanceNumber)
                })
            except Exception as e:
                print(f"  Error reading {f.name}: {e}")

    # Sort by slice location
    dicom_files.sort(key=lambda x: x['slice_location'])
    return dicom_files

print("=" * 60)
print("CT-XRF DATA INTEGRATION")
print("=" * 60)

# Read SE1 (TOP)
print("\nReading SE1 (TOP) series...")
se1_path = CT_PATH / "SE1"
se1_files = read_dicom_series(se1_path)
print(f"  Loaded {len(se1_files)} slices")

# Read SE2 (BOTTOM)
print("\nReading SE2 (BOTTOM) series...")
se2_path = CT_PATH / "SE2"
se2_files = read_dicom_series(se2_path)
print(f"  Loaded {len(se2_files)} slices")

# ==============================================================================
# 2. EXTRACT PIXEL DATA AND CREATE DENSITY PROFILES
# ==============================================================================

def extract_pixel_array(dcm):
    """Extract pixel array with proper windowing."""
    try:
        pixel_array = dcm.pixel_array
        # Apply VOI LUT (windowing) if available
        if hasattr(dcm, 'WindowCenter') and hasattr(dcm, 'WindowWidth'):
            pixel_array = apply_voi_lut(pixel_array, dcm)
        return pixel_array
    except Exception as e:
        print(f"  Error extracting pixels: {e}")
        return None

def extract_density_profile(dicom_files, profile_width=20):
    """
    Extract density profile along core axis.

    Parameters:
    -----------
    dicom_files : list of dict with 'dcm' and 'slice_location' keys
    profile_width : int, number of columns to average (centered)

    Returns:
    --------
    DataFrame with slice_location and density values
    """
    profiles = []

    for item in dicom_files:
        dcm = item['dcm']
        pixel_array = extract_pixel_array(dcm)

        if pixel_array is not None:
            # Get center columns
            mid_col = pixel_array.shape[1] // 2
            start_col = mid_col - profile_width // 2
            end_col = mid_col + profile_width // 2

            # Extract central strip
            central_strip = pixel_array[:, start_col:end_col]

            # Calculate statistics
            mean_density = np.mean(central_strip)
            std_density = np.std(central_strip)
            min_density = np.min(central_strip)
            max_density = np.max(central_strip)

            # Row-wise profile (depth profile within slice)
            row_profile = np.mean(central_strip, axis=1)

            profiles.append({
                'slice_location': item['slice_location'],
                'instance_number': item['instance_number'],
                'mean_density': mean_density,
                'std_density': std_density,
                'min_density': min_density,
                'max_density': max_density,
                'row_profile': row_profile,
                'n_rows': pixel_array.shape[0],
                'n_cols': pixel_array.shape[1]
            })

    return pd.DataFrame(profiles)

print("\n" + "=" * 60)
print("EXTRACTING DENSITY PROFILES")
print("=" * 60)

print("\nProcessing SE1 (TOP)...")
se1_profiles = extract_density_profile(se1_files)
print(f"  Extracted profiles for {len(se1_profiles)} slices")
print(f"  Depth range: {se1_profiles['slice_location'].min():.1f} to {se1_profiles['slice_location'].max():.1f} mm")
print(f"  Mean density range: {se1_profiles['mean_density'].min():.0f} to {se1_profiles['mean_density'].max():.0f}")

print("\nProcessing SE2 (BOTTOM)...")
se2_profiles = extract_density_profile(se2_files)
print(f"  Extracted profiles for {len(se2_profiles)} slices")
print(f"  Depth range: {se2_profiles['slice_location'].min():.1f} to {se2_profiles['slice_location'].max():.1f} mm")
print(f"  Mean density range: {se2_profiles['mean_density'].min():.0f} to {se2_profiles['mean_density'].max():.0f}")

# Add series identifier
se1_profiles['series'] = 'SE1_TOP'
se2_profiles['series'] = 'SE2_BTTM'

# Combine
ct_profiles = pd.concat([se1_profiles, se2_profiles], ignore_index=True)

# ==============================================================================
# 3. CREATE STACKED CT IMAGE (CORONAL VIEW)
# ==============================================================================

def create_stacked_image(dicom_files, output_path, series_name):
    """Create a stacked image showing all slices side by side."""
    images = []

    for item in dicom_files:
        pixel_array = extract_pixel_array(item['dcm'])
        if pixel_array is not None:
            # Normalize to 0-255
            img_norm = ((pixel_array - pixel_array.min()) /
                       (pixel_array.max() - pixel_array.min()) * 255).astype(np.uint8)
            images.append(img_norm)

    if images:
        # Stack horizontally (each slice is a column in the final image)
        # Take central strip from each image
        strips = []
        strip_width = 50  # pixels

        for img in images:
            mid = img.shape[1] // 2
            strip = img[:, mid - strip_width//2 : mid + strip_width//2]
            strips.append(strip)

        stacked = np.hstack(strips)

        # Save as PNG using PIL
        from PIL import Image
        img_pil = Image.fromarray(stacked)
        img_pil.save(output_path)
        print(f"  Saved: {output_path.name}")

        return stacked
    return None

print("\n" + "=" * 60)
print("CREATING CT VISUALIZATIONS")
print("=" * 60)

print("\nCreating stacked CT image for SE1...")
se1_stacked = create_stacked_image(se1_files, FIG_PATH / "ct_stacked_SE1_TOP.png", "SE1_TOP")

print("\nCreating stacked CT image for SE2...")
se2_stacked = create_stacked_image(se2_files, FIG_PATH / "ct_stacked_SE2_BTTM.png", "SE2_BTTM")

# ==============================================================================
# 4. SAVE CT DENSITY DATA
# ==============================================================================

print("\n" + "=" * 60)
print("SAVING CT DENSITY DATA")
print("=" * 60)

# Save summary profiles (without row_profile array for CSV)
ct_summary = ct_profiles.drop(columns=['row_profile']).copy()
ct_summary.to_csv(OUTPUT_PATH / "tables" / "ct_density_profiles.csv", index=False)
print(f"\nSaved: ct_density_profiles.csv ({len(ct_summary)} rows)")

# ==============================================================================
# 5. LOAD XRF DATA FOR COMPARISON
# ==============================================================================

print("\n" + "=" * 60)
print("LOADING XRF DATA FOR ALIGNMENT")
print("=" * 60)

xrf_data = pd.read_csv(OUTPUT_PATH / "tables" / "xrf_data_stacked.csv")
print(f"Loaded {len(xrf_data)} XRF measurements")

# XRF depth summary
xrf_summary = xrf_data.groupby('section').agg({
    'position_mm': ['min', 'max', 'count'],
    'Ca': 'median',
    'Ti': 'median'
}).reset_index()
xrf_summary.columns = ['section', 'min_mm', 'max_mm', 'n_points', 'Ca_median', 'Ti_median']
xrf_summary['range_mm'] = xrf_summary['max_mm'] - xrf_summary['min_mm']

# CT depth ranges
ct_se1_range = (se1_profiles['slice_location'].min(), se1_profiles['slice_location'].max())
ct_se2_range = (se2_profiles['slice_location'].min(), se2_profiles['slice_location'].max())

print(f"\nCT depth ranges:")
print(f"  SE1 (TOP): {ct_se1_range[0]:.1f} to {ct_se1_range[1]:.1f} mm ({ct_se1_range[1]-ct_se1_range[0]:.0f} mm)")
print(f"  SE2 (BTTM): {ct_se2_range[0]:.1f} to {ct_se2_range[1]:.1f} mm ({ct_se2_range[1]-ct_se2_range[0]:.0f} mm)")

# Find potential matching XRF sections based on depth range
ct_total_range = max(ct_se1_range[1], ct_se2_range[1]) - min(ct_se1_range[0], ct_se2_range[0])
print(f"\nTotal CT coverage: ~{ct_total_range:.0f} mm")

print("\nXRF sections with similar depth ranges (50-150 mm):")
similar = xrf_summary[(xrf_summary['range_mm'] >= 50) & (xrf_summary['range_mm'] <= 150)]
print(similar[['section', 'min_mm', 'max_mm', 'range_mm', 'n_points']].to_string(index=False))

# ==============================================================================
# 6. CREATE DEPTH-ALIGNED COMPARISON
# ==============================================================================

print("\n" + "=" * 60)
print("CREATING DEPTH-ALIGNED CT-XRF COMPARISON")
print("=" * 60)

# For demonstration, use TAM-5AB-6-7-A which has a similar depth range
# and starts near 0 (38mm) which may correspond to CT origin

test_section = "TAM-5AB-6-7-A"
xrf_test = xrf_data[xrf_data['section'] == test_section].copy()

if len(xrf_test) > 0:
    print(f"\nUsing {test_section} for alignment test:")
    print(f"  XRF range: {xrf_test['position_mm'].min():.0f} to {xrf_test['position_mm'].max():.0f} mm")
    print(f"  XRF measurements: {len(xrf_test)}")

    # Normalize depths for comparison
    # Assume CT slice_location 0 corresponds to XRF position at start of section
    xrf_test['depth_normalized'] = xrf_test['position_mm'] - xrf_test['position_mm'].min()

    # CT depth from SE1
    ct_test = se1_profiles.copy()
    ct_test['depth_normalized'] = ct_test['slice_location'] - ct_test['slice_location'].min()

    # Save aligned data
    xrf_test[['depth_normalized', 'position_mm', 'Ca', 'Ti', 'Fe', 'Mn', 'Ca_Ti', 'Fe_Mn']].to_csv(
        OUTPUT_PATH / "tables" / "xrf_for_ct_alignment.csv", index=False)
    ct_test[['depth_normalized', 'slice_location', 'mean_density', 'std_density']].to_csv(
        OUTPUT_PATH / "tables" / "ct_for_xrf_alignment.csv", index=False)

    print("\nSaved alignment data:")
    print("  - xrf_for_ct_alignment.csv")
    print("  - ct_for_xrf_alignment.csv")

# ==============================================================================
# 7. CREATE MATPLOTLIB VISUALIZATION
# ==============================================================================

try:
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    print("\n" + "=" * 60)
    print("CREATING INTEGRATED VISUALIZATION")
    print("=" * 60)

    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 3, height_ratios=[1, 1], width_ratios=[1, 1, 1])

    # Panel A: CT density profile SE1
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(se1_profiles['slice_location'], se1_profiles['mean_density'], 'b-', linewidth=1.5)
    ax1.fill_between(se1_profiles['slice_location'],
                     se1_profiles['mean_density'] - se1_profiles['std_density'],
                     se1_profiles['mean_density'] + se1_profiles['std_density'],
                     alpha=0.3)
    ax1.set_xlabel('Slice Location (mm)')
    ax1.set_ylabel('CT Density (mean intensity)')
    ax1.set_title('A. CT Density Profile - SE1 (TOP)')
    ax1.grid(True, alpha=0.3)

    # Panel B: CT density profile SE2
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(se2_profiles['slice_location'], se2_profiles['mean_density'], 'r-', linewidth=1.5)
    ax2.fill_between(se2_profiles['slice_location'],
                     se2_profiles['mean_density'] - se2_profiles['std_density'],
                     se2_profiles['mean_density'] + se2_profiles['std_density'],
                     alpha=0.3, color='red')
    ax2.set_xlabel('Slice Location (mm)')
    ax2.set_ylabel('CT Density (mean intensity)')
    ax2.set_title('B. CT Density Profile - SE2 (BTTM)')
    ax2.grid(True, alpha=0.3)

    # Panel C: CT stacked image
    ax3 = fig.add_subplot(gs[0, 2])
    if se1_stacked is not None:
        ax3.imshow(se1_stacked, cmap='gray', aspect='auto')
        ax3.set_xlabel('Slice Number')
        ax3.set_ylabel('Row (pixels)')
        ax3.set_title('C. CT Stacked Image - SE1')

    # Panel D: XRF Ca profile (test section)
    ax4 = fig.add_subplot(gs[1, 0])
    if len(xrf_test) > 0:
        ax4.plot(xrf_test['depth_normalized'], xrf_test['Ca'], 'g-', linewidth=1.5)
        ax4.set_xlabel('Depth (mm, normalized)')
        ax4.set_ylabel('Ca (counts)')
        ax4.set_title(f'D. XRF Ca Profile - {test_section}')
        ax4.grid(True, alpha=0.3)

    # Panel E: CT vs depth overlay
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.plot(ct_test['depth_normalized'], ct_test['mean_density'], 'b-', linewidth=1.5, label='CT Density')
    ax5.set_xlabel('Depth (mm, normalized)')
    ax5.set_ylabel('CT Density', color='b')
    ax5.tick_params(axis='y', labelcolor='b')
    ax5.set_title('E. CT-XRF Depth Alignment')

    # Overlay XRF Ca on secondary axis
    if len(xrf_test) > 0:
        ax5b = ax5.twinx()
        ax5b.plot(xrf_test['depth_normalized'], xrf_test['Ca'], 'g-', linewidth=1.5, alpha=0.7, label='XRF Ca')
        ax5b.set_ylabel('XRF Ca (counts)', color='g')
        ax5b.tick_params(axis='y', labelcolor='g')
    ax5.grid(True, alpha=0.3)

    # Panel F: Correlation scatter
    ax6 = fig.add_subplot(gs[1, 2])
    # Interpolate CT to XRF depths for correlation
    if len(xrf_test) > 0 and len(ct_test) > 0:
        from scipy import interpolate

        # Only use overlapping depth range
        depth_min = max(ct_test['depth_normalized'].min(), xrf_test['depth_normalized'].min())
        depth_max = min(ct_test['depth_normalized'].max(), xrf_test['depth_normalized'].max())

        if depth_max > depth_min:
            # Interpolate CT to XRF depths
            ct_interp = interpolate.interp1d(ct_test['depth_normalized'], ct_test['mean_density'],
                                             kind='linear', fill_value='extrapolate')

            xrf_overlap = xrf_test[(xrf_test['depth_normalized'] >= depth_min) &
                                   (xrf_test['depth_normalized'] <= depth_max)].copy()
            xrf_overlap['ct_density'] = ct_interp(xrf_overlap['depth_normalized'])

            ax6.scatter(xrf_overlap['Ca'], xrf_overlap['ct_density'], alpha=0.6, s=30)
            ax6.set_xlabel('XRF Ca (counts)')
            ax6.set_ylabel('CT Density (interpolated)')
            ax6.set_title('F. CT Density vs XRF Ca')

            # Add correlation coefficient
            if len(xrf_overlap) > 5:
                corr = xrf_overlap['Ca'].corr(xrf_overlap['ct_density'])
                ax6.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax6.transAxes,
                        fontsize=12, verticalalignment='top')
            ax6.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIG_PATH / "ct_xrf_integration_summary.png", dpi=300, bbox_inches='tight')
    plt.savefig(FIG_PATH / "ct_xrf_integration_summary.pdf", bbox_inches='tight')
    print("\nSaved: ct_xrf_integration_summary.png/pdf")

except ImportError as e:
    print(f"\nMatplotlib not available: {e}")
    print("Skipping visualization...")

# ==============================================================================
# 8. SUMMARY
# ==============================================================================

print("\n" + "=" * 60)
print("CT-XRF INTEGRATION SUMMARY")
print("=" * 60)

print(f"""
CT Data:
  - SE1 (TOP): {len(se1_profiles)} slices, {ct_se1_range[0]:.1f} to {ct_se1_range[1]:.1f} mm
  - SE2 (BTTM): {len(se2_profiles)} slices, {ct_se2_range[0]:.1f} to {ct_se2_range[1]:.1f} mm
  - Total coverage: ~{ct_total_range:.0f} mm

Output Files:
  - ct_density_profiles.csv: CT density data for all slices
  - ct_stacked_SE1_TOP.png: Stacked CT image (SE1)
  - ct_stacked_SE2_BTTM.png: Stacked CT image (SE2)
  - ct_xrf_integration_summary.png/pdf: Integrated visualization
  - xrf_for_ct_alignment.csv: XRF data formatted for alignment
  - ct_for_xrf_alignment.csv: CT data formatted for alignment

Next Steps:
  1. Verify CT-XRF core correspondence using sample records
  2. Refine depth alignment with visual feature matching
  3. Calculate formal CT-Ca correlation coefficient
  4. Integrate CT density into stratigraphic figures
""")
