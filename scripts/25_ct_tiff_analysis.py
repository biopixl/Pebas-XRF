#!/usr/bin/env python3
"""
CT TIFF Data Analysis for Pebas Formation Cores
================================================
Analyzes converted CT TIFF images to extract density profiles,
identify core regions, and create visualizations.
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ==============================================================================
# CONFIGURATION
# ==============================================================================

BASE_PATH = Path("/Users/isaac/Documents/GitHub/PhD_Projects/Pebas-XRF")
TIFF_PATH = BASE_PATH / "Pebas-CT" / "TIFF"
OUTPUT_PATH = BASE_PATH / "output"
FIG_PATH = OUTPUT_PATH / "figures" / "ct_tiff_analysis"
TABLE_PATH = OUTPUT_PATH / "tables"

FIG_PATH.mkdir(parents=True, exist_ok=True)

# Threshold for identifying valid core material (normalized 0-1 scale)
# Air/background will be near 0, core material will be higher
CORE_THRESHOLD = 0.05  # 5% of normalized range

# ==============================================================================
# 1. LOAD TIFF FILES
# ==============================================================================

def load_tiff_series(series_dir, series_name):
    """Load all TIFF files from a series directory."""
    tiff_files = sorted(series_dir.glob("*.tif"))

    print(f"\n{series_name}: Loading {len(tiff_files)} TIFF files...")

    slices = []
    for tiff_file in tiff_files:
        # Parse filename for metadata
        # Format: SE1_TOP_slice002_loc-039.52mm.tif
        parts = tiff_file.stem.split("_")
        slice_num = int(parts[2].replace("slice", ""))
        loc_str = parts[3].replace("loc", "").replace("mm", "")
        slice_loc = float(loc_str)

        # Load image
        img = np.array(Image.open(tiff_file))

        slices.append({
            'file': tiff_file.name,
            'slice_num': slice_num,
            'slice_location': slice_loc,
            'image': img,
            'shape': img.shape,
            'dtype': img.dtype
        })

    # Sort by slice location
    slices.sort(key=lambda x: x['slice_location'])

    print(f"  Loaded {len(slices)} slices")
    print(f"  Image shape: {slices[0]['shape']}")
    print(f"  Slice location range: {slices[0]['slice_location']:.2f} to {slices[-1]['slice_location']:.2f} mm")

    return slices

# ==============================================================================
# 2. EXTRACT DENSITY PROFILES
# ==============================================================================

def extract_density_profile(slices, series_name):
    """Extract density statistics for each slice."""
    profiles = []

    for s in slices:
        img = s['image'].astype(np.float64)

        # Normalize to 0-1 if 16-bit
        if img.max() > 1:
            img = img / 65535.0

        # Calculate statistics
        mean_density = np.mean(img)
        std_density = np.std(img)
        min_density = np.min(img)
        max_density = np.max(img)

        # Identify core region (central portion of image)
        center_y, center_x = img.shape[0] // 2, img.shape[1] // 2
        radius = min(img.shape) // 4

        # Create circular mask for core region
        y, x = np.ogrid[:img.shape[0], :img.shape[1]]
        mask = ((x - center_x)**2 + (y - center_y)**2) <= radius**2

        core_mean = np.mean(img[mask])
        core_std = np.std(img[mask])

        # Determine if this slice contains valid core material
        is_core = core_mean > CORE_THRESHOLD

        profiles.append({
            'series': series_name,
            'slice_num': s['slice_num'],
            'slice_location': s['slice_location'],
            'mean_density': mean_density,
            'std_density': std_density,
            'min_density': min_density,
            'max_density': max_density,
            'core_mean': core_mean,
            'core_std': core_std,
            'is_core': is_core
        })

    return pd.DataFrame(profiles)

# ==============================================================================
# 3. ANALYZE CORE STRUCTURE
# ==============================================================================

def analyze_core_structure(profile_df, series_name):
    """Identify core segments and gaps."""
    print(f"\n{series_name} Core Structure Analysis:")

    # Identify contiguous core regions
    profile_df['core_change'] = profile_df['is_core'].diff().fillna(0).abs()

    segments = []
    current_segment = None

    for idx, row in profile_df.iterrows():
        if row['is_core']:
            if current_segment is None:
                current_segment = {
                    'start_loc': row['slice_location'],
                    'start_idx': idx,
                    'densities': [row['core_mean']]
                }
            else:
                current_segment['densities'].append(row['core_mean'])
        else:
            if current_segment is not None:
                current_segment['end_loc'] = profile_df.loc[idx-1, 'slice_location']
                current_segment['end_idx'] = idx - 1
                current_segment['length'] = current_segment['end_loc'] - current_segment['start_loc']
                current_segment['mean_density'] = np.mean(current_segment['densities'])
                segments.append(current_segment)
                current_segment = None

    # Handle last segment
    if current_segment is not None:
        current_segment['end_loc'] = profile_df.iloc[-1]['slice_location']
        current_segment['end_idx'] = len(profile_df) - 1
        current_segment['length'] = current_segment['end_loc'] - current_segment['start_loc']
        current_segment['mean_density'] = np.mean(current_segment['densities'])
        segments.append(current_segment)

    # Report
    total_core = sum(s['length'] for s in segments)
    total_range = profile_df['slice_location'].max() - profile_df['slice_location'].min()

    print(f"  Total slices: {len(profile_df)}")
    print(f"  Valid core slices: {profile_df['is_core'].sum()}")
    print(f"  Core segments: {len(segments)}")
    print(f"  Total core coverage: {total_core:.1f} mm")
    print(f"  Total range: {total_range:.1f} mm")

    for i, seg in enumerate(segments):
        print(f"    Segment {i+1}: {seg['start_loc']:.1f} to {seg['end_loc']:.1f} mm "
              f"({seg['length']:.1f} mm, mean density: {seg['mean_density']:.3f})")

    return segments

# ==============================================================================
# 4. CREATE VISUALIZATIONS
# ==============================================================================

def create_density_plot(se1_profile, se2_profile, segments_se1, segments_se2):
    """Create density profile visualization."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # SE1_TOP density profile
    ax1 = axes[0, 0]
    colors = ['#2ecc71' if x else '#e74c3c' for x in se1_profile['is_core']]
    ax1.bar(se1_profile['slice_location'], se1_profile['core_mean'],
            width=2, color=colors, edgecolor='none', alpha=0.7)
    ax1.axhline(y=CORE_THRESHOLD, color='gray', linestyle='--', label=f'Core threshold ({CORE_THRESHOLD})')
    ax1.set_xlabel('Slice Location (mm)')
    ax1.set_ylabel('Normalized Density')
    ax1.set_title('SE1_TOP Density Profile')
    ax1.legend()

    # SE2_BTTM density profile
    ax2 = axes[0, 1]
    colors = ['#2ecc71' if x else '#e74c3c' for x in se2_profile['is_core']]
    ax2.bar(se2_profile['slice_location'], se2_profile['core_mean'],
            width=2, color=colors, edgecolor='none', alpha=0.7)
    ax2.axhline(y=CORE_THRESHOLD, color='gray', linestyle='--', label=f'Core threshold ({CORE_THRESHOLD})')
    ax2.set_xlabel('Slice Location (mm)')
    ax2.set_ylabel('Normalized Density')
    ax2.set_title('SE2_BTTM Density Profile')
    ax2.legend()

    # SE1_TOP segment visualization
    ax3 = axes[1, 0]
    for i, seg in enumerate(segments_se1):
        ax3.barh(0, seg['length'], left=seg['start_loc'], height=0.5,
                color=plt.cm.viridis(seg['mean_density']), edgecolor='black')
        ax3.text(seg['start_loc'] + seg['length']/2, 0, f"{seg['length']:.0f}mm",
                ha='center', va='center', fontsize=9, color='white', fontweight='bold')
    ax3.set_xlim(se1_profile['slice_location'].min() - 5, se1_profile['slice_location'].max() + 5)
    ax3.set_ylim(-0.5, 0.5)
    ax3.set_xlabel('Slice Location (mm)')
    ax3.set_title(f'SE1_TOP Core Segments ({len(segments_se1)} segments)')
    ax3.set_yticks([])

    # SE2_BTTM segment visualization
    ax4 = axes[1, 1]
    for i, seg in enumerate(segments_se2):
        ax4.barh(0, seg['length'], left=seg['start_loc'], height=0.5,
                color=plt.cm.viridis(seg['mean_density']), edgecolor='black')
        ax4.text(seg['start_loc'] + seg['length']/2, 0, f"{seg['length']:.0f}mm",
                ha='center', va='center', fontsize=9, color='white', fontweight='bold')
    ax4.set_xlim(se2_profile['slice_location'].min() - 5, se2_profile['slice_location'].max() + 5)
    ax4.set_ylim(-0.5, 0.5)
    ax4.set_xlabel('Slice Location (mm)')
    ax4.set_title(f'SE2_BTTM Core Segments ({len(segments_se2)} segments)')
    ax4.set_yticks([])

    plt.tight_layout()
    plt.savefig(FIG_PATH / 'ct_density_profiles.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {FIG_PATH / 'ct_density_profiles.png'}")

def create_slice_montage(slices, series_name, n_cols=8):
    """Create a montage of CT slices."""
    n_slices = len(slices)
    n_rows = (n_slices + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 2*n_rows))
    axes = axes.flatten()

    for i, s in enumerate(slices):
        img = s['image'].astype(np.float64)
        if img.max() > 1:
            img = img / 65535.0

        axes[i].imshow(img, cmap='gray', vmin=0, vmax=0.3)
        axes[i].set_title(f"{s['slice_location']:.1f}mm", fontsize=8)
        axes[i].axis('off')

    # Hide unused axes
    for i in range(n_slices, len(axes)):
        axes[i].axis('off')

    plt.suptitle(f'{series_name} CT Slices', fontsize=14)
    plt.tight_layout()
    plt.savefig(FIG_PATH / f'ct_montage_{series_name}.png', dpi=100, bbox_inches='tight')
    plt.close()
    print(f"Saved: {FIG_PATH / f'ct_montage_{series_name}.png'}")

def create_3d_reconstruction(slices, series_name):
    """Create a simple 3D maximum intensity projection."""
    # Stack images
    stack = np.array([s['image'] for s in slices])

    if stack.max() > 1:
        stack = stack / 65535.0

    # Maximum intensity projections
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Axial MIP (along slice axis)
    mip_axial = np.max(stack, axis=0)
    axes[0].imshow(mip_axial, cmap='gray', vmin=0, vmax=0.3)
    axes[0].set_title('Axial MIP (top view)')
    axes[0].axis('off')

    # Sagittal MIP
    mip_sagittal = np.max(stack, axis=2)
    axes[1].imshow(mip_sagittal.T, cmap='gray', aspect='auto', vmin=0, vmax=0.3)
    axes[1].set_title('Sagittal MIP (side view)')
    axes[1].set_xlabel('Slice')
    axes[1].set_ylabel('Row')

    # Coronal MIP
    mip_coronal = np.max(stack, axis=1)
    axes[2].imshow(mip_coronal.T, cmap='gray', aspect='auto', vmin=0, vmax=0.3)
    axes[2].set_title('Coronal MIP (front view)')
    axes[2].set_xlabel('Slice')
    axes[2].set_ylabel('Column')

    plt.suptitle(f'{series_name} Maximum Intensity Projections', fontsize=14)
    plt.tight_layout()
    plt.savefig(FIG_PATH / f'ct_mip_{series_name}.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {FIG_PATH / f'ct_mip_{series_name}.png'}")

# ==============================================================================
# 5. MAIN ANALYSIS
# ==============================================================================

def main():
    print("=" * 70)
    print("CT TIFF DATA ANALYSIS")
    print("=" * 70)

    # Load TIFF files
    se1_slices = load_tiff_series(TIFF_PATH / "SE1_TOP", "SE1_TOP")
    se2_slices = load_tiff_series(TIFF_PATH / "SE2_BTTM", "SE2_BTTM")

    # Extract density profiles
    print("\n" + "=" * 70)
    print("EXTRACTING DENSITY PROFILES")
    print("=" * 70)

    se1_profile = extract_density_profile(se1_slices, "SE1_TOP")
    se2_profile = extract_density_profile(se2_slices, "SE2_BTTM")

    # Analyze core structure
    print("\n" + "=" * 70)
    print("CORE STRUCTURE ANALYSIS")
    print("=" * 70)

    segments_se1 = analyze_core_structure(se1_profile, "SE1_TOP")
    segments_se2 = analyze_core_structure(se2_profile, "SE2_BTTM")

    # Save profiles to CSV
    all_profiles = pd.concat([se1_profile, se2_profile], ignore_index=True)
    all_profiles.to_csv(TABLE_PATH / 'ct_tiff_density_profiles.csv', index=False)
    print(f"\nSaved: {TABLE_PATH / 'ct_tiff_density_profiles.csv'}")

    # Create visualizations
    print("\n" + "=" * 70)
    print("CREATING VISUALIZATIONS")
    print("=" * 70)

    create_density_plot(se1_profile, se2_profile, segments_se1, segments_se2)
    create_slice_montage(se1_slices, "SE1_TOP")
    create_slice_montage(se2_slices, "SE2_BTTM")
    create_3d_reconstruction(se1_slices, "SE1_TOP")
    create_3d_reconstruction(se2_slices, "SE2_BTTM")

    # Summary
    print("\n" + "=" * 70)
    print("ANALYSIS SUMMARY")
    print("=" * 70)

    print(f"\nSE1_TOP:")
    print(f"  Slices: {len(se1_slices)}")
    print(f"  Core segments: {len(segments_se1)}")
    print(f"  Total core: {sum(s['length'] for s in segments_se1):.1f} mm")

    print(f"\nSE2_BTTM:")
    print(f"  Slices: {len(se2_slices)}")
    print(f"  Core segments: {len(segments_se2)}")
    print(f"  Total core: {sum(s['length'] for s in segments_se2):.1f} mm")

    print(f"\nOutput files:")
    print(f"  {TABLE_PATH / 'ct_tiff_density_profiles.csv'}")
    print(f"  {FIG_PATH / 'ct_density_profiles.png'}")
    print(f"  {FIG_PATH / 'ct_montage_SE1_TOP.png'}")
    print(f"  {FIG_PATH / 'ct_montage_SE2_BTTM.png'}")
    print(f"  {FIG_PATH / 'ct_mip_SE1_TOP.png'}")
    print(f"  {FIG_PATH / 'ct_mip_SE2_BTTM.png'}")

if __name__ == "__main__":
    main()
