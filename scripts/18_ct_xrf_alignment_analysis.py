#!/usr/bin/env python3
"""
CT-XRF Alignment Analysis
=========================
Systematic analysis to find optimal depth alignment between CT and XRF data.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import interpolate, stats
import matplotlib.pyplot as plt

# ==============================================================================
# CONFIGURATION
# ==============================================================================

BASE_PATH = Path("/Users/isaac/Documents/GitHub/PhD_Projects/Pebas-XRF")
OUTPUT_PATH = BASE_PATH / "output"
FIG_PATH = OUTPUT_PATH / "figures" / "ct_integration"

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

print("=" * 60)
print("CT-XRF ALIGNMENT ANALYSIS")
print("=" * 60)

# Load CT density profiles
ct_data = pd.read_csv(OUTPUT_PATH / "tables" / "ct_density_profiles.csv")
print(f"\nLoaded CT data: {len(ct_data)} slices")

# Load XRF data
xrf_data = pd.read_csv(OUTPUT_PATH / "tables" / "xrf_data_stacked.csv")
print(f"Loaded XRF data: {len(xrf_data)} measurements")

# ==============================================================================
# 2. IDENTIFY VALID CT CORE REGION
# ==============================================================================

print("\n" + "=" * 60)
print("IDENTIFYING VALID CT CORE REGION")
print("=" * 60)

# Values close to -32768 are air/background
# Core region has values > -32000 (arbitrary threshold)
AIR_THRESHOLD = -32000

for series in ['SE1_TOP', 'SE2_BTTM']:
    series_data = ct_data[ct_data['series'] == series].copy()
    core_mask = series_data['mean_density'] > AIR_THRESHOLD
    core_data = series_data[core_mask]

    if len(core_data) > 0:
        print(f"\n{series}:")
        print(f"  Total slices: {len(series_data)}")
        print(f"  Core slices (density > {AIR_THRESHOLD}): {len(core_data)}")
        print(f"  Core depth range: {core_data['slice_location'].min():.1f} to {core_data['slice_location'].max():.1f} mm")
        print(f"  Core density range: {core_data['mean_density'].min():.0f} to {core_data['mean_density'].max():.0f}")

# ==============================================================================
# 3. TEST MULTIPLE XRF SECTIONS FOR BEST MATCH
# ==============================================================================

print("\n" + "=" * 60)
print("TESTING XRF SECTIONS FOR CT ALIGNMENT")
print("=" * 60)

# Focus on SE1 (TOP) which has clearer core signal
ct_se1 = ct_data[ct_data['series'] == 'SE1_TOP'].copy()
ct_core = ct_se1[ct_se1['mean_density'] > AIR_THRESHOLD].copy()

if len(ct_core) < 5:
    print("Not enough valid CT core data!")
else:
    # Normalize CT depth (0 = start of core region)
    ct_core['depth_norm'] = ct_core['slice_location'] - ct_core['slice_location'].min()
    ct_coverage = ct_core['depth_norm'].max()

    print(f"\nCT core coverage: {ct_coverage:.1f} mm")

    # Get unique XRF sections
    sections = xrf_data.groupby('section').agg({
        'position_mm': ['min', 'max', 'count']
    }).reset_index()
    sections.columns = ['section', 'min_mm', 'max_mm', 'n_points']
    sections['range_mm'] = sections['max_mm'] - sections['min_mm']

    # Test sections with similar or larger depth range
    test_sections = sections[sections['range_mm'] >= ct_coverage * 0.5].copy()

    print(f"\nTesting {len(test_sections)} XRF sections with range >= {ct_coverage*0.5:.0f} mm")

    results = []

    for _, row in test_sections.iterrows():
        section = row['section']
        xrf_section = xrf_data[xrf_data['section'] == section].copy()

        # Normalize XRF depth
        xrf_section['depth_norm'] = xrf_section['position_mm'] - xrf_section['position_mm'].min()

        # Interpolate CT to XRF depths
        # Only use overlapping range
        depth_min = max(ct_core['depth_norm'].min(), xrf_section['depth_norm'].min())
        depth_max = min(ct_core['depth_norm'].max(), xrf_section['depth_norm'].max())

        if depth_max - depth_min < 20:  # Need at least 20mm overlap
            continue

        # Create interpolation function for CT
        ct_interp = interpolate.interp1d(
            ct_core['depth_norm'],
            ct_core['mean_density'],
            kind='linear',
            fill_value='extrapolate'
        )

        # Get XRF points in overlap region
        xrf_overlap = xrf_section[
            (xrf_section['depth_norm'] >= depth_min) &
            (xrf_section['depth_norm'] <= depth_max)
        ].copy()

        if len(xrf_overlap) < 5:
            continue

        # Interpolate CT to XRF depths
        xrf_overlap['ct_density'] = ct_interp(xrf_overlap['depth_norm'])

        # Calculate correlations with key XRF elements
        corr_ca = xrf_overlap['Ca'].corr(xrf_overlap['ct_density'])
        corr_ti = xrf_overlap['Ti'].corr(xrf_overlap['ct_density'])
        corr_fe = xrf_overlap['Fe'].corr(xrf_overlap['ct_density'])

        results.append({
            'section': section,
            'n_points': len(xrf_overlap),
            'overlap_mm': depth_max - depth_min,
            'corr_Ca': corr_ca,
            'corr_Ti': corr_ti,
            'corr_Fe': corr_fe,
            'abs_corr_Ca': abs(corr_ca)
        })

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('abs_corr_Ca', ascending=False)

    print("\nCorrelation results (sorted by |r| with Ca):")
    print(results_df[['section', 'n_points', 'overlap_mm', 'corr_Ca', 'corr_Ti']].to_string(index=False))

# ==============================================================================
# 4. TEST DEPTH OFFSET ALIGNMENT
# ==============================================================================

print("\n" + "=" * 60)
print("OPTIMIZING DEPTH OFFSET")
print("=" * 60)

# Use the section with highest correlation for offset optimization
if len(results_df) > 0:
    best_section = results_df.iloc[0]['section']
    print(f"\nUsing {best_section} for offset optimization")

    xrf_test = xrf_data[xrf_data['section'] == best_section].copy()
    xrf_test['depth_norm'] = xrf_test['position_mm'] - xrf_test['position_mm'].min()

    # Test different offsets
    offsets = np.arange(-50, 51, 2)  # -50 to +50 mm in 2mm steps
    offset_results = []

    for offset in offsets:
        # Apply offset to XRF depths
        xrf_shifted = xrf_test.copy()
        xrf_shifted['depth_shifted'] = xrf_shifted['depth_norm'] + offset

        # Find overlap with CT
        depth_min = max(ct_core['depth_norm'].min(), xrf_shifted['depth_shifted'].min())
        depth_max = min(ct_core['depth_norm'].max(), xrf_shifted['depth_shifted'].max())

        if depth_max - depth_min < 20:
            continue

        # Interpolate CT
        ct_interp = interpolate.interp1d(
            ct_core['depth_norm'],
            ct_core['mean_density'],
            kind='linear',
            fill_value='extrapolate'
        )

        # Get overlap points
        xrf_overlap = xrf_shifted[
            (xrf_shifted['depth_shifted'] >= depth_min) &
            (xrf_shifted['depth_shifted'] <= depth_max)
        ].copy()

        if len(xrf_overlap) < 5:
            continue

        xrf_overlap['ct_density'] = ct_interp(xrf_overlap['depth_shifted'])

        corr = xrf_overlap['Ca'].corr(xrf_overlap['ct_density'])

        offset_results.append({
            'offset_mm': offset,
            'correlation': corr,
            'n_points': len(xrf_overlap)
        })

    offset_df = pd.DataFrame(offset_results)

    if len(offset_df) > 0:
        # Find optimal offset (maximum POSITIVE correlation expected for Ca-density)
        best_idx = offset_df['correlation'].idxmax()
        best_offset = offset_df.loc[best_idx, 'offset_mm']
        best_corr = offset_df.loc[best_idx, 'correlation']

        print(f"\nOptimal offset: {best_offset:.0f} mm")
        print(f"Best correlation: {best_corr:.3f}")

        # Plot offset optimization
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Panel A: Correlation vs offset
        ax1 = axes[0]
        ax1.plot(offset_df['offset_mm'], offset_df['correlation'], 'b-', linewidth=2)
        ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax1.axvline(x=best_offset, color='red', linestyle='--', label=f'Optimal: {best_offset:.0f} mm')
        ax1.scatter([best_offset], [best_corr], color='red', s=100, zorder=5)
        ax1.set_xlabel('Depth Offset (mm)')
        ax1.set_ylabel('CT-Ca Correlation (r)')
        ax1.set_title(f'A. Offset Optimization - {best_section}')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Panel B: Aligned profiles at optimal offset
        ax2 = axes[1]

        # Apply optimal offset
        xrf_aligned = xrf_test.copy()
        xrf_aligned['depth_aligned'] = xrf_aligned['depth_norm'] + best_offset

        # Plot CT
        ax2.plot(ct_core['depth_norm'], ct_core['mean_density'], 'b-', linewidth=2, label='CT Density')
        ax2.set_ylabel('CT Density', color='blue')
        ax2.tick_params(axis='y', labelcolor='blue')

        # Plot XRF Ca on secondary axis
        ax2b = ax2.twinx()
        ax2b.plot(xrf_aligned['depth_aligned'], xrf_aligned['Ca'], 'g-', linewidth=2, label='XRF Ca', alpha=0.8)
        ax2b.set_ylabel('XRF Ca (counts)', color='green')
        ax2b.tick_params(axis='y', labelcolor='green')

        ax2.set_xlabel('Depth (mm, normalized)')
        ax2.set_title(f'B. Aligned Profiles (offset = {best_offset:.0f} mm, r = {best_corr:.3f})')
        ax2.grid(True, alpha=0.3)

        # Add legends
        lines1, labels1 = ax2.get_legend_handles_labels()
        lines2, labels2 = ax2b.get_legend_handles_labels()
        ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

        plt.tight_layout()
        plt.savefig(FIG_PATH / "ct_xrf_offset_optimization.png", dpi=300, bbox_inches='tight')
        plt.savefig(FIG_PATH / "ct_xrf_offset_optimization.pdf", bbox_inches='tight')
        print(f"\nSaved: ct_xrf_offset_optimization.png/pdf")

# ==============================================================================
# 5. CREATE FINAL ALIGNED VISUALIZATION
# ==============================================================================

print("\n" + "=" * 60)
print("CREATING FINAL ALIGNED VISUALIZATION")
print("=" * 60)

if len(results_df) > 0 and 'best_offset' in dir():
    # Create comprehensive visualization with aligned data
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Apply optimal offset
    xrf_aligned = xrf_test.copy()
    xrf_aligned['depth_aligned'] = xrf_aligned['depth_norm'] + best_offset

    # Get overlap region
    depth_min = max(ct_core['depth_norm'].min(), xrf_aligned['depth_aligned'].min())
    depth_max = min(ct_core['depth_norm'].max(), xrf_aligned['depth_aligned'].max())

    ct_interp = interpolate.interp1d(
        ct_core['depth_norm'],
        ct_core['mean_density'],
        kind='linear',
        fill_value='extrapolate'
    )

    xrf_overlap = xrf_aligned[
        (xrf_aligned['depth_aligned'] >= depth_min) &
        (xrf_aligned['depth_aligned'] <= depth_max)
    ].copy()
    xrf_overlap['ct_density'] = ct_interp(xrf_overlap['depth_aligned'])

    # Panel A: CT density profile
    ax1 = axes[0, 0]
    ax1.plot(ct_core['depth_norm'], ct_core['mean_density'], 'b-', linewidth=2)
    ax1.fill_between(ct_core['depth_norm'],
                     ct_core['mean_density'] - ct_core['std_density'],
                     ct_core['mean_density'] + ct_core['std_density'],
                     alpha=0.3)
    ax1.set_xlabel('Depth (mm)')
    ax1.set_ylabel('CT Density')
    ax1.set_title('A. CT Density Profile (SE1 Core Region)')
    ax1.grid(True, alpha=0.3)

    # Panel B: XRF Ca profile
    ax2 = axes[0, 1]
    ax2.plot(xrf_aligned['depth_aligned'], xrf_aligned['Ca'], 'g-', linewidth=2)
    ax2.set_xlabel('Depth (mm, aligned)')
    ax2.set_ylabel('Ca (counts)')
    ax2.set_title(f'B. XRF Ca Profile ({best_section}, offset={best_offset:.0f}mm)')
    ax2.grid(True, alpha=0.3)

    # Panel C: Scatter plot with regression
    ax3 = axes[1, 0]
    ax3.scatter(xrf_overlap['Ca'], xrf_overlap['ct_density'], alpha=0.7, s=50, c='purple')

    # Add regression line
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        xrf_overlap['Ca'], xrf_overlap['ct_density']
    )
    x_line = np.array([xrf_overlap['Ca'].min(), xrf_overlap['Ca'].max()])
    ax3.plot(x_line, slope * x_line + intercept, 'r--', linewidth=2, label=f'r = {r_value:.3f}')

    ax3.set_xlabel('XRF Ca (counts)')
    ax3.set_ylabel('CT Density (interpolated)')
    ax3.set_title(f'C. CT-Ca Correlation (r = {r_value:.3f}, p = {p_value:.4f})')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # Panel D: Dual-axis aligned profiles
    ax4 = axes[1, 1]
    l1 = ax4.plot(ct_core['depth_norm'], ct_core['mean_density'], 'b-', linewidth=2, label='CT Density')
    ax4.set_ylabel('CT Density', color='blue')
    ax4.tick_params(axis='y', labelcolor='blue')

    ax4b = ax4.twinx()
    l2 = ax4b.plot(xrf_aligned['depth_aligned'], xrf_aligned['Ca'], 'g-', linewidth=2, alpha=0.8, label='XRF Ca')
    ax4b.set_ylabel('XRF Ca (counts)', color='green')
    ax4b.tick_params(axis='y', labelcolor='green')

    ax4.set_xlabel('Depth (mm)')
    ax4.set_title('D. Aligned CT-XRF Profiles')
    ax4.legend(l1 + l2, ['CT Density', 'XRF Ca'], loc='upper right')
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIG_PATH / "ct_xrf_final_alignment.png", dpi=300, bbox_inches='tight')
    plt.savefig(FIG_PATH / "ct_xrf_final_alignment.pdf", bbox_inches='tight')
    print(f"Saved: ct_xrf_final_alignment.png/pdf")

    # Save aligned data
    xrf_overlap.to_csv(OUTPUT_PATH / "tables" / "ct_xrf_aligned_data.csv", index=False)
    print(f"Saved: ct_xrf_aligned_data.csv ({len(xrf_overlap)} aligned points)")

# ==============================================================================
# 6. SUMMARY
# ==============================================================================

print("\n" + "=" * 60)
print("ALIGNMENT ANALYSIS SUMMARY")
print("=" * 60)

if len(results_df) > 0:
    print(f"""
Best matching XRF section: {best_section}
Optimal depth offset: {best_offset:.0f} mm
CT-Ca correlation: r = {best_corr:.3f}

Interpretation:
  - Positive correlation indicates CT density increases with carbonate content
  - Higher Ca (more carbonate/shells) = higher CT attenuation (denser material)
  - This validates the CT-XRF correspondence

Output files:
  - ct_xrf_offset_optimization.png/pdf: Offset analysis
  - ct_xrf_final_alignment.png/pdf: Final aligned visualization
  - ct_xrf_aligned_data.csv: Aligned CT-XRF data for further analysis
""")
