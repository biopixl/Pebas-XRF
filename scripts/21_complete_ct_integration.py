#!/usr/bin/env python3
"""
Complete CT-XRF Integration Analysis
=====================================
Comprehensive alignment and correlation analysis using ALL CT data
(both SE1_TOP and SE2_BTTM series) across all XRF sections.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import interpolate, stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# CONFIGURATION
# ==============================================================================

BASE_PATH = Path("/Users/isaac/Documents/GitHub/PhD_Projects/Pebas-XRF")
OUTPUT_PATH = BASE_PATH / "output"
FIG_PATH = OUTPUT_PATH / "figures" / "ct_integration"
FIG_PATH.mkdir(parents=True, exist_ok=True)

# Air/background threshold for identifying valid core regions
AIR_THRESHOLD = -32000

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

print("=" * 70)
print("COMPLETE CT-XRF INTEGRATION ANALYSIS")
print("=" * 70)

# Load CT density profiles
ct_data = pd.read_csv(OUTPUT_PATH / "tables" / "ct_density_profiles.csv")
print(f"\nLoaded CT data: {len(ct_data)} total slices")

# Separate by series
ct_se1 = ct_data[ct_data['series'] == 'SE1_TOP'].copy()
ct_se2 = ct_data[ct_data['series'] == 'SE2_BTTM'].copy()
print(f"  SE1_TOP: {len(ct_se1)} slices")
print(f"  SE2_BTTM: {len(ct_se2)} slices")

# Load XRF data
xrf_data = pd.read_csv(OUTPUT_PATH / "tables" / "xrf_data_stacked.csv")
print(f"\nLoaded XRF data: {len(xrf_data)} measurements")

# ==============================================================================
# 2. IDENTIFY VALID CORE REGIONS IN BOTH CT SERIES
# ==============================================================================

print("\n" + "=" * 70)
print("IDENTIFYING VALID CORE REGIONS")
print("=" * 70)

def get_valid_core_region(ct_series, series_name):
    """Identify valid core material (excluding air/background)."""
    core_mask = ct_series['mean_density'] > AIR_THRESHOLD
    core_data = ct_series[core_mask].copy()

    if len(core_data) > 0:
        core_data['depth_norm'] = core_data['slice_location'] - core_data['slice_location'].min()
        coverage = core_data['depth_norm'].max()
        print(f"\n{series_name}:")
        print(f"  Total slices: {len(ct_series)}")
        print(f"  Valid core slices: {len(core_data)}")
        print(f"  Depth range: {core_data['slice_location'].min():.1f} to {core_data['slice_location'].max():.1f} mm")
        print(f"  Core coverage: {coverage:.1f} mm")
        print(f"  Density range: {core_data['mean_density'].min():.0f} to {core_data['mean_density'].max():.0f}")
        return core_data
    return pd.DataFrame()

ct_se1_core = get_valid_core_region(ct_se1, "SE1_TOP")
ct_se2_core = get_valid_core_region(ct_se2, "SE2_BTTM")

# ==============================================================================
# 3. TEST ALL XRF SECTIONS AGAINST BOTH CT SERIES
# ==============================================================================

print("\n" + "=" * 70)
print("TESTING ALL XRF SECTIONS AGAINST BOTH CT SERIES")
print("=" * 70)

def test_ct_xrf_alignment(ct_core, xrf_section_data, offset_range=(-100, 101, 2)):
    """
    Test alignment between CT core and XRF section data.

    Returns best correlation, offset, and aligned data.
    """
    if len(ct_core) < 5 or len(xrf_section_data) < 5:
        return None

    # Normalize XRF depth
    xrf_data_norm = xrf_section_data.copy()
    xrf_data_norm['depth_norm'] = xrf_data_norm['position_mm'] - xrf_data_norm['position_mm'].min()

    ct_coverage = ct_core['depth_norm'].max()
    xrf_coverage = xrf_data_norm['depth_norm'].max()

    # Test different offsets
    offsets = np.arange(*offset_range)
    best_corr = -2
    best_offset = 0
    best_data = None

    for offset in offsets:
        xrf_shifted = xrf_data_norm.copy()
        xrf_shifted['depth_shifted'] = xrf_shifted['depth_norm'] + offset

        # Find overlap
        depth_min = max(ct_core['depth_norm'].min(), xrf_shifted['depth_shifted'].min())
        depth_max = min(ct_core['depth_norm'].max(), xrf_shifted['depth_shifted'].max())

        if depth_max - depth_min < 15:  # Need at least 15mm overlap
            continue

        # Interpolate CT to XRF depths
        try:
            ct_interp = interpolate.interp1d(
                ct_core['depth_norm'],
                ct_core['mean_density'],
                kind='linear',
                fill_value='extrapolate'
            )

            xrf_overlap = xrf_shifted[
                (xrf_shifted['depth_shifted'] >= depth_min) &
                (xrf_shifted['depth_shifted'] <= depth_max)
            ].copy()

            if len(xrf_overlap) < 5:
                continue

            xrf_overlap['ct_density'] = ct_interp(xrf_overlap['depth_shifted'])

            # Calculate correlation with Ca
            corr = xrf_overlap['Ca'].corr(xrf_overlap['ct_density'])

            if not np.isnan(corr) and corr > best_corr:
                best_corr = corr
                best_offset = offset
                best_data = xrf_overlap.copy()
        except:
            continue

    if best_data is not None and len(best_data) >= 5:
        return {
            'correlation': best_corr,
            'offset': best_offset,
            'n_points': len(best_data),
            'overlap_mm': best_data['depth_shifted'].max() - best_data['depth_shifted'].min(),
            'aligned_data': best_data
        }
    return None

# Get unique sections
sections = xrf_data['section'].unique()

# Test all combinations
results = []

for section in sections:
    xrf_section = xrf_data[xrf_data['section'] == section].copy()

    # Skip excluded or very short sections
    xrf_valid = xrf_section[~xrf_section['excluded']]
    if len(xrf_valid) < 10:
        continue

    # Test against SE1_TOP
    if len(ct_se1_core) > 0:
        result = test_ct_xrf_alignment(ct_se1_core, xrf_valid)
        if result:
            results.append({
                'section': section,
                'ct_series': 'SE1_TOP',
                **result
            })

    # Test against SE2_BTTM
    if len(ct_se2_core) > 0:
        result = test_ct_xrf_alignment(ct_se2_core, xrf_valid)
        if result:
            results.append({
                'section': section,
                'ct_series': 'SE2_BTTM',
                **result
            })

# Create results dataframe
if results:
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('correlation', ascending=False)

    print("\nTop 10 CT-XRF correlations (all combinations):")
    print("-" * 80)
    display_cols = ['section', 'ct_series', 'correlation', 'offset', 'n_points', 'overlap_mm']
    print(results_df[display_cols].head(10).to_string(index=False))

    # Save results
    results_df[display_cols].to_csv(OUTPUT_PATH / "tables" / "ct_xrf_all_correlations.csv", index=False)
    print(f"\nSaved: ct_xrf_all_correlations.csv ({len(results_df)} combinations tested)")

# ==============================================================================
# 4. DETAILED ANALYSIS OF BEST MATCHES
# ==============================================================================

print("\n" + "=" * 70)
print("DETAILED ANALYSIS OF BEST MATCHES")
print("=" * 70)

# Get best match for each CT series
best_se1 = results_df[results_df['ct_series'] == 'SE1_TOP'].iloc[0] if len(results_df[results_df['ct_series'] == 'SE1_TOP']) > 0 else None
best_se2 = results_df[results_df['ct_series'] == 'SE2_BTTM'].iloc[0] if len(results_df[results_df['ct_series'] == 'SE2_BTTM']) > 0 else None

best_matches = []
if best_se1 is not None:
    best_matches.append(('SE1_TOP', best_se1, ct_se1_core))
if best_se2 is not None:
    best_matches.append(('SE2_BTTM', best_se2, ct_se2_core))

# Calculate detailed statistics for best matches
detailed_results = []

for series_name, match, ct_core in best_matches:
    section = match['section']
    offset = match['offset']
    aligned_data = match['aligned_data']

    print(f"\n{series_name} Best Match: {section}")
    print(f"  Offset: {offset} mm")
    print(f"  N points: {match['n_points']}")
    print(f"  Overlap: {match['overlap_mm']:.1f} mm")
    print(f"  Correlation: r = {match['correlation']:.4f}")

    # Calculate p-value
    if len(aligned_data) >= 3:
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            aligned_data['Ca'], aligned_data['ct_density']
        )
        print(f"  P-value: {p_value:.6f}")
        print(f"  Significant: {'Yes' if p_value < 0.05 else 'No'} (p < 0.05)")

        detailed_results.append({
            'ct_series': series_name,
            'section': section,
            'offset_mm': offset,
            'n_points': match['n_points'],
            'overlap_mm': match['overlap_mm'],
            'correlation': match['correlation'],
            'p_value': p_value,
            'significant': p_value < 0.05,
            'slope': slope,
            'intercept': intercept
        })

# Save detailed results
if detailed_results:
    detailed_df = pd.DataFrame(detailed_results)
    detailed_df.to_csv(OUTPUT_PATH / "tables" / "ct_xrf_best_matches.csv", index=False)
    print(f"\nSaved: ct_xrf_best_matches.csv")

# ==============================================================================
# 5. CREATE COMPREHENSIVE VISUALIZATION
# ==============================================================================

print("\n" + "=" * 70)
print("CREATING COMPREHENSIVE VISUALIZATION")
print("=" * 70)

fig = plt.figure(figsize=(16, 12))
gs = gridspec.GridSpec(3, 3, height_ratios=[1, 1, 1])

# Panel A: SE1 CT density profile
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(ct_se1_core['depth_norm'], ct_se1_core['mean_density'], 'b-', linewidth=2)
ax1.fill_between(ct_se1_core['depth_norm'],
                 ct_se1_core['mean_density'] - ct_se1_core['std_density'],
                 ct_se1_core['mean_density'] + ct_se1_core['std_density'],
                 alpha=0.3)
ax1.set_xlabel('Depth (mm, normalized)')
ax1.set_ylabel('CT Density (raw units)')
ax1.set_title('A. CT Density - SE1_TOP')
ax1.grid(True, alpha=0.3)

# Panel B: SE2 CT density profile
ax2 = fig.add_subplot(gs[0, 1])
if len(ct_se2_core) > 0:
    ax2.plot(ct_se2_core['depth_norm'], ct_se2_core['mean_density'], 'r-', linewidth=2)
    ax2.fill_between(ct_se2_core['depth_norm'],
                     ct_se2_core['mean_density'] - ct_se2_core['std_density'],
                     ct_se2_core['mean_density'] + ct_se2_core['std_density'],
                     alpha=0.3, color='red')
ax2.set_xlabel('Depth (mm, normalized)')
ax2.set_ylabel('CT Density (raw units)')
ax2.set_title('B. CT Density - SE2_BTTM')
ax2.grid(True, alpha=0.3)

# Panel C: Correlation ranking
ax3 = fig.add_subplot(gs[0, 2])
top_results = results_df.head(15)
colors = ['blue' if s == 'SE1_TOP' else 'red' for s in top_results['ct_series']]
bars = ax3.barh(range(len(top_results)), top_results['correlation'], color=colors, alpha=0.7)
ax3.set_yticks(range(len(top_results)))
ax3.set_yticklabels([f"{r['section'][:15]}" for _, r in top_results.iterrows()], fontsize=8)
ax3.set_xlabel('Correlation (r)')
ax3.set_title('C. Top 15 CT-XRF Correlations')
ax3.axvline(x=0, color='gray', linestyle='--')
ax3.legend([plt.Rectangle((0,0),1,1,fc='blue',alpha=0.7),
            plt.Rectangle((0,0),1,1,fc='red',alpha=0.7)],
           ['SE1_TOP', 'SE2_BTTM'], loc='lower right')
ax3.invert_yaxis()

# Panels D-F: Best SE1 match details
if best_se1 is not None:
    aligned_data = best_se1['aligned_data']
    section = best_se1['section']
    offset = best_se1['offset']

    # Panel D: Aligned profiles
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.plot(ct_se1_core['depth_norm'], ct_se1_core['mean_density'], 'b-', linewidth=2, label='CT Density')
    ax4.set_ylabel('CT Density', color='blue')
    ax4.tick_params(axis='y', labelcolor='blue')

    ax4b = ax4.twinx()
    xrf_section = xrf_data[xrf_data['section'] == section].copy()
    xrf_section['depth_aligned'] = xrf_section['position_mm'] - xrf_section['position_mm'].min() + offset
    ax4b.plot(xrf_section['depth_aligned'], xrf_section['Ca'], 'g-', linewidth=1.5, alpha=0.8, label='XRF Ca')
    ax4b.set_ylabel('XRF Ca (counts)', color='green')
    ax4b.tick_params(axis='y', labelcolor='green')

    ax4.set_xlabel('Depth (mm)')
    ax4.set_title(f'D. SE1 Aligned: {section}\n(offset={offset}mm, r={best_se1["correlation"]:.3f})')
    ax4.grid(True, alpha=0.3)

    # Panel E: Scatter plot
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.scatter(aligned_data['Ca'], aligned_data['ct_density'], alpha=0.7, s=50, c='purple')

    # Regression line
    if len(aligned_data) >= 3:
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            aligned_data['Ca'], aligned_data['ct_density']
        )
        x_line = np.array([aligned_data['Ca'].min(), aligned_data['Ca'].max()])
        ax5.plot(x_line, slope * x_line + intercept, 'r--', linewidth=2)
        ax5.text(0.05, 0.95, f'r = {r_value:.3f}\np = {p_value:.4f}\nn = {len(aligned_data)}',
                transform=ax5.transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax5.set_xlabel('XRF Ca (counts)')
    ax5.set_ylabel('CT Density (interpolated)')
    ax5.set_title('E. SE1 CT-Ca Correlation')
    ax5.grid(True, alpha=0.3)

# Panel F: Multiple element correlations for SE1
ax6 = fig.add_subplot(gs[1, 2])
if best_se1 is not None:
    aligned_data = best_se1['aligned_data']
    elements = ['Ca', 'Ti', 'Fe', 'Sr']
    colors_elem = ['green', 'brown', 'orange', 'purple']
    correlations = []

    for elem in elements:
        if elem in aligned_data.columns:
            corr = aligned_data[elem].corr(aligned_data['ct_density'])
            correlations.append(corr)
        else:
            correlations.append(0)

    ax6.bar(elements, correlations, color=colors_elem, alpha=0.7)
    ax6.axhline(y=0, color='gray', linestyle='--')
    ax6.set_ylabel('Correlation with CT Density')
    ax6.set_title('F. Element-CT Correlations (SE1 Match)')
    ax6.set_ylim(-1, 1)

# Panels G-I: Best SE2 match details
if best_se2 is not None:
    aligned_data = best_se2['aligned_data']
    section = best_se2['section']
    offset = best_se2['offset']

    # Panel G: Aligned profiles
    ax7 = fig.add_subplot(gs[2, 0])
    ax7.plot(ct_se2_core['depth_norm'], ct_se2_core['mean_density'], 'r-', linewidth=2, label='CT Density')
    ax7.set_ylabel('CT Density', color='red')
    ax7.tick_params(axis='y', labelcolor='red')

    ax7b = ax7.twinx()
    xrf_section = xrf_data[xrf_data['section'] == section].copy()
    xrf_section['depth_aligned'] = xrf_section['position_mm'] - xrf_section['position_mm'].min() + offset
    ax7b.plot(xrf_section['depth_aligned'], xrf_section['Ca'], 'g-', linewidth=1.5, alpha=0.8, label='XRF Ca')
    ax7b.set_ylabel('XRF Ca (counts)', color='green')
    ax7b.tick_params(axis='y', labelcolor='green')

    ax7.set_xlabel('Depth (mm)')
    ax7.set_title(f'G. SE2 Aligned: {section}\n(offset={offset}mm, r={best_se2["correlation"]:.3f})')
    ax7.grid(True, alpha=0.3)

    # Panel H: Scatter plot
    ax8 = fig.add_subplot(gs[2, 1])
    ax8.scatter(aligned_data['Ca'], aligned_data['ct_density'], alpha=0.7, s=50, c='darkred')

    if len(aligned_data) >= 3:
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            aligned_data['Ca'], aligned_data['ct_density']
        )
        x_line = np.array([aligned_data['Ca'].min(), aligned_data['Ca'].max()])
        ax8.plot(x_line, slope * x_line + intercept, 'r--', linewidth=2)
        ax8.text(0.05, 0.95, f'r = {r_value:.3f}\np = {p_value:.4f}\nn = {len(aligned_data)}',
                transform=ax8.transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax8.set_xlabel('XRF Ca (counts)')
    ax8.set_ylabel('CT Density (interpolated)')
    ax8.set_title('H. SE2 CT-Ca Correlation')
    ax8.grid(True, alpha=0.3)

# Panel I: Summary statistics
ax9 = fig.add_subplot(gs[2, 2])
ax9.axis('off')

summary_text = "COMPLETE CT-XRF INTEGRATION SUMMARY\n" + "=" * 40 + "\n\n"
summary_text += f"CT Data Coverage:\n"
summary_text += f"  SE1_TOP: {len(ct_se1_core)} valid slices, {ct_se1_core['depth_norm'].max():.0f}mm\n"
summary_text += f"  SE2_BTTM: {len(ct_se2_core)} valid slices, {ct_se2_core['depth_norm'].max():.0f}mm\n\n"

summary_text += f"XRF Sections Tested: {len(sections)}\n"
summary_text += f"Total Combinations: {len(results_df)}\n\n"

summary_text += f"Best Matches:\n"
if best_se1 is not None:
    summary_text += f"  SE1: {best_se1['section']}\n"
    summary_text += f"       r={best_se1['correlation']:.3f}, n={best_se1['n_points']}\n"
if best_se2 is not None:
    summary_text += f"  SE2: {best_se2['section']}\n"
    summary_text += f"       r={best_se2['correlation']:.3f}, n={best_se2['n_points']}\n"

summary_text += f"\nSignificant Correlations (p<0.05):\n"
sig_count = len([r for r in detailed_results if r['significant']])
summary_text += f"  {sig_count} of {len(detailed_results)} best matches"

ax9.text(0.1, 0.9, summary_text, transform=ax9.transAxes, fontsize=10,
         verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.3))

plt.tight_layout()
plt.savefig(FIG_PATH / "ct_xrf_complete_integration.png", dpi=300, bbox_inches='tight')
plt.savefig(FIG_PATH / "ct_xrf_complete_integration.pdf", bbox_inches='tight')
print("Saved: ct_xrf_complete_integration.png/pdf")

# ==============================================================================
# 6. SAVE ALL ALIGNED DATA
# ==============================================================================

print("\n" + "=" * 70)
print("SAVING ALIGNED DATA")
print("=" * 70)

# Save aligned data for best matches
all_aligned = []

for series_name, match, ct_core in best_matches:
    aligned = match['aligned_data'].copy()
    aligned['ct_series'] = series_name
    aligned['best_match_section'] = match['section']
    aligned['offset_mm'] = match['offset']
    all_aligned.append(aligned)

if all_aligned:
    combined_aligned = pd.concat(all_aligned, ignore_index=True)
    # Select key columns
    key_cols = ['ct_series', 'best_match_section', 'offset_mm', 'position_mm',
                'depth_shifted', 'Ca', 'Ti', 'Fe', 'ct_density']
    key_cols = [c for c in key_cols if c in combined_aligned.columns]
    combined_aligned[key_cols].to_csv(OUTPUT_PATH / "tables" / "ct_xrf_complete_aligned.csv", index=False)
    print(f"Saved: ct_xrf_complete_aligned.csv ({len(combined_aligned)} aligned points)")

# ==============================================================================
# 7. FINAL SUMMARY
# ==============================================================================

print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

print(f"""
CT Data:
  SE1_TOP: {len(ct_se1_core)} valid slices, {ct_se1_core['depth_norm'].max():.0f} mm coverage
  SE2_BTTM: {len(ct_se2_core)} valid slices, {ct_se2_core['depth_norm'].max():.0f} mm coverage

XRF Alignment Results:
  Sections tested: {len(sections)}
  Combinations analyzed: {len(results_df)}

Best Matches:
""")

for dr in detailed_results:
    sig_str = "SIGNIFICANT" if dr['significant'] else "not significant"
    print(f"  {dr['ct_series']}: {dr['section']}")
    print(f"    Correlation: r = {dr['correlation']:.4f}")
    print(f"    P-value: {dr['p_value']:.6f} ({sig_str})")
    print(f"    N points: {dr['n_points']}, Overlap: {dr['overlap_mm']:.1f} mm")
    print()

print("""
Output Files:
  - ct_xrf_all_correlations.csv: All tested combinations
  - ct_xrf_best_matches.csv: Best match for each CT series
  - ct_xrf_complete_aligned.csv: Aligned data for best matches
  - ct_xrf_complete_integration.png/pdf: Comprehensive visualization
""")
