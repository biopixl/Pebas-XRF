#!/usr/bin/env python3
"""
CT-XRF Section Matching Analysis
================================
Attempts to identify which XRF section corresponds to the CT scan
by comparing density profiles with elemental chemistry.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from PIL import Image
import matplotlib.pyplot as plt
from scipy import stats
from scipy.interpolate import interp1d

# ==============================================================================
# CONFIGURATION
# ==============================================================================

BASE_PATH = Path("/Users/isaac/Documents/GitHub/PhD_Projects/Pebas-XRF")
TIFF_PATH = BASE_PATH / "Pebas-CT" / "TIFF"
OUTPUT_PATH = BASE_PATH / "output"
FIG_PATH = OUTPUT_PATH / "figures" / "ct_xrf_matching"

FIG_PATH.mkdir(parents=True, exist_ok=True)

# Candidate sections (length 60-150mm, similar to CT coverage)
CANDIDATE_SECTIONS = [
    'TAM-1-2-3B-C',      # 84 mm - closest to CT
    'SC-5-6-7ABC-E',     # 90 mm
    'SC-5-6-7ABC-C',     # 99 mm
    'SC-5-6-7ABC-D',     # 105 mm
    'TAM-5AB-6-7-A',     # 114 mm
    'SC-3AB-4ABCD-C',    # 123 mm
    'SC-3AB-4ABCD-B',    # 129 mm
]

# ==============================================================================
# 1. LOAD CT DENSITY DATA
# ==============================================================================

def load_ct_profiles():
    """Load CT density profiles from TIFF analysis."""
    ct_df = pd.read_csv(OUTPUT_PATH / "tables" / "ct_tiff_density_profiles.csv")

    profiles = {}
    for series in ['SE1_TOP', 'SE2_BTTM']:
        df = ct_df[ct_df['series'] == series].copy()
        df = df.sort_values('slice_location')

        # Normalize depth to start at 0
        df['depth_norm'] = df['slice_location'] - df['slice_location'].min()

        profiles[series] = df

    return profiles

# ==============================================================================
# 2. LOAD XRF DATA FOR CANDIDATES
# ==============================================================================

def load_xrf_candidates():
    """Load XRF data for candidate sections."""
    xrf_df = pd.read_csv(OUTPUT_PATH / "tables" / "xrf_data_stacked.csv")

    candidates = {}
    for section in CANDIDATE_SECTIONS:
        df = xrf_df[xrf_df['section'] == section].copy()
        if len(df) > 0:
            df = df.sort_values('position_mm')

            # Normalize depth to start at 0
            df['depth_norm'] = df['position_mm'] - df['position_mm'].min()

            candidates[section] = df

    return candidates

# ==============================================================================
# 3. EXTRACT OPTICAL LINE PROFILES
# ==============================================================================

def extract_optical_profile(section_name, xrf_df):
    """Extract optical intensity profile from core image for XRF-covered region."""
    # Find optical image
    optical_paths = list(BASE_PATH.glob(f"TAM-SC-IsaacA/*/{section_name.split('-')[0]}*/{section_name}/optical.tif"))
    if not optical_paths:
        optical_paths = list(BASE_PATH.glob(f"TAM-SC-IsaacA/*/*/{section_name}/optical.tif"))

    if not optical_paths:
        return None

    # Use first match (skip 'copied' directories)
    optical_path = None
    for p in optical_paths:
        if 'copied' not in str(p):
            optical_path = p
            break

    if optical_path is None:
        return None

    try:
        img = Image.open(optical_path)
        img_array = np.array(img)

        # Convert to grayscale if RGB
        if len(img_array.shape) == 3:
            img_gray = np.mean(img_array, axis=2)
        else:
            img_gray = img_array

        # Get XRF position range
        start_mm = xrf_df['position_mm'].min()
        end_mm = xrf_df['position_mm'].max()

        # Convert to pixel coordinates (0.2mm per pixel for Itrax, along width)
        # Note: Itrax images have core along width (x-axis)
        start_px = int(start_mm / 0.2)
        end_px = int(end_mm / 0.2)

        # Ensure within bounds
        start_px = max(0, start_px)
        end_px = min(img_gray.shape[1], end_px)

        # Extract mean intensity along core axis
        profile = np.mean(img_gray[:, start_px:end_px], axis=0)

        # Create position array
        positions = np.linspace(start_mm, end_mm, len(profile))

        return pd.DataFrame({
            'position_mm': positions,
            'optical_intensity': profile,
            'depth_norm': positions - positions.min()
        })

    except Exception as e:
        print(f"  Error loading optical for {section_name}: {e}")
        return None

# ==============================================================================
# 4. COMPUTE CORRELATIONS
# ==============================================================================

def compute_correlation(ct_profile, xrf_profile, element='Ca'):
    """
    Compute correlation between CT density and XRF element.
    Interpolates to common depth scale.
    """
    # Get depth ranges
    ct_depth = ct_profile['depth_norm'].values
    ct_density = ct_profile['core_mean'].values

    xrf_depth = xrf_profile['depth_norm'].values
    xrf_values = xrf_profile[element].values

    # Find overlapping range
    max_depth = min(ct_depth.max(), xrf_depth.max())

    if max_depth < 10:  # Need at least 10mm overlap
        return None, None, 0

    # Create common depth grid
    n_points = min(len(ct_depth), len(xrf_depth))
    common_depth = np.linspace(0, max_depth, n_points)

    # Interpolate both to common grid
    try:
        ct_interp = interp1d(ct_depth, ct_density, kind='linear', fill_value='extrapolate')
        xrf_interp = interp1d(xrf_depth, xrf_values, kind='linear', fill_value='extrapolate')

        ct_common = ct_interp(common_depth)
        xrf_common = xrf_interp(common_depth)

        # Compute correlation
        r, p = stats.pearsonr(ct_common, xrf_common)

        return r, p, max_depth

    except Exception as e:
        return None, None, 0

# ==============================================================================
# 5. VISUALIZATION
# ==============================================================================

def plot_comparison(ct_series, ct_profile, section_name, xrf_profile, optical_profile=None):
    """Create comparison plot between CT and XRF/optical."""
    n_panels = 3 if optical_profile is not None else 2
    fig, axes = plt.subplots(n_panels, 1, figsize=(12, 3*n_panels), sharex=True)

    # CT density profile
    ax1 = axes[0]
    ax1.plot(ct_profile['depth_norm'], ct_profile['core_mean'], 'b-', linewidth=1.5, label='CT Density')
    ax1.fill_between(ct_profile['depth_norm'], ct_profile['core_mean'], alpha=0.3)
    ax1.set_ylabel('CT Density\n(normalized)')
    ax1.set_title(f'{ct_series} vs {section_name}')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)

    # XRF Ca profile
    ax2 = axes[1]
    ax2.plot(xrf_profile['depth_norm'], xrf_profile['Ca'] / 1000, 'g-', linewidth=1.5, label='Ca (kcps)')
    ax2.fill_between(xrf_profile['depth_norm'], xrf_profile['Ca'] / 1000, alpha=0.3, color='green')
    ax2.set_ylabel('Ca\n(kcps)')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)

    # Optical profile if available
    if optical_profile is not None:
        ax3 = axes[2]
        ax3.plot(optical_profile['depth_norm'], optical_profile['optical_intensity'], 'gray', linewidth=1, label='Optical')
        ax3.set_ylabel('Optical\nIntensity')
        ax3.legend(loc='upper right')
        ax3.grid(True, alpha=0.3)

    axes[-1].set_xlabel('Depth (mm, normalized)')

    plt.tight_layout()
    return fig

def create_correlation_summary(results):
    """Create summary plot of all correlations."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Separate by CT series
    for i, ct_series in enumerate(['SE1_TOP', 'SE2_BTTM']):
        ax = axes[i]

        series_results = [r for r in results if r['ct_series'] == ct_series]
        if not series_results:
            continue

        sections = [r['section'] for r in series_results]
        correlations = [r['correlation'] for r in series_results]
        overlaps = [r['overlap'] for r in series_results]

        # Color by significance
        colors = ['green' if r['p_value'] and r['p_value'] < 0.05 else 'gray' for r in series_results]

        bars = ax.barh(sections, correlations, color=colors, edgecolor='black')

        # Add overlap annotation
        for j, (bar, overlap) in enumerate(zip(bars, overlaps)):
            ax.text(bar.get_width() + 0.02, bar.get_y() + bar.get_height()/2,
                   f'{overlap:.0f}mm', va='center', fontsize=9)

        ax.set_xlabel('Correlation with Ca')
        ax.set_title(f'{ct_series} Correlations')
        ax.axvline(x=0, color='black', linewidth=0.5)
        ax.set_xlim(-1, 1.2)
        ax.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()
    return fig

# ==============================================================================
# 6. MAIN ANALYSIS
# ==============================================================================

def main():
    print("=" * 70)
    print("CT-XRF SECTION MATCHING ANALYSIS")
    print("=" * 70)

    # Load data
    print("\nLoading CT density profiles...")
    ct_profiles = load_ct_profiles()

    print("Loading XRF candidate sections...")
    xrf_candidates = load_xrf_candidates()

    print(f"  Loaded {len(xrf_candidates)} candidate sections")

    # Compute correlations
    print("\n" + "=" * 70)
    print("CORRELATION ANALYSIS")
    print("=" * 70)

    results = []

    for ct_series, ct_profile in ct_profiles.items():
        print(f"\n{ct_series}:")
        ct_length = ct_profile['depth_norm'].max()
        print(f"  CT coverage: {ct_length:.1f} mm")

        for section, xrf_df in xrf_candidates.items():
            xrf_length = xrf_df['depth_norm'].max()

            # Compute correlation with Ca
            r, p, overlap = compute_correlation(ct_profile, xrf_df, 'Ca')

            if r is not None:
                sig = "*" if p < 0.05 else ""
                print(f"  {section}: r={r:+.3f}{sig} (p={p:.4f}, overlap={overlap:.0f}mm)")

                results.append({
                    'ct_series': ct_series,
                    'section': section,
                    'correlation': r,
                    'p_value': p,
                    'overlap': overlap,
                    'xrf_length': xrf_length
                })

                # Create comparison plot for significant correlations
                if abs(r) > 0.3:
                    optical = extract_optical_profile(section, xrf_df)
                    fig = plot_comparison(ct_series, ct_profile, section, xrf_df, optical)
                    fig.savefig(FIG_PATH / f'comparison_{ct_series}_{section}.png', dpi=150, bbox_inches='tight')
                    plt.close(fig)

    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(OUTPUT_PATH / "tables" / "ct_xrf_matching_results.csv", index=False)
    print(f"\nSaved: {OUTPUT_PATH / 'tables' / 'ct_xrf_matching_results.csv'}")

    # Create summary plot
    fig = create_correlation_summary(results)
    fig.savefig(FIG_PATH / 'correlation_summary.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {FIG_PATH / 'correlation_summary.png'}")

    # Report best matches
    print("\n" + "=" * 70)
    print("BEST MATCHES")
    print("=" * 70)

    for ct_series in ['SE1_TOP', 'SE2_BTTM']:
        series_results = [r for r in results if r['ct_series'] == ct_series]
        if series_results:
            best = max(series_results, key=lambda x: abs(x['correlation']))
            print(f"\n{ct_series} best match: {best['section']}")
            print(f"  Correlation: r = {best['correlation']:.3f}")
            print(f"  P-value: {best['p_value']:.4f}")
            print(f"  Overlap: {best['overlap']:.0f} mm")

if __name__ == "__main__":
    main()
