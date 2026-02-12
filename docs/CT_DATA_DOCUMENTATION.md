# Pebas Formation CT Core Scanning: Complete Data Documentation

## Overview

This document archives all computed tomography (CT) data and analysis results for sediment core material from the Pebas Formation, Peru. The CT dataset is documented separately from the main XRF manuscript because **the specific core section identity could not be confirmed**, precluding direct integration with XRF measurements.

---

## 1. Data Acquisition

### Scanner Information
- **Scanner**: Toshiba Aquilion clinical CT
- **Reconstruction kernel**: FC30 (Bone 2.0)
- **Orientation**: Coronal
- **Scan date**: January 23, 2026 (inferred from file timestamps)

### DICOM Metadata
| Field | SE1_TOP | SE2_BTTM |
|-------|---------|----------|
| PatientID | 1 23 26 | 1 23 26 |
| StudyID | 22166 | 22166 |
| SeriesDescription | Bone 2.0 Coron TOP.Ref | Bone 2.0 Coron BTTM.Ref |
| Slices | 38 | 42 |
| Slice Thickness | 2.0 mm | 2.0 mm |
| Image Matrix | 1624 × 1624 | 1304 × 1304 |
| Pixel Spacing | 0.616 mm | 0.468 mm |
| SliceLocation Range | -39.5 to +34.5 mm | -48.4 to +33.6 mm |

### Sample Identification Issue
The PatientID "1 23 26" and StudyID "22166" do not correspond to any known XRF core section identifier. Attempts to match CT density profiles to XRF chemical profiles were unsuccessful (see Section 5).

---

## 2. Data Processing Pipeline

### 2.1 DICOM to TIFF Conversion
**Script**: `scripts/24_convert_dicom_to_tiff.R`

Converted 80 DICOM files to 16-bit TIFF format:
- SE1_TOP: 38 TIFF files (~5 MB each, 1624×1624 px)
- SE2_BTTM: 42 TIFF files (~3 MB each, 1304×1304 px)

**Output location**: `Pebas-CT/TIFF/`

**Naming convention**: `{series}_slice{###}_loc{±##.##}mm.tif`

### 2.2 Density Profile Extraction
**Script**: `scripts/25_ct_tiff_analysis.py`

Extracted per-slice statistics:
- Mean density (full image)
- Core mean density (central circular ROI)
- Standard deviation
- Min/max values

**Output**: `output/tables/ct_tiff_density_profiles.csv`

### 2.3 CT-XRF Matching Analysis
**Script**: `scripts/26_ct_xrf_matching.py`

Attempted to identify CT sample by correlating density profiles with XRF Ca profiles from candidate sections.

**Output**: `output/tables/ct_xrf_matching_results.csv`

---

## 3. CT Data Summary

### SE1_TOP Series
| Parameter | Value |
|-----------|-------|
| Total slices | 38 |
| Coverage | 74 mm |
| Mean density | 0.221 (normalized) |
| Density range | 0.10 - 0.42 |
| Pattern | Bimodal - peaks at -5mm and +25mm |

### SE2_BTTM Series
| Parameter | Value |
|-----------|-------|
| Total slices | 42 |
| Coverage | 82 mm |
| Mean density | 0.301 (normalized) |
| Density range | 0.12 - 0.76 |
| Pattern | Layered - clear banding visible |

### Density Profile Characteristics
Both CT series show:
1. Continuous core material (no air gaps in normalized data)
2. Density variations suggesting lithological layering
3. Overlapping SliceLocation ranges, indicating same physical location scanned with different parameters

---

## 4. Visualization Products

### Generated Figures
| Figure | Description | Location |
|--------|-------------|----------|
| ct_density_profiles.png | Density bar charts for both series | output/figures/ct_tiff_analysis/ |
| ct_montage_SE1_TOP.png | All 38 SE1 slices | output/figures/ct_tiff_analysis/ |
| ct_montage_SE2_BTTM.png | All 42 SE2 slices | output/figures/ct_tiff_analysis/ |
| ct_mip_SE1_TOP.png | Maximum intensity projections SE1 | output/figures/ct_tiff_analysis/ |
| ct_mip_SE2_BTTM.png | Maximum intensity projections SE2 | output/figures/ct_tiff_analysis/ |
| correlation_summary.png | CT-XRF correlation results | output/figures/ct_xrf_matching/ |

### Slice Montage Description
The montage images show:
- Core material as bright vertical/elongated features
- Black regions = air/background outside core
- Density banding visible as horizontal striations
- Two-segment structure in some views

### MIP (Maximum Intensity Projection) Description
- **Axial MIP**: Top-down view showing core cross-section shape
- **Sagittal MIP**: Side view revealing density variation along slices
- **Coronal MIP**: Front view showing horizontal layering

---

## 5. CT-XRF Matching Analysis Results

### Candidate Sections Tested
Sections with XRF coverage length similar to CT (60-150 mm):
1. TAM-1-2-3B-C (84 mm)
2. SC-5-6-7ABC-E (90 mm)
3. SC-5-6-7ABC-C (99 mm)
4. SC-5-6-7ABC-D (105 mm)
5. TAM-5AB-6-7-A (114 mm)
6. SC-3AB-4ABCD-C (123 mm)
7. SC-3AB-4ABCD-B (129 mm)

### Correlation Results

#### SE1_TOP Correlations with XRF Ca
| Section | r | p-value | Overlap |
|---------|---|---------|---------|
| TAM-5AB-6-7-A | -0.430 | 0.022* | 74 mm |
| TAM-1-2-3B-C | -0.312 | 0.129 | 74 mm |
| SC-3AB-4ABCD-C | +0.284 | 0.084 | 74 mm |
| SC-5-6-7ABC-C | +0.153 | 0.438 | 74 mm |
| SC-5-6-7ABC-E | +0.125 | 0.621 | 74 mm |
| SC-5-6-7ABC-D | +0.099 | 0.661 | 74 mm |
| SC-3AB-4ABCD-B | +0.060 | 0.790 | 74 mm |

#### SE2_BTTM Correlations with XRF Ca
| Section | r | p-value | Overlap |
|---------|---|---------|---------|
| SC-3AB-4ABCD-B | -0.377 | 0.084 | 82 mm |
| TAM-1-2-3B-C | -0.294 | 0.154 | 82 mm |
| TAM-5AB-6-7-A | -0.267 | 0.170 | 82 mm |
| SC-5-6-7ABC-C | +0.185 | 0.346 | 82 mm |
| SC-5-6-7ABC-E | -0.169 | 0.501 | 82 mm |
| SC-3AB-4ABCD-C | +0.158 | 0.345 | 82 mm |
| SC-5-6-7ABC-D | -0.086 | 0.703 | 82 mm |

### Key Finding
**No valid match identified.** The only statistically significant correlation (SE1_TOP vs TAM-5AB-6-7-A, r = -0.430, p = 0.022) is **negative**, which contradicts the expected positive relationship between CT density and carbonate content.

### Interpretation
1. The CT sample is likely from a section NOT included in the XRF dataset
2. Alternatively, the depth coordinate systems are incompatible
3. CT density may be influenced by factors other than carbonate content

---

## 6. File Inventory

### Raw Data
```
Pebas-CT/
├── ST1/
│   ├── SE1/          # 39 DICOM files (IM1-IM39)
│   └── SE2/          # 43 DICOM files (IM1-IM43)
├── ST1.zip           # Original archive
└── TIFF/
    ├── SE1_TOP/      # 38 converted TIFF files
    └── SE2_BTTM/     # 42 converted TIFF files
```

### Processed Data
```
output/tables/
├── ct_dicom_metadata.csv       # DICOM header information
├── ct_series_summary.csv       # Series-level summary
├── ct_density_profiles.csv     # Original DICOM density analysis
├── ct_tiff_density_profiles.csv # TIFF-based density profiles
└── ct_xrf_matching_results.csv # Correlation analysis results
```

### Analysis Scripts
```
scripts/
├── 16_ct_integration.R           # Initial DICOM metadata extraction
├── 24_convert_dicom_to_tiff.R    # DICOM to TIFF conversion
├── 25_ct_tiff_analysis.py        # TIFF density analysis
└── 26_ct_xrf_matching.py         # CT-XRF correlation analysis
```

### Figures
```
output/figures/
├── ct_tiff_analysis/
│   ├── ct_density_profiles.png
│   ├── ct_montage_SE1_TOP.png
│   ├── ct_montage_SE2_BTTM.png
│   ├── ct_mip_SE1_TOP.png
│   └── ct_mip_SE2_BTTM.png
└── ct_xrf_matching/
    ├── correlation_summary.png
    ├── comparison_SE1_TOP_TAM-5AB-6-7-A.png
    ├── comparison_SE1_TOP_TAM-1-2-3B-C.png
    └── comparison_SE2_BTTM_SC-3AB-4ABCD-B.png
```

---

## 7. Recommendations for Future Work

### To Enable CT-XRF Integration
1. **Contact CT facility**: Obtain records for StudyID 22166 to identify sample provenance
2. **Decode PatientID**: Determine if "1 23 26" encodes sample information
3. **Physical inspection**: Compare CT slice images with core photographs for visual matching
4. **New CT acquisitions**: If integration is essential, CT scan identified core sections with clear labeling

### Alternative Uses of CT Data
Even without XRF integration, the CT data provides:
- 3D density structure of ~75-80 mm of Pebas Formation sediment
- Internal layering and fabric information
- Potential for porosity/compaction analysis
- Reference dataset for future CT-sediment studies

---

## 8. Data Availability

The CT DICOM data, converted TIFF files, and all analysis outputs are archived in the project repository. Raw DICOM files are available upon request for researchers who wish to apply alternative analysis methods.

**Repository**: https://github.com/biopixl/Pebas-XRF

**Contact**: [Author contact information]

---

## Document History

| Date | Version | Changes |
|------|---------|---------|
| 2026-02-12 | 1.0 | Initial documentation created |

---

*This documentation was generated as part of the Pebas-XRF project to archive CT data that could not be integrated with XRF measurements due to unknown sample provenance.*
