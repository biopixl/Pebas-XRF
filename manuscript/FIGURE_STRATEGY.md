# Manuscript Figure Strategy and Synthesis Recommendations

## Current Figure Inventory

### Main Text Figures (10 total)
| Figure | File | Description | Status |
|--------|------|-------------|--------|
| 1 | integrated_GROUP1_pub.png | TAM GROUP1 stratigraphy | OK |
| 2 | integrated_GROUP2_pub.png | TAM GROUP2 stratigraphy | OK |
| 3 | integrated_GROUP3_pub.png | TAM GROUP3 stratigraphy | OK |
| 4 | integrated_GROUP4_pub.png | SC GROUP4 stratigraphy | OK |
| 5 | integrated_GROUP5_pub.png | SC GROUP5 stratigraphy | OK |
| 6 | integrated_GROUP6_pub.png | SC GROUP6 stratigraphy | OK |
| 7 | integrated_GROUP7_pub.png | SC GROUP7 stratigraphy | OK |
| 8 | stacked_TAM.png | TAM composite column | OK |
| 9 | stacked_SC.png | SC composite column | OK |
| 10 | ct_xrf_complete_integration.png | CT-XRF validation | OK |

### Proxy Evaluation Section Figures (6 total)
| Figure | File | Status |
|--------|------|--------|
| S1 | fig_S1_redundancy_analysis.png | OK |
| S2 | fig_S2_element_variability.png | OK |
| S3 | fig_S3_correlation_matrix.png | OK |
| S4 | fig_S4_crossplots.png | OK |
| S5 | fig_S5_distributions.png | OK |
| Main | fig_main_stratigraphic_example.png | OK |

### Filter Verification Section Figures (7 total)
| Figure | File | Status |
|--------|------|--------|
| V1 | fig_V1_foam_spectral_signatures.png | **MISSING** |
| V2 | fig_V2_inc_coh_foam_detection.png | **MISSING** |
| V3 | fig_V3_facies_threshold_validation.png | **MISSING** |
| V4 | fig_V4_redox_threshold_validation.png | **MISSING** |
| V5 | fig_V5_window_size_optimization.png | **MISSING** |
| V6 | fig_V6_qc_filter_effectiveness.png | **MISSING** |
| V7 | fig_V7_alignment_verification.png | **MISSING** |

---

## Critical Issue: Missing Filter Validation Figures

The `filter_verification.tex` section references 7 figures in `../figures/filter_validation/` that do not exist. These must be generated before the manuscript will compile correctly.

**Required figures:**
1. Foam spectral signatures (boxplots comparing sediment vs foam zones)
2. Inc/Coh foam detection criterion (density plot with threshold lines)
3. Facies threshold validation (Ca/Ti distribution with threshold markers)
4. Redox threshold validation (Fe/Mn distribution and cross-plot)
5. Window size optimization (smoothing comparison panel)
6. QC filter effectiveness (MSE, CPS, surface distance distributions)
7. Alignment verification (overlay of Ca/Ti and Fe/Mn profiles)

---

## Strategic Recommendations for Synthesis Communication

### 1. Add Study Area Map (High Priority)
**Purpose:** Orient readers to Pebas Formation context and locality positions
**Content:**
- Regional map showing Western Amazonia
- Inset with TAM and SC locality positions
- Geological context (Pebas Formation extent)
- Coordinate markers

**Literature support:** Standard for sediment core studies (Croudace & Rothwell 2015)

### 2. Add Conceptual Synthesis Diagram (High Priority)
**Purpose:** Communicate multi-proxy integration visually
**Suggested content:**
- Schematic cross-section of Pebas mega-wetland
- Proxy indicators mapped to environmental signals
- TAM vs SC depositional setting comparison
- CT-XRF validation pathway illustrated

**Placement:** Before or within Discussion section

### 3. Add Site Comparison Summary Figure (Medium Priority)
**Purpose:** Synthesize TAM vs SC differences in single visualization
**Content:**
- Side-by-side box plots of key proxies (Ca, Ti, Fe/Mn, Zr/Rb)
- Facies pie charts for each site
- Redox condition comparison
- Statistical significance indicators

**Available data:** Can be generated from existing summary statistics

### 4. Consider PCA Biplot Inclusion (Medium Priority)
**Available files:**
- pca_biplot.png
- pca_scree.png
- pca_analysis.png

**Value:** Demonstrates geochemical associations between elements and validates proxy suite selection. Shows carbonate vs detrital element clustering.

**Placement:** Could enhance proxy_evaluation section

### 5. Add Stratigraphic Correlation Panel (Low Priority)
**Purpose:** Show inter-core correlations and stratigraphic equivalence
**Content:**
- TAM and SC columns side-by-side
- Potential tie-lines between equivalent horizons
- Shell-rich interval markers
- CT-validated intervals highlighted

---

## Recommended Figure Order for Publication

### Main Text (suggested 12-14 figures)
1. **Study Area Map** (NEW - to create)
2. Integrated GROUP1 (TAM)
3. Integrated GROUP2 (TAM)
4. Integrated GROUP3 (TAM)
5. Integrated GROUP4 (SC)
6. Integrated GROUP5 (SC)
7. Integrated GROUP6 (SC)
8. Integrated GROUP7 (SC)
9. Stacked TAM composite
10. Stacked SC composite
11. CT-XRF integration (9-panel)
12. **Site Comparison Summary** (NEW - to create)
13. **Conceptual Synthesis Diagram** (NEW - to create)

### Supplementary/Supporting Figures
- S1-S5: Proxy evaluation figures (existing)
- S6: PCA biplot (existing, add to supplement)
- V1-V7: Filter validation figures (TO CREATE)

---

## Literature-Aligned Figure Practices

Based on references.bib citations:

### From Croudace & Rothwell (2015) - XRF Core Scanning
- Multi-proxy stratigraphic panels with optical images
- Element ratio profiles with threshold annotations
- QC documentation with filter effectiveness metrics

### From Wesselingh et al. (2002) - Pebas Formation
- Paleogeographic reconstruction maps
- Facies distribution diagrams
- Mollusk assemblage zones

### From Vonhof et al. (2003) - Isotope Studies
- Cross-plots for proxy validation (Ca-Sr, Fe-Mn)
- Site comparison statistics

### From Davies et al. (2015) - Palaeolimnology
- Composite stratigraphic columns with multiple proxies
- PCA/ordination diagrams for geochemical data

---

## Action Items

### Immediate (Required for Compilation)
1. [ ] Create filter_validation figure directory
2. [ ] Generate 7 missing filter validation figures (V1-V7)

### Short-term (Enhance Synthesis)
3. [ ] Create study area/locality map
4. [ ] Create site comparison summary figure
5. [ ] Consider adding PCA biplot to supplement

### Optional (Further Enhancement)
6. [ ] Create conceptual synthesis diagram
7. [ ] Create stratigraphic correlation panel
8. [ ] Update figure captions with consistent formatting

---

## Technical Notes

### Figure Paths in LaTeX
- Main figures: `figures/*.png`
- Proxy evaluation: `figures/proxy_evaluation/*.png`
- Filter validation: `../figures/filter_validation/*.png` (relative path from sections/)

### Recommended Image Specifications
- Resolution: 300 DPI minimum
- Width: Full page (textwidth) for composite figures
- Format: PNG or PDF for vector graphics
- Color scheme: Consistent palette (Set2 for groups, blue/orange for TAM/SC)

---

*Generated: 2026-02-06*
*Based on manuscript review and literature analysis*
