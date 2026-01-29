# XRF Proxy Evaluation for Pebas Formation Sediments

## Miocene Western Amazonian Lacustrine-Wetland System

---

## Executive Summary

This report evaluates XRF elemental proxies for paleoenvironmental reconstruction of Middle-Late Miocene (c. 17-10 Ma) sediments from the Pebas Mega-Wetland System. Evaluation criteria include: (1) detection quality with Mo tube, (2) relevance to tropical freshwater setting, (3) empirical utility in this dataset.

**Key Finding:** Several commonly-used ratios (especially Ca/Ti) are statistically redundant with raw element counts in this dataset. We recommend a streamlined proxy suite based on empirical analysis.

---

## 1. Paleoenvironmental Context

| Parameter | Pebas Formation Characteristics |
|-----------|--------------------------------|
| Age | Middle-Late Miocene (c. 17-10 Ma) |
| Setting | Freshwater mega-wetland, >1 million km² |
| Salinity | Predominantly freshwater (⁸⁷Sr/⁸⁶Sr: 0.7081-0.7085); rare oligohaline episodes |
| Climate | Tropical humid; intense chemical weathering (CIA 80-100) |
| Drainage | Dominantly Andean-sourced |

---

## 2. Empirical Proxy Evaluation

### Redundancy Analysis

Ratios where r(element, ratio) > 0.85 provide no additional information beyond the raw element.

| Ratio | Correlation with Numerator | Assessment |
|-------|---------------------------|------------|
| Ca/Ti | r = 0.89 | **REDUNDANT** - Use raw Ca instead |
| Ba/Ti | r = 0.75 | Borderline - consider raw Ba |
| Zr/Rb | r = 0.60 | **USEFUL** - ratio adds information |
| Si/Ti | r = 0.56 | **USEFUL** - ratio adds information |
| Fe/Ti | r = 0.48 | **USEFUL** - ratio adds information |
| K/Ti | r = 0.39 | **USEFUL** - ratio adds information |
| Fe/Mn | r = 0.25 | **USEFUL** - ratio adds information |
| Rb/Sr | r = -0.17 | **USEFUL** - ratio adds information |

### Element Variability

| Element | CV (%) | Detection | Notes |
|---------|--------|-----------|-------|
| Ca | 81% | Excellent | Highest variance; carbonate signal |
| Fe | 46% | Excellent | Moderate variance; detrital + redox |
| Zr | 43% | Excellent | Heavy mineral indicator |
| K | 37% | Good | Clay/weathering indicator |
| Ti | 36% | Excellent | Stable detrital normalizer |
| Al | 36% | Marginal | Mo tube detection limited |
| Rb | 29% | Excellent | Clay mineral indicator |

---

## 3. Recommended Proxy Suite for Manuscript

### Tier 1: Primary Proxies (Plot in Main Figures)

| Proxy | Type | Interpretation | Justification |
|-------|------|----------------|---------------|
| **Ca** | Element | Authigenic/biogenic carbonate | High variance (CV=81%); Ca/Ti redundant |
| **Ti** | Element | Terrigenous detrital flux | Stable conservative element; excellent detection |
| **Fe/Mn** | Ratio | Redox conditions | r=0.25 with Fe; adds unique information |
| **Zr/Rb** | Ratio | Grain size / hydrodynamic energy | r=0.60; robust grain size proxy |

### Tier 2: Supporting Proxies (Plot in Supplementary or Select Figures)

| Proxy | Type | Interpretation | Justification |
|-------|------|----------------|---------------|
| **Fe** | Element | Lateritic input / catchment weathering | Independent of Ti normalization issues |
| **Sr** | Element | Carbonate mineralogy / salinity indicator | Track aragonite vs calcite; incursion events |
| **K/Ti** | Ratio | Relative weathering intensity | r=0.39; useful despite limited tropical relevance |
| **Rb/Sr** | Ratio | Carbonate influence (inverse) | r=-0.17; tracks Sr enrichment in carbonates |

### Tier 3: Contextual (Report but Don't Feature)

| Proxy | Notes |
|-------|-------|
| Inc/Coh | Matrix composition; organic matter indicator |
| Ba | Potential salinity indicator (not productivity in freshwater) |
| Mn | Context for Fe/Mn interpretation |

### Not Recommended

| Proxy | Reason |
|-------|--------|
| Ca/Ti | Statistically redundant with Ca (r=0.89) |
| Si/Al | Poor Mo tube detection for both elements |
| CIA | Requires reliable Al; tropical setting at weathering ceiling |

---

## 4. Proxy Mechanisms for Pebas Formation

### Ca (Carbonate Signal)

**Mechanism:** Authigenic carbonate precipitation in shallow lakes during evaporative/dry periods. Ostracod and mollusc shell material. Possibly detrital carbonate from catchment.

**Interpretation:**
- High Ca → Carbonate-dominated sedimentation; reduced clastic input; drier/evaporative conditions
- Low Ca → Terrigenous-dominated sedimentation; enhanced runoff; wetter conditions

**Caveats:** Cannot distinguish authigenic vs detrital carbonate without Sr/Ca or isotopic analysis.

---

### Ti (Terrigenous Flux)

**Mechanism:** Ti is hosted in resistant heavy minerals (ilmenite, rutile, titanite) derived from Andean and cratonic sources. Immobile during weathering and diagenesis.

**Interpretation:**
- High Ti → Enhanced terrigenous sediment delivery; higher runoff/discharge
- Low Ti → Reduced clastic input; carbonate or organic-dominated intervals

**Caveats:** Ti concentration affected by dilution from carbonates or organics.

---

### Fe/Mn (Redox Proxy)

**Mechanism:** Under reducing (anoxic) conditions at the sediment-water interface, Mn²⁺ is more soluble than Fe²⁺ and diffuses upward out of sediments.

**Interpretation:**
- High Fe/Mn → More reducing bottom water conditions
- Low Fe/Mn → More oxic conditions; Mn retained in sediments

**Critical Caveats for Pebas:**
1. Detrital Fe from lateritic catchment can obscure redox signal
2. Under permanent anoxia, Mn may be trapped as authigenic carbonates
3. Lake-specific; do not over-interpret small variations
4. Best used to identify major redox transitions, not subtle changes

---

### Zr/Rb (Grain Size Proxy)

**Mechanism:** Zr concentrates in heavy minerals (zircon) in coarse fraction; Rb substitutes for K in clay minerals in fine fraction.

**Interpretation:**
- High Zr/Rb → Coarser sediment; higher energy depositional environment
- Low Zr/Rb → Finer sediment; lower energy; suspended load dominance

**Caveats:**
1. Provenance effects possible (Andean vs cratonic sources have different Zr/Rb)
2. Heavy mineral concentration (placer effects) can elevate Zr independently
3. Use ln(Zr/Rb) for improved linearity with grain size

---

### Sr (Salinity/Carbonate Mineralogy)

**Mechanism:** Sr substitutes readily into aragonite (high Sr/Ca) but less into calcite (low Sr/Ca). Marine/brackish waters have higher Sr.

**Interpretation:**
- High Sr → Aragonite-rich intervals (molluscs, ostracods); possible brackish influence
- Low Sr → Calcite-dominated; freshwater authigenic carbonate

**Application:** Use Sr alongside Ca to discriminate carbonate types and detect oligohaline incursion events documented in Pebas Formation.

---

### K/Ti (Weathering/Clay Composition)

**Mechanism:** K hosted in illite, muscovite, K-feldspar; depleted during intense chemical weathering. Ti is conservative.

**Interpretation:**
- High K/Ti → Less weathered material; illite-rich clays; possibly different provenance
- Low K/Ti → Intensely weathered material; kaolinite-dominated

**Limited Relevance:** Miocene Amazonian catchment was already at weathering ceiling (CIA 80-100). K/Ti may show minimal dynamic range. Most useful for detecting anomalous K-enriched layers (volcanic ash, feldspar-rich inputs).

---

## 5. Recommended Figure Structure for Manuscript

### Main Text Figure: Composite Stratigraphic Profile

**Panel Layout (left to right):**

1. **Core Image** - Optical scan with depth scale
2. **Ca** (log scale) - Carbonate/authigenic signal
3. **Ti** - Terrigenous flux
4. **Fe/Mn** (log scale) - Redox indicator
5. **Zr/Rb** (log scale) - Grain size proxy
6. **Facies Interpretation** - Summary column

**Rationale:** Four geochemical tracks plus image and interpretation. Each track provides independent information (confirmed by redundancy analysis).

### Supplementary Figure: Extended Element Suite

**Additional panels:**
- Fe, Sr, K, Mn, Rb (raw elements)
- K/Ti, Rb/Sr (secondary ratios)
- Inc/Coh (matrix indicator)

### Cross-Plot Figures

1. **Ca vs Ti** - Carbonate-detrital mixing
2. **Fe vs Mn** - Redox systematics
3. **Zr vs Rb** - Grain size relationship
4. **Sr vs Ca** - Carbonate mineralogy

---

## 6. Statistical Recommendations

### Data Presentation

- Use **log₁₀ transformation** for ratios (Fe/Mn, Zr/Rb) to improve visualization
- Present elements in **counts per second (cps)** or calibrated units
- Report **median ± IQR** rather than mean ± SD for skewed distributions

### Avoid

- Plotting redundant ratio alongside its numerator element
- Over-interpreting small Fe/Mn variations without independent redox evidence
- Applying marine proxy interpretations (Ba/Ti productivity) to freshwater setting

---

## 7. Key References

1. Croudace, I.W. & Rothwell, R.G. (2015). *Micro-XRF Studies of Sediment Cores*. Springer.
2. Wesselingh, F.P. et al. (2002). Lake Pebas: a palaeoecological reconstruction. *Cainozoic Research* 1, 35-81.
3. Vonhof, H.B. et al. (1998). Reconstruction of the Miocene western Amazonian aquatic system. *Palaeogeogr. Palaeoclimatol. Palaeoecol.* 141, 85-93.
4. Weltje, G.J. & Tjallingii, R. (2008). Calibration of XRF core scanners. *Earth Planet. Sci. Lett.* 274, 423-438.
5. Wu et al. (2020). Evaluating Zr/Rb as grain-size indicator. *Geochem. Geophys. Geosyst.* 21, e2020GC009350.

---

*Report version: 2.0 | Updated with empirical redundancy analysis*
