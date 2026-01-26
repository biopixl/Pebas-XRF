# Pebas-XRF

XRF core scanner analysis of sediment cores from the Pebas Formation, Peru.

## Core Sites
- **TAM** - Tamshiyacu
- **SC** - Santa Corina

## Data
Itrax XRF core scanner data (Mo tube, 30 kV, 10s exposure, 3mm step).

## Analysis Pipeline

```r
source("scripts/00_run_all.R")
```

### Scripts
| Script | Description |
|--------|-------------|
| `01_data_import.R` | Import Itrax result files, apply QC filters |
| `02_element_ratios.R` | Calculate paleoenvironmental element ratios |
| `03_stratigraphic_plots.R` | Multi-panel depth profiles |
| `04_pca_analysis.R` | CLR-transformed PCA analysis |

### Key Element Ratios
| Ratio | Interpretation |
|-------|---------------|
| Ca/Ti | Carbonate vs terrigenous input |
| Fe/Mn | Redox conditions |
| K/Ti | Weathering intensity |
| Zr/Rb | Grain size proxy |

## Requirements
- R 4.0+
- tidyverse, compositions, patchwork, here, fs
