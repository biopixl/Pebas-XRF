# ==============================================================================
# Pebas-XRF: Element Ratio Calculations
# ==============================================================================
# Calculates paleoenvironmentally-relevant element ratios from XRF data
# Optimized for Miocene Pebas Formation tropical lacustrine setting
# See docs/proxy_relevance_report.md for detailed proxy evaluation
# ==============================================================================

library(tidyverse)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

base_path <- here::here()
output_path <- file.path(base_path, "output")

# Load QC-filtered data
xrf_qc <- read_csv(file.path(output_path, "tables", "xrf_data_qc.csv"),
                   show_col_types = FALSE)

# Filter to QC-passing measurements (keep excluded for optional masking)
xrf <- xrf_qc %>% filter(qc_pass)

# Ensure 'excluded' column exists (may be missing in older data)
if (!"excluded" %in% names(xrf)) {
  xrf <- xrf %>% mutate(excluded = FALSE)
}

message(sprintf("Loaded %d QC-passing measurements (%d in exclusion zones)",
                nrow(xrf), sum(xrf$excluded, na.rm = TRUE)))

# ==============================================================================
# 2. PROXY DEFINITIONS - Based on Empirical Redundancy Analysis
# ==============================================================================
#
# TIER 1 (Primary - Main Figures):
#   Ca      - Carbonate/authigenic signal (Ca/Ti redundant, r=0.89)
#   Ti      - Terrigenous detrital flux
#   Fe/Mn   - Redox indicator (r=0.25 with Fe; adds information)
#   Zr/Rb   - Grain size proxy (r=0.60 with Zr; adds information)
#
# TIER 2 (Supporting - Supplementary):
#   Fe      - Lateritic input
#   Sr      - Carbonate mineralogy / salinity
#   K/Ti    - Weathering indicator (r=0.39; useful)
#   Rb/Sr   - Inverse carbonate indicator (r=-0.17; useful)
#
# NOT RECOMMENDED:
#   Ca/Ti   - Redundant with Ca (r=0.89)
#   Si/Al   - Poor Mo tube detection
#   Ba/Ti   - Marine productivity interpretation invalid for freshwater
# ==============================================================================

#' Calculate element ratios for Pebas Formation analysis
#'
#' @param data XRF data tibble with element columns
#' @return Data with added ratio columns
calculate_element_ratios <- function(data) {

  data %>%
    mutate(
      # =================================================================
      # TIER 1: PRIMARY PROXIES (Main Figures)
      # =================================================================

      # Fe/Mn: Redox conditions at sediment-water interface
      # r = 0.25 with Fe - ratio adds substantial information
      # Higher = more reducing (anoxic) conditions
      # CAVEAT: Detrital Fe can obscure signal; use for major transitions only
      Fe_Mn = Fe / Mn,
      log_Fe_Mn = log10(Fe / Mn),

      # Zr/Rb: Grain size / hydrodynamic energy proxy
      # r = 0.60 with Zr - ratio adds information
      # Higher = coarser sediment, higher energy environment
      Zr_Rb = Zr / Rb,
      log_Zr_Rb = log10(Zr / Rb),

      # =================================================================
      # TIER 2: SUPPORTING PROXIES (Supplementary Figures)
      # =================================================================

      # K/Ti: Clay composition / weathering indicator
      # r = 0.39 with K - ratio adds information
      # Limited relevance in tropical setting (already at weathering ceiling)
      # Useful for detecting K-enriched anomalies (ash layers?)
      K_Ti = K / Ti,
      log_K_Ti = log10(K / Ti),

      # Rb/Sr: Inverse carbonate indicator
      # r = -0.17 with Rb - strongly differentiated from raw Rb
      # Lower values = more Sr-rich carbonate intervals
      Rb_Sr = Rb / Sr,
      log_Rb_Sr = log10(Rb / Sr),

      # Sr/Ca: Carbonate mineralogy discriminator
      # Aragonite (molluscs, ostracods) has higher Sr/Ca than calcite
      # Useful for detecting brackish incursion events
      Sr_Ca = Sr / Ca,

      # =================================================================
      # TIER 3: CONTEXTUAL (Report but don't feature)
      # =================================================================

      # Fe/Ti: Detrital composition
      # r = 0.48 with Fe - provides some additional information
      Fe_Ti = Fe / Ti,

      # Si/Ti: Biogenic silica indicator (if diatoms present)
      # r = 0.56 with Si - useful but Si detection marginal
      Si_Ti = Si / Ti,

      # Ba/Ti: Salinity indicator (NOT productivity in freshwater)
      # r = 0.75 with Ba - borderline redundant
      Ba_Ti = Ba / Ti,

      # Mn/Ca: Alternative redox-carbonate coupling
      Mn_Ca = Mn / Ca,

      # Inc/Coh: Matrix composition / organic matter indicator
      inc_coh = Mo_inc / Mo_coh,

      # =================================================================
      # DEPRECATED - Retained for backwards compatibility only
      # =================================================================

      # Ca/Ti: REDUNDANT with Ca (r = 0.89)
      # Ti variance too low relative to Ca; ratio dominated by Ca
      # USE RAW Ca INSTEAD for carbonate signal
      Ca_Ti = Ca / Ti,
      log_Ca_Ti = log10(Ca / Ti),

      # Si/Al: NOT RECOMMENDED - poor Mo tube detection for both
      Si_Al = Si / Al,

      # K/Rb: Limited diagnostic value in tropical setting
      K_Rb = K / Rb
    )
}

# ==============================================================================
# 3. CALCULATE RATIOS
# ==============================================================================

message("Calculating element ratios...")

xrf_ratios <- xrf %>%
  calculate_element_ratios()

# Define ratio columns for QC checking
ratio_cols <- c(
  # Tier 1
  "Fe_Mn", "Zr_Rb",
  # Tier 2
  "K_Ti", "Rb_Sr", "Sr_Ca",
  # Tier 3
  "Fe_Ti", "Si_Ti", "Ba_Ti", "Mn_Ca", "inc_coh",
  # Deprecated
  "Ca_Ti", "Si_Al", "K_Rb"
)

# Check for infinite values (division by zero)
for (col in ratio_cols) {
  n_inf <- sum(is.infinite(xrf_ratios[[col]]), na.rm = TRUE)
  n_na <- sum(is.na(xrf_ratios[[col]]))
  if (n_inf > 0 || n_na > 0) {
    message(sprintf("  %s: %d Inf, %d NA values", col, n_inf, n_na))
  }
}

# Replace Inf with NA
xrf_ratios <- xrf_ratios %>%
  mutate(across(all_of(ratio_cols), ~ifelse(is.infinite(.), NA_real_, .)))

# ==============================================================================
# 4. SUMMARY STATISTICS - Primary Proxies
# ==============================================================================

message("\n=== PRIMARY PROXY SUMMARY (Tier 1) ===")

# Summary by core series - focus on recommended proxies
ratio_summary <- xrf_ratios %>%
  filter(!excluded) %>%
  group_by(core_series) %>%
  summarise(
    n = n(),
    # Tier 1 elements
    Ca_median = median(Ca, na.rm = TRUE),
    Ti_median = median(Ti, na.rm = TRUE),
    # Tier 1 ratios
    Fe_Mn_median = median(Fe_Mn, na.rm = TRUE),
    Zr_Rb_median = median(Zr_Rb, na.rm = TRUE),
    # Tier 2
    Sr_median = median(Sr, na.rm = TRUE),
    Fe_median = median(Fe, na.rm = TRUE)
  )

print(ratio_summary)

# ==============================================================================
# 5. REDUNDANCY CHECK - Confirm Ca/Ti is redundant
# ==============================================================================

message("\n=== REDUNDANCY ANALYSIS ===")
message("Ratios with r > 0.85 vs numerator are redundant:\n")

redundancy_check <- tibble(
  ratio = c("Ca/Ti", "Fe/Mn", "Zr/Rb", "K/Ti", "Rb/Sr", "Fe/Ti", "Ba/Ti"),
  numerator = c("Ca", "Fe", "Zr", "K", "Rb", "Fe", "Ba"),
  r = c(
    cor(xrf_ratios$Ca, xrf_ratios$Ca_Ti, use = "complete.obs"),
    cor(xrf_ratios$Fe, xrf_ratios$Fe_Mn, use = "complete.obs"),
    cor(xrf_ratios$Zr, xrf_ratios$Zr_Rb, use = "complete.obs"),
    cor(xrf_ratios$K, xrf_ratios$K_Ti, use = "complete.obs"),
    cor(xrf_ratios$Rb, xrf_ratios$Rb_Sr, use = "complete.obs"),
    cor(xrf_ratios$Fe, xrf_ratios$Fe_Ti, use = "complete.obs"),
    cor(xrf_ratios$Ba, xrf_ratios$Ba_Ti, use = "complete.obs")
  )
) %>%
  mutate(
    status = ifelse(abs(r) > 0.85, "REDUNDANT", "USEFUL"),
    r = round(r, 3)
  ) %>%
  arrange(desc(abs(r)))

print(redundancy_check)

# ==============================================================================
# 6. CORRELATION MATRIX - Key Elements
# ==============================================================================

message("\n=== ELEMENT CORRELATION MATRIX ===")

element_cors <- xrf_ratios %>%
  filter(!excluded) %>%
  select(all_of(c("Ca", "Ti", "Fe", "Mn", "K", "Rb", "Sr", "Zr"))) %>%
  cor(use = "pairwise.complete.obs")

write.csv(element_cors, file.path(output_path, "tables", "element_correlations.csv"))
print(round(element_cors, 2))

# ==============================================================================
# 7. SAVE RESULTS
# ==============================================================================

write_csv(xrf_ratios, file.path(output_path, "tables", "xrf_data_ratios.csv"))
message(sprintf("\nRatio data saved to: %s",
                file.path(output_path, "tables", "xrf_data_ratios.csv")))

# Save redundancy analysis
write_csv(redundancy_check, file.path(output_path, "tables", "proxy_redundancy_analysis.csv"))

message("\n=== RECOMMENDED PROXIES FOR MANUSCRIPT ===")
message("MAIN FIGURES: Ca, Ti, Fe/Mn, Zr/Rb")
message("SUPPLEMENTARY: Fe, Sr, K/Ti, Rb/Sr")
message("AVOID: Ca/Ti (redundant), Si/Al (poor detection)")
message("\nRatio calculation complete! Next step: run 03_stratigraphic_plots.R")
