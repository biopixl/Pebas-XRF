# ==============================================================================
# Pebas-XRF: Element Ratio Calculations
# ==============================================================================
# Calculates paleoenvironmentally-relevant element ratios from XRF data
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
# 2. ELEMENT RATIO DEFINITIONS
# ==============================================================================

# Key paleoenvironmental ratios for lacustrine/fluvial sediments
# Adapted for Miocene Pebas Formation context (tropical freshwater mega-wetland)
# Based on literature review: Croudace & Rothwell 2015; Davies et al. 2015;
# Wesselingh et al. 2002; Vonhof et al. 1998, 2003
# See docs/proxy_relevance_report.md for detailed proxy evaluation

#' Calculate standard element ratios
#'
#' @param data XRF data tibble with element columns
#' @return Data with added ratio columns
calculate_element_ratios <- function(data) {

  data %>%
    mutate(
      # ===================================================================
      # PRIMARY PROXIES - High confidence for Pebas Formation
      # ===================================================================

      # -------------------------------------------------------------------
      # TERRIGENOUS vs BIOGENIC INPUT (HIGH RELEVANCE)
      # -------------------------------------------------------------------
      # Ca/Ti: Carbonate (authigenic/biogenic) vs detrital input
      # In Pebas context: Higher = authigenic carbonate precipitation
      # during drier/evaporative periods; Lower = terrigenous dominance
      Ca_Ti = Ca / Ti,
      log_Ca_Ti = log10(Ca / Ti),

      # -------------------------------------------------------------------
      # CARBONATE MINERALOGY / SALINITY (HIGH RELEVANCE - NEW)
      # -------------------------------------------------------------------
      # Sr/Ca: Discriminates carbonate type and salinity influence
      # Higher Sr/Ca = aragonite (ostracods, molluscs) or brackish influence
      # Lower Sr/Ca = freshwater calcite precipitation
      # Useful for detecting oligohaline incursion events in Pebas
      Sr_Ca = Sr / Ca,
      log_Sr_Ca = log10(Sr / Ca),

      # -------------------------------------------------------------------
      # DETRITAL COMPOSITION / WEATHERING (HIGH RELEVANCE)
      # -------------------------------------------------------------------
      # Fe/Ti: Lateritic weathering indicator, runoff intensity
      # Higher = enhanced transport of weathered material from catchment
      # Excellent detection with Mo tube; less redox-affected than Fe/Mn
      Fe_Ti = Fe / Ti,
      log_Fe_Ti = log10(Fe / Ti),

      # -------------------------------------------------------------------
      # GRAIN SIZE / HYDRODYNAMIC ENERGY (MODERATE-HIGH RELEVANCE)
      # -------------------------------------------------------------------
      # Zr/Rb: Coarse (heavy minerals) vs fine (clays) fraction
      # Higher = coarser sediment, higher energy environment
      # Note: Some provenance effects possible (Andean vs cratonic sources)
      Zr_Rb = Zr / Rb,
      log_Zr_Rb = log10(Zr / Rb),

      # ===================================================================
      # SECONDARY PROXIES - Use with caveats
      # ===================================================================

      # -------------------------------------------------------------------
      # REDOX CONDITIONS (MODERATE RELEVANCE - USE CAUTIOUSLY)
      # -------------------------------------------------------------------
      # Fe/Mn: Redox indicator at sediment-water interface
      # CAUTION: Interpretation complicated by:
      # - Detrital Fe input from lateritic catchment
      # - Diagenetic Mn trapping under permanent anoxia
      # - Lake-specific calibration required
      Fe_Mn = Fe / Mn,
      log_Fe_Mn = log10(Fe / Mn),

      # Mn/Ca: Alternative redox-carbonate coupling indicator
      # May track bottom water oxygenation during carbonate precipitation
      Mn_Ca = Mn / Ca,

      # -------------------------------------------------------------------
      # BIOGENIC SILICA (MODERATE RELEVANCE - DETECTION LIMITED)
      # -------------------------------------------------------------------
      # Si/Ti: Biogenic silica (diatoms) vs detrital silica
      # CAUTION: Si detection marginal with Mo tube
      # Needs microscopic validation of diatom presence
      Si_Ti = Si / Ti,
      log_Si_Ti = log10(Si / Ti),

      # -------------------------------------------------------------------
      # SALINITY INDICATOR (MODERATE - REINTERPRETED)
      # -------------------------------------------------------------------
      # Ba/Ti: In freshwater settings, tracks salinity NOT productivity
      # Higher Ba = possible oligohaline conditions (barite precipitation)
      # Marine Ba/Ti productivity interpretation NOT valid for Pebas
      Ba_Ti = Ba / Ti,

      # -------------------------------------------------------------------
      # Rb/Sr: REINTERPRETED as carbonate/salinity indicator
      # NOT weathering proxy in Pebas context (Sr dominated by carbonates)
      # Lower Rb/Sr = carbonate-rich or brackish intervals
      # -------------------------------------------------------------------
      Rb_Sr = Rb / Sr,
      log_Rb_Sr = log10(Rb / Sr),

      # ===================================================================
      # LIMITED RELEVANCE PROXIES - Retained for completeness
      # ===================================================================

      # K/Ti: LIMITED RELEVANCE in tropical setting
      # Catchment already intensely weathered (CIA 80-100); K depleted
      # May only detect unusual K-enriched layers (volcanic ash?)
      K_Ti = K / Ti,
      log_K_Ti = log10(K / Ti),

      # Si/Al: NOT RECOMMENDED - poor Al and Si detection with Mo tube
      # Use Zr/Rb for grain size instead
      Si_Al = Si / Al,

      # K/Rb: LIMITED - both elements depleted in tropical weathering
      K_Rb = K / Rb,

      # -------------------------------------------------------------------
      # MATRIX CORRECTION (MODERATE RELEVANCE)
      # -------------------------------------------------------------------
      # Incoherent/Coherent scatter ratio (Compton correction)
      # Related to average atomic number of matrix
      # May track organic-rich vs mineral-rich intervals
      inc_coh = Mo_inc / Mo_coh
    )
}

# ==============================================================================
# 3. CALCULATE RATIOS
# ==============================================================================

message("Calculating element ratios...")

xrf_ratios <- xrf %>%
  calculate_element_ratios()

# Add Br_Ti if Br column exists (not always in calibrated data)
if ("Br" %in% names(xrf_ratios)) {
  xrf_ratios <- xrf_ratios %>%
    mutate(Br_Ti = Br / Ti)
  message("  Br_Ti ratio calculated (Br available in data)")
} else {
  message("  Br_Ti ratio skipped (Br not available in calibrated data)")
}

# Check for infinite values (division by zero)
# Updated ratio list including new proxies for Pebas Formation context
ratio_cols <- c(
  # Primary proxies (high confidence)
  "Ca_Ti", "Sr_Ca", "Fe_Ti", "Zr_Rb",
  # Secondary proxies (use with caveats)
  "Fe_Mn", "Mn_Ca", "Si_Ti", "Ba_Ti", "Rb_Sr",
  # Limited relevance (retained for completeness)
  "K_Ti", "Si_Al", "K_Rb",
  # Matrix correction
  "inc_coh"
)
if ("Br_Ti" %in% names(xrf_ratios)) ratio_cols <- c(ratio_cols, "Br_Ti")

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
# 4. SUMMARY STATISTICS
# ==============================================================================

# Summary by core series - focus on primary proxies for Pebas context
ratio_summary <- xrf_ratios %>%
  group_by(core_series) %>%
  summarise(
    n = n(),
    across(
      all_of(c("Ca_Ti", "Sr_Ca", "Fe_Ti", "Zr_Rb", "Fe_Mn")),
      list(
        median = ~median(., na.rm = TRUE),
        iqr = ~IQR(., na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  )

print(ratio_summary)

# ==============================================================================
# 5. CORRELATION ANALYSIS
# ==============================================================================

# Correlation matrix for key elements
element_cors <- xrf_ratios %>%
  select(all_of(c("Al", "Si", "K", "Ca", "Ti", "Fe", "Mn", "Rb", "Sr", "Zr"))) %>%
  cor(use = "pairwise.complete.obs")

# Save correlation matrix
write.csv(element_cors, file.path(output_path, "tables", "element_correlations.csv"))

message("\nElement correlation matrix:")
print(round(element_cors, 2))

# ==============================================================================
# 6. SAVE RESULTS
# ==============================================================================

write_csv(xrf_ratios, file.path(output_path, "tables", "xrf_data_ratios.csv"))
message(sprintf("\nRatio data saved to: %s",
                file.path(output_path, "tables", "xrf_data_ratios.csv")))

message("\nRatio calculation complete! Next step: run 03_stratigraphic_plots.R")
