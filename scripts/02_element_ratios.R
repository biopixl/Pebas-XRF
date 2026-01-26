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
# Based on literature review (Croudace & Rothwell 2015; Davies et al. 2015)

#' Calculate standard element ratios
#'
#' @param data XRF data tibble with element columns
#' @return Data with added ratio columns
calculate_element_ratios <- function(data) {

  data %>%
    mutate(
      # -------------------------------------------------------------------
      # TERRIGENOUS vs BIOGENIC INPUT
      # -------------------------------------------------------------------
      # Ca/Ti: Carbonate (biogenic/authigenic) vs detrital input
      # Higher = more carbonate precipitation or less clastic input
      Ca_Ti = Ca / Ti,
      log_Ca_Ti = log10(Ca / Ti),

      # Si/Ti: Biogenic silica (diatoms) vs detrital silica
      # Higher = more diatom productivity (in lacustrine settings)
      Si_Ti = Si / Ti,
      log_Si_Ti = log10(Si / Ti),

      # -------------------------------------------------------------------
      # REDOX CONDITIONS
      # -------------------------------------------------------------------
      # Fe/Mn: Redox indicator at sediment-water interface
      # Higher = more reducing (anoxic) conditions
      # Mn is preferentially mobilized under reducing conditions
      Fe_Mn = Fe / Mn,
      log_Fe_Mn = log10(Fe / Mn),

      # -------------------------------------------------------------------
      # CHEMICAL WEATHERING INTENSITY
      # -------------------------------------------------------------------
      # K/Ti: Clay composition, weathering intensity
      # Higher K = more illite/muscovite, less weathered
      K_Ti = K / Ti,
      log_K_Ti = log10(K / Ti),

      # Rb/Sr: Silicate weathering intensity
      # Higher = more intense chemical weathering
      Rb_Sr = Rb / Sr,
      log_Rb_Sr = log10(Rb / Sr),

      # -------------------------------------------------------------------
      # GRAIN SIZE / ENERGY PROXIES
      # -------------------------------------------------------------------
      # Zr/Rb: Coarse (heavy minerals) vs fine (clays) fraction
      # Higher = coarser sediment, higher energy environment
      Zr_Rb = Zr / Rb,
      log_Zr_Rb = log10(Zr / Rb),

      # Si/Al: Quartz (coarse) vs clay (fine) content
      # Higher = coarser sediment
      Si_Al = Si / Al,

      # -------------------------------------------------------------------
      # PRODUCTIVITY / ORGANIC MATTER
      # -------------------------------------------------------------------
      # Ba/Ti: Paleoproductivity (Ba scavenged by organic particles)
      Ba_Ti = Ba / Ti,

      # -------------------------------------------------------------------
      # DETRITAL COMPOSITION
      # -------------------------------------------------------------------
      # Fe/Ti: Fe-bearing minerals relative to Ti
      Fe_Ti = Fe / Ti,

      # K/Rb: Potassium feldspar vs clay
      K_Rb = K / Rb,

      # -------------------------------------------------------------------
      # MATRIX CORRECTION
      # -------------------------------------------------------------------
      # Incoherent/Coherent scatter ratio (Compton correction)
      # Related to average atomic number of matrix
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
ratio_cols <- c("Ca_Ti", "Si_Ti", "Fe_Mn", "K_Ti", "Rb_Sr", "Zr_Rb", "Si_Al",
                "Ba_Ti", "Fe_Ti", "K_Rb", "inc_coh")
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

# Summary by core series
ratio_summary <- xrf_ratios %>%
  group_by(core_series) %>%
  summarise(
    n = n(),
    across(
      all_of(c("Ca_Ti", "Fe_Mn", "K_Ti", "Zr_Rb")),
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
