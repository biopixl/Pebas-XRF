# ==============================================================================
# Pebas-XRF: QC Review and Section Cropping
# ==============================================================================
# Interactive review of XRF data to identify and exclude:
# - Foam fills
# - Core gaps/breaks
# - Disturbed sections
# - Edge effects
# ==============================================================================

library(tidyverse)
library(patchwork)

# ==============================================================================
# 1. LOAD RAW DATA (before QC filtering)
# ==============================================================================

base_path <- here::here()
output_path <- file.path(base_path, "output")

# Re-import to get all data including invalid measurements
source(file.path(base_path, "scripts", "01_data_import.R"))

# Use xrf_qc which has all measurements with QC flags
xrf_all <- xrf_qc

message(sprintf("Loaded %d total measurements", nrow(xrf_all)))

# ==============================================================================
# 2. DIAGNOSTIC INDICATORS FOR FOAM/GAPS
# ==============================================================================

# Foam and gaps typically show:
# - Low total counts (cps) - X-rays not absorbed by dense material
# - High sample_surface - detector far from surface
# - Low Fe, Ti, Al - absence of mineral matter
# - Anomalous Mo scatter ratios

xrf_diagnostics <- xrf_all %>%
  group_by(section) %>%
  mutate(
    # Z-scores within each section
    cps_zscore = (cps - mean(cps, na.rm = TRUE)) / sd(cps, na.rm = TRUE),
    Fe_zscore = (Fe - mean(Fe, na.rm = TRUE)) / sd(Fe, na.rm = TRUE),
    Ti_zscore = (Ti - mean(Ti, na.rm = TRUE)) / sd(Ti, na.rm = TRUE),

    # Foam detection flags
    flag_low_cps = cps < quantile(cps, 0.05, na.rm = TRUE),
    flag_high_surface = sample_surface > 7.5,
    flag_low_Fe = Fe < quantile(Fe, 0.1, na.rm = TRUE),
    flag_low_Ti = Ti < quantile(Ti, 0.1, na.rm = TRUE),

    # Combined foam probability score
    foam_score = flag_low_cps + flag_high_surface + flag_low_Fe + flag_low_Ti,
    likely_foam = foam_score >= 3
  ) %>%
  ungroup()

# Summary of potential foam
foam_summary <- xrf_diagnostics %>%
  group_by(section) %>%
  summarise(
    n_total = n(),
    n_likely_foam = sum(likely_foam, na.rm = TRUE),
    pct_foam = round(100 * n_likely_foam / n_total, 1),
    position_range = paste0(min(position_mm), "-", max(position_mm), " mm")
  ) %>%
  arrange(desc(pct_foam))

message("\n=== FOAM/GAP DETECTION SUMMARY ===")
print(foam_summary, n = 30)

# ==============================================================================
# 3. GENERATE QC DIAGNOSTIC PLOTS
# ==============================================================================

#' Create QC diagnostic plot for a section
#'
#' Shows multiple indicators to identify foam/gaps
plot_section_qc <- function(data, section_name) {

  sect_data <- data %>% filter(section == section_name)

  if (nrow(sect_data) == 0) return(NULL)

  # Panel 1: Total counts (cps)
  p1 <- ggplot(sect_data, aes(x = cps, y = position_mm)) +
    geom_path(color = "gray50") +
    geom_point(aes(color = likely_foam), size = 1.5) +
    scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red")) +
    scale_y_reverse() +
    labs(x = "Total counts (cps)", color = "Likely\nFoam") +
    theme_minimal() +
    theme(legend.position = "none")

  # Panel 2: Sample surface distance
  p2 <- ggplot(sect_data, aes(x = sample_surface, y = position_mm)) +
    geom_path(color = "gray50") +
    geom_point(aes(color = likely_foam), size = 1.5) +
    scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red")) +
    geom_vline(xintercept = 7.5, linetype = "dashed", color = "red") +
    scale_y_reverse() +
    labs(x = "Surface dist (mm)") +
    theme_minimal() +
    theme(legend.position = "none")

  # Panel 3: Fe counts
  p3 <- ggplot(sect_data, aes(x = Fe, y = position_mm)) +
    geom_path(color = "gray50") +
    geom_point(aes(color = likely_foam), size = 1.5) +
    scale_color_manual(values = c("FALSE" = "darkorange", "TRUE" = "red")) +
    scale_x_log10() +
    scale_y_reverse() +
    labs(x = "Fe (cps, log)") +
    theme_minimal() +
    theme(legend.position = "none")

  # Panel 4: Ti counts
  p4 <- ggplot(sect_data, aes(x = Ti, y = position_mm)) +
    geom_path(color = "gray50") +
    geom_point(aes(color = likely_foam), size = 1.5) +
    scale_color_manual(values = c("FALSE" = "forestgreen", "TRUE" = "red")) +
    scale_x_log10() +
    scale_y_reverse() +
    labs(x = "Ti (cps, log)") +
    theme_minimal() +
    theme(legend.position = "none")

  # Panel 5: MSE (spectral fit quality)
  p5 <- ggplot(sect_data, aes(x = MSE, y = position_mm)) +
    geom_path(color = "gray50") +
    geom_point(aes(color = MSE > 10), size = 1.5) +
    scale_color_manual(values = c("FALSE" = "purple", "TRUE" = "red")) +
    geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
    scale_y_reverse() +
    labs(x = "MSE") +
    theme_minimal() +
    theme(legend.position = "none")

  # Panel 6: Foam score
  p6 <- ggplot(sect_data, aes(x = foam_score, y = position_mm)) +
    geom_point(aes(color = factor(foam_score)), size = 2) +
    scale_color_manual(values = c("0" = "green", "1" = "yellow",
                                   "2" = "orange", "3" = "red", "4" = "darkred")) +
    scale_y_reverse() +
    labs(x = "Foam Score (0-4)", color = "Score") +
    theme_minimal()

  # Combine
  combined <- (p1 | p2 | p3 | p4 | p5 | p6) +
    plot_annotation(
      title = sprintf("QC Review: %s", section_name),
      subtitle = sprintf("Red points = likely foam/gap (n=%d of %d)",
                        sum(sect_data$likely_foam), nrow(sect_data))
    )

  return(combined)
}

# Generate QC plots for all sections
message("\nGenerating QC diagnostic plots...")

qc_dir <- file.path(output_path, "figures", "qc_review")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

sections <- unique(xrf_diagnostics$section)

for (sect in sections) {
  p <- plot_section_qc(xrf_diagnostics, sect)
  if (!is.null(p)) {
    filename <- paste0("qc_", gsub("[^A-Za-z0-9_-]", "_", sect), ".png")
    ggsave(file.path(qc_dir, filename), p, width = 16, height = 8, dpi = 150, bg = "white")
  }
}

message(sprintf("QC plots saved to: %s", qc_dir))

# ==============================================================================
# 4. EXCLUSION ZONES CONFIGURATION
# ==============================================================================

# Create a template for manual exclusion zones
# Edit this file to define position ranges to exclude for each section

exclusion_template <- tibble(
  section = sections,
  exclude_start_mm = NA_real_,
  exclude_end_mm = NA_real_,
  notes = ""
)

# Pre-populate with detected foam zones
foam_zones <- xrf_diagnostics %>%
  filter(likely_foam) %>%
  group_by(section) %>%
  summarise(
    exclude_start_mm = min(position_mm),
    exclude_end_mm = max(position_mm),
    notes = sprintf("Auto-detected foam (%d points)", n())
  )

# Merge
exclusion_config <- exclusion_template %>%
  rows_update(foam_zones, by = "section") %>%
  arrange(section)

# Save configuration file
config_path <- file.path(base_path, "data", "exclusion_zones.csv")
dir.create(dirname(config_path), showWarnings = FALSE)
write_csv(exclusion_config, config_path)

message(sprintf("\nExclusion zone template saved to: %s", config_path))
message("Edit this file to define/adjust exclusion zones, then re-run analysis.")

# ==============================================================================
# 5. SUMMARY REPORT
# ==============================================================================

message("\n")
message("=" |> rep(60) |> paste(collapse = ""))
message("QC REVIEW SUMMARY")
message("=" |> rep(60) |> paste(collapse = ""))

message(sprintf("\nTotal measurements: %d", nrow(xrf_all)))
message(sprintf("Likely foam/gap: %d (%.1f%%)",
                sum(xrf_diagnostics$likely_foam, na.rm = TRUE),
                100 * mean(xrf_diagnostics$likely_foam, na.rm = TRUE)))

message("\nSections with >10% potential foam:")
foam_summary %>%
  filter(pct_foam > 10) %>%
  print()

message("\n--- NEXT STEPS ---")
message("1. Review QC plots in: output/figures/qc_review/")
message("2. Edit exclusion zones in: data/exclusion_zones.csv")
message("3. Re-run 01_data_import.R to apply exclusions")

# ==============================================================================
# 6. APPLY EXCLUSIONS (if config exists)
# ==============================================================================

apply_exclusions <- function(data, config_path) {

  if (!file.exists(config_path)) {
    message("No exclusion config found. Skipping.")
    return(data)
  }

  exclusions <- read_csv(config_path, show_col_types = FALSE) %>%
    filter(!is.na(exclude_start_mm) & !is.na(exclude_end_mm))

  if (nrow(exclusions) == 0) {
    message("No exclusion zones defined. Skipping.")
    return(data)
  }

  message(sprintf("Applying %d exclusion zones...", nrow(exclusions)))

  # Mark excluded points
  data_excluded <- data %>%
    left_join(exclusions %>% select(section, exclude_start_mm, exclude_end_mm),
              by = "section") %>%
    mutate(
      excluded = !is.na(exclude_start_mm) &
                 position_mm >= exclude_start_mm &
                 position_mm <= exclude_end_mm
    ) %>%
    select(-exclude_start_mm, -exclude_end_mm)

  n_excluded <- sum(data_excluded$excluded, na.rm = TRUE)
  message(sprintf("Excluded %d measurements (%.1f%%)",
                  n_excluded, 100 * n_excluded / nrow(data)))

  # Return filtered data
  data_excluded %>% filter(!excluded) %>% select(-excluded)
}

# Check if exclusions should be applied
if (file.exists(config_path)) {
  exclusions <- read_csv(config_path, show_col_types = FALSE)
  if (any(!is.na(exclusions$exclude_start_mm))) {
    message("\nExclusion zones detected. Applying...")
    xrf_cleaned <- apply_exclusions(xrf_diagnostics, config_path)

    # Save cleaned data
    write_csv(xrf_cleaned, file.path(output_path, "tables", "xrf_data_cleaned.csv"))
    message(sprintf("Cleaned data saved: %d measurements", nrow(xrf_cleaned)))
  }
}
