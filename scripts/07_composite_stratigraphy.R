# ==============================================================================
# Pebas-XRF: Composite Stratigraphic Plots
# ==============================================================================
# Properly align sections within cores for continuous paleoenvironmental records
# ==============================================================================

library(tidyverse)
library(patchwork)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

base_path <- here::here()
output_path <- file.path(base_path, "output")
fig_path <- file.path(output_path, "figures", "stratigraphy")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# LOAD DATA
# ==============================================================================

xrf_data <- read_csv(file.path(output_path, "tables", "xrf_data_ratios.csv"),
                     show_col_types = FALSE) %>%
  filter(qc_pass == TRUE)

message(sprintf("Loaded %d QC-passed measurements", nrow(xrf_data)))

# ==============================================================================
# DEFINE CORE COMPOSITES
# ==============================================================================

# Each GROUP represents a continuous core - sections are consecutive segments
# Position (mm) represents actual depth along the core

core_info <- xrf_data %>%
  group_by(core_series, group) %>%
  summarise(
    sections = paste(unique(section), collapse = ", "),
    n_sections = n_distinct(section),
    depth_min = min(position_mm),
    depth_max = max(position_mm),
    n_measurements = n(),
    .groups = "drop"
  ) %>%
  arrange(core_series, group)

message("\n=== CORE STRUCTURE ===")
print(core_info)

# ==============================================================================
# CREATE COMPOSITE PLOT FUNCTION
# ==============================================================================

create_composite_plot <- function(data, core_name, title_prefix = "") {

  # Sort by depth
  data <- data %>% arrange(position_mm)

  depth_range <- range(data$position_mm)

  # Create individual element/ratio panels

  # Panel 1: Fe (log scale) - terrigenous indicator
  p_fe <- ggplot(data, aes(x = Fe, y = position_mm)) +
    geom_path(color = "brown", linewidth = 0.5) +
    geom_point(size = 0.8, color = "brown", alpha = 0.5) +
    scale_x_log10() +
    scale_y_reverse() +
    labs(x = "Fe (cps)", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.title.x = element_text(size = 8))

  # Panel 2: Ca (log scale) - carbonate
  p_ca <- ggplot(data, aes(x = Ca, y = position_mm)) +
    geom_path(color = "steelblue", linewidth = 0.5) +
    geom_point(size = 0.8, color = "steelblue", alpha = 0.5) +
    scale_x_log10() +
    scale_y_reverse() +
    labs(x = "Ca (cps)", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.title.x = element_text(size = 8))

  # Panel 3: Ca/Ti ratio - carbonate vs terrigenous
  p_cati <- ggplot(data, aes(x = Ca_Ti, y = position_mm)) +
    geom_path(color = "darkgreen", linewidth = 0.5) +
    geom_point(size = 0.8, color = "darkgreen", alpha = 0.5) +
    geom_vline(xintercept = median(data$Ca_Ti, na.rm = TRUE),
               linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_y_reverse() +
    labs(x = "Ca/Ti", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.title.x = element_text(size = 8))

  # Panel 4: Fe/Mn ratio - redox proxy
  p_femn <- ggplot(data, aes(x = Fe_Mn, y = position_mm)) +
    geom_path(color = "purple", linewidth = 0.5) +
    geom_point(size = 0.8, color = "purple", alpha = 0.5) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "red", linewidth = 0.3) +
    annotate("text", x = 55, y = depth_range[1], label = "oxic/anoxic",
             color = "red", size = 2, hjust = 0) +
    scale_y_reverse() +
    labs(x = "Fe/Mn", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.title.x = element_text(size = 8))

  # Panel 5: K/Ti ratio - weathering intensity
  p_kti <- ggplot(data, aes(x = K_Ti, y = position_mm)) +
    geom_path(color = "orange", linewidth = 0.5) +
    geom_point(size = 0.8, color = "orange", alpha = 0.5) +
    scale_y_reverse() +
    labs(x = "K/Ti", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.title.x = element_text(size = 8))

  # Panel 6: Zr/Rb ratio - grain size
  p_zrrb <- ggplot(data, aes(x = Zr_Rb, y = position_mm)) +
    geom_path(color = "darkred", linewidth = 0.5) +
    geom_point(size = 0.8, color = "darkred", alpha = 0.5) +
    scale_y_reverse() +
    labs(x = "Zr/Rb", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.title.x = element_text(size = 8))

  # Add depth axis to first panel only
  p_fe <- p_fe + labs(y = "Depth (mm)")

  # Combine panels
  combined <- p_fe + p_ca + p_cati + p_femn + p_kti + p_zrrb +
    plot_layout(nrow = 1) +
    plot_annotation(
      title = paste0(title_prefix, core_name),
      subtitle = sprintf("Depth: %.0f - %.0f mm | N = %d measurements",
                         depth_range[1], depth_range[2], nrow(data)),
      theme = theme(
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40")
      )
    )

  return(combined)
}

# ==============================================================================
# GENERATE COMPOSITE PLOTS BY CORE (GROUP)
# ==============================================================================

message("\n=== GENERATING COMPOSITE PLOTS ===")

# Process each core (GROUP)
cores <- unique(xrf_data$group)

for (core in cores) {
  core_data <- xrf_data %>% filter(group == core)
  series <- unique(core_data$core_series)

  message(sprintf("Processing: %s (%s)", core, series))

  p <- create_composite_plot(core_data, core, paste0(series, " - "))

  filename <- paste0("composite_", core, ".png")
  ggsave(file.path(fig_path, filename), p,
         width = 14, height = 8, dpi = 150, bg = "white")
}

# ==============================================================================
# CREATE SITE-LEVEL COMPOSITES (ALL CORES PER SITE)
# ==============================================================================

message("\n=== GENERATING SITE COMPOSITES ===")

for (series in c("TAM", "SC")) {

  site_data <- xrf_data %>% filter(core_series == series)
  cores_in_site <- unique(site_data$group)

  # Create a multi-panel figure with one row per core
  plot_list <- list()

  for (i in seq_along(cores_in_site)) {
    core <- cores_in_site[i]
    core_data <- site_data %>% filter(group == core)

    p <- create_composite_plot(core_data, core, "")
    plot_list[[i]] <- p
  }

  # Stack vertically
  combined_site <- wrap_plots(plot_list, ncol = 1) +
    plot_annotation(
      title = sprintf("%s Cores - Composite Stratigraphy",
                      ifelse(series == "TAM", "Tamshiyacu", "Santa Corina")),
      subtitle = "Fe, Ca, Ca/Ti (carbonate), Fe/Mn (redox), K/Ti (weathering), Zr/Rb (grain size)",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray40")
      )
    )

  filename <- paste0("site_composite_", series, ".png")
  n_cores <- length(cores_in_site)
  ggsave(file.path(fig_path, filename), combined_site,
         width = 14, height = 6 * n_cores, dpi = 150, bg = "white")

  message(sprintf("Saved: %s (%d cores)", filename, n_cores))
}

# ==============================================================================
# CREATE CORRELATION PANEL FOR CYCLIC ANALYSIS
# ==============================================================================

message("\n=== GENERATING CORRELATION PLOTS ===")

# For each core, show element correlations to identify consistent signals

for (series in c("TAM", "SC")) {

  site_data <- xrf_data %>% filter(core_series == series)

  # Rolling averages to smooth noise and reveal trends
  window_size <- 5  # 5-point rolling average (15mm at 3mm step)

  site_smoothed <- site_data %>%
    group_by(group) %>%
    arrange(position_mm) %>%
    mutate(
      Fe_smooth = zoo::rollmean(Fe, window_size, fill = NA),
      Ca_smooth = zoo::rollmean(Ca, window_size, fill = NA),
      Ca_Ti_smooth = zoo::rollmean(Ca_Ti, window_size, fill = NA),
      Fe_Mn_smooth = zoo::rollmean(Fe_Mn, window_size, fill = NA)
    ) %>%
    ungroup()

  # Cross-correlation between Ca/Ti and Fe/Mn (should be related to depositional conditions)
  p_corr <- ggplot(site_smoothed, aes(x = Ca_Ti_smooth, y = Fe_Mn_smooth, color = group)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    scale_x_log10() + scale_y_log10() +
    labs(
      title = sprintf("%s: Ca/Ti vs Fe/Mn Correlation by Core",
                      ifelse(series == "TAM", "Tamshiyacu", "Santa Corina")),
      subtitle = "Relationship between carbonate accumulation and redox conditions",
      x = "Ca/Ti (carbonate proxy, smoothed)",
      y = "Fe/Mn (redox proxy, smoothed)",
      color = "Core"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")

  ggsave(file.path(fig_path, paste0("correlation_", series, ".png")), p_corr,
         width = 10, height = 7, dpi = 150, bg = "white")
}

message("\n=== COMPLETE ===")
message(sprintf("Figures saved to: %s", fig_path))
