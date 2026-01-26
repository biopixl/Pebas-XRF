# ==============================================================================
# Pebas-XRF: Stacked Stratigraphic Columns
# ==============================================================================
# Stack cores in correct stratigraphic order for continuous depth profiles
# Based on IMG_4825.jpeg core order documentation
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
# STRATIGRAPHIC ORDER (from IMG_4825.jpeg)
# ==============================================================================

# TAM cores (top to bottom): Group1 -> Group2 -> Group3
# SC cores (top to bottom): Group4 -> Group5 -> Group6 -> Group7

tam_order <- c("GROUP1", "GROUP2", "GROUP3")
sc_order <- c("GROUP4", "GROUP5", "GROUP6", "GROUP7")

# ==============================================================================
# LOAD DATA
# ==============================================================================

xrf_data <- read_csv(file.path(output_path, "tables", "xrf_data_ratios.csv"),
                     show_col_types = FALSE) %>%
  filter(qc_pass == TRUE)

message(sprintf("Loaded %d QC-passed measurements", nrow(xrf_data)))

# ==============================================================================
# CALCULATE CUMULATIVE DEPTH
# ==============================================================================

calculate_cumulative_depth <- function(data, group_order) {

  cumulative_offset <- 0
  result <- tibble()

  for (grp in group_order) {
    grp_data <- data %>% filter(group == grp)

    if (nrow(grp_data) == 0) next

    # Get the depth range for this group
    min_pos <- min(grp_data$position_mm, na.rm = TRUE)
    max_pos <- max(grp_data$position_mm, na.rm = TRUE)

    # Normalize to start at 0, then add cumulative offset
    grp_data <- grp_data %>%
      mutate(
        relative_depth = position_mm - min_pos,
        cumulative_depth = relative_depth + cumulative_offset,
        strat_order = which(group_order == grp)
      )

    result <- bind_rows(result, grp_data)

    # Update offset for next group (add this group's thickness + small gap)
    cumulative_offset <- cumulative_offset + (max_pos - min_pos) + 50  # 50mm gap between cores
  }

  return(result)
}

# Apply to each site
tam_stacked <- xrf_data %>%
  filter(core_series == "TAM") %>%
  calculate_cumulative_depth(tam_order)

sc_stacked <- xrf_data %>%
  filter(core_series == "SC") %>%
  calculate_cumulative_depth(sc_order)

message(sprintf("TAM: %.0f mm total depth", max(tam_stacked$cumulative_depth)))
message(sprintf("SC: %.0f mm total depth", max(sc_stacked$cumulative_depth)))

# ==============================================================================
# CREATE STACKED STRATIGRAPHIC PLOT FUNCTION
# ==============================================================================

create_stacked_plot <- function(data, site_name, group_order) {

  # Get group boundaries for annotation
  group_bounds <- data %>%
    group_by(group) %>%
    summarise(
      depth_top = min(cumulative_depth),
      depth_bottom = max(cumulative_depth),
      .groups = "drop"
    ) %>%
    arrange(depth_top)

  depth_max <- max(data$cumulative_depth, na.rm = TRUE)

  # Common theme
  common_theme <- theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 9),
      plot.margin = margin(5, 5, 5, 5)
    )

  # Panel 1: Fe (terrigenous)
  p_fe <- ggplot(data, aes(x = Fe, y = cumulative_depth, color = group)) +
    geom_path(linewidth = 0.4, alpha = 0.8) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_x_log10() +
    scale_y_reverse(limits = c(depth_max, 0)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Fe (cps)", y = "Depth (mm)") +
    common_theme +
    theme(legend.position = "none")

  # Panel 2: Ca (carbonate)
  p_ca <- ggplot(data, aes(x = Ca, y = cumulative_depth, color = group)) +
    geom_path(linewidth = 0.4, alpha = 0.8) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_x_log10() +
    scale_y_reverse(limits = c(depth_max, 0)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Ca (cps)", y = NULL) +
    common_theme +
    theme(legend.position = "none", axis.text.y = element_blank())

  # Panel 3: Ca/Ti (carbonate vs terrigenous)
  p_cati <- ggplot(data, aes(x = Ca_Ti, y = cumulative_depth, color = group)) +
    geom_path(linewidth = 0.4, alpha = 0.8) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_vline(xintercept = median(data$Ca_Ti, na.rm = TRUE),
               linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_y_reverse(limits = c(depth_max, 0)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Ca/Ti", y = NULL) +
    common_theme +
    theme(legend.position = "none", axis.text.y = element_blank())

  # Panel 4: Fe/Mn (redox)
  p_femn <- ggplot(data, aes(x = Fe_Mn, y = cumulative_depth, color = group)) +
    geom_path(linewidth = 0.4, alpha = 0.8) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "red", linewidth = 0.4) +
    scale_y_reverse(limits = c(depth_max, 0)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Fe/Mn", y = NULL) +
    common_theme +
    theme(legend.position = "none", axis.text.y = element_blank())

  # Panel 5: K/Ti (weathering)
  p_kti <- ggplot(data, aes(x = K_Ti, y = cumulative_depth, color = group)) +
    geom_path(linewidth = 0.4, alpha = 0.8) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_y_reverse(limits = c(depth_max, 0)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "K/Ti", y = NULL) +
    common_theme +
    theme(legend.position = "none", axis.text.y = element_blank())

  # Panel 6: Zr/Rb (grain size)
  p_zrrb <- ggplot(data, aes(x = Zr_Rb, y = cumulative_depth, color = group)) +
    geom_path(linewidth = 0.4, alpha = 0.8) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_y_reverse(limits = c(depth_max, 0)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Zr/Rb", y = NULL) +
    common_theme +
    theme(legend.position = "none", axis.text.y = element_blank())

  # Create core labels panel
  p_labels <- ggplot(group_bounds) +
    geom_rect(aes(xmin = 0, xmax = 1,
                  ymin = depth_top, ymax = depth_bottom,
                  fill = group), alpha = 0.3) +
    geom_text(aes(x = 0.5, y = (depth_top + depth_bottom) / 2, label = group),
              size = 2.5, fontface = "bold") +
    scale_y_reverse(limits = c(depth_max, 0)) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(legend.position = "none")

  # Combine all panels
  combined <- p_labels + p_fe + p_ca + p_cati + p_femn + p_kti + p_zrrb +
    plot_layout(widths = c(0.8, 1, 1, 1, 1, 1, 1)) +
    plot_annotation(
      title = sprintf("%s - Composite Stratigraphic Column", site_name),
      subtitle = sprintf("Cores stacked in stratigraphic order (top = youngest) | Total depth: %.0f mm | N = %d",
                         depth_max, nrow(data)),
      caption = "Proxies: Fe/Ca (elements), Ca/Ti (carbonate), Fe/Mn (redox, red line = oxic/anoxic threshold), K/Ti (weathering), Zr/Rb (grain size)",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        plot.caption = element_text(size = 8, color = "gray50")
      )
    )

  return(combined)
}

# ==============================================================================
# GENERATE STACKED PLOTS
# ==============================================================================

message("\n=== GENERATING STACKED STRATIGRAPHIC PLOTS ===")

# TAM
p_tam <- create_stacked_plot(tam_stacked, "Tamshiyacu (TAM)", tam_order)
ggsave(file.path(fig_path, "stacked_TAM.png"), p_tam,
       width = 14, height = 12, dpi = 150, bg = "white")
message("Saved: stacked_TAM.png")

# SC
p_sc <- create_stacked_plot(sc_stacked, "Santa Corina (SC)", sc_order)
ggsave(file.path(fig_path, "stacked_SC.png"), p_sc,
       width = 14, height = 14, dpi = 150, bg = "white")
message("Saved: stacked_SC.png")

# ==============================================================================
# SAVE STACKED DATA FOR FURTHER ANALYSIS
# ==============================================================================

stacked_data <- bind_rows(
  tam_stacked %>% mutate(site = "TAM"),
  sc_stacked %>% mutate(site = "SC")
)

write_csv(stacked_data, file.path(output_path, "tables", "xrf_data_stacked.csv"))
message(sprintf("\nStacked data saved: %d measurements", nrow(stacked_data)))

# ==============================================================================
# SUMMARY
# ==============================================================================

message("\n=== STRATIGRAPHIC SUMMARY ===")

stacked_data %>%
  group_by(site, group) %>%
  summarise(
    depth_top = min(cumulative_depth),
    depth_bottom = max(cumulative_depth),
    thickness = max(cumulative_depth) - min(cumulative_depth),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(site, depth_top) %>%
  print(n = 20)
