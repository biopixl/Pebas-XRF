#!/usr/bin/env Rscript
# ==============================================================================
# Pebas-XRF: Stacked Stratigraphic Columns with Section Labels
# ==============================================================================
# Generate stacked stratigraphic figures with both GROUP and section labels
# ==============================================================================

library(tidyverse)
library(patchwork)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

base_path <- here::here()
output_path <- file.path(base_path, "output")
fig_path <- file.path(base_path, "manuscript", "figures")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

# Stratigraphic order
tam_order <- c("GROUP1", "GROUP2", "GROUP3")
sc_order <- c("GROUP4", "GROUP5", "GROUP6", "GROUP7")

# Section short labels
section_labels <- tribble(
  ~section, ~short_label,
  "TAM-1-2-3B-A", "A", "TAM-1-2-3B-B", "B", "TAM-1-2-3B-C", "C",
  "TAM-3A-4-5CDE-A", "A", "TAM-3A-4-5CDE-B", "B",
  "TAM-3A-4-5CDE-RUN2-C", "C", "TAM-3A-4-5CDE-RUN2-D", "D", "TAM-3A-4-5CDE-RUN2-E", "E",
  "TAM-5AB-6-7-A", "A", "TAM-5AB-6-7-B-RUN2", "B", "TAM-5AB-6-7-C", "C",
  "SC-1ABC-2-3C-A-RUN1", "A", "SC-1ABC-2-3C-RUN2-B", "B", "SC-1ABC-2-3C-RUN2-C", "C",
  "SC-3AB-4ABCD-A", "A", "SC-3AB-4ABCD-B", "B", "SC-3AB-4ABCD-C", "C",
  "SC-3AB-4ABCD-D", "D", "SC-3AB-4ABCD-RUN2-D", "D2",
  "SC-3AB-4ABCD-RUN2-E", "E", "SC-3AB-4ABCD-RUN2-F", "F",
  "SC-5-6-7ABC-A", "A", "SC-5-6-7ABC-B", "B", "SC-5-6-7ABC-C", "C",
  "SC-5-6-7ABC-D", "D", "SC-5-6-7ABC-E", "E",
  "SC8-A", "A"
)

# ==============================================================================
# LOAD DATA
# ==============================================================================

xrf_data <- read_csv(file.path(output_path, "tables", "xrf_data_stacked.csv"),
                     show_col_types = FALSE)

message(sprintf("Loaded %d measurements", nrow(xrf_data)))

# ==============================================================================
# CALCULATE CUMULATIVE DEPTH
# ==============================================================================

calculate_cumulative_depth <- function(data, group_order) {
  cumulative_offset <- 0
  result <- tibble()

  for (grp in group_order) {
    grp_data <- data %>% filter(group == grp)
    if (nrow(grp_data) == 0) next

    min_pos <- min(grp_data$position_mm, na.rm = TRUE)
    max_pos <- max(grp_data$position_mm, na.rm = TRUE)

    grp_data <- grp_data %>%
      mutate(
        relative_depth = position_mm - min_pos,
        cumulative_depth = relative_depth + cumulative_offset,
        strat_order = which(group_order == grp)
      )

    result <- bind_rows(result, grp_data)
    cumulative_offset <- cumulative_offset + (max_pos - min_pos) + 50
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

# ==============================================================================
# CREATE STACKED PLOT WITH SECTION LABELS
# ==============================================================================

create_stacked_plot <- function(data, site_name, group_order) {

  # Get group boundaries
  group_bounds <- data %>%
    group_by(group) %>%
    summarise(
      depth_top = min(cumulative_depth),
      depth_bottom = max(cumulative_depth),
      .groups = "drop"
    ) %>%
    arrange(depth_top)

  # Get section boundaries with labels
  section_bounds <- data %>%
    left_join(section_labels, by = "section") %>%
    group_by(group, section, short_label) %>%
    summarise(
      depth_top = min(cumulative_depth),
      depth_bottom = max(cumulative_depth),
      mid_depth = (min(cumulative_depth) + max(cumulative_depth)) / 2,
      .groups = "drop"
    ) %>%
    arrange(depth_top)

  depth_max <- max(data$cumulative_depth, na.rm = TRUE)

  common_theme <- theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 9),
      plot.margin = margin(5, 5, 5, 5)
    )

  # Panel: Group labels
  p_groups <- ggplot(group_bounds) +
    geom_rect(aes(xmin = 0, xmax = 1,
                  ymin = depth_top, ymax = depth_bottom,
                  fill = group), alpha = 0.3) +
    geom_text(aes(x = 0.5, y = (depth_top + depth_bottom) / 2, label = group),
              size = 2.5, fontface = "bold") +
    scale_y_reverse(limits = c(depth_max, 0)) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = "Depth (mm)", title = "Group") +
    theme_void() +
    theme(legend.position = "none",
          axis.title.y = element_text(size = 9, angle = 90),
          axis.text.y = element_text(size = 8))

  # Panel: Section labels
  p_sections <- ggplot(section_bounds) +
    geom_segment(aes(x = 0, xend = 1, y = depth_top, yend = depth_top),
                 color = "gray80", linewidth = 0.3) +
    geom_text(aes(x = 0.5, y = mid_depth, label = short_label),
              size = 2.5, fontface = "bold") +
    scale_y_reverse(limits = c(depth_max, 0)) +
    labs(x = NULL, y = NULL, title = "Sec.") +
    theme_void() +
    theme(plot.title = element_text(size = 9, hjust = 0.5))

  # Panel: Ca (carbonate)
  p_ca <- ggplot(data, aes(x = Ca, y = cumulative_depth, color = group)) +
    geom_path(linewidth = 0.4, alpha = 0.8) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_x_log10() +
    scale_y_reverse(limits = c(depth_max, 0)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Ca (cps)", y = NULL, title = "Carbonate") +
    common_theme +
    theme(legend.position = "none", axis.text.y = element_blank())

  # Panel: Ti (terrigenous)
  p_ti <- ggplot(data, aes(x = Ti, y = cumulative_depth, color = group)) +
    geom_path(linewidth = 0.4, alpha = 0.8) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_x_log10() +
    scale_y_reverse(limits = c(depth_max, 0)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Ti (cps)", y = NULL, title = "Terrigenous") +
    common_theme +
    theme(legend.position = "none", axis.text.y = element_blank())

  # Panel: Ca/Ti
  p_cati <- ggplot(data, aes(x = Ca_Ti, y = cumulative_depth, color = group)) +
    geom_path(linewidth = 0.4, alpha = 0.8) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_vline(xintercept = c(2, 5, 10), linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_y_reverse(limits = c(depth_max, 0)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Ca/Ti", y = NULL, title = "Ca/Ti") +
    common_theme +
    theme(legend.position = "none", axis.text.y = element_blank())

  # Panel: Fe/Mn (redox)
  p_femn <- ggplot(data, aes(x = Fe_Mn, y = cumulative_depth, color = group)) +
    geom_path(linewidth = 0.4, alpha = 0.8) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "red", linewidth = 0.4) +
    scale_y_reverse(limits = c(depth_max, 0)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Fe/Mn", y = NULL, title = "Redox") +
    common_theme +
    theme(legend.position = "none", axis.text.y = element_blank())

  # Panel: Zr/Rb (grain size)
  p_zrrb <- ggplot(data, aes(x = Zr_Rb, y = cumulative_depth, color = group)) +
    geom_path(linewidth = 0.4, alpha = 0.8) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_y_reverse(limits = c(depth_max, 0)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Zr/Rb", y = NULL, title = "Grain Size") +
    common_theme +
    theme(legend.position = "none", axis.text.y = element_blank())

  # Combine all panels
  combined <- p_groups + p_sections + p_ca + p_ti + p_cati + p_femn + p_zrrb +
    plot_layout(widths = c(0.8, 0.4, 1, 1, 1, 1, 1)) +
    plot_annotation(
      title = sprintf("%s - Composite Stratigraphic Column with Section Labels", site_name),
      subtitle = sprintf("Cores stacked in stratigraphic order (top = youngest) | Total depth: %.0f mm | N = %d",
                         depth_max, nrow(data)),
      caption = "Ca/Ti thresholds: 2, 5, 10 | Fe/Mn = 50 (oxic/reducing boundary)",
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
       width = 14, height = 12, dpi = 300, bg = "white")
message("Saved: stacked_TAM.png")

# SC
p_sc <- create_stacked_plot(sc_stacked, "Santa Corina (SC)", sc_order)
ggsave(file.path(fig_path, "stacked_SC.png"), p_sc,
       width = 14, height = 14, dpi = 300, bg = "white")
message("Saved: stacked_SC.png")

message("\n=== STACKED FIGURES COMPLETE ===")
