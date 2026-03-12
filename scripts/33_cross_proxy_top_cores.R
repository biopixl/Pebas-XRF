#!/usr/bin/env Rscript
# ==============================================================================
# Cross-Proxy Comparison: Ca/Ti, Fe/Ti, Mn/Ti at Top of Cores
# ==============================================================================

library(tidyverse)
library(patchwork)
library(magick)
library(zoo)

# Paths
base_path <- here::here()
output_dir <- file.path(base_path, "output/figures/cross_proxy")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
xrf_data <- read_csv(file.path(base_path, "output/tables/xrf_data_stacked.csv"),
                     show_col_types = FALSE) %>%
  mutate(
    Ca_Ti = Ca / Ti,
    Fe_Ti = Fe / Ti,
    Mn_Ti = Mn / Ti,
    site = if_else(str_detect(section, "^TAM"), "TAM", "SC")
  )

xrf_valid <- xrf_data %>% filter(qc_pass, !excluded)

# Section configuration
section_config <- tribble(
  ~section, ~optical_path, ~group, ~site, ~strat_order,
  "TAM-5AB-6-7-C", "GROUP3/TAM-5AB-6-7/TAM-5AB-6-7-C", "GROUP3", "TAM", 11,
  "TAM-5AB-6-7-B-RUN2", "GROUP3/TAM-5AB-6-7/TAM-5AB-6-7-B-RUN2", "GROUP3", "TAM", 10,
  "TAM-5AB-6-7-A", "GROUP3/TAM-5AB-6-7/TAM-5AB-6-7-A", "GROUP3", "TAM", 9,
  "SC8-A", "GROUP7/SC8-A/SC8-A", "GROUP7", "SC", 27
)

data_path <- file.path(base_path, "TAM-SC-IsaacA")

# Image processing functions
read_document_params <- function(doc_path) {
  if (!file.exists(doc_path)) return(NULL)
  lines <- readLines(doc_path, warn = FALSE)
  xrf_line <- grep("Start coordinate", lines, value = TRUE)[1]
  if (is.na(xrf_line)) return(NULL)
  parts <- strsplit(xrf_line, "\t")[[1]]
  xrf_start <- as.numeric(parts[2])
  xrf_end <- as.numeric(parts[4])
  optical_line <- grep("Optical Start", lines, value = TRUE)[1]
  if (!is.na(optical_line)) {
    parts <- strsplit(optical_line, "\t")[[1]]
    optical_start <- as.numeric(parts[2])
    optical_end <- as.numeric(parts[4])
  } else {
    optical_start <- 0
    optical_end <- 1222
  }
  list(xrf_start = xrf_start, xrf_end = xrf_end,
       optical_start = optical_start, optical_end = optical_end)
}

process_core_image <- function(section_name, optical_path, data_start, data_end) {
  section_dir <- file.path(data_path, optical_path)
  img_path <- file.path(section_dir, "optical.tif")
  doc_path <- file.path(section_dir, "document.txt")

  if (!file.exists(img_path)) return(NULL)
  params <- read_document_params(doc_path)
  if (is.null(params)) return(NULL)

  img <- image_read(img_path)
  img_info <- image_info(img)
  optical_range <- params$optical_end - params$optical_start
  pixels_per_mm <- img_info$width / optical_range

  xrf_start_px <- max(0, round((data_start - params$optical_start) * pixels_per_mm))
  xrf_end_px <- min(img_info$width, round((data_end - params$optical_start) * pixels_per_mm))
  crop_width <- xrf_end_px - xrf_start_px
  if (crop_width <= 0) return(NULL)

  crop_geom <- sprintf("%dx%d+%d+0", crop_width, img_info$height, xrf_start_px)
  img_cropped <- image_crop(img, crop_geom)

  # Enhance
  img_enhanced <- img_cropped %>%
    image_level(black_point = 2, white_point = 60, mid_point = 0.5) %>%
    image_normalize() %>%
    image_modulate(brightness = 180) %>%
    image_contrast(sharpen = 25)

  list(image = img_enhanced, section = section_name,
       xrf_start = data_start, xrf_end = data_end)
}

# Color scheme for Ti-normalized ratios
proxy_colors <- c(
  "Ca_Ti" = "#2166ac",  # Blue for carbonate
  "Fe_Ti" = "#b2182b",  # Red for iron
  "Mn_Ti" = "#1a9850"   # Green for manganese
)

proxy_labels <- c(
  "Ca_Ti" = "Ca/Ti",
  "Fe_Ti" = "Fe/Ti",
  "Mn_Ti" = "Mn/Ti"
)

# Site colors
site_colors <- c("TAM" = "#D55E00", "SC" = "#0072B2")

# Theme
theme_strat <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      legend.position = "none",
      plot.margin = margin(5, 10, 5, 5)
    )
}

# ==============================================================================
# Generate single section figure with Ca/Ti, Fe/Ti, Mn/Ti
# ==============================================================================

generate_cross_proxy_figure <- function(section_name, title_suffix = "") {
  sect_info <- section_config %>% filter(section == section_name)
  if (nrow(sect_info) == 0) return(NULL)

  sect_data <- xrf_valid %>%
    filter(section == section_name) %>%
    arrange(position_mm)

  if (nrow(sect_data) < 10) return(NULL)

  data_start <- min(sect_data$position_mm)
  data_end <- max(sect_data$position_mm)
  depth_cm <- (data_end - data_start) / 10

  # Apply smoothing
  for (proxy in c("Ca_Ti", "Fe_Ti", "Mn_Ti")) {
    smooth_col <- paste0(proxy, "_smooth")
    sect_data[[smooth_col]] <- rollmean(sect_data[[proxy]], 5, fill = NA, align = "center")
  }

  # Position in cm from start
  sect_data <- sect_data %>%
    mutate(depth_cm = (position_mm - data_start) / 10)

  plot_list <- list()

  # Core image
  img_result <- tryCatch(
    process_core_image(section_name, sect_info$optical_path, data_start, data_end),
    error = function(e) NULL
  )

  if (!is.null(img_result)) {
    img_rotated <- image_rotate(img_result$image, 90)
    target_width <- 300
    target_height <- round(depth_cm * 10)
    img_scaled <- image_resize(img_rotated, sprintf("%dx%d!", target_width, target_height))
    core_raster <- as.raster(img_scaled)

    site_label <- if (sect_info$site == "TAM") "TAM" else "SC"

    core_plot <- ggplot() +
      annotation_raster(core_raster, xmin = 0, xmax = 1, ymin = 0, ymax = depth_cm) +
      scale_y_reverse(limits = c(depth_cm, 0), expand = c(0, 0)) +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = NULL, y = "Depth (cm)", title = "Core") +
      theme_strat() +
      theme(axis.text.x = element_blank(),
            panel.background = element_rect(fill = "gray30", color = NA))

    plot_list$core <- core_plot
  }

  # Proxy panels
  for (proxy in c("Ca_Ti", "Fe_Ti", "Mn_Ti")) {
    smooth_col <- paste0(proxy, "_smooth")
    color <- proxy_colors[proxy]
    label <- proxy_labels[proxy]

    p <- ggplot(sect_data, aes(y = depth_cm)) +
      geom_point(aes(x = .data[[proxy]]), color = color, alpha = 0.4, size = 1.2) +
      geom_path(aes(x = .data[[smooth_col]]), color = color, linewidth = 1.2, na.rm = TRUE) +
      scale_x_log10() +
      scale_y_reverse(limits = c(depth_cm, 0), expand = c(0.01, 0)) +
      labs(x = label, y = NULL, title = label) +
      theme_strat() +
      theme(axis.text.y = element_blank())

    plot_list[[proxy]] <- p
  }

  # Combine
  n_panels <- length(plot_list)
  widths <- c(0.5, rep(1, n_panels - 1))

  combined <- wrap_plots(plot_list, nrow = 1, widths = widths) +
    plot_annotation(
      title = sprintf("%s%s", section_name, title_suffix),
      subtitle = "Ti-normalized elemental ratios: Ca/Ti (carbonate), Fe/Ti (iron), Mn/Ti (manganese)"
    )

  combined
}

# ==============================================================================
# Generate side-by-side comparison figure
# ==============================================================================

generate_comparison_figure <- function(tam_sections, sc_sections, title) {

  # Process TAM sections
  tam_data <- xrf_valid %>%
    filter(section %in% tam_sections) %>%
    arrange(strat_order, position_mm)

  # Process SC sections
  sc_data <- xrf_valid %>%
    filter(section %in% sc_sections) %>%
    arrange(strat_order, position_mm)

  if (nrow(tam_data) < 10 || nrow(sc_data) < 10) return(NULL)

  # Apply smoothing
  for (proxy in c("Ca_Ti", "Fe_Ti", "Mn_Ti")) {
    smooth_col <- paste0(proxy, "_smooth")
    tam_data[[smooth_col]] <- rollmean(tam_data[[proxy]], 5, fill = NA, align = "center")
    sc_data[[smooth_col]] <- rollmean(sc_data[[proxy]], 5, fill = NA, align = "center")
  }

  # Combine with site labels
  tam_data$site <- "TAM"
  sc_data$site <- "SC"

  # Calculate cumulative depth per site
  tam_data <- tam_data %>%
    mutate(depth_cm = (cumulative_depth - min(cumulative_depth)) / 10)
  sc_data <- sc_data %>%
    mutate(depth_cm = (cumulative_depth - min(cumulative_depth)) / 10)

  plot_data <- bind_rows(tam_data, sc_data)

  # Create faceted comparison plot
  plots <- list()

  for (proxy in c("Ca_Ti", "Fe_Ti", "Mn_Ti")) {
    smooth_col <- paste0(proxy, "_smooth")
    color <- proxy_colors[proxy]
    label <- proxy_labels[proxy]

    p <- ggplot(plot_data, aes(x = depth_cm, color = site)) +
      geom_point(aes(y = .data[[proxy]]), alpha = 0.3, size = 0.8) +
      geom_line(aes(y = .data[[smooth_col]]), linewidth = 1, na.rm = TRUE) +
      facet_wrap(~site, ncol = 1, scales = "free_x") +
      scale_color_manual(values = site_colors) +
      scale_y_log10() +
      labs(x = "Depth (cm)", y = label, title = label) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(face = "bold", size = 10)
      )

    plots[[proxy]] <- p
  }

  combined <- wrap_plots(plots, ncol = 3) +
    plot_annotation(title = title)

  combined
}

# ==============================================================================
# Generate density comparison figure
# ==============================================================================

generate_density_comparison <- function() {
  # Get top sections
  tam_top <- xrf_valid %>%
    filter(section %in% c("TAM-5AB-6-7-C", "TAM-5AB-6-7-B-RUN2", "TAM-5AB-6-7-A")) %>%
    mutate(site = "TAM")

  sc_top <- xrf_valid %>%
    filter(section == "SC8-A") %>%
    mutate(site = "SC")

  plot_data <- bind_rows(tam_top, sc_top)

  plots <- list()

  for (proxy in c("Ca_Ti", "Fe_Ti", "Mn_Ti")) {
    label <- proxy_labels[proxy]

    p <- ggplot(plot_data, aes(x = .data[[proxy]], fill = site, color = site)) +
      geom_density(alpha = 0.4, linewidth = 0.8) +
      geom_rug(alpha = 0.3, length = unit(0.02, "npc"), sides = "b") +
      scale_fill_manual(values = site_colors, name = "Site") +
      scale_color_manual(values = site_colors, name = "Site") +
      scale_x_log10() +
      labs(x = label, y = "Density", title = label) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = if (proxy == "Mn_Ti") c(0.85, 0.85) else "none",
        legend.background = element_rect(fill = "white", color = "gray80")
      )

    plots[[proxy]] <- p
  }

  combined <- wrap_plots(plots, ncol = 3) +
    plot_annotation(
      title = "Cross-Proxy Distribution: Top of Cores Comparison",
      subtitle = "TAM GROUP3 (youngest TAM, ~12.9 Ma) vs SC8-A (youngest SC, ~13.3 Ma)"
    )

  combined
}

# ==============================================================================
# Generate crossplot figure (Mn/Ti vs Fe/Ti, colored by Ca/Ti)
# ==============================================================================

generate_crossplot <- function() {
  tam_top <- xrf_valid %>%
    filter(section %in% c("TAM-5AB-6-7-C", "TAM-5AB-6-7-B-RUN2", "TAM-5AB-6-7-A")) %>%
    mutate(site = "TAM")

  sc_top <- xrf_valid %>%
    filter(section == "SC8-A") %>%
    mutate(site = "SC")

  plot_data <- bind_rows(tam_top, sc_top)

  # Mn/Ti vs Fe/Ti crossplot
  p1 <- ggplot(plot_data, aes(x = Mn_Ti, y = Fe_Ti, color = site)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = site_colors, name = "Site") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Mn/Ti", y = "Fe/Ti",
         title = "Redox Crossplot: Mn/Ti vs Fe/Ti") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = c(0.85, 0.15),
      legend.background = element_rect(fill = "white", color = "gray80")
    )

  # Ca/Ti vs Mn/Ti crossplot
  p2 <- ggplot(plot_data, aes(x = Mn_Ti, y = Ca_Ti, color = site)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = site_colors, name = "Site") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Mn/Ti", y = "Ca/Ti",
         title = "Carbonate-Redox Crossplot: Ca/Ti vs Mn/Ti") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )

  # Fe/Ti vs Ca/Ti crossplot
  p3 <- ggplot(plot_data, aes(x = Ca_Ti, y = Fe_Ti, color = site)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = site_colors, name = "Site") +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Ca/Ti", y = "Fe/Ti",
         title = "Iron-Carbonate Crossplot: Fe/Ti vs Ca/Ti") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )

  combined <- (p1 | p2 | p3) +
    plot_annotation(
      title = "Cross-Proxy Relationships: Top of Cores",
      subtitle = "TAM (orange) shows lower Mn/Ti, higher Fe/Ti indicating more reducing conditions than SC (blue)"
    )

  combined
}

# ==============================================================================
# Generate all figures
# ==============================================================================

message("Generating cross-proxy figures...")

# Individual sections
fig1 <- generate_cross_proxy_figure("TAM-5AB-6-7-C", " (Youngest TAM, ~12.9 Ma)")
ggsave(file.path(output_dir, "tam_top_cross_proxy.png"), fig1,
       width = 14, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "tam_top_cross_proxy.pdf"), fig1,
       width = 14, height = 8, bg = "white")

fig2 <- generate_cross_proxy_figure("SC8-A", " (Youngest SC, ~13.3 Ma)")
ggsave(file.path(output_dir, "sc_top_cross_proxy.png"), fig2,
       width = 14, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "sc_top_cross_proxy.pdf"), fig2,
       width = 14, height = 8, bg = "white")

# Density comparison
fig3 <- generate_density_comparison()
ggsave(file.path(output_dir, "density_comparison.png"), fig3,
       width = 14, height = 5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "density_comparison.pdf"), fig3,
       width = 14, height = 5, bg = "white")

# Crossplots
fig4 <- generate_crossplot()
ggsave(file.path(output_dir, "crossplots.png"), fig4,
       width = 16, height = 5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "crossplots.pdf"), fig4,
       width = 16, height = 5, bg = "white")

# Statistics summary
tam_top <- xrf_valid %>%
  filter(section %in% c("TAM-5AB-6-7-C", "TAM-5AB-6-7-B-RUN2", "TAM-5AB-6-7-A"))
sc_top <- xrf_valid %>%
  filter(section == "SC8-A")

stats <- tribble(
  ~Proxy, ~TAM_median, ~TAM_mean, ~SC_median, ~SC_mean, ~Ratio,
  "Ca/Ti", median(tam_top$Ca_Ti, na.rm=TRUE), mean(tam_top$Ca_Ti, na.rm=TRUE),
          median(sc_top$Ca_Ti, na.rm=TRUE), mean(sc_top$Ca_Ti, na.rm=TRUE),
          median(tam_top$Ca_Ti, na.rm=TRUE) / median(sc_top$Ca_Ti, na.rm=TRUE),
  "Fe/Ti", median(tam_top$Fe_Ti, na.rm=TRUE), mean(tam_top$Fe_Ti, na.rm=TRUE),
          median(sc_top$Fe_Ti, na.rm=TRUE), mean(sc_top$Fe_Ti, na.rm=TRUE),
          median(tam_top$Fe_Ti, na.rm=TRUE) / median(sc_top$Fe_Ti, na.rm=TRUE),
  "Mn/Ti", median(tam_top$Mn_Ti, na.rm=TRUE), mean(tam_top$Mn_Ti, na.rm=TRUE),
          median(sc_top$Mn_Ti, na.rm=TRUE), mean(sc_top$Mn_Ti, na.rm=TRUE),
          median(tam_top$Mn_Ti, na.rm=TRUE) / median(sc_top$Mn_Ti, na.rm=TRUE)
)

write_csv(stats, file.path(output_dir, "cross_proxy_stats.csv"))

message("Done! Files saved to: ", output_dir)
message("\nSummary statistics:")
print(stats)
