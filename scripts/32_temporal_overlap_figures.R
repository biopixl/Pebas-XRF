#!/usr/bin/env Rscript
# ==============================================================================
# Pebas-XRF: Temporal Overlap Figure Generator
# ==============================================================================
# Generates Shiny-style figures (core image + XRF traces) for the temporal
# overlap period (13.275-13.446 Ma) comparing TAM and SC sites
# ==============================================================================

library(tidyverse)
library(magick)
library(patchwork)
library(here)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

base_path <- here::here()
data_path <- file.path(base_path, "TAM-SC-IsaacA")
output_path <- file.path(base_path, "output")
fig_path <- file.path(output_path, "figures", "temporal_overlap")

# Create output directory
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

# Image enhancement settings (from Shiny app)
BRIGHTNESS_BOOST <- 80
CONTRAST_BOOST <- 25
GAMMA_CORRECTION <- 2.0

# Temporal framework
# TAM: 12.935-13.446 Ma (511 kyr, 398 cm, SR = 778 cm/Ma)
# SC: 13.275-14.298 Ma (1023 kyr, 424 cm, SR = 414 cm/Ma)
# Overlap: 13.275-13.446 Ma (171 kyr)

# Section configuration (all sections for figures)
section_config <- tribble(
  ~section, ~optical_path, ~group, ~site, ~strat_order,
  # GROUP1 - TAM (oldest)
  "TAM-1-2-3B-A", "GROUP1/TAM-1-2-3B/TAM-1-2-3B-A", "GROUP1", "TAM", 1,
  "TAM-1-2-3B-B", "GROUP1/TAM-1-2-3B/TAM-1-2-3B-B", "GROUP1", "TAM", 2,
  "TAM-1-2-3B-C", "GROUP1/TAM-1-2-3B/TAM-1-2-3B-C", "GROUP1", "TAM", 3,
  # GROUP2 - TAM
  "TAM-3A-4-5CDE-A", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-A", "GROUP2", "TAM", 4,
  "TAM-3A-4-5CDE-B", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-B", "GROUP2", "TAM", 5,
  "TAM-3A-4-5CDE-RUN2-C", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-RUN2-C", "GROUP2", "TAM", 6,
  "TAM-3A-4-5CDE-RUN2-D", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-RUN2-D", "GROUP2", "TAM", 7,
  "TAM-3A-4-5CDE-RUN2-E", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-RUN2-E", "GROUP2", "TAM", 8,
  # GROUP3 - TAM (youngest TAM)
  "TAM-5AB-6-7-A", "GROUP3/TAM-5AB-6-7/TAM-5AB-6-7-A", "GROUP3", "TAM", 9,
  "TAM-5AB-6-7-B-RUN2", "GROUP3/TAM-5AB-6-7/TAM-5AB-6-7-B-RUN2", "GROUP3", "TAM", 10,
  "TAM-5AB-6-7-C", "GROUP3/TAM-5AB-6-7/TAM-5AB-6-7-C", "GROUP3", "TAM", 11,
  # GROUP4 - SC (oldest SC)
  "SC-1ABC-2-3C-A-RUN1", "GROUP4/SC-1ABC-2-3C-RUN1/SC-1ABC-2-3C-A-RUN1", "GROUP4", "SC", 12,
  "SC-1ABC-2-3C-RUN2-B", "GROUP4/SC-1ABC-2-3C-RUN2/SC-1ABC-2-3C-RUN2-B", "GROUP4", "SC", 13,
  "SC-1ABC-2-3C-RUN2-C", "GROUP4/SC-1ABC-2-3C-RUN2/SC-1ABC-2-3C-RUN2-C", "GROUP4", "SC", 14,
  # GROUP5 - SC
  "SC-3AB-4ABCD-A", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-A", "GROUP5", "SC", 15,
  "SC-3AB-4ABCD-B", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-B", "GROUP5", "SC", 16,
  "SC-3AB-4ABCD-C", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-C", "GROUP5", "SC", 17,
  "SC-3AB-4ABCD-D", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-D", "GROUP5", "SC", 18,
  "SC-3AB-4ABCD-RUN2-D", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-RUN2-D", "GROUP5", "SC", 19,
  "SC-3AB-4ABCD-RUN2-E", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-RUN2-E", "GROUP5", "SC", 20,
  "SC-3AB-4ABCD-RUN2-F", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-RUN2-F", "GROUP5", "SC", 21,
  # GROUP6 - SC
  "SC-5-6-7ABC-A", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-A", "GROUP6", "SC", 22,
  "SC-5-6-7ABC-B", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-B", "GROUP6", "SC", 23,
  "SC-5-6-7ABC-C", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-C", "GROUP6", "SC", 24,
  "SC-5-6-7ABC-D", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-D", "GROUP6", "SC", 25,
  "SC-5-6-7ABC-E", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-E", "GROUP6", "SC", 26,
  # GROUP7 - SC (youngest)
  "SC8-A", "GROUP7/SC8-A/SC8-A", "GROUP7", "SC", 27
)

# Proxy colors (from Shiny app)
proxy_colors <- c(
  "Ca_Ti" = "#2166ac", "Fe_Mn" = "#7b3294", "Mn_Ti" = "#637939",
  "K_Ti" = "#e08214", "Zr_Rb" = "#1a9850", "Sr_Ca" = "#d95f02"
)

proxy_labels <- c(
  "Ca_Ti" = "Ca/Ti", "Fe_Mn" = "Fe/Mn", "Mn_Ti" = "Mn/Ti",
  "K_Ti" = "K/Ti", "Zr_Rb" = "Zr/Rb", "Sr_Ca" = "Sr/Ca"
)

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("Loading XRF data...")
xrf_data <- read_csv(file.path(output_path, "tables", "xrf_data_stacked.csv"),
                     show_col_types = FALSE) %>%
  mutate(
    Sr_Ca = Sr / Ca,
    Mn_Ti = Mn / Ti
  )

xrf_valid <- xrf_data %>% filter(qc_pass, !excluded)

exclusion_zones <- read_csv(file.path(base_path, "data", "exclusion_zones.csv"),
                            show_col_types = FALSE) %>%
  filter(!is.na(exclude_start_mm))

message(sprintf("Loaded %d valid measurements", nrow(xrf_valid)))

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

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

brighten_image_dynamic <- function(img, brightness = 80, contrast = 25, gamma = 2.0) {
  img %>%
    image_level(black_point = 2, white_point = 60, mid_point = 1/gamma) %>%
    image_normalize() %>%
    image_modulate(brightness = 100 + brightness) %>%
    image_contrast(sharpen = contrast)
}

process_core_image <- function(section_name, optical_path, data_start, data_end,
                               brightness = 80, contrast = 25, gamma = 2.0) {
  section_dir <- file.path(data_path, optical_path)
  img_path <- file.path(section_dir, "optical.tif")
  doc_path <- file.path(section_dir, "document.txt")

  if (!file.exists(img_path)) {
    message(sprintf("Image not found: %s", img_path))
    return(NULL)
  }
  params <- read_document_params(doc_path)
  if (is.null(params)) {
    message(sprintf("Document params not found: %s", doc_path))
    return(NULL)
  }

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
  img_bright <- brighten_image_dynamic(img_cropped, brightness, contrast, gamma)

  target_width <- max(crop_width * 4, 2000)
  img_hires <- image_resize(img_bright, sprintf("%dx", target_width))

  list(image = img_hires, section = section_name,
       xrf_start = data_start, xrf_end = data_end)
}

# ==============================================================================
# GENERATE SINGLE SECTION FIGURE
# ==============================================================================

generate_section_figure <- function(section_name, proxies = c("Mn_Ti", "Ca_Ti", "Fe_Mn"),
                                    smooth_window = 5, show_exclusions = TRUE) {

  sect_info <- section_config %>% filter(section == section_name)
  if (nrow(sect_info) == 0) {
    message(sprintf("Section not found: %s", section_name))
    return(NULL)
  }

  # Get data for this section
  sect_data <- xrf_valid %>%
    filter(section == section_name) %>%
    arrange(position_mm)

  if (nrow(sect_data) == 0) {
    message(sprintf("No valid data for section: %s", section_name))
    return(NULL)
  }

  # Get full data range
  data_start <- min(sect_data$position_mm)
  data_end <- max(sect_data$position_mm)

  # Convert position to depth in cm (relative to section start)
  sect_data <- sect_data %>%
    mutate(depth_cm = (position_mm - data_start) / 10)

  depth_max <- max(sect_data$depth_cm, na.rm = TRUE)
  depth_min <- 0

  # Apply smoothing
  for (proxy in proxies) {
    if (proxy %in% names(sect_data)) {
      smooth_col <- paste0(proxy, "_smooth")
      sect_data[[smooth_col]] <- zoo::rollmean(sect_data[[proxy]],
                                                smooth_window,
                                                fill = NA, align = "center")
    }
  }

  # Get exclusion zones for this section
  sect_excl <- exclusion_zones %>%
    filter(section == section_name) %>%
    mutate(
      start_cm = (exclude_start_mm - data_start) / 10,
      end_cm = (exclude_end_mm - data_start) / 10
    ) %>%
    filter(end_cm > 0, start_cm < depth_max)

  # Theme
  theme_strat <- function(base_size = 11) {
    theme_minimal(base_size = base_size) +
      theme(
        plot.title = element_text(size = base_size + 1, face = "bold"),
        axis.title = element_text(size = base_size),
        axis.text = element_text(size = base_size - 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        legend.position = "none",
        plot.margin = margin(5, 10, 5, 5)
      )
  }

  plot_list <- list()

  # Process core image
  img_result <- tryCatch(
    process_core_image(section_name, sect_info$optical_path, data_start, data_end),
    error = function(e) {
      message(sprintf("Error processing image: %s", e$message))
      NULL
    }
  )

  if (!is.null(img_result)) {
    # Rotate and convert to raster
    img_rotated <- image_rotate(img_result$image, 90)
    core_raster <- as.raster(img_rotated)

    core_plot <- ggplot() +
      annotation_raster(core_raster,
                        xmin = 0, xmax = 1,
                        ymin = depth_min, ymax = depth_max) +
      scale_y_reverse(limits = c(depth_max, depth_min), expand = c(0, 0)) +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = NULL, y = "Depth (cm)", title = "Core") +
      theme_strat() +
      theme(axis.text.x = element_blank(),
            panel.background = element_rect(fill = "gray30", color = NA))

    # Add exclusion shading
    if (show_exclusions && nrow(sect_excl) > 0) {
      for (i in seq_len(nrow(sect_excl))) {
        core_plot <- core_plot + annotate("rect",
                                           xmin = -Inf, xmax = Inf,
                                           ymin = sect_excl$start_cm[i],
                                           ymax = sect_excl$end_cm[i],
                                           fill = "gray50", alpha = 0.6)
      }
    }

    plot_list$core <- core_plot
  }

  # Create proxy panels
  for (i in seq_along(proxies)) {
    proxy <- proxies[i]
    smooth_col <- paste0(proxy, "_smooth")

    if (!proxy %in% names(sect_data)) {
      message(sprintf("Proxy not found: %s", proxy))
      next
    }

    color <- proxy_colors[proxy]
    if (is.na(color)) color <- "#333333"
    label <- proxy_labels[proxy]
    if (is.na(label)) label <- proxy

    p <- ggplot(sect_data, aes(y = depth_cm)) +
      geom_point(aes_string(x = proxy), color = color, alpha = 0.4, size = 1) +
      scale_x_log10() +
      scale_y_reverse(limits = c(depth_max, depth_min), expand = c(0.01, 0))

    if (smooth_col %in% names(sect_data)) {
      p <- p + geom_path(aes_string(x = smooth_col), color = color, linewidth = 1.2, na.rm = TRUE)
    }

    p <- p +
      labs(x = label,
           y = if (i == 1 && is.null(plot_list$core)) "Depth (cm)" else NULL,
           title = label) +
      theme_strat()

    if (i > 1 || !is.null(plot_list$core)) {
      p <- p + theme(axis.text.y = element_blank())
    }

    plot_list[[proxy]] <- p
  }

  # Combine panels
  if (length(plot_list) == 0) return(NULL)

  n_panels <- length(plot_list)
  if ("core" %in% names(plot_list)) {
    widths <- c(0.6, rep(1, n_panels - 1))
  } else {
    widths <- rep(1, n_panels)
  }

  combined <- wrap_plots(plot_list, nrow = 1, widths = widths) +
    plot_annotation(
      title = sprintf("%s (%s)", section_name, sect_info$site),
      subtitle = sprintf("%d measurements | %.1f cm",
                         nrow(sect_data), depth_max),
      caption = "Gray shading = excluded zones"
    )

  combined
}

# ==============================================================================
# GENERATE SIDE-BY-SIDE COMPARISON FIGURE
# ==============================================================================

generate_comparison_figure <- function(tam_sections, sc_sections,
                                       proxies = c("Mn_Ti", "Ca_Ti", "Fe_Mn"),
                                       smooth_window = 5) {

  # Process TAM data
  tam_data <- xrf_valid %>%
    filter(section %in% tam_sections) %>%
    arrange(section, position_mm)

  # Process SC data
  sc_data <- xrf_valid %>%
    filter(section %in% sc_sections) %>%
    arrange(section, position_mm)

  # Get TAM section info
  tam_config <- section_config %>%
    filter(section %in% tam_sections) %>%
    arrange(strat_order)

  # Get SC section info
  sc_config <- section_config %>%
    filter(section %in% sc_sections) %>%
    arrange(strat_order)

  # Calculate cumulative depth for TAM
  tam_ranges <- tam_data %>%
    group_by(section) %>%
    summarise(data_start = min(position_mm), data_end = max(position_mm), .groups = "drop")

  cumulative_offset <- 0
  tam_processed <- tibble()

  for (sect in tam_config$section) {
    sect_data <- tam_data %>% filter(section == sect)
    if (nrow(sect_data) == 0) next

    range_info <- tam_ranges %>% filter(section == sect)
    data_start <- range_info$data_start

    sect_data <- sect_data %>%
      mutate(
        relative_depth = position_mm - data_start,
        cumulative_depth = relative_depth + cumulative_offset,
        cumulative_depth_cm = cumulative_depth / 10
      )

    cumulative_offset <- cumulative_offset + (range_info$data_end - data_start) + 30
    tam_processed <- bind_rows(tam_processed, sect_data)
  }

  # Calculate cumulative depth for SC
  sc_ranges <- sc_data %>%
    group_by(section) %>%
    summarise(data_start = min(position_mm), data_end = max(position_mm), .groups = "drop")

  cumulative_offset <- 0
  sc_processed <- tibble()

  for (sect in sc_config$section) {
    sect_data <- sc_data %>% filter(section == sect)
    if (nrow(sect_data) == 0) next

    range_info <- sc_ranges %>% filter(section == sect)
    data_start <- range_info$data_start

    sect_data <- sect_data %>%
      mutate(
        relative_depth = position_mm - data_start,
        cumulative_depth = relative_depth + cumulative_offset,
        cumulative_depth_cm = cumulative_depth / 10
      )

    cumulative_offset <- cumulative_offset + (range_info$data_end - data_start) + 30
    sc_processed <- bind_rows(sc_processed, sect_data)
  }

  # Apply smoothing
  for (proxy in proxies) {
    if (proxy %in% names(tam_processed)) {
      smooth_col <- paste0(proxy, "_smooth")
      tam_processed[[smooth_col]] <- zoo::rollmean(tam_processed[[proxy]],
                                                    smooth_window,
                                                    fill = NA, align = "center")
      sc_processed[[smooth_col]] <- zoo::rollmean(sc_processed[[proxy]],
                                                   smooth_window,
                                                   fill = NA, align = "center")
    }
  }

  # Theme
  theme_strat <- function(base_size = 10) {
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

  # Generate TAM plots
  tam_depth_max <- max(tam_processed$cumulative_depth_cm, na.rm = TRUE)
  tam_plots <- list()

  for (i in seq_along(proxies)) {
    proxy <- proxies[i]
    smooth_col <- paste0(proxy, "_smooth")

    if (!proxy %in% names(tam_processed)) next

    color <- proxy_colors[proxy]
    if (is.na(color)) color <- "#333333"
    label <- proxy_labels[proxy]
    if (is.na(label)) label <- proxy

    p <- ggplot(tam_processed, aes(y = cumulative_depth_cm)) +
      geom_point(aes_string(x = proxy), color = color, alpha = 0.4, size = 0.8) +
      scale_x_log10() +
      scale_y_reverse(limits = c(tam_depth_max, 0), expand = c(0.01, 0))

    if (smooth_col %in% names(tam_processed)) {
      p <- p + geom_path(aes_string(x = smooth_col), color = color, linewidth = 1, na.rm = TRUE)
    }

    p <- p +
      labs(x = label,
           y = if (i == 1) "Depth (cm)" else NULL,
           title = if (i == 1) "TAM (Tamshiyacu)" else NULL) +
      theme_strat()

    if (i > 1) {
      p <- p + theme(axis.text.y = element_blank())
    }

    tam_plots[[proxy]] <- p
  }

  # Generate SC plots
  sc_depth_max <- max(sc_processed$cumulative_depth_cm, na.rm = TRUE)
  sc_plots <- list()

  for (i in seq_along(proxies)) {
    proxy <- proxies[i]
    smooth_col <- paste0(proxy, "_smooth")

    if (!proxy %in% names(sc_processed)) next

    color <- proxy_colors[proxy]
    if (is.na(color)) color <- "#333333"
    label <- proxy_labels[proxy]
    if (is.na(label)) label <- proxy

    p <- ggplot(sc_processed, aes(y = cumulative_depth_cm)) +
      geom_point(aes_string(x = proxy), color = color, alpha = 0.4, size = 0.8) +
      scale_x_log10() +
      scale_y_reverse(limits = c(sc_depth_max, 0), expand = c(0.01, 0))

    if (smooth_col %in% names(sc_processed)) {
      p <- p + geom_path(aes_string(x = smooth_col), color = color, linewidth = 1, na.rm = TRUE)
    }

    p <- p +
      labs(x = label,
           y = if (i == 1) "Depth (cm)" else NULL,
           title = if (i == 1) "SC (Santa Corina)" else NULL) +
      theme_strat()

    if (i > 1) {
      p <- p + theme(axis.text.y = element_blank())
    }

    sc_plots[[proxy]] <- p
  }

  # Combine TAM plots
  tam_combined <- wrap_plots(tam_plots, nrow = 1, widths = rep(1, length(tam_plots)))

  # Combine SC plots
  sc_combined <- wrap_plots(sc_plots, nrow = 1, widths = rep(1, length(sc_plots)))

  # Stack TAM and SC
  final <- tam_combined / sc_combined +
    plot_annotation(
      title = "Temporal Overlap Period (13.275–13.446 Ma)",
      subtitle = sprintf("TAM: %s | SC: %s",
                         paste(tam_sections, collapse = ", "),
                         paste(sc_sections, collapse = ", ")),
      caption = "Both sites recording same ~171 kyr interval"
    )

  final
}

# ==============================================================================
# GENERATE FIGURES
# ==============================================================================

message("\n=== Generating Temporal Overlap Figures ===\n")

# Figure 1: TAM-1-2-3B-A (oldest TAM overlap section)
message("Generating Figure 1: TAM-1-2-3B-A...")
fig1 <- generate_section_figure("TAM-1-2-3B-A",
                                 proxies = c("Mn_Ti", "Ca_Ti", "Fe_Mn"))
if (!is.null(fig1)) {
  ggsave(file.path(fig_path, "fig1_TAM-1-2-3B-A.png"), fig1,
         width = 14, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(fig_path, "fig1_TAM-1-2-3B-A.pdf"), fig1,
         width = 14, height = 8, bg = "white")
  message("  Saved fig1_TAM-1-2-3B-A")
}

# Figure 2: TAM-1-2-3B-B
message("Generating Figure 2: TAM-1-2-3B-B...")
fig2 <- generate_section_figure("TAM-1-2-3B-B",
                                 proxies = c("Mn_Ti", "Ca_Ti", "Fe_Mn"))
if (!is.null(fig2)) {
  ggsave(file.path(fig_path, "fig2_TAM-1-2-3B-B.png"), fig2,
         width = 14, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(fig_path, "fig2_TAM-1-2-3B-B.pdf"), fig2,
         width = 14, height = 8, bg = "white")
  message("  Saved fig2_TAM-1-2-3B-B")
}

# Figure 3: TAM-1-2-3B-C
message("Generating Figure 3: TAM-1-2-3B-C...")
fig3 <- generate_section_figure("TAM-1-2-3B-C",
                                 proxies = c("Mn_Ti", "Ca_Ti", "Fe_Mn"))
if (!is.null(fig3)) {
  ggsave(file.path(fig_path, "fig3_TAM-1-2-3B-C.png"), fig3,
         width = 14, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(fig_path, "fig3_TAM-1-2-3B-C.pdf"), fig3,
         width = 14, height = 8, bg = "white")
  message("  Saved fig3_TAM-1-2-3B-C")
}

# Figure 4: SC8-A (youngest SC overlap section)
message("Generating Figure 4: SC8-A...")
fig4 <- generate_section_figure("SC8-A",
                                 proxies = c("Mn_Ti", "Ca_Ti", "Fe_Mn"))
if (!is.null(fig4)) {
  ggsave(file.path(fig_path, "fig4_SC8-A.png"), fig4,
         width = 14, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(fig_path, "fig4_SC8-A.pdf"), fig4,
         width = 14, height = 8, bg = "white")
  message("  Saved fig4_SC8-A")
}

# Figure 5: Side-by-side comparison (all overlap sections)
message("Generating Figure 5: Side-by-side comparison...")
fig5 <- generate_comparison_figure(
  tam_sections = c("TAM-1-2-3B-A", "TAM-1-2-3B-B", "TAM-1-2-3B-C"),
  sc_sections = c("SC8-A"),
  proxies = c("Mn_Ti", "Ca_Ti", "Fe_Mn")
)
if (!is.null(fig5)) {
  ggsave(file.path(fig_path, "fig5_overlap_comparison.png"), fig5,
         width = 14, height = 12, dpi = 300, bg = "white")
  ggsave(file.path(fig_path, "fig5_overlap_comparison.pdf"), fig5,
         width = 14, height = 12, bg = "white")
  message("  Saved fig5_overlap_comparison")
}

# ==============================================================================
# CASE STUDY FIGURES
# ==============================================================================

message("\n=== Generating Case Study Figures ===\n")

# Case Study 1: TAM-3A-4-5CDE-A (high Mn/Ti example - anoxia)
message("Generating Case Study 1: High Mn/Ti (Anoxic conditions)...")
case1 <- generate_section_figure("TAM-3A-4-5CDE-A",
                                  proxies = c("Mn_Ti", "Fe_Mn", "Ca_Ti"))
if (!is.null(case1)) {
  ggsave(file.path(fig_path, "case1_anoxic_TAM-3A-4-5CDE-A.png"), case1,
         width = 14, height = 8, dpi = 300, bg = "white")
  message("  Saved case1_anoxic_TAM-3A-4-5CDE-A")
}

# Case Study 2: SC-3AB-4ABCD-A (carbonate-rich example)
message("Generating Case Study 2: High Ca/Ti (Carbonate input)...")
case2 <- generate_section_figure("SC-3AB-4ABCD-A",
                                  proxies = c("Ca_Ti", "Mn_Ti", "Fe_Mn"))
if (!is.null(case2)) {
  ggsave(file.path(fig_path, "case2_carbonate_SC-3AB-4ABCD-A.png"), case2,
         width = 14, height = 8, dpi = 300, bg = "white")
  message("  Saved case2_carbonate_SC-3AB-4ABCD-A")
}

# Case Study 3: TAM-5AB-6-7-A (environmental transition)
message("Generating Case Study 3: Environmental transition...")
case3 <- generate_section_figure("TAM-5AB-6-7-A",
                                  proxies = c("Mn_Ti", "Ca_Ti", "Fe_Mn", "K_Ti"))
if (!is.null(case3)) {
  ggsave(file.path(fig_path, "case3_transition_TAM-5AB-6-7-A.png"), case3,
         width = 16, height = 8, dpi = 300, bg = "white")
  message("  Saved case3_transition_TAM-5AB-6-7-A")
}

# Case Study 4: SC-5-6-7ABC-A (detailed redox variations)
message("Generating Case Study 4: Redox variations...")
case4 <- generate_section_figure("SC-5-6-7ABC-A",
                                  proxies = c("Fe_Mn", "Mn_Ti", "Ca_Ti"))
if (!is.null(case4)) {
  ggsave(file.path(fig_path, "case4_redox_SC-5-6-7ABC-A.png"), case4,
         width = 14, height = 8, dpi = 300, bg = "white")
  message("  Saved case4_redox_SC-5-6-7ABC-A")
}

# ==============================================================================
# TOP OF CORES FIGURES (Youngest sections from each site)
# ==============================================================================

message("\n=== Generating Top of Cores Figures ===\n")

# TAM youngest section (top of TAM core, ~12.9 Ma)
message("Generating TAM top section: TAM-5AB-6-7-C...")
fig_tam_top <- generate_section_figure("TAM-5AB-6-7-C",
                                        proxies = c("Mn_Ti", "Ca_Ti", "Fe_Mn"))
if (!is.null(fig_tam_top)) {
  ggsave(file.path(fig_path, "top_TAM-5AB-6-7-C.png"), fig_tam_top,
         width = 14, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(fig_path, "top_TAM-5AB-6-7-C.pdf"), fig_tam_top,
         width = 14, height = 8, bg = "white")
  message("  Saved top_TAM-5AB-6-7-C")
}

# TAM-5AB-6-7-B-RUN2
message("Generating TAM section: TAM-5AB-6-7-B-RUN2...")
fig_tam_b <- generate_section_figure("TAM-5AB-6-7-B-RUN2",
                                      proxies = c("Mn_Ti", "Ca_Ti", "Fe_Mn"))
if (!is.null(fig_tam_b)) {
  ggsave(file.path(fig_path, "top_TAM-5AB-6-7-B-RUN2.png"), fig_tam_b,
         width = 14, height = 8, dpi = 300, bg = "white")
  message("  Saved top_TAM-5AB-6-7-B-RUN2")
}

# SC youngest section (top of SC core, ~13.3 Ma) - already generated as fig4_SC8-A

# Additional good SC sections for case studies
message("Generating SC section: SC-1ABC-2-3C-A-RUN1 (oldest SC)...")
fig_sc_oldest <- generate_section_figure("SC-1ABC-2-3C-A-RUN1",
                                          proxies = c("Mn_Ti", "Ca_Ti", "Fe_Mn"))
if (!is.null(fig_sc_oldest)) {
  ggsave(file.path(fig_path, "bottom_SC-1ABC-2-3C-A-RUN1.png"), fig_sc_oldest,
         width = 14, height = 8, dpi = 300, bg = "white")
  message("  Saved bottom_SC-1ABC-2-3C-A-RUN1")
}

# Full TAM stacked comparison
message("Generating full TAM GROUP3 comparison (youngest TAM)...")
fig_tam_group3 <- generate_comparison_figure(
  tam_sections = c("TAM-5AB-6-7-A", "TAM-5AB-6-7-B-RUN2", "TAM-5AB-6-7-C"),
  sc_sections = c("SC8-A"),
  proxies = c("Mn_Ti", "Ca_Ti", "Fe_Mn")
)
if (!is.null(fig_tam_group3)) {
  ggsave(file.path(fig_path, "top_cores_comparison.png"), fig_tam_group3,
         width = 14, height = 12, dpi = 300, bg = "white")
  ggsave(file.path(fig_path, "top_cores_comparison.pdf"), fig_tam_group3,
         width = 14, height = 12, bg = "white")
  message("  Saved top_cores_comparison")
}

message("\n=== Figure Generation Complete ===")
message(sprintf("Figures saved to: %s", fig_path))
