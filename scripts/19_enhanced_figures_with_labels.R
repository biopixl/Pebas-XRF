#!/usr/bin/env Rscript
# ==============================================================================
# Pebas-XRF: Enhanced Figures with Section Labels
# ==============================================================================
# Generate publication-quality figures with clear section labels for core
# visualization in the manuscript.
# ==============================================================================

library(tidyverse)
library(magick)
library(zoo)
library(patchwork)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

base_path <- here::here()
data_path <- file.path(base_path, "TAM-SC-IsaacA")
output_path <- file.path(base_path, "output")
fig_path <- file.path(base_path, "manuscript", "figures")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

# Image settings - SIGNIFICANTLY INCREASED for visibility
BRIGHTNESS_BOOST <- 80      # Increased from 40
CONTRAST_BOOST <- 30        # Increased from 10
GAMMA_CORRECTION <- 2.0     # Increased from 1.3 (brightens midtones significantly)
WINDOW_SIZE <- 5

# Consistent color for ALL excluded/gap zones
EXCLUSION_GRAY <- "#c8c8c8"

# Load exclusion zones
exclusion_zones <- read_csv(file.path(base_path, "data", "exclusion_zones.csv"),
                            show_col_types = FALSE) %>%
  filter(!is.na(exclude_start_mm) & notes == "foam")

# ==============================================================================
# SECTION TO PATH MAPPING
# ==============================================================================

section_paths <- tribble(
  ~section, ~optical_path, ~short_label,
  # GROUP1 - TAM
  "TAM-1-2-3B-A", "GROUP1/TAM-1-2-3B/TAM-1-2-3B-A", "A",
  "TAM-1-2-3B-B", "GROUP1/TAM-1-2-3B/TAM-1-2-3B-B", "B",
  "TAM-1-2-3B-C", "GROUP1/TAM-1-2-3B/TAM-1-2-3B-C", "C",

  # GROUP2 - TAM
  "TAM-3A-4-5CDE-A", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-A", "A",
  "TAM-3A-4-5CDE-B", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-B", "B",
  "TAM-3A-4-5CDE-RUN2-C", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-RUN2-C", "C",
  "TAM-3A-4-5CDE-RUN2-D", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-RUN2-D", "D",
  "TAM-3A-4-5CDE-RUN2-E", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-RUN2-E", "E",

  # GROUP3 - TAM
  "TAM-5AB-6-7-A", "GROUP3/TAM-5AB-6-7/TAM-5AB-6-7-A", "A",
  "TAM-5AB-6-7-B-RUN2", "GROUP3/TAM-5AB-6-7/TAM-5AB-6-7-B-RUN2", "B",
  "TAM-5AB-6-7-C", "GROUP3/TAM-5AB-6-7/TAM-5AB-6-7-C", "C",

  # GROUP4 - SC
  "SC-1ABC-2-3C-A-RUN1", "GROUP4/SC-1ABC-2-3C-RUN1/SC-1ABC-2-3C-A-RUN1", "A",
  "SC-1ABC-2-3C-RUN2-B", "GROUP4/SC-1ABC-2-3C-RUN2/SC-1ABC-2-3C-RUN2-B", "B",
  "SC-1ABC-2-3C-RUN2-C", "GROUP4/SC-1ABC-2-3C-RUN2/SC-1ABC-2-3C-RUN2-C", "C",

  # GROUP5 - SC
  "SC-3AB-4ABCD-A", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-A", "A",
  "SC-3AB-4ABCD-B", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-B", "B",
  "SC-3AB-4ABCD-C", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-C", "C",
  "SC-3AB-4ABCD-D", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-D", "D",
  "SC-3AB-4ABCD-RUN2-D", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-RUN2-D", "D2",
  "SC-3AB-4ABCD-RUN2-E", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-RUN2-E", "E",
  "SC-3AB-4ABCD-RUN2-F", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-RUN2-F", "F",

  # GROUP6 - SC
  "SC-5-6-7ABC-A", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-A", "A",
  "SC-5-6-7ABC-B", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-B", "B",
  "SC-5-6-7ABC-C", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-C", "C",
  "SC-5-6-7ABC-D", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-D", "D",
  "SC-5-6-7ABC-E", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-E", "E",

  # GROUP7 - SC
  "SC8-A", "GROUP7/SC8-A/SC8-A", "A"
)

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

  return(list(
    xrf_start = xrf_start,
    xrf_end = xrf_end,
    optical_start = optical_start,
    optical_end = optical_end
  ))
}

brighten_image <- function(img, brightness = BRIGHTNESS_BOOST, contrast = CONTRAST_BOOST, gamma = GAMMA_CORRECTION) {
  # Apply histogram equalization to maximize contrast
  img <- image_equalize(img)

  # Aggressive gamma correction to lift dark sediment
  img <- image_level(img, black_point = 0, white_point = 100, mid_point = 1/5.0)

  # Strong brightness boost
  img <- image_modulate(img, brightness = 250)

  # Enhance contrast to restore definition
  img <- image_contrast(img, sharpen = 40)

  # Another gamma pass to lift remaining shadows
  img <- image_level(img, black_point = 0, white_point = 100, mid_point = 0.3)

  # Final brightness push
  img <- image_modulate(img, brightness = 150)

  # One more gamma for stubborn dark areas
  img <- image_level(img, black_point = 0, white_point = 100, mid_point = 0.5)

  # Replace remaining very dark pixels with uniform gray
  # Scanner background stays black even after brightening (fuzz=8 catches near-black)
  img <- image_transparent(img, color = "black", fuzz = 8)
  img <- image_background(img, color = EXCLUSION_GRAY)

  return(img)
}

mask_foam_sections <- function(img, section_name, xrf_start, xrf_end, pixels_per_mm) {
  foam_zones <- exclusion_zones %>%
    filter(section == section_name)

  if (nrow(foam_zones) == 0) return(img)

  img_info <- image_info(img)

  for (i in seq_len(nrow(foam_zones))) {
    zone <- foam_zones[i, ]

    foam_start_px <- round((zone$exclude_start_mm - xrf_start) * pixels_per_mm)
    foam_end_px <- round((zone$exclude_end_mm - xrf_start) * pixels_per_mm)

    if (foam_end_px < 0 || foam_start_px > img_info$width) next

    foam_start_px <- max(0, foam_start_px)
    foam_end_px <- min(img_info$width, foam_end_px)
    foam_width <- foam_end_px - foam_start_px

    if (foam_width > 0) {
      # Use consistent light gray for all excluded zones
      foam_mask <- image_blank(foam_width, img_info$height, color = EXCLUSION_GRAY)
      img <- image_composite(img, foam_mask, offset = sprintf("+%d+0", foam_start_px))
    }
  }

  return(img)
}

process_section <- function(section_name, optical_path, mask_foam = TRUE,
                            data_start = NULL, data_end = NULL) {
  section_dir <- file.path(data_path, optical_path)
  img_path <- file.path(section_dir, "optical.tif")
  doc_path <- file.path(section_dir, "document.txt")

  if (!file.exists(img_path)) {
    return(NULL)
  }

  params <- read_document_params(doc_path)
  if (is.null(params)) return(NULL)

  img <- image_read(img_path)
  img_info <- image_info(img)

  optical_range <- params$optical_end - params$optical_start
  pixels_per_mm <- img_info$width / optical_range

  xrf_start <- if (!is.null(data_start)) data_start else params$xrf_start
  xrf_end <- if (!is.null(data_end)) data_end else params$xrf_end

  xrf_start_px <- max(0, round((xrf_start - params$optical_start) * pixels_per_mm))
  xrf_end_px <- min(img_info$width, round((xrf_end - params$optical_start) * pixels_per_mm))
  crop_width <- xrf_end_px - xrf_start_px

  if (crop_width <= 0) return(NULL)

  # Crop horizontally (along core length)
  crop_geom <- sprintf("%dx%d+%d+0", crop_width, img_info$height, xrf_start_px)
  img_cropped <- image_crop(img, crop_geom)

  # Crop vertically to remove scanner background edges (keep central 60% of height)
  # The core is typically centered with black scanner background on top/bottom
  cropped_info <- image_info(img_cropped)
  margin_pct <- 0.20  # Remove 20% from top and bottom
  top_margin <- round(cropped_info$height * margin_pct)
  core_height <- round(cropped_info$height * (1 - 2 * margin_pct))
  vertical_crop <- sprintf("%dx%d+0+%d", cropped_info$width, core_height, top_margin)
  img_cropped <- image_crop(img_cropped, vertical_crop)

  cropped_pixels_per_mm <- crop_width / (xrf_end - xrf_start)

  if (mask_foam) {
    img_cropped <- mask_foam_sections(img_cropped, section_name,
                                       xrf_start, xrf_end,
                                       cropped_pixels_per_mm)
  }

  img_bright <- brighten_image(img_cropped)

  return(list(
    image = img_bright,
    section = section_name,
    xrf_start = xrf_start,
    xrf_end = xrf_end
  ))
}

#' Create composite strip (without embedded labels - labels added via ggplot)
create_composite_strip <- function(section_images, target_width = 150,
                                    total_depth_range = NULL,
                                    pixels_per_mm = 0.5) {
  if (length(section_images) == 0) return(NULL)

  # Sort sections by xrf_start position
  order_idx <- order(sapply(section_images, function(s) s$xrf_start))
  section_images <- section_images[order_idx]

  if (is.null(total_depth_range)) {
    total_depth_range <- max(sapply(section_images, function(s) s$xrf_end)) -
                         min(sapply(section_images, function(s) s$xrf_start))
  }

  # Calculate total height
  total_height <- round(total_depth_range * pixels_per_mm)

  # Create canvas - consistent light gray for all gaps (same as foam zones)
  composite <- image_blank(target_width, total_height, color = EXCLUSION_GRAY)

  # Get the minimum position
  min_pos <- min(sapply(section_images, function(s) s$xrf_start))

  # Place each section at correct position
  for (i in seq_along(section_images)) {
    s <- section_images[[i]]

    # Rotate and resize section image
    img_rotated <- image_rotate(s$image, 90)
    img_resized <- image_resize(img_rotated, sprintf("%dx", target_width))

    # Calculate vertical position
    section_depth <- s$xrf_end - s$xrf_start
    section_height <- round(section_depth * pixels_per_mm)
    y_offset <- round((s$xrf_start - min_pos) * pixels_per_mm)

    # Resize to exact section height
    img_scaled <- image_resize(img_resized, sprintf("%dx%d!", target_width, section_height))

    # Composite onto canvas
    composite <- image_composite(composite, img_scaled,
                                  offset = sprintf("+0+%d", y_offset))

    # Add section boundary line (except for last section)
    if (i < length(section_images)) {
      next_start <- section_images[[i + 1]]$xrf_start
      boundary_y <- round((s$xrf_end - min_pos) * pixels_per_mm)

      # Only draw thin subtle line if there's a gap (slightly darker than background)
      if (next_start > s$xrf_end + 5) {
        line <- image_blank(target_width, 2, color = "#a0a0a0")
        composite <- image_composite(composite, line,
                                      offset = sprintf("+0+%d", boundary_y))
      }
    }
  }

  return(composite)
}

# ==============================================================================
# LOAD XRF DATA
# ==============================================================================

message("Loading XRF data...")
xrf_data <- read_csv(file.path(output_path, "tables", "xrf_data_stacked.csv"),
                     show_col_types = FALSE)

groups <- unique(xrf_data$group)

# ==============================================================================
# GENERATE ENHANCED FIGURES WITH SECTION LABELS
# ==============================================================================

facies_colors <- c("Shell-rich" = "#2166AC", "Carbonate" = "#67A9CF",
                   "Mixed" = "#D1E5F0", "Clastic" = "#B2182B")

message("\n=== GENERATING ENHANCED FIGURES WITH SECTION LABELS ===\n")

for (grp in groups) {
  message(sprintf("Processing %s...", grp))

  # Get sections for this group
  group_data <- xrf_data %>% filter(group == grp)
  group_sections <- unique(group_data$section)

  # Calculate actual XRF data position ranges per section
  section_ranges <- group_data %>%
    group_by(section) %>%
    summarise(
      data_start = min(position_mm),
      data_end = max(position_mm),
      .groups = "drop"
    )

  # Create section labels lookup
  section_labels <- section_paths %>%
    filter(section %in% group_sections) %>%
    select(section, short_label) %>%
    deframe()

  # Process optical images
  section_images <- list()

  for (sec in group_sections) {
    path_info <- section_paths %>% filter(section == sec)

    if (nrow(path_info) == 0) next

    sec_range <- section_ranges %>% filter(section == sec)
    data_start <- sec_range$data_start
    data_end <- sec_range$data_end

    result <- process_section(sec, path_info$optical_path[1],
                               data_start = data_start, data_end = data_end)

    if (!is.null(result)) {
      section_images[[sec]] <- result
    }
  }

  # Calculate depth range
  depth_min_full <- min(group_data$position_mm)
  depth_max_full <- max(group_data$position_mm)
  total_depth_range <- depth_max_full - depth_min_full

  # Create composite (labels added via ggplot panel)
  composite <- NULL
  if (length(section_images) > 0) {
    pixels_per_mm <- 0.5
    composite <- create_composite_strip(
      section_images,
      target_width = 150,
      total_depth_range = total_depth_range,
      pixels_per_mm = pixels_per_mm
    )
  }

  # Prepare data for plotting
  group_data_all <- group_data %>%
    arrange(position_mm) %>%
    mutate(
      depth_cm = position_mm / 10,
      facies = case_when(
        Ca_Ti > 10 ~ "Shell-rich",
        Ca_Ti > 5 ~ "Carbonate",
        Ca_Ti > 2 ~ "Mixed",
        TRUE ~ "Clastic"
      ),
      facies = factor(facies, levels = c("Shell-rich", "Carbonate", "Mixed", "Clastic"))
    )

  group_data_filt <- group_data_all %>%
    filter(!excluded) %>%
    mutate(
      Ca_Ti_filt = zoo::rollmean(Ca_Ti, WINDOW_SIZE, fill = NA, align = "center"),
      Fe_Mn_filt = zoo::rollmean(Fe_Mn, WINDOW_SIZE, fill = NA, align = "center")
    )

  # Calculate statistics
  site <- ifelse(grepl("^TAM", group_sections[1]), "TAM", "SC")
  site_name <- ifelse(site == "TAM", "Tamshiyacu", "Santa Corina")

  stats <- group_data_filt %>%
    summarise(
      n = n(),
      n_excluded = sum(group_data_all$excluded),
      pct_reducing = mean(Fe_Mn > 50, na.rm = TRUE) * 100
    )

  # Define y-axis limits
  depth_min <- depth_min_full
  depth_max <- depth_max_full
  y_limits <- c(depth_max/10, depth_min/10)

  # Build section boundaries data for annotation
  section_bounds <- section_ranges %>%
    left_join(section_paths %>% select(section, short_label), by = "section") %>%
    mutate(
      mid_depth = (data_start + data_end) / 2 / 10
    )

  # Create figure panels
  plot_list <- list()

  # Panel 1: Core image with section labels
  if (!is.null(composite)) {
    core_raster <- as.raster(composite)

    plot_list$core <- ggplot() +
      annotation_raster(core_raster,
                        xmin = 0, xmax = 1,
                        ymin = depth_min/10, ymax = depth_max/10) +
      scale_y_reverse(limits = c(depth_max/10, depth_min/10)) +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = NULL, y = "Depth (cm)", title = "Core") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            panel.grid = element_blank())
  }

  # Panel 2: Section labels as separate panel
  plot_list$sections <- ggplot(section_bounds) +
    geom_segment(aes(x = 0, xend = 1, y = data_start/10, yend = data_start/10),
                 color = "gray70", linewidth = 0.3) +
    geom_text(aes(x = 0.5, y = mid_depth, label = short_label),
              size = 3, fontface = "bold") +
    scale_y_reverse(limits = y_limits) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = "Sec.") +
    theme_minimal(base_size = 10) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())

  # Panel 3: Facies
  facies_data <- group_data_all %>%
    mutate(
      facies_display = ifelse(excluded, "Excluded", as.character(facies)),
      facies_display = factor(facies_display,
                               levels = c("Shell-rich", "Carbonate", "Mixed", "Clastic", "Excluded"))
    )

  facies_colors_ext <- c(facies_colors, "Excluded" = EXCLUSION_GRAY)  # Consistent light gray

  plot_list$facies <- ggplot(facies_data, aes(y = depth_cm)) +
    geom_tile(aes(x = 0.5, fill = facies_display), width = 1, height = 0.4) +
    scale_fill_manual(values = facies_colors_ext, name = "Facies", drop = FALSE) +
    scale_y_reverse(limits = y_limits) +
    labs(x = NULL, y = NULL, title = "Facies") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank(), legend.position = "none",
          axis.text.y = element_blank())

  # Panel 4: Ca/Ti
  plot_list$cati <- ggplot(group_data_filt, aes(y = depth_cm)) +
    geom_point(aes(x = Ca_Ti), alpha = 0.3, size = 0.8, color = "darkgreen") +
    geom_path(aes(x = Ca_Ti_filt), color = "darkgreen", linewidth = 0.6, na.rm = TRUE) +
    geom_vline(xintercept = c(2, 5, 10), linetype = "dashed", color = "gray60", linewidth = 0.3) +
    scale_y_reverse(limits = y_limits) +
    labs(x = "Ca/Ti", y = NULL, title = "Carbonate") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank())

  # Panel 5: Fe/Mn
  plot_list$femn <- ggplot(group_data_filt, aes(y = depth_cm)) +
    geom_ribbon(aes(xmin = 0, xmax = Fe_Mn_filt,
                    fill = ifelse(Fe_Mn_filt > 50, "Reducing", "Oxic")),
                alpha = 0.4, na.rm = TRUE) +
    geom_path(aes(x = Fe_Mn_filt), color = "purple", linewidth = 0.6, na.rm = TRUE) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "red", linewidth = 0.5) +
    scale_fill_manual(values = c("Reducing" = "purple", "Oxic" = "orange"), guide = "none") +
    scale_y_reverse(limits = y_limits) +
    labs(x = "Fe/Mn", y = NULL, title = "Redox") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank())

  # Combine panels
  if (!is.null(composite)) {
    fig <- plot_list$core + plot_list$sections + plot_list$facies +
      plot_list$cati + plot_list$femn +
      plot_layout(widths = c(0.5, 0.15, 0.2, 1, 1))
  } else {
    fig <- plot_list$sections + plot_list$facies +
      plot_list$cati + plot_list$femn +
      plot_layout(widths = c(0.2, 0.3, 1, 1))
  }

  # Get section list for subtitle
  section_list <- paste(section_bounds$short_label, collapse = ", ")

  fig <- fig + plot_annotation(
    title = sprintf("%s: Integrated Stratigraphy with Section Labels", grp),
    subtitle = sprintf("%s | Sections: %s | n=%d (%d excluded) | %.0f%% reducing",
                       site_name, section_list, stats$n, stats$n_excluded, stats$pct_reducing),
    caption = "Facies: Ca/Ti thresholds (2, 5, 10) | Fe/Mn = 50 (oxic/reducing)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      plot.caption = element_text(size = 8, color = "gray50")
    )
  )

  # Save publication figure
  fig_file <- sprintf("integrated_%s_pub.png", grp)
  ggsave(file.path(fig_path, fig_file), fig,
         width = 12, height = 10, dpi = 300, bg = "white")
  message(sprintf("  Saved: %s", fig_file))
}

message("\n=== ENHANCED FIGURES COMPLETE ===")
message("Figures saved to: manuscript/figures/")
