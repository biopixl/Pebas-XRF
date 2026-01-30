# ==============================================================================
# Pebas-XRF: Case Study Reports for All Groups
# ==============================================================================
# Generate standardized case study figures with brightened optical images
# for all 7 core groups (TAM: GROUP1-3, SC: GROUP4-7)
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
fig_path <- file.path(output_path, "figures", "case_study")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

# Image brightness adjustment
BRIGHTNESS_BOOST <- 40
CONTRAST_BOOST <- 10
GAMMA_CORRECTION <- 1.3  # >1 brightens midtones
WINDOW_SIZE <- 5

# Load exclusion zones for foam masking
exclusion_zones <- read_csv(file.path(base_path, "data", "exclusion_zones.csv"),
                            show_col_types = FALSE) %>%
  filter(!is.na(exclude_start_mm) & notes == "foam")

# ==============================================================================
# SECTION TO PATH MAPPING (based on actual directory structure)
# ==============================================================================

section_paths <- tribble(
  ~section, ~optical_path,
  # GROUP1 - TAM
  "TAM-1-2-3B-A", "GROUP1/TAM-1-2-3B/TAM-1-2-3B-A",
  "TAM-1-2-3B-B", "GROUP1/TAM-1-2-3B/TAM-1-2-3B-B",
  "TAM-1-2-3B-C", "GROUP1/TAM-1-2-3B/TAM-1-2-3B-C",

  # GROUP2 - TAM
  "TAM-3A-4-5CDE-A", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-A",
  "TAM-3A-4-5CDE-B", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-B",
  "TAM-3A-4-5CDE-RUN2-C", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-RUN2-C",
  "TAM-3A-4-5CDE-RUN2-D", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-RUN2-D",
  "TAM-3A-4-5CDE-RUN2-E", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-RUN2-E",

  # GROUP3 - TAM
  "TAM-5AB-6-7-A", "GROUP3/TAM-5AB-6-7/TAM-5AB-6-7-A",
  "TAM-5AB-6-7-B-RUN2", "GROUP3/TAM-5AB-6-7/TAM-5AB-6-7-B-RUN2",
  "TAM-5AB-6-7-C", "GROUP3/TAM-5AB-6-7/TAM-5AB-6-7-C",

  # GROUP4 - SC
  "SC-1ABC-2-3C-A-RUN1", "GROUP4/SC-1ABC-2-3C-RUN1/SC-1ABC-2-3C-A-RUN1",
  "SC-1ABC-2-3C-RUN2-B", "GROUP4/SC-1ABC-2-3C-RUN2/SC-1ABC-2-3C-RUN2-B",
  "SC-1ABC-2-3C-RUN2-C", "GROUP4/SC-1ABC-2-3C-RUN2/SC-1ABC-2-3C-RUN2-C",

  # GROUP5 - SC
  "SC-3AB-4ABCD-A", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-A",
  "SC-3AB-4ABCD-B", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-B",
  "SC-3AB-4ABCD-C", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-C",
  "SC-3AB-4ABCD-D", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-D",
  "SC-3AB-4ABCD-RUN2-D", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-RUN2-D",
  "SC-3AB-4ABCD-RUN2-E", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-RUN2-E",
  "SC-3AB-4ABCD-RUN2-F", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-RUN2-F",

  # GROUP6 - SC
  "SC-5-6-7ABC-A", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-A",
  "SC-5-6-7ABC-B", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-B",
  "SC-5-6-7ABC-C", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-C",
  "SC-5-6-7ABC-D", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-D",
  "SC-5-6-7ABC-E", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-E",

  # GROUP7 - SC
  "SC8-A", "GROUP7/SC8-A/SC8-A"
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
  # Apply gamma correction first (brightens midtones without clipping)
  if (gamma != 1.0) {
    img <- image_level(img, black_point = 0, white_point = 100, mid_point = 1/gamma)
  }

  # Then brightness and contrast
  img %>%
    image_modulate(brightness = 100 + brightness) %>%
    image_contrast(sharpen = contrast)
}

mask_foam_sections <- function(img, section_name, xrf_start, xrf_end, pixels_per_mm) {
  # Get foam zones for this section
  foam_zones <- exclusion_zones %>%
    filter(section == section_name)

  if (nrow(foam_zones) == 0) return(img)

  img_info <- image_info(img)

  for (i in seq_len(nrow(foam_zones))) {
    zone <- foam_zones[i, ]

    # Convert foam mm positions to pixels relative to XRF range
    foam_start_px <- round((zone$exclude_start_mm - xrf_start) * pixels_per_mm)
    foam_end_px <- round((zone$exclude_end_mm - xrf_start) * pixels_per_mm)

    # Skip if outside visible range
    if (foam_end_px < 0 || foam_start_px > img_info$width) next

    # Clamp to image bounds
    foam_start_px <- max(0, foam_start_px)
    foam_end_px <- min(img_info$width, foam_end_px)
    foam_width <- foam_end_px - foam_start_px

    if (foam_width > 0) {
      # Create gray mask for foam region
      foam_mask <- image_blank(foam_width, img_info$height, color = "#808080")
      img <- image_composite(img, foam_mask, offset = sprintf("+%d+0", foam_start_px))
    }
  }

  return(img)
}

process_section <- function(section_name, optical_path, mask_foam = TRUE) {
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

  xrf_start_px <- max(0, round((params$xrf_start - params$optical_start) * pixels_per_mm))
  xrf_end_px <- min(img_info$width, round((params$xrf_end - params$optical_start) * pixels_per_mm))
  crop_width <- xrf_end_px - xrf_start_px

  if (crop_width <= 0) return(NULL)

  crop_geom <- sprintf("%dx%d+%d+0", crop_width, img_info$height, xrf_start_px)
  img_cropped <- image_crop(img, crop_geom)

  # Calculate pixels_per_mm for the cropped image (same as original)
  cropped_pixels_per_mm <- crop_width / (params$xrf_end - params$xrf_start)

  # Apply foam masking before brightening
  if (mask_foam) {
    img_cropped <- mask_foam_sections(img_cropped, section_name,
                                       params$xrf_start, params$xrf_end,
                                       cropped_pixels_per_mm)
  }

  # Apply brightness and gamma correction
  img_bright <- brighten_image(img_cropped)

  return(list(
    image = img_bright,
    section = section_name,
    xrf_start = params$xrf_start,
    xrf_end = params$xrf_end
  ))
}

#' Create composite strip with proper gap handling
#'
#' @param section_images List of processed section images with xrf_start/xrf_end
#' @param target_width Width of output strip in pixels
#' @param total_depth_range Total depth range from XRF data (for scaling)
#' @param pixels_per_mm Scaling factor for image height
#' @return Composite image with gaps properly represented
create_composite_strip <- function(section_images, target_width = 150,
                                    total_depth_range = NULL, pixels_per_mm = 1.0) {
  if (length(section_images) == 0) return(NULL)


  # Sort sections by xrf_start position
  section_images <- section_images[order(sapply(section_images, function(s) s$xrf_start))]

  # If no total range provided, calculate from sections (old behavior)
  if (is.null(total_depth_range)) {
    processed <- map(section_images, function(s) {
      img_rotated <- image_rotate(s$image, 90)
      image_resize(img_rotated, sprintf("%dx", target_width))
    })
    composite <- image_append(do.call(c, processed), stack = TRUE)
    return(composite)
  }

  # Calculate total height based on actual depth range
  total_height <- round(total_depth_range * pixels_per_mm)

  # Start with blank canvas (gray background for gaps)
  composite <- image_blank(target_width, total_height, color = "#404040")

  # Get the minimum position (start of depth scale)
  min_pos <- min(sapply(section_images, function(s) s$xrf_start))

  # Place each section at correct position
  for (s in section_images) {
    # Rotate and resize section image
    img_rotated <- image_rotate(s$image, 90)
    img_resized <- image_resize(img_rotated, sprintf("%dx", target_width))

    # Calculate vertical position for this section
    section_depth <- s$xrf_end - s$xrf_start
    section_height <- round(section_depth * pixels_per_mm)
    y_offset <- round((s$xrf_start - min_pos) * pixels_per_mm)

    # Resize to exact section height to match depth scale
    img_scaled <- image_resize(img_resized, sprintf("%dx%d!", target_width, section_height))

    # Composite onto canvas at correct position
    composite <- image_composite(composite, img_scaled,
                                  offset = sprintf("+0+%d", y_offset))
  }

  return(composite)
}

# ==============================================================================
# LOAD XRF DATA
# ==============================================================================

message("Loading XRF data...")
xrf_data <- read_csv(file.path(output_path, "tables", "xrf_data_stacked.csv"),
                     show_col_types = FALSE)

# Get unique groups
groups <- unique(xrf_data$group)

# ==============================================================================
# PROCESS ALL GROUPS
# ==============================================================================

facies_colors <- c("Shell-rich" = "#2166AC", "Carbonate" = "#67A9CF",
                   "Mixed" = "#D1E5F0", "Clastic" = "#B2182B")

all_group_stats <- list()

for (grp in groups) {
  message(sprintf("\n=== Processing %s ===", grp))

  # Get sections for this group
  group_data <- xrf_data %>% filter(group == grp)
  group_sections <- unique(group_data$section)

  # Process optical images for each section
  section_images <- list()

  for (sec in group_sections) {
    path_info <- section_paths %>% filter(section == sec)

    if (nrow(path_info) == 0) {
      message(sprintf("  No path mapping for %s", sec))
      next
    }

    result <- process_section(sec, path_info$optical_path[1])

    if (!is.null(result)) {
      section_images[[sec]] <- result
      message(sprintf("  %s: %.0f-%.0f mm (brightened)",
                      sec, result$xrf_start, result$xrf_end))
    } else {
      message(sprintf("  %s: image not found or invalid", sec))
    }
  }

  if (length(section_images) == 0) {
    message(sprintf("  No valid images for %s, creating plot without core image...", grp))
  }

  # Calculate actual depth range from XRF data for proper image scaling
  depth_min_full <- min(group_data$position_mm)
  depth_max_full <- max(group_data$position_mm)
  total_depth_range <- depth_max_full - depth_min_full

  # Create composite if we have images - now with proper gap handling
  composite <- NULL
  if (length(section_images) > 0) {
    # Use 0.5 pixels per mm for reasonable image size
    pixels_per_mm <- 0.5
    composite <- create_composite_strip(section_images, target_width = 150,
                                         total_depth_range = total_depth_range,
                                         pixels_per_mm = pixels_per_mm)
    composite_path <- file.path(fig_path, sprintf("core_optical_%s.png", grp))
    image_write(composite, composite_path)
    message(sprintf("  Saved: core_optical_%s.png (with gap alignment)", grp))
  }

  # Prepare group data for plotting - FILTER EXCLUDED POINTS for spectra
  # Keep all points for facies column but filter for curves
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

  # Filter to non-excluded points for spectra curves
  group_data <- group_data_all %>%
    filter(!excluded) %>%
    mutate(
      Ca_Ti_filt = zoo::rollmean(Ca_Ti, WINDOW_SIZE, fill = NA, align = "center"),
      Fe_Mn_filt = zoo::rollmean(Fe_Mn, WINDOW_SIZE, fill = NA, align = "center"),
      K_Ti_filt = zoo::rollmean(K_Ti, WINDOW_SIZE, fill = NA, align = "center"),
      Zr_Rb_filt = zoo::rollmean(Zr_Rb, WINDOW_SIZE, fill = NA, align = "center")
    )

  # Calculate statistics (from non-excluded data only)
  site <- ifelse(grepl("^TAM", group_sections[1]), "TAM", "SC")
  site_name <- ifelse(site == "TAM", "Tamshiyacu", "Santa Corina")

  stats <- group_data %>%
    summarise(
      group = grp,
      site = site,
      n = n(),
      n_excluded = sum(group_data_all$excluded),
      depth_range_mm = depth_max_full - depth_min_full,
      Ca_Ti_mean = mean(Ca_Ti, na.rm = TRUE),
      Ca_Ti_sd = sd(Ca_Ti, na.rm = TRUE),
      Fe_Mn_mean = mean(Fe_Mn, na.rm = TRUE),
      Fe_Mn_sd = sd(Fe_Mn, na.rm = TRUE),
      pct_reducing = mean(Fe_Mn > 50, na.rm = TRUE) * 100,
      pct_shell_rich = mean(Ca_Ti > 10, na.rm = TRUE) * 100,
      pct_clastic = mean(Ca_Ti < 2, na.rm = TRUE) * 100
    )
  all_group_stats[[grp]] <- stats

  # Use FULL depth range for image mapping (consistent with composite)
  depth_min <- depth_min_full
  depth_max <- depth_max_full

  # Create figure panels
  plot_list <- list()

  # Panel 1: Core image (if available)
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

  # Panel: Facies (use ALL data including excluded for complete depth coverage)
  facies_y_label <- if (is.null(composite)) "Depth (cm)" else ""

  # Create facies data with foam indicator
  facies_data <- group_data_all %>%
    mutate(
      facies_display = ifelse(excluded, "Excluded", as.character(facies)),
      facies_display = factor(facies_display,
                               levels = c("Shell-rich", "Carbonate", "Mixed", "Clastic", "Excluded"))
    )

  facies_colors_ext <- c(facies_colors, "Excluded" = "#808080")

  plot_list$facies <- ggplot(facies_data, aes(y = depth_cm)) +
    geom_tile(aes(x = 0.5, fill = facies_display), width = 1, height = 0.4) +
    scale_fill_manual(values = facies_colors_ext, name = "Facies", drop = FALSE) +
    scale_y_reverse() +
    labs(x = NULL, y = facies_y_label, title = "Facies") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank(), legend.position = "none")

  if (!is.null(composite)) {
    plot_list$facies <- plot_list$facies + theme(axis.text.y = element_blank())
  }

  # Panel: Ca/Ti
  plot_list$cati <- ggplot(group_data, aes(y = depth_cm)) +
    geom_point(aes(x = Ca_Ti), alpha = 0.3, size = 0.8, color = "darkgreen") +
    geom_path(aes(x = Ca_Ti_filt), color = "darkgreen", linewidth = 0.6, na.rm = TRUE) +
    geom_vline(xintercept = c(2, 5, 10), linetype = "dashed", color = "gray60", linewidth = 0.3) +
    scale_y_reverse() +
    labs(x = "Ca/Ti", y = NULL, title = "Carbonate") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank())

  # Panel: Fe/Mn
  plot_list$femn <- ggplot(group_data, aes(y = depth_cm)) +
    geom_ribbon(aes(xmin = 0, xmax = Fe_Mn_filt,
                    fill = ifelse(Fe_Mn_filt > 50, "Reducing", "Oxic")),
                alpha = 0.4, na.rm = TRUE) +
    geom_path(aes(x = Fe_Mn_filt), color = "purple", linewidth = 0.6, na.rm = TRUE) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "red", linewidth = 0.5) +
    scale_fill_manual(values = c("Reducing" = "purple", "Oxic" = "orange"), guide = "none") +
    scale_y_reverse() +
    labs(x = "Fe/Mn", y = NULL, title = "Redox") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank())

  # Panel: K/Ti
  plot_list$kti <- ggplot(group_data, aes(y = depth_cm)) +
    geom_point(aes(x = K_Ti), alpha = 0.3, size = 0.8, color = "darkorange") +
    geom_path(aes(x = K_Ti_filt), color = "darkorange", linewidth = 0.6, na.rm = TRUE) +
    scale_y_reverse() +
    labs(x = "K/Ti", y = NULL, title = "Weathering") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank())

  # Panel: Zr/Rb
  plot_list$zrrb <- ggplot(group_data, aes(y = depth_cm)) +
    geom_point(aes(x = Zr_Rb), alpha = 0.3, size = 0.8, color = "darkred") +
    geom_path(aes(x = Zr_Rb_filt), color = "darkred", linewidth = 0.6, na.rm = TRUE) +
    scale_y_reverse() +
    labs(x = "Zr/Rb", y = NULL, title = "Grain Size") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank())

  # Combine panels
  if (!is.null(composite)) {
    fig <- plot_list$core + plot_list$facies + plot_list$cati +
      plot_list$femn + plot_list$kti + plot_list$zrrb +
      plot_layout(widths = c(0.6, 0.3, 1, 1, 1, 1))
  } else {
    fig <- plot_list$facies + plot_list$cati +
      plot_list$femn + plot_list$kti + plot_list$zrrb +
      plot_layout(widths = c(0.4, 1, 1, 1, 1))
  }

  fig <- fig + plot_annotation(
    title = sprintf("%s %s: Core Image with Geochemical Stratigraphy", site, grp),
    subtitle = sprintf("%s | n=%d measurements (%d excluded) | %.1f cm depth | %.0f%% reducing",
                       site_name, nrow(group_data), stats$n_excluded,
                       (depth_max - depth_min)/10, stats$pct_reducing),
    caption = "Facies: Ca/Ti thresholds (2, 5, 10) | Fe/Mn = 50 (oxic/reducing) | Gray = excluded foam/gaps",
    theme = theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      plot.caption = element_text(size = 8, color = "gray50")
    )
  )

  # Save
  fig_file <- sprintf("integrated_%s.png", grp)
  ggsave(file.path(fig_path, fig_file), fig,
         width = 14, height = 10, dpi = 300, bg = "white")
  message(sprintf("  Saved: %s", fig_file))
}

# ==============================================================================
# SUMMARY TABLE
# ==============================================================================

summary_table <- bind_rows(all_group_stats) %>%
  distinct() %>%  # Remove duplicate rows from grouped summarise
  arrange(site, group)

write_csv(summary_table, file.path(output_path, "tables", "all_groups_summary.csv"))

message("\n=== ALL GROUPS SUMMARY ===")
print(summary_table)

message("\n=== CASE STUDY GENERATION COMPLETE ===")
