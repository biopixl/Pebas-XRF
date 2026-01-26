# ==============================================================================
# Pebas-XRF: Core Image Composite for Case Study
# ==============================================================================
# Create composite optical image aligned with XRF geochemical data
# GROUP3: TAM-5AB-6-7 sections A, B-RUN2, C
# ==============================================================================

library(tidyverse)
library(magick)
library(patchwork)
library(grid)
library(png)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

base_path <- here::here()
data_path <- file.path(base_path, "TAM-SC-IsaacA", "GROUP3", "TAM-5AB-6-7")
output_path <- file.path(base_path, "output")
fig_path <- file.path(output_path, "figures", "case_study")

# Section parameters from document.txt files
sections <- tribble(
  ~section,              ~xrf_start, ~xrf_end, ~optical_start, ~optical_end,
  "TAM-5AB-6-7-A",       38,         158,      -0.7,           1222.3,
  "TAM-5AB-6-7-B-RUN2",  204,        659,      -0.7,           1222.3,
  "TAM-5AB-6-7-C",       678,        1175,     -0.7,           1222.3
)

# ==============================================================================
# PROCESS OPTICAL IMAGES
# ==============================================================================

message("Processing optical images for GROUP3...")

process_section_image <- function(section_name, xrf_start, xrf_end, optical_start, optical_end) {

  # Load optical image
  img_path <- file.path(data_path, section_name, "optical.tif")

  if (!file.exists(img_path)) {
    message(sprintf("  Warning: %s not found", img_path))
    return(NULL)
  }

  img <- image_read(img_path)
  img_info <- image_info(img)

  # Calculate pixel positions for XRF range
  # Optical image spans optical_start to optical_end (mm)
  optical_range <- optical_end - optical_start
  pixels_per_mm <- img_info$width / optical_range

  # Convert XRF positions to pixel coordinates
  xrf_start_px <- round((xrf_start - optical_start) * pixels_per_mm)
  xrf_end_px <- round((xrf_end - optical_start) * pixels_per_mm)
  crop_width <- xrf_end_px - xrf_start_px

  # Crop image to XRF scan range
  # geometry format: WxH+X+Y
  crop_geom <- sprintf("%dx%d+%d+0", crop_width, img_info$height, xrf_start_px)
  img_cropped <- image_crop(img, crop_geom)

  message(sprintf("  %s: XRF %.0f-%.0f mm, cropped to %d pixels",
                  section_name, xrf_start, xrf_end, crop_width))

  return(list(
    image = img_cropped,
    section = section_name,
    xrf_start = xrf_start,
    xrf_end = xrf_end,
    depth_range = xrf_end - xrf_start
  ))
}

# Process each section
section_images <- pmap(sections, function(section, xrf_start, xrf_end, optical_start, optical_end) {
  process_section_image(section, xrf_start, xrf_end, optical_start, optical_end)
})

# Filter out any NULL results
section_images <- compact(section_images)

# ==============================================================================
# CREATE COMPOSITE STRIP
# ==============================================================================

message("\nCreating composite core strip...")

# Stack images vertically (convert to same width first)
target_width <- 200  # pixels for the strip

composite_images <- map(section_images, function(s) {
  # Resize to consistent width, maintaining aspect ratio for height
  img_info <- image_info(s$image)
  scale_factor <- target_width / img_info$width
  new_height <- round(img_info$height * scale_factor)

  # For core images, we want to show the length, so rotate 90 degrees
  # Original: width = along core, height = across core
  # We want: width = across core (narrow), height = along core (long)
  img_rotated <- image_rotate(s$image, 90)
  img_resized <- image_resize(img_rotated, sprintf("%dx", target_width))

  return(img_resized)
})

# Stack vertically
composite <- image_append(do.call(c, composite_images), stack = TRUE)

# Save composite strip
image_write(composite, file.path(fig_path, "core_optical_GROUP3.png"))
message(sprintf("  Saved: core_optical_GROUP3.png (%s)",
                paste(image_info(composite)$width, "x", image_info(composite)$height)))

# ==============================================================================
# CREATE FIGURE WITH OPTICAL + GEOCHEMISTRY
# ==============================================================================

message("\nGenerating integrated figure with optical image...")

# Load XRF data
xrf_data <- read_csv(file.path(output_path, "tables", "xrf_data_stacked.csv"),
                     show_col_types = FALSE)

g3 <- xrf_data %>%
  filter(group == "GROUP3") %>%
  arrange(position_mm) %>%
  mutate(
    depth_cm = position_mm / 10,
    Ca_Ti_filt = zoo::rollmean(Ca_Ti, 5, fill = NA, align = "center"),
    Fe_Mn_filt = zoo::rollmean(Fe_Mn, 5, fill = NA, align = "center"),
    facies = case_when(
      Ca_Ti > 10 ~ "Shell-rich",
      Ca_Ti > 5 ~ "Carbonate",
      Ca_Ti > 2 ~ "Mixed",
      TRUE ~ "Clastic"
    ),
    facies = factor(facies, levels = c("Shell-rich", "Carbonate", "Mixed", "Clastic"))
  )

# Read the composite image for plotting
core_img <- image_read(file.path(fig_path, "core_optical_GROUP3.png"))
core_info <- image_info(core_img)

# Convert to raster for ggplot
core_raster <- as.raster(core_img)

# Calculate depth mapping
# Total XRF depth range
depth_min <- min(g3$position_mm)
depth_max <- max(g3$position_mm)

# Create the integrated plot
facies_colors <- c("Shell-rich" = "#2166AC", "Carbonate" = "#67A9CF",
                   "Mixed" = "#D1E5F0", "Clastic" = "#B2182B")

# Panel 1: Core image
p_core <- ggplot() +
  annotation_raster(core_raster,
                    xmin = 0, xmax = 1,
                    ymin = depth_min/10, ymax = depth_max/10) +
  scale_y_reverse(limits = c(depth_max/10, depth_min/10)) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Depth (cm)", title = "Core") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )

# Panel 2: Facies interpretation
p_facies <- ggplot(g3, aes(y = depth_cm)) +
  geom_tile(aes(x = 0.5, fill = facies), width = 1, height = 0.4) +
  scale_fill_manual(values = facies_colors, name = "Facies") +
  scale_y_reverse() +
  labs(x = NULL, y = NULL, title = "Facies") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Panel 3: Ca/Ti
p_cati <- ggplot(g3, aes(y = depth_cm)) +
  geom_point(aes(x = Ca_Ti), alpha = 0.3, size = 0.8, color = "darkgreen") +
  geom_path(aes(x = Ca_Ti_filt), color = "darkgreen", linewidth = 0.6) +
  geom_vline(xintercept = c(2, 5, 10), linetype = "dashed", color = "gray60", linewidth = 0.3) +
  scale_y_reverse() +
  labs(x = "Ca/Ti", y = NULL, title = "Carbonate") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_blank())

# Panel 4: Fe/Mn
p_femn <- ggplot(g3, aes(y = depth_cm)) +
  geom_ribbon(aes(xmin = 0, xmax = Fe_Mn_filt,
                  fill = ifelse(Fe_Mn_filt > 50, "Reducing", "Oxic")),
              alpha = 0.4, na.rm = TRUE) +
  geom_path(aes(x = Fe_Mn_filt), color = "purple", linewidth = 0.6) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_fill_manual(values = c("Reducing" = "purple", "Oxic" = "orange"), guide = "none") +
  scale_y_reverse() +
  labs(x = "Fe/Mn", y = NULL, title = "Redox") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_blank())

# Combine with optical image
fig_integrated <- p_core + p_facies + p_cati + p_femn +
  plot_layout(widths = c(0.8, 0.4, 1, 1)) +
  plot_annotation(
    title = "TAM GROUP3: Core Image with Geochemical Stratigraphy",
    subtitle = "Optical image aligned with XRF proxy data | Tamshiyacu, Pebas Formation",
    caption = "Facies: Shell-rich (Ca/Ti>10), Carbonate (5-10), Mixed (2-5), Clastic (<2) | Red line: Fe/Mn=50 oxic/reducing threshold",
    theme = theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      plot.caption = element_text(size = 8, color = "gray50")
    )
  )

ggsave(file.path(fig_path, "fig5_core_integrated_GROUP3.png"), fig_integrated,
       width = 12, height = 12, dpi = 300, bg = "white")
ggsave(file.path(fig_path, "fig5_core_integrated_GROUP3.pdf"), fig_integrated,
       width = 12, height = 12, bg = "white")

message("  Saved: fig5_core_integrated_GROUP3.png/pdf")

# ==============================================================================
# INDIVIDUAL SECTION IMAGES (FULL RESOLUTION)
# ==============================================================================

message("\nSaving individual section images...")

for (s in section_images) {
  # Save cropped section at full resolution
  section_filename <- sprintf("optical_%s.png", gsub("-", "_", s$section))

  # Rotate for vertical display
  img_rotated <- image_rotate(s$image, 90)

  image_write(img_rotated, file.path(fig_path, section_filename))
  message(sprintf("  Saved: %s", section_filename))
}

message("\n=== CORE IMAGE PROCESSING COMPLETE ===")
