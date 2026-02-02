# ==============================================================================
# Pebas-XRF: Image Brightness Enhancement for Publication
# ==============================================================================
# Enhances core optical images with improved brightness/contrast for publication
# ==============================================================================

library(magick)
library(tidyverse)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

base_path <- here::here()
input_path <- file.path(base_path, "output", "figures", "case_study")
output_path <- file.path(base_path, "output", "figures", "enhanced")
manuscript_path <- file.path(base_path, "manuscript", "figures")

# Create output directory
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# Enhancement parameters - aggressive for dark core images
BRIGHTNESS_BOOST <- 40    # Percentage increase in brightness (0-100)
CONTRAST_BOOST <- 25      # Percentage increase in contrast (0-100)
GAMMA <- 1.8              # Gamma correction (>1 = brighter midtones)

# ==============================================================================
# ENHANCE CORE OPTICAL IMAGES
# ==============================================================================

message("=== Enhancing core optical images for publication ===\n")

# Find all core optical images
core_images <- list.files(input_path, pattern = "core_optical_GROUP.*\\.png$", full.names = TRUE)

message(sprintf("Found %d core optical images to enhance", length(core_images)))

for (img_path in core_images) {
  img_name <- basename(img_path)
  message(sprintf("  Processing: %s", img_name))

  # Read image
  img <- image_read(img_path)

  # Apply enhancements
  img_enhanced <- img %>%
    # Adjust brightness and contrast
    image_modulate(brightness = 100 + BRIGHTNESS_BOOST) %>%
    # Adjust contrast using sigmoidal
    image_contrast(sharpen = CONTRAST_BOOST) %>%
    # Apply gamma correction for midtone adjustment
    image_modulate(brightness = 100) %>%
    image_normalize()  # Stretch histogram for better dynamic range

  # Aggressive level adjustment for dark core images
  # Equivalent to Photoshop levels: stretch histogram significantly
  img_enhanced <- img %>%
    image_level(black_point = 5, white_point = 70, mid_point = 2.0) %>%
    image_normalize()  # Stretch to full dynamic range

  # Save enhanced image
  output_file <- file.path(output_path, gsub("\\.png$", "_enhanced.png", img_name))
  image_write(img_enhanced, output_file, format = "png", quality = 100)
  message(sprintf("    Saved: %s", basename(output_file)))
}

# ==============================================================================
# ENHANCE INTEGRATED FIGURES
# ==============================================================================

message("\n=== Enhancing integrated stratigraphic figures ===\n")

integrated_images <- list.files(manuscript_path, pattern = "integrated_GROUP.*\\.png$", full.names = TRUE)

message(sprintf("Found %d integrated figures to enhance", length(integrated_images)))

for (img_path in integrated_images) {
  img_name <- basename(img_path)
  message(sprintf("  Processing: %s", img_name))

  # Read image
  img <- image_read(img_path)

  # For integrated figures, apply moderate enhancement
  # Balance between brightening core image and preserving spectra colors
  img_enhanced <- img %>%
    image_level(black_point = 2, white_point = 85, mid_point = 1.4)

  # Save enhanced version
  output_file <- file.path(output_path, gsub("\\.png$", "_enhanced.png", img_name))
  image_write(img_enhanced, output_file, format = "png", quality = 100)
  message(sprintf("    Saved: %s", basename(output_file)))

  # Also save to manuscript folder as _pub version
  pub_file <- file.path(manuscript_path, gsub("\\.png$", "_pub.png", img_name))
  image_write(img_enhanced, pub_file, format = "png", quality = 100)
  message(sprintf("    Publication copy: %s", basename(pub_file)))
}

# ==============================================================================
# UPDATE REPORT
# ==============================================================================

message("\n=== Image Enhancement Complete ===")
message(sprintf("Enhanced images saved to: %s", output_path))
message(sprintf("Publication versions saved to: %s", manuscript_path))
message("\nEnhancement settings applied:")
message(sprintf("  - Level adjustment: black=0, white=95-97%%, gamma=1.1-1.2"))
message("\nTo use in LaTeX, update figure paths to *_pub.png versions")
message("Or use original images with LaTeX adjustbox brightness adjustment")
