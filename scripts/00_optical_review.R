# ==============================================================================
# Pebas-XRF: Optical Image Review for QC
# ==============================================================================
# Brighten optical images and add position scale for visual QC
# ==============================================================================

library(tidyverse)
library(tiff)
library(grid)
library(gridExtra)

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

base_path <- here::here()
data_path <- file.path(base_path, "TAM-SC-IsaacA")
output_path <- file.path(base_path, "output", "figures", "optical_review")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 2. FIND ALL OPTICAL IMAGES
# ==============================================================================

optical_files <- list.files(
  data_path,
  pattern = "optical\\.tif$",
  recursive = TRUE,
  full.names = TRUE
)

# Exclude 'copied' directories
optical_files <- optical_files[!grepl("/copied/", optical_files)]

message(sprintf("Found %d optical images", length(optical_files)))

# ==============================================================================
# 3. GET POSITION INFO FROM DOCUMENT.TXT
# ==============================================================================

get_scan_positions <- function(optical_path) {
  # Find corresponding document.txt
  doc_path <- file.path(dirname(optical_path), "document.txt")

  if (!file.exists(doc_path)) {
    return(list(start = NA, stop = NA, step = NA))
  }

  doc <- readLines(doc_path, warn = FALSE)

  # Parse start/stop coordinates
  coord_line <- grep("Start coordinate", doc, value = TRUE)
  if (length(coord_line) > 0) {
    parts <- strsplit(coord_line, "\t")[[1]]
    start <- as.numeric(parts[2])
    stop <- as.numeric(parts[4])
  } else {
    start <- NA
    stop <- NA
  }

  # Parse step size
  step_line <- grep("^Step size", doc, value = TRUE)[1]
  if (length(step_line) > 0 && !is.na(step_line)) {
    step <- as.numeric(strsplit(step_line, "\t")[[1]][2]) / 1000  # microns to mm
  } else {
    step <- 3  # default 3mm
  }

  return(list(start = start, stop = stop, step = step))
}

# ==============================================================================
# 4. PROCESS AND BRIGHTEN IMAGES
# ==============================================================================

process_optical_image <- function(optical_path, brightness_factor = 1.5) {

  section_name <- basename(dirname(optical_path))
  message(sprintf("Processing: %s", section_name))

  # Read TIF
  img <- tryCatch(
    readTIFF(optical_path, native = FALSE),
    error = function(e) {
      warning(sprintf("Failed to read: %s", optical_path))
      return(NULL)
    }
  )

  if (is.null(img)) return(NULL)

  # Get position info
  positions <- get_scan_positions(optical_path)

  # Brighten image
  img_bright <- pmin(img * brightness_factor, 1)

  # Get dimensions
  dims <- dim(img_bright)
  height <- dims[1]
  width <- dims[2]

  # Create output filename
  out_file <- file.path(output_path, paste0("optical_", section_name, ".png"))

  # Save as PNG with position scale annotation
  png(out_file, width = width, height = height + 100, res = 150)

  # Set up layout
  par(mar = c(0, 0, 2, 0))

  # Plot image
  if (length(dims) == 3) {
    # RGB image
    plot(0, 0, type = "n", xlim = c(0, width), ylim = c(0, height),
         asp = 1, xlab = "", ylab = "", axes = FALSE)
    rasterImage(img_bright, 0, 0, width, height)
  } else {
    # Grayscale
    image(t(img_bright[nrow(img_bright):1, ]), col = gray.colors(256),
          axes = FALSE, asp = height/width)
  }

  # Add title with position info
  if (!is.na(positions$start)) {
    title(sprintf("%s\nPosition: %.0f - %.0f mm",
                  section_name, positions$start, positions$stop),
          cex.main = 1.2)
  } else {
    title(section_name, cex.main = 1.2)
  }

  # Add scale bar
  if (!is.na(positions$start) && !is.na(positions$stop)) {
    # Calculate mm per pixel
    mm_range <- positions$stop - positions$start
    mm_per_px <- mm_range / width

    # Add scale markers every 100mm
    scale_interval <- 100  # mm
    px_interval <- scale_interval / mm_per_px

    # Draw scale
    for (pos_mm in seq(0, mm_range, by = scale_interval)) {
      px <- pos_mm / mm_per_px
      if (px < width) {
        segments(px, height * 0.98, px, height, col = "red", lwd = 2)
        text(px, height * 0.96, sprintf("%.0f", positions$start + pos_mm),
             col = "red", cex = 0.7, adj = c(0.5, 1))
      }
    }
  }

  dev.off()

  message(sprintf("  Saved: %s", basename(out_file)))

  return(out_file)
}

# ==============================================================================
# 5. PROCESS ALL IMAGES
# ==============================================================================

message("\nProcessing optical images (brightened)...")

processed_files <- map_chr(optical_files, function(f) {
  tryCatch(
    process_optical_image(f, brightness_factor = 1.8),
    error = function(e) {
      warning(sprintf("Error processing %s: %s", f, e$message))
      return(NA_character_)
    }
  )
})

message(sprintf("\n%d images processed", sum(!is.na(processed_files))))
message(sprintf("Output directory: %s", output_path))

# ==============================================================================
# 6. CREATE SUMMARY TABLE
# ==============================================================================

# Load result data to get position ranges
result_files <- list.files(
  data_path,
  pattern = "result\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)
result_files <- result_files[!grepl("/copied/", result_files)]

position_summary <- map_dfr(result_files, function(f) {
  section <- basename(dirname(f))

  # Read first few lines to get position range
  data <- tryCatch(
    read_tsv(f, skip = 2, n_max = 1000, show_col_types = FALSE),
    error = function(e) return(NULL)
  )

  if (is.null(data) || !"position (mm)" %in% names(data)) {
    return(tibble(section = section, start_mm = NA, end_mm = NA, n_points = NA))
  }

  tibble(
    section = section,
    start_mm = min(data$`position (mm)`, na.rm = TRUE),
    end_mm = max(data$`position (mm)`, na.rm = TRUE),
    n_points = nrow(data)
  )
})

# Print summary for reference
message("\n=== SECTION POSITION SUMMARY ===")
print(position_summary, n = 30)

# Save for reference
write_csv(position_summary, file.path(output_path, "section_positions.csv"))

message("\n--- NEXT STEPS ---")
message("1. Review brightened optical images in: output/figures/optical_review/")
message("2. Note position ranges (mm) where foam/gaps appear")
message("3. Edit data/exclusion_zones.csv with those ranges")
message("4. Re-run analysis pipeline")
