# ==============================================================================
# Convert DICOM CT data to TIFF format
# ==============================================================================

library(oro.dicom)
library(tiff)

# Paths
base_path <- here::here()
ct_path <- file.path(base_path, "Pebas-CT", "ST1")
output_path <- file.path(base_path, "Pebas-CT", "TIFF")

# Create output directories
dir.create(file.path(output_path, "SE1_TOP"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_path, "SE2_BTTM"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Function to convert a single DICOM file to TIFF
# ==============================================================================

convert_dicom_file <- function(dicom_file, output_dir, series_name) {
  tryCatch({
    # Read DICOM with pixel data
    dcm <- readDICOMFile(dicom_file, pixelData = TRUE)

    # Extract metadata
    hdr <- dcm$hdr

    # Get slice location and instance number
    slice_loc_idx <- which(hdr$name == "SliceLocation")
    instance_idx <- which(hdr$name == "InstanceNumber")

    if (length(slice_loc_idx) == 0 || length(instance_idx) == 0) {
      message(sprintf("  Skipping %s: missing metadata", basename(dicom_file)))
      return(FALSE)
    }

    slice_loc <- as.numeric(hdr$value[slice_loc_idx[1]])
    instance <- as.numeric(hdr$value[instance_idx[1]])

    # Get pixel data
    img <- dcm$img

    if (is.null(img) || length(img) == 0) {
      message(sprintf("  Skipping %s: no pixel data", basename(dicom_file)))
      return(FALSE)
    }

    # Normalize to 0-1 range for TIFF
    img_min <- min(img, na.rm = TRUE)
    img_max <- max(img, na.rm = TRUE)

    if (img_max > img_min) {
      img_norm <- (img - img_min) / (img_max - img_min)
    } else {
      img_norm <- matrix(0, nrow = nrow(img), ncol = ncol(img))
    }

    # Create filename with slice info
    tiff_name <- sprintf("%s_slice%03d_loc%+07.2fmm.tif", series_name, instance, slice_loc)
    tiff_path <- file.path(output_dir, tiff_name)

    # Write 16-bit TIFF
    writeTIFF(img_norm, tiff_path, bits.per.sample = 16)

    message(sprintf("  %s -> %s", basename(dicom_file), tiff_name))
    return(TRUE)

  }, error = function(e) {
    message(sprintf("  Error with %s: %s", basename(dicom_file), e$message))
    return(FALSE)
  })
}

# ==============================================================================
# Convert series
# ==============================================================================

convert_series <- function(series_dir, output_subdir, series_name) {
  input_dir <- file.path(ct_path, series_dir)
  output_dir <- file.path(output_path, output_subdir)

  # List DICOM files (exclude IM1 which is often just metadata)
  files <- list.files(input_dir, pattern = "^IM[0-9]+$", full.names = TRUE)
  files <- files[!grepl("IM1$", files)]  # Exclude IM1

  message(sprintf("\n%s: Converting %d DICOM files...", series_name, length(files)))

  success_count <- 0
  for (f in files) {
    if (convert_dicom_file(f, output_dir, series_name)) {
      success_count <- success_count + 1
    }
  }

  message(sprintf("%s: %d files converted successfully", series_name, success_count))
  return(success_count)
}

# ==============================================================================
# Main conversion
# ==============================================================================

message("="
, strrep("=", 59))
message("DICOM to TIFF Conversion")
message(strrep("=", 60))

se1_count <- convert_series("SE1", "SE1_TOP", "SE1_TOP")
se2_count <- convert_series("SE2", "SE2_BTTM", "SE2_BTTM")

message("\n", strrep("=", 60))
message("CONVERSION COMPLETE")
message(strrep("=", 60))
message(sprintf("Output directory: %s", output_path))
message(sprintf("SE1_TOP: %d TIFF files", length(list.files(file.path(output_path, "SE1_TOP"), pattern = "\\.tif$"))))
message(sprintf("SE2_BTTM: %d TIFF files", length(list.files(file.path(output_path, "SE2_BTTM"), pattern = "\\.tif$"))))
