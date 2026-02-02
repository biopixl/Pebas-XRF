# ==============================================================================
# Pebas-XRF: CT Data Integration
# ==============================================================================
# Reads DICOM CT data and prepares for integration with XRF core scanner data
# CT scans provide density information that complements XRF elemental data
# ==============================================================================

library(tidyverse)
library(oro.dicom)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

base_path <- here::here()
ct_path <- file.path(base_path, "Pebas-CT", "ST1")
output_path <- file.path(base_path, "output")
fig_path <- file.path(output_path, "figures", "ct_integration")

dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. READ DICOM HEADERS (Metadata only - fast)
# ==============================================================================

message("=== Reading CT DICOM Metadata ===\n")

# Function to read DICOM header only
read_dicom_header <- function(file_path) {
  tryCatch({
    dcm <- readDICOMFile(file_path, pixelData = FALSE)

    # Extract key metadata
    header_list <- list(
      file = basename(file_path),
      series = basename(dirname(file_path)),
      PatientID = dcm$hdr$value[dcm$hdr$name == "PatientID"],
      StudyID = dcm$hdr$value[dcm$hdr$name == "StudyID"],
      SeriesDescription = dcm$hdr$value[dcm$hdr$name == "SeriesDescription"],
      InstanceNumber = as.numeric(dcm$hdr$value[dcm$hdr$name == "InstanceNumber"]),
      SliceLocation = as.numeric(dcm$hdr$value[dcm$hdr$name == "SliceLocation"]),
      SliceThickness = as.numeric(dcm$hdr$value[dcm$hdr$name == "SliceThickness"]),
      ImagePositionPatient = dcm$hdr$value[dcm$hdr$name == "ImagePositionPatient"],
      PixelSpacing = dcm$hdr$value[dcm$hdr$name == "PixelSpacing"],
      Rows = as.numeric(dcm$hdr$value[dcm$hdr$name == "Rows"]),
      Columns = as.numeric(dcm$hdr$value[dcm$hdr$name == "Columns"])
    )

    return(as_tibble(header_list))
  }, error = function(e) {
    message(sprintf("  Error reading %s: %s", basename(file_path), e$message))
    return(NULL)
  })
}

# Read all DICOM headers
se1_files <- list.files(file.path(ct_path, "SE1"), pattern = "^IM", full.names = TRUE)
se2_files <- list.files(file.path(ct_path, "SE2"), pattern = "^IM", full.names = TRUE)

message(sprintf("Found %d SE1 files, %d SE2 files", length(se1_files), length(se2_files)))

# Read headers
message("\nReading SE1 (TOP) headers...")
se1_headers <- map_dfr(se1_files, read_dicom_header) %>%
  arrange(InstanceNumber)

message("\nReading SE2 (BOTTOM) headers...")
se2_headers <- map_dfr(se2_files, read_dicom_header) %>%
  arrange(InstanceNumber)

# Combine
ct_metadata <- bind_rows(
  se1_headers %>% mutate(series_name = "SE1_TOP"),
  se2_headers %>% mutate(series_name = "SE2_BTTM")
)

message(sprintf("\nTotal slices: %d", nrow(ct_metadata)))

# ==============================================================================
# 2. SUMMARIZE CT METADATA
# ==============================================================================

message("\n=== CT Data Summary ===\n")

ct_summary <- ct_metadata %>%
  group_by(series_name, SeriesDescription) %>%
  summarise(
    n_slices = n(),
    min_depth = min(SliceLocation, na.rm = TRUE),
    max_depth = max(SliceLocation, na.rm = TRUE),
    depth_range = max_depth - min_depth,
    slice_thickness = first(SliceThickness),
    pixel_spacing = first(PixelSpacing),
    .groups = "drop"
  )

print(ct_summary)

# Save metadata
write_csv(ct_metadata, file.path(output_path, "tables", "ct_dicom_metadata.csv"))
write_csv(ct_summary, file.path(output_path, "tables", "ct_series_summary.csv"))

# ==============================================================================
# 3. CREATE DEPTH PROFILE TABLE FOR XRF ALIGNMENT
# ==============================================================================

message("\n=== Creating Depth Profile for XRF Alignment ===\n")

ct_depth_profile <- ct_metadata %>%
  select(series_name, InstanceNumber, SliceLocation, SliceThickness) %>%
  arrange(series_name, SliceLocation) %>%
  mutate(
    # Create depth bins that can match XRF 3mm resolution
    depth_bin_3mm = round(SliceLocation / 3) * 3
  )

message("CT depth profile created with following ranges:")
message(sprintf("  SE1 (TOP): %.1f to %.1f mm (%d slices)",
                min(ct_depth_profile$SliceLocation[ct_depth_profile$series_name == "SE1_TOP"]),
                max(ct_depth_profile$SliceLocation[ct_depth_profile$series_name == "SE1_TOP"]),
                sum(ct_depth_profile$series_name == "SE1_TOP")))
message(sprintf("  SE2 (BTTM): %.1f to %.1f mm (%d slices)",
                min(ct_depth_profile$SliceLocation[ct_depth_profile$series_name == "SE2_BTTM"]),
                max(ct_depth_profile$SliceLocation[ct_depth_profile$series_name == "SE2_BTTM"]),
                sum(ct_depth_profile$series_name == "SE2_BTTM")))

# ==============================================================================
# 4. READ PIXEL DATA FOR REPRESENTATIVE SLICES
# ==============================================================================

message("\n=== Reading Representative CT Slices ===\n")

# Read full DICOM data for a few slices to demonstrate density variation
read_dicom_with_pixels <- function(file_path) {
  tryCatch({
    dcm <- readDICOMFile(file_path, pixelData = TRUE)

    # Extract image and convert to Hounsfield Units if possible
    img <- dcm$img

    # Get window center/width for display
    wc <- as.numeric(dcm$hdr$value[dcm$hdr$name == "WindowCenter"][1])
    ww <- as.numeric(dcm$hdr$value[dcm$hdr$name == "WindowWidth"][1])

    list(
      image = img,
      slice_location = as.numeric(dcm$hdr$value[dcm$hdr$name == "SliceLocation"]),
      window_center = wc,
      window_width = ww
    )
  }, error = function(e) {
    message(sprintf("  Error reading pixels from %s: %s", basename(file_path), e$message))
    return(NULL)
  })
}

# Read middle slices from each series
mid_se1 <- se1_files[length(se1_files) %/% 2]
mid_se2 <- se2_files[length(se2_files) %/% 2]

message(sprintf("Reading middle slice from SE1: %s", basename(mid_se1)))
se1_sample <- read_dicom_with_pixels(mid_se1)

message(sprintf("Reading middle slice from SE2: %s", basename(mid_se2)))
se2_sample <- read_dicom_with_pixels(mid_se2)

# ==============================================================================
# 5. EXTRACT DENSITY PROFILES
# ==============================================================================

message("\n=== Extracting CT Density Profiles ===\n")

# Extract central column profile (approximate core axis)
extract_central_profile <- function(dcm_data, label) {
  if (is.null(dcm_data) || is.null(dcm_data$image)) return(NULL)

  img <- dcm_data$image
  mid_col <- ncol(img) %/% 2

  # Average across 10 central columns for noise reduction
  cols <- (mid_col - 5):(mid_col + 4)
  profile <- rowMeans(img[, cols])

  tibble(
    row = 1:length(profile),
    density = profile,
    series = label,
    slice_location = dcm_data$slice_location
  )
}

se1_profile <- extract_central_profile(se1_sample, "SE1_TOP")
se2_profile <- extract_central_profile(se2_sample, "SE2_BTTM")

ct_profiles <- bind_rows(se1_profile, se2_profile)

# ==============================================================================
# 6. CREATE VISUALIZATION
# ==============================================================================

message("\n=== Creating CT Visualization ===\n")

# CT slice images
if (!is.null(se1_sample$image)) {
  png(file.path(fig_path, "ct_slice_SE1_sample.png"), width = 800, height = 800)
  image(se1_sample$image, col = gray.colors(256), axes = FALSE,
        main = sprintf("CT Slice SE1 (TOP) at %.1f mm", se1_sample$slice_location))
  dev.off()
  message("  Saved: ct_slice_SE1_sample.png")
}

if (!is.null(se2_sample$image)) {
  png(file.path(fig_path, "ct_slice_SE2_sample.png"), width = 800, height = 800)
  image(se2_sample$image, col = gray.colors(256), axes = FALSE,
        main = sprintf("CT Slice SE2 (BTTM) at %.1f mm", se2_sample$slice_location))
  dev.off()
  message("  Saved: ct_slice_SE2_sample.png")
}

# Density profile plot
if (nrow(ct_profiles) > 0) {
  p_profile <- ggplot(ct_profiles, aes(x = row, y = density, color = series)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~series, ncol = 1, scales = "free_y") +
    labs(
      title = "CT Density Profiles (Central Axis)",
      subtitle = sprintf("SE1 at %.1f mm | SE2 at %.1f mm",
                        se1_sample$slice_location, se2_sample$slice_location),
      x = "Row (pixels)", y = "CT Intensity",
      color = "Series"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(file.path(fig_path, "ct_density_profiles.png"), p_profile,
         width = 10, height = 8, dpi = 300, bg = "white")
  message("  Saved: ct_density_profiles.png")
}

# ==============================================================================
# 7. DEPTH COMPARISON WITH XRF
# ==============================================================================

message("\n=== Comparing CT and XRF Depth Coverage ===\n")

# Load XRF data for comparison
xrf_data <- read_csv(file.path(output_path, "tables", "xrf_data_stacked.csv"),
                     show_col_types = FALSE)

# XRF depth summary by section
xrf_depths <- xrf_data %>%
  group_by(section) %>%
  summarise(
    xrf_min = min(position_mm),
    xrf_max = max(position_mm),
    xrf_range = xrf_max - xrf_min,
    n_measurements = n(),
    .groups = "drop"
  )

# Find potential matching sections (based on depth range similarity)
ct_se1_range <- diff(range(ct_depth_profile$SliceLocation[ct_depth_profile$series_name == "SE1_TOP"]))
ct_se2_range <- diff(range(ct_depth_profile$SliceLocation[ct_depth_profile$series_name == "SE2_BTTM"]))

message("CT depth ranges:")
message(sprintf("  SE1: %.1f mm", ct_se1_range))
message(sprintf("  SE2: %.1f mm", ct_se2_range))

message("\nXRF sections with similar depth ranges (50-120 mm):")
similar_sections <- xrf_depths %>%
  filter(xrf_range >= 50 & xrf_range <= 120) %>%
  arrange(xrf_range)
print(similar_sections)

# ==============================================================================
# 8. SUMMARY
# ==============================================================================

message("\n=== CT Integration Summary ===")
message(sprintf("CT data location: %s", ct_path))
message(sprintf("Total CT slices: %d", nrow(ct_metadata)))
message(sprintf("CT series: SE1 (TOP, %d slices), SE2 (BTTM, %d slices)",
                sum(ct_metadata$series_name == "SE1_TOP"),
                sum(ct_metadata$series_name == "SE2_BTTM")))
message(sprintf("\nOutput saved to: %s", fig_path))
message("\nNote: CT-XRF alignment requires additional metadata to identify which")
message("core section was scanned. StudyID 22166 should be matched to core records.")
message("\nNext steps:")
message("  1. Identify CT sample from core records using StudyID/PatientID")
message("  2. Create depth-matched CT density profile")
message("  3. Correlate CT density with XRF Ca signal (carbonate content)")
