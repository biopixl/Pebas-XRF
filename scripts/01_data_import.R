# ==============================================================================
# Pebas-XRF: Itrax Core Scanner Data Import and Processing
# ==============================================================================
# This script imports and processes XRF data from Itrax core scanner results
# for TAM (Tamshiyacu) and SC (Santa Cruz) sediment cores
# ==============================================================================

# Load required packages
library(tidyverse)
library(fs)

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

# Define paths
base_path <- here::here()
data_path <- file.path(base_path, "TAM-SC-IsaacA")
output_path <- file.path(base_path, "output")

# Define key elements for analysis (detected well with Mo tube)
major_elements <- c("Al", "Si", "K", "Ca", "Ti", "Fe")
minor_elements <- c("P", "S", "Mn", "Rb", "Sr", "Zr", "Ba")
trace_elements <- c("V", "Cr", "Co", "Ni", "Cu", "Zn", "As", "Br", "Pb")

all_elements <- c(major_elements, minor_elements, trace_elements)

# ==============================================================================
# 2. DATA IMPORT FUNCTIONS
# ==============================================================================

#' Import a single Itrax Results.txt file (post-processing calibrated data)
#'
#' @param file_path Path to the Results.txt file
#' @return A tibble with XRF data and metadata
import_itrax_result <- function(file_path) {

  # Read the file
  raw_data <- read_tsv(
    file_path,
    skip = 2,  # Skip header lines
    col_types = cols(.default = col_double(), filename = col_character()),
    show_col_types = FALSE
  )

  # Extract metadata from path
  path_parts <- str_split(file_path, "/")[[1]]

  # Find GROUP and section info
  group_idx <- which(str_detect(path_parts, "^GROUP"))
  group <- if (length(group_idx) > 0) path_parts[group_idx[1]] else NA_character_

  # Get section name (parent directory of Results.txt)
  section <- basename(dirname(file_path))

  # Determine core series (TAM or SC)
  core_series <- case_when(
    str_detect(section, "^TAM") ~ "TAM",
    str_detect(section, "^SC") ~ "SC",
    TRUE ~ NA_character_
  )

  # Clean up data
  result <- raw_data %>%
    rename(
      position_mm = `position (mm)`,
      sample_surface = `sample surface`,
      Mo_inc = `Mo inc`,
      Mo_coh = `Mo coh`
    ) %>%
    mutate(
      group = group,
      section = section,
      core_series = core_series,
      file_path = file_path,
      # Convert position to depth (relative to section start)
      depth_mm = position_mm - min(position_mm, na.rm = TRUE)
    ) %>%
    # Filter out invalid measurements
    filter(validity == 1 | is.na(validity)) %>%
    select(
      # Metadata
      group, core_series, section, position_mm, depth_mm,
      sample_surface, validity, cps, MSE,
      # Major elements
      all_of(intersect(names(.), major_elements)),
      # Minor elements
      all_of(intersect(names(.), minor_elements)),
      # Trace elements
      all_of(intersect(names(.), trace_elements)),
      # Scattering peaks for normalization
      Mo_inc, Mo_coh,
      # Keep other columns
      everything()
    )

  return(result)
}

#' Import all Itrax result files from a directory
#'
#' @param base_dir Base directory containing GROUP folders
#' @param exclude_copied Whether to exclude files in 'copied' directories
#' @return A tibble with all XRF data combined
import_all_itrax <- function(base_dir, exclude_copied = TRUE) {


  # Find all Results.txt files (capital R = post-processing calibrated data)
  result_files <- dir_ls(
    base_dir,
    recurse = TRUE,
    regexp = "Results\\.txt$"
  )

  # Exclude copied directories if requested

if (exclude_copied) {
    result_files <- result_files[!str_detect(result_files, "/copied/")]
  }

  # Exclude stopped runs (incomplete scans)
  result_files <- result_files[!str_detect(result_files, "-stopped/")]

  message(sprintf("Found %d result files to import", length(result_files)))

  # Import all files
  all_data <- map_dfr(result_files, function(f) {
    message(sprintf("  Importing: %s", basename(dirname(f))))
    tryCatch(
      import_itrax_result(f),
      error = function(e) {
        warning(sprintf("Failed to import %s: %s", f, e$message))
        return(NULL)
      }
    )
  })

  return(all_data)
}

# ==============================================================================
# 3. QUALITY CONTROL FUNCTIONS
# ==============================================================================

#' Apply quality control filters to XRF data
#'
#' @param data XRF data tibble
#' @param mse_threshold Maximum MSE (mean squared error) for spectral fit
#' @param cps_min Minimum counts per second
#' @return Filtered data with QC flags
apply_qc_filters <- function(data, mse_threshold = 10, cps_min = 20000) {

  data %>%
    mutate(
      # QC flags
      qc_mse = MSE <= mse_threshold,
      qc_cps = cps >= cps_min,
      qc_surface = sample_surface < 8,  # Not too far from detector
      qc_pass = qc_mse & qc_cps & qc_surface
    ) %>%
    # Report QC statistics
    {
      qc_summary <- summarise(.,
        n_total = n(),
        n_pass_mse = sum(qc_mse, na.rm = TRUE),
        n_pass_cps = sum(qc_cps, na.rm = TRUE),
        n_pass_surface = sum(qc_surface, na.rm = TRUE),
        n_pass_all = sum(qc_pass, na.rm = TRUE)
      )
      message(sprintf("QC Summary: %d/%d measurements pass all filters (%.1f%%)",
                      qc_summary$n_pass_all, qc_summary$n_total,
                      100 * qc_summary$n_pass_all / qc_summary$n_total))
      .
    }
}

# ==============================================================================
# 4. RUN IMPORT
# ==============================================================================

message("Starting Itrax XRF data import...")
message(sprintf("Data path: %s", data_path))

# Import all data
xrf_raw <- import_all_itrax(data_path, exclude_copied = TRUE)

message(sprintf("\nImported %d measurements from %d sections",
                nrow(xrf_raw),
                n_distinct(xrf_raw$section)))

# Apply QC filters
xrf_qc <- apply_qc_filters(xrf_raw, mse_threshold = 10, cps_min = 20000)

# ==============================================================================
# Apply user-defined exclusion zones (foam, gaps, etc.)
# ==============================================================================

exclusion_file <- file.path(base_path, "data", "exclusion_zones.csv")

if (file.exists(exclusion_file)) {
  exclusions <- read_csv(exclusion_file, show_col_types = FALSE) %>%
    filter(!is.na(exclude_start_mm))

  if (nrow(exclusions) > 0) {
    message(sprintf("\nApplying %d exclusion zones...", nrow(exclusions)))

    # Mark excluded points
    xrf_qc <- xrf_qc %>%
      mutate(excluded = FALSE)

    for (i in 1:nrow(exclusions)) {
      sect <- exclusions$section[i]
      start_mm <- exclusions$exclude_start_mm[i]
      end_mm <- exclusions$exclude_end_mm[i]

      n_excluded <- sum(xrf_qc$section == sect &
                        xrf_qc$position_mm >= start_mm &
                        xrf_qc$position_mm <= end_mm)

      xrf_qc <- xrf_qc %>%
        mutate(excluded = excluded |
                 (section == sect & position_mm >= start_mm & position_mm <= end_mm))

      message(sprintf("  %s: %.0f-%.0f mm (%s) - %d points excluded",
                      sect, start_mm, end_mm, exclusions$notes[i], n_excluded))
    }

    # Note: Do NOT filter out excluded zones - mark them for optional masking
    # This keeps data available for visualization with toggle option
    message(sprintf("Total marked as excluded: %d points", sum(xrf_qc$excluded)))
    message(sprintf("QC-pass points (excluding foam/gaps): %d", sum(xrf_qc$qc_pass & !xrf_qc$excluded)))
    message(sprintf("QC-pass points (including foam/gaps): %d", sum(xrf_qc$qc_pass)))
  }
}

# Summary by core series
xrf_qc %>%
  group_by(core_series) %>%
  summarise(
    n_sections = n_distinct(section),
    n_measurements = n(),
    n_qc_pass = sum(qc_pass, na.rm = TRUE),
    depth_range_mm = paste0(round(min(depth_mm)), "-", round(max(depth_mm)))
  ) %>%
  print()

# Save processed data
write_csv(xrf_qc, file.path(output_path, "tables", "xrf_data_qc.csv"))
message(sprintf("\nData saved to: %s", file.path(output_path, "tables", "xrf_data_qc.csv")))

# ==============================================================================
# 5. QUICK DIAGNOSTIC PLOTS
# ==============================================================================

# Element distribution by core series
p_fe_dist <- ggplot(xrf_qc %>% filter(qc_pass), aes(x = Fe, fill = core_series)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_x_log10() +
  labs(title = "Fe Distribution by Core Series",
       x = "Fe (cps)", y = "Count") +
  theme_minimal()

ggsave(file.path(output_path, "figures", "diagnostic_fe_distribution.png"),
       p_fe_dist, width = 8, height = 5)

message("\nImport complete! Next step: run 02_element_ratios.R")
