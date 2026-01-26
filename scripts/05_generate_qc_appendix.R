# ==============================================================================
# Pebas-XRF: Generate QC Appendix Figures
# ==============================================================================
# Creates publication-ready figures showing excluded sections for manuscript
# ==============================================================================

library(tidyverse)
library(patchwork)
library(png)
library(grid)
library(gridExtra)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

base_path <- here::here()
output_path <- file.path(base_path, "output")
optical_path <- file.path(output_path, "figures", "optical_review")
appendix_path <- file.path(output_path, "figures", "appendix")
dir.create(appendix_path, showWarnings = FALSE, recursive = TRUE)

exclusion_file <- file.path(base_path, "data", "exclusion_zones.csv")

# ==============================================================================
# LOAD DATA
# ==============================================================================

# Load XRF data
xrf_data <- read_csv(file.path(output_path, "tables", "xrf_data_qc.csv"),
                     show_col_types = FALSE)

# Load exclusions
exclusions <- read_csv(exclusion_file, show_col_types = FALSE) %>%
  filter(!is.na(exclude_start_mm))

message(sprintf("Loaded %d exclusion zones", nrow(exclusions)))

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Get optical image info from document.txt
get_optical_info <- function(section_name, data_path) {
  # Find the section folder
  section_dirs <- list.dirs(data_path, recursive = TRUE)
  section_dir <- section_dirs[basename(section_dirs) == section_name &
                               !grepl("/copied/", section_dirs)]

  if (length(section_dir) == 0) return(NULL)
  section_dir <- section_dir[1]

  doc_file <- file.path(section_dir, "document.txt")
  if (!file.exists(doc_file)) return(NULL)

  doc <- readLines(doc_file, warn = FALSE, encoding = "latin1")

  # Parse XRF scan coordinates
  xrf_line <- grep("^Start coordinate", doc, value = TRUE)
  if (length(xrf_line) > 0) {
    parts <- strsplit(xrf_line, "\t")[[1]]
    xrf_start <- as.numeric(parts[2])
    xrf_stop <- as.numeric(parts[4])
  } else {
    return(NULL)
  }

  # Parse optical image coordinates
  optical_line <- grep("^Optical Start", doc, value = TRUE)
  if (length(optical_line) > 0) {
    parts <- strsplit(optical_line, "\t")[[1]]
    optical_start <- as.numeric(parts[2])
    optical_end <- as.numeric(parts[4])
  } else {
    optical_start <- 0
    optical_end <- 1500
  }

  list(
    optical_start = optical_start,
    optical_end = optical_end,
    xrf_start = xrf_start,
    xrf_stop = xrf_stop
  )
}

# ==============================================================================
# GENERATE INDIVIDUAL SECTION FIGURES
# ==============================================================================

#' Create QC figure for a single section
#'
#' Shows optical image aligned with XRF profile on shared x-axis
create_section_qc_figure <- function(section_name, xrf_data, exclusions,
                                      optical_path, data_path, element = "Fe") {

  # Get section data
  sect_data <- xrf_data %>% filter(section == section_name)

  if (nrow(sect_data) == 0) {
    message(sprintf("  No data for %s, skipping", section_name))
    return(NULL)
  }

  # Get exclusions for this section
  sect_excl <- exclusions %>% filter(section == section_name)

  # Get optical image info for alignment
  info <- get_optical_info(section_name, data_path)

  # Position range from XRF data
  pos_min <- min(sect_data$position_mm, na.rm = TRUE)
  pos_max <- max(sect_data$position_mm, na.rm = TRUE)

  # Try to load optical image
  optical_file <- file.path(optical_path, paste0("optical_", section_name, ".png"))

  # Create the figure using base R for proper alignment
  fig_file <- tempfile(fileext = ".png")

  png(fig_file, width = 14, height = 8, units = "in", res = 150)

  # Set up layout: optical on top, XRF below
  par(mfrow = c(2, 1), mar = c(1, 5, 3, 2), oma = c(3, 0, 2, 0))

  # X-axis range with small padding
  x_range <- pos_max - pos_min
  x_min <- pos_min - x_range * 0.02
  x_max <- pos_max + x_range * 0.02

  # === PANEL 1: Optical Image ===
  if (file.exists(optical_file) && !is.null(info)) {
    img <- readPNG(optical_file)
    img_width <- ncol(img)
    img_height <- nrow(img)

    # Calculate pixel positions for XRF range
    mm_per_pixel <- (info$optical_end - info$optical_start) / img_width

    px_start <- (x_min - info$optical_start) / mm_per_pixel
    px_end <- (x_max - info$optical_start) / mm_per_pixel

    # Crop image to XRF range (with bounds checking)
    px_start <- max(1, floor(px_start))
    px_end <- min(img_width, ceiling(px_end))

    if (px_end > px_start) {
      img_crop <- img[, px_start:px_end, , drop = FALSE]

      # Plot cropped image aligned with XRF x-axis
      plot(0, 0, type = "n", xlim = c(x_min, x_max), ylim = c(0, 1),
           xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           main = "")

      rasterImage(img_crop, x_min, 0, x_max, 1)

      # Add exclusion zones as semi-transparent overlays
      if (nrow(sect_excl) > 0) {
        for (i in 1:nrow(sect_excl)) {
          rect(sect_excl$exclude_start_mm[i], 0,
               sect_excl$exclude_end_mm[i], 1,
               col = rgb(1, 0, 0, 0.4), border = "red", lwd = 3)
        }
      }

      # Add position scale ticks
      axis(1, at = pretty(c(x_min, x_max)), labels = FALSE, tck = -0.03)
      mtext("Optical Image (aligned to XRF scan range)", side = 3, line = 0.5, cex = 1)
    }
  } else {
    plot.new()
    text(0.5, 0.5, "Optical image not available", cex = 1.5)
  }

  # === PANEL 2: XRF Profile ===
  y_vals <- sect_data[[element]]
  y_vals_log <- log10(pmax(y_vals, 1))

  plot(sect_data$position_mm, y_vals_log, type = "l", col = "steelblue", lwd = 2,
       xlim = c(x_min, x_max),
       xlab = "", ylab = paste0("log10(", element, " cps)"),
       main = paste(element, "Profile"), cex.main = 1.2, cex.lab = 1.1)
  points(sect_data$position_mm, y_vals_log, pch = 16, cex = 0.6, col = "steelblue")

  # Add exclusion zones
  if (nrow(sect_excl) > 0) {
    for (i in 1:nrow(sect_excl)) {
      rect(sect_excl$exclude_start_mm[i], par("usr")[3],
           sect_excl$exclude_end_mm[i], par("usr")[4],
           col = rgb(1, 0, 0, 0.3), border = NA)
      abline(v = c(sect_excl$exclude_start_mm[i], sect_excl$exclude_end_mm[i]),
             col = "red", lty = 2, lwd = 2)

      # Label the exclusion
      mid_x <- mean(c(sect_excl$exclude_start_mm[i], sect_excl$exclude_end_mm[i]))
      text(mid_x, par("usr")[4] - 0.1 * diff(par("usr")[3:4]),
           sect_excl$notes[i], col = "red", cex = 1, font = 2)
    }
  }

  grid(nx = NA, ny = NULL, col = "gray80", lty = 1)

  # Shared x-axis label
  mtext("Position (mm)", side = 1, outer = TRUE, line = 1.5, cex = 1.2)

  # Main title
  n_excl <- nrow(sect_excl)
  excl_text <- ifelse(n_excl > 0, sprintf(" | %d exclusion zone(s)", n_excl), "")
  mtext(sprintf("%s | Position: %.0f - %.0f mm | N = %d%s",
                section_name, pos_min, pos_max, nrow(sect_data), excl_text),
        side = 3, outer = TRUE, line = 0.5, cex = 1.3, font = 2)

  dev.off()

  return(fig_file)
}

# ==============================================================================
# GENERATE FIGURES ONLY FOR SECTIONS WITH EXCLUSIONS
# ==============================================================================

message("\nGenerating appendix figures for excluded sections only...")

# Data path for reading document.txt files
data_path <- file.path(base_path, "TAM-SC-IsaacA")

# Only process sections that have exclusions
sections_with_exclusions <- exclusions %>%
  filter(!is.na(exclude_start_mm)) %>%
  pull(section) %>%
  unique() %>%
  sort()

message(sprintf("Generating figures for %d sections with exclusions", length(sections_with_exclusions)))

for (sect in sections_with_exclusions) {
  message(sprintf("Processing: %s", sect))

  fig_temp <- create_section_qc_figure(sect, xrf_data, exclusions, optical_path, data_path)

  if (!is.null(fig_temp)) {
    filename <- paste0("appendix_qc_", gsub("[^A-Za-z0-9_-]", "_", sect), ".png")
    file.copy(fig_temp, file.path(appendix_path, filename), overwrite = TRUE)
    unlink(fig_temp)
  }
}

# ==============================================================================
# CREATE COMBINED SUMMARY FIGURE
# ==============================================================================

message("\nCreating QC exclusion summary figure...")

# Calculate exclusion statistics
exclusion_stats <- exclusions %>%
  filter(!is.na(exclude_start_mm)) %>%
  mutate(excluded_mm = exclude_end_mm - exclude_start_mm) %>%
  group_by(section) %>%
  summarise(
    n_zones = n(),
    total_excluded_mm = sum(excluded_mm),
    .groups = "drop"
  ) %>%
  left_join(
    xrf_data %>%
      group_by(section) %>%
      summarise(
        total_mm = max(position_mm) - min(position_mm),
        n_points = n(),
        .groups = "drop"
      ),
    by = "section"
  ) %>%
  mutate(
    pct_excluded = 100 * total_excluded_mm / total_mm,
    core_series = ifelse(grepl("^TAM", section), "TAM", "SC")
  ) %>%
  arrange(core_series, section)

# Create summary bar plot
png(file.path(appendix_path, "qc_exclusion_summary.png"),
    width = 12, height = 8, units = "in", res = 150)

par(mar = c(8, 5, 3, 2))

barplot(
  exclusion_stats$pct_excluded,
  names.arg = exclusion_stats$section,
  las = 2,
  col = ifelse(exclusion_stats$core_series == "TAM", "steelblue", "darkorange"),
  ylab = "% of scan length excluded",
  main = "QC Exclusions by Section (foam fills removed)",
  cex.names = 0.7,
  ylim = c(0, max(exclusion_stats$pct_excluded) * 1.1)
)

# Add legend
legend("topright", legend = c("TAM (Tamshiyacu)", "SC (Santa Corina)"),
       fill = c("steelblue", "darkorange"), bty = "n")

# Add count labels
text(
  x = seq(0.7, by = 1.2, length.out = nrow(exclusion_stats)),
  y = exclusion_stats$pct_excluded + 1,
  labels = paste0(exclusion_stats$n_zones, " zone(s)"),
  cex = 0.6, col = "gray30"
)

dev.off()

message("Summary figure saved: qc_exclusion_summary.png")

# ==============================================================================
# CREATE SUMMARY TABLE FOR MANUSCRIPT
# ==============================================================================

# Summarize exclusions by section
exclusion_summary <- xrf_data %>%
  group_by(section, core_series) %>%
  summarise(
    n_total = n(),
    pos_start = min(position_mm, na.rm = TRUE),
    pos_end = max(position_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    exclusions %>%
      group_by(section) %>%
      summarise(
        n_exclusions = n(),
        excluded_ranges = paste(
          sprintf("%.0f-%.0f", exclude_start_mm, exclude_end_mm),
          collapse = "; "
        ),
        exclusion_notes = paste(notes, collapse = "; "),
        .groups = "drop"
      ),
    by = "section"
  ) %>%
  mutate(
    n_exclusions = replace_na(n_exclusions, 0),
    excluded_ranges = replace_na(excluded_ranges, "None"),
    exclusion_notes = replace_na(exclusion_notes, "")
  )

# Save summary table
write_csv(exclusion_summary, file.path(output_path, "tables", "exclusion_summary.csv"))

# Create LaTeX table for manuscript
latex_table <- exclusion_summary %>%
  select(
    Section = section,
    Series = core_series,
    `N Total` = n_total,
    `Position (mm)` = pos_start,
    `Exclusions` = n_exclusions,
    `Excluded Ranges` = excluded_ranges
  ) %>%
  mutate(
    `Position (mm)` = sprintf("%.0f-%.0f",
                               exclusion_summary$pos_start,
                               exclusion_summary$pos_end)
  )

# Write LaTeX table
sink(file.path(appendix_path, "exclusion_table.tex"))
cat("\\begin{longtable}{llrlrl}\n")
cat("\\caption{Quality control exclusion zones for XRF core scanner data.}\n")
cat("\\label{tab:exclusions}\\\\\n")
cat("\\toprule\n")
cat("Section & Series & N & Position (mm) & Excl. & Excluded Ranges \\\\\n")
cat("\\midrule\n")
cat("\\endfirsthead\n")
cat("\\multicolumn{6}{c}{\\tablename\\ \\thetable\\ -- continued}\\\\\n")
cat("\\toprule\n")
cat("Section & Series & N & Position (mm) & Excl. & Excluded Ranges \\\\\n")
cat("\\midrule\n")
cat("\\endhead\n")
cat("\\bottomrule\n")
cat("\\endfoot\n")

for (i in 1:nrow(latex_table)) {
  row <- latex_table[i, ]
  cat(sprintf("%s & %s & %d & %s & %d & %s \\\\\n",
              gsub("_", "\\\\_", row$Section),
              row$Series,
              row$`N Total`,
              row$`Position (mm)`,
              row$Exclusions,
              gsub(";", ",", row$`Excluded Ranges`)))
}

cat("\\end{longtable}\n")
sink()

message(sprintf("\nAppendix figures saved to: %s", appendix_path))
message(sprintf("Exclusion summary saved to: %s", file.path(output_path, "tables", "exclusion_summary.csv")))
message(sprintf("LaTeX table saved to: %s", file.path(appendix_path, "exclusion_table.tex")))

# ==============================================================================
# SUMMARY
# ==============================================================================

message("\n=== EXCLUSION SUMMARY ===")
message(sprintf("Total sections: %d", n_distinct(xrf_data$section)))
message(sprintf("Sections with exclusions: %d", length(sections_with_exclusions)))
message(sprintf("Total exclusion zones: %d", sum(exclusion_summary$n_exclusions)))

# Print sections with exclusions
if (any(exclusion_summary$n_exclusions > 0)) {
  message("\nSections with exclusion zones:")
  exclusion_summary %>%
    filter(n_exclusions > 0) %>%
    select(section, n_exclusions, excluded_ranges) %>%
    print(n = 50)
}
