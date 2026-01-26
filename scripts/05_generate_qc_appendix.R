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
# GENERATE INDIVIDUAL SECTION FIGURES
# ==============================================================================

#' Create QC figure for a single section
#'
#' Shows optical image with XRF profile and exclusion zones marked
create_section_qc_figure <- function(section_name, xrf_data, exclusions,
                                      optical_path, element = "Fe") {

  # Get section data
  sect_data <- xrf_data %>% filter(section == section_name)

  if (nrow(sect_data) == 0) {
    message(sprintf("  No data for %s, skipping", section_name))
    return(NULL)
  }

  # Get exclusions for this section
  sect_excl <- exclusions %>% filter(section == section_name)

  # Position range
  pos_min <- min(sect_data$position_mm, na.rm = TRUE)
  pos_max <- max(sect_data$position_mm, na.rm = TRUE)

  # Create XRF profile plot
  p_xrf <- ggplot(sect_data, aes(x = position_mm, y = .data[[element]])) +
    geom_line(color = "steelblue", linewidth = 0.6) +
    scale_y_log10() +
    labs(x = "Position (mm)", y = paste0(element, " (cps)")) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 9)
    )

  # Add exclusion zones
  if (nrow(sect_excl) > 0) {
    for (i in 1:nrow(sect_excl)) {
      p_xrf <- p_xrf +
        annotate("rect",
                 xmin = sect_excl$exclude_start_mm[i],
                 xmax = sect_excl$exclude_end_mm[i],
                 ymin = -Inf, ymax = Inf,
                 fill = "red", alpha = 0.3) +
        annotate("segment",
                 x = sect_excl$exclude_start_mm[i],
                 xend = sect_excl$exclude_start_mm[i],
                 y = -Inf, yend = Inf,
                 color = "red", linetype = "dashed", linewidth = 0.5) +
        annotate("segment",
                 x = sect_excl$exclude_end_mm[i],
                 xend = sect_excl$exclude_end_mm[i],
                 y = -Inf, yend = Inf,
                 color = "red", linetype = "dashed", linewidth = 0.5)
    }
  }

  # Try to load optical image
  optical_file <- file.path(optical_path, paste0("optical_", section_name, ".png"))

  if (file.exists(optical_file)) {
    img <- readPNG(optical_file)

    # Create optical plot using annotation_custom
    p_optical <- ggplot() +
      annotation_custom(
        rasterGrob(img, interpolate = TRUE),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      ) +
      theme_void() +
      coord_fixed(ratio = nrow(img) / ncol(img) * 10)

    # Combine plots
    combined <- p_optical / p_xrf +
      plot_layout(heights = c(1, 2)) +
      plot_annotation(
        title = section_name,
        subtitle = sprintf("Position: %.0f - %.0f mm | %d measurements | %d exclusion zones",
                           pos_min, pos_max, nrow(sect_data), nrow(sect_excl)),
        theme = theme(
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "gray40")
        )
      )
  } else {
    # No optical image, just XRF plot
    combined <- p_xrf +
      ggtitle(
        section_name,
        subtitle = sprintf("Position: %.0f - %.0f mm | %d measurements | %d exclusion zones",
                           pos_min, pos_max, nrow(sect_data), nrow(sect_excl))
      )
  }

  return(combined)
}

# ==============================================================================
# GENERATE ALL FIGURES
# ==============================================================================

message("\nGenerating appendix figures...")

sections <- sort(unique(xrf_data$section))

for (sect in sections) {
  message(sprintf("Processing: %s", sect))

  fig <- create_section_qc_figure(sect, xrf_data, exclusions, optical_path)

  if (!is.null(fig)) {
    filename <- paste0("appendix_qc_", gsub("[^A-Za-z0-9_-]", "_", sect), ".png")
    ggsave(
      file.path(appendix_path, filename),
      fig,
      width = 10, height = 6, dpi = 150, bg = "white"
    )
  }
}

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
message(sprintf("Total sections: %d", length(sections)))
message(sprintf("Sections with exclusions: %d", sum(exclusion_summary$n_exclusions > 0)))
message(sprintf("Total exclusion zones: %d", sum(exclusion_summary$n_exclusions)))

# Print sections with exclusions
if (any(exclusion_summary$n_exclusions > 0)) {
  message("\nSections with exclusion zones:")
  exclusion_summary %>%
    filter(n_exclusions > 0) %>%
    select(section, n_exclusions, excluded_ranges) %>%
    print(n = 50)
}
