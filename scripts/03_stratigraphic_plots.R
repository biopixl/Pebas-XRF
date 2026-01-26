# ==============================================================================
# Pebas-XRF: Stratigraphic Plotting Functions
# ==============================================================================
# Multi-panel stratigraphic plots for XRF core data
# TAM = Tamshiyacu, SC = Santa Corina
# ==============================================================================

library(tidyverse)
library(patchwork)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

base_path <- here::here()
output_path <- file.path(base_path, "output")

xrf <- read_csv(file.path(output_path, "tables", "xrf_data_ratios.csv"),
                show_col_types = FALSE)

message(sprintf("Loaded %d measurements", nrow(xrf)))

# ==============================================================================
# 2. PLOTTING THEME
# ==============================================================================

# Custom theme for stratigraphic plots
theme_strat <- function() {
  theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 8),
      strip.text = element_text(face = "bold", size = 9),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "none"
    )
}

# Color palette for elements
element_colors <- c(
  # Major elements
  "Al" = "#4575b4", "Si" = "#74add1", "K" = "#abd9e9",
  "Ca" = "#fee090", "Ti" = "#fdae61", "Fe" = "#d73027",
  # Minor/Trace
  "Mn" = "#f46d43", "Rb" = "#762a83", "Sr" = "#9970ab",
  "Zr" = "#5aae61", "Ba" = "#1b7837"
)

# Colors for ratios
ratio_colors <- c(
  "Ca_Ti" = "#2166ac",     # Carbonate/terrigenous

  "Fe_Mn" = "#b2182b",     # Redox
  "K_Ti" = "#4393c3",      # Weathering
  "Zr_Rb" = "#762a83",     # Grain size
  "Si_Ti" = "#1b7837",     # Biogenic silica
  "Rb_Sr" = "#e08214"      # Weathering intensity
)

# ==============================================================================
# 3. PLOTTING FUNCTIONS
# ==============================================================================

#' Create a single element profile panel
#'
#' @param data XRF data for one section
#' @param element Element column name
#' @param depth_col Depth column name
#' @param color Fill color
#' @param log_scale Use log10 scale for x-axis
#' @return ggplot object
plot_element_profile <- function(data, element, depth_col = "position_mm",
                                  color = "steelblue", log_scale = FALSE) {

  p <- ggplot(data, aes(x = .data[[element]], y = .data[[depth_col]])) +
    geom_path(color = color, linewidth = 0.5) +
    geom_area(alpha = 0.3, fill = color) +
    scale_y_reverse() +
    labs(x = element) +
    theme_strat()

  if (log_scale) {
    p <- p + scale_x_log10()
  }

  return(p)
}

#' Create multi-panel stratigraphic plot for one section
#'
#' @param data XRF data for one section
#' @param elements Vector of element names to plot
#' @param depth_col Depth column name
#' @param title Plot title
#' @param log_scale Use log10 scale
#' @return Combined patchwork plot
plot_section_stratigraphy <- function(data, elements, depth_col = "position_mm",
                                       title = NULL, log_scale = TRUE) {

  plots <- map(elements, function(elem) {
    col <- if (elem %in% names(element_colors)) element_colors[[elem]] else "steelblue"
    plot_element_profile(data, elem, depth_col, color = col, log_scale = log_scale)
  })

  combined <- wrap_plots(plots, nrow = 1) +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(face = "bold", size = 12))
    )

  return(combined)
}

#' Create multi-panel ratio plot for one section
#'
#' @param data XRF data with ratios
#' @param ratios Vector of ratio column names
#' @param depth_col Depth column name
#' @param title Plot title
#' @return Combined patchwork plot
plot_section_ratios <- function(data, ratios, depth_col = "position_mm",
                                 title = NULL) {

  plots <- map(ratios, function(ratio) {
    # Use log version if available
    log_ratio <- paste0("log_", ratio)
    use_col <- if (log_ratio %in% names(data)) log_ratio else ratio

    col <- if (ratio %in% names(ratio_colors)) ratio_colors[[ratio]] else "darkorange"

    ggplot(data, aes(x = .data[[use_col]], y = .data[[depth_col]])) +
      geom_path(color = col, linewidth = 0.6) +
      scale_y_reverse() +
      labs(x = gsub("_", "/", ratio)) +
      theme_strat()
  })

  combined <- wrap_plots(plots, nrow = 1) +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(face = "bold", size = 12))
    )

  return(combined)
}

#' Create complete stratigraphic figure for a section
#'
#' @param data XRF data with ratios
#' @param section_name Section to plot
#' @param elements Elements to include
#' @param ratios Ratios to include
#' @return Combined figure
plot_complete_section <- function(data, section_name,
                                   elements = c("Ti", "K", "Ca", "Fe", "Mn"),
                                   ratios = c("Ca_Ti", "Fe_Mn", "K_Ti", "Zr_Rb")) {

  section_data <- data %>% filter(section == section_name)

  if (nrow(section_data) == 0) {
    warning(sprintf("No data for section: %s", section_name))
    return(NULL)
  }

  p_elements <- plot_section_stratigraphy(
    section_data, elements,
    title = paste0(section_name, " - Major Elements (log cps)")
  )

  p_ratios <- plot_section_ratios(
    section_data, ratios,
    title = paste0(section_name, " - Element Ratios (log)")
  )

  combined <- p_elements / p_ratios +
    plot_layout(heights = c(1, 1))

  return(combined)
}

# ==============================================================================
# 4. GENERATE PLOTS FOR ALL SECTIONS
# ==============================================================================

message("Generating stratigraphic plots...")

# Get unique sections
sections <- unique(xrf$section)
message(sprintf("Found %d sections to plot", length(sections)))

# Key elements and ratios
elements_to_plot <- c("Ti", "K", "Ca", "Fe", "Mn", "Sr")
ratios_to_plot <- c("Ca_Ti", "Fe_Mn", "K_Ti", "Zr_Rb")

# Generate plots for each section
for (sect in sections) {

  p <- plot_complete_section(xrf, sect, elements_to_plot, ratios_to_plot)

  if (!is.null(p)) {
    # Clean filename
    filename <- paste0("strat_", gsub("[^A-Za-z0-9_-]", "_", sect), ".png")
    filepath <- file.path(output_path, "figures", filename)

    ggsave(filepath, p, width = 12, height = 8, dpi = 150, bg = "white")
    message(sprintf("  Saved: %s", filename))
  }
}

# ==============================================================================
# 5. COMPOSITE PLOTS BY CORE SERIES
# ==============================================================================

#' Create composite plot comparing all sections in a core series
#'
#' @param data XRF data
#' @param series Core series name (TAM or SC)
#' @param element Element to plot
#' @return ggplot faceted by section
plot_series_composite <- function(data, series, element = "Fe") {

  series_data <- data %>%
    filter(core_series == series) %>%
    arrange(section, position_mm)

  ggplot(series_data, aes(x = .data[[element]], y = position_mm)) +
    geom_path(color = element_colors[[element]], linewidth = 0.4) +
    geom_area(alpha = 0.2, fill = element_colors[[element]]) +
    scale_y_reverse() +
    scale_x_log10() +
    facet_wrap(~section, scales = "free", ncol = 4) +
    labs(
      title = sprintf("%s Cores - %s Profiles",
                      ifelse(series == "TAM", "Tamshiyacu", "Santa Corina"),
                      element),
      x = paste0(element, " (cps)"),
      y = "Position (mm)"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      panel.grid.minor = element_blank()
    )
}

# TAM composite
p_tam <- plot_series_composite(xrf, "TAM", "Fe")
ggsave(file.path(output_path, "figures", "composite_TAM_Fe.png"),
       p_tam, width = 14, height = 10, dpi = 150, bg = "white")

# SC composite
p_sc <- plot_series_composite(xrf, "SC", "Fe")
ggsave(file.path(output_path, "figures", "composite_SC_Fe.png"),
       p_sc, width = 14, height = 10, dpi = 150, bg = "white")

message("\nComposite plots saved!")

# ==============================================================================
# 6. RATIO COMPARISON PLOT
# ==============================================================================

# Compare key ratios between TAM and SC
p_ratio_compare <- xrf %>%
  select(core_series, section, position_mm, Ca_Ti, Fe_Mn, K_Ti, Zr_Rb) %>%
  pivot_longer(cols = c(Ca_Ti, Fe_Mn, K_Ti, Zr_Rb),
               names_to = "ratio", values_to = "value") %>%
  ggplot(aes(x = log10(value), fill = core_series)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ratio, scales = "free", ncol = 2) +
  scale_fill_manual(
    values = c("TAM" = "#2166ac", "SC" = "#b2182b"),
    labels = c("TAM" = "Tamshiyacu", "SC" = "Santa Corina")
  ) +
  labs(
    title = "Element Ratio Distributions by Core Series",
    x = "log10(Ratio)", y = "Density",
    fill = "Core Series"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(output_path, "figures", "ratio_distributions.png"),
       p_ratio_compare, width = 10, height = 8, dpi = 150, bg = "white")

message("\nStratigraphic plotting complete! Next step: run 04_pca_analysis.R")
