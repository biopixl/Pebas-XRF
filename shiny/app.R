#!/usr/bin/env Rscript
# ==============================================================================
# Pebas-XRF: Interactive Stratigraphic Figure Builder
# ==============================================================================
# Shiny app for building custom stratigraphic figures with core images
# and selectable proxies/elements
# ==============================================================================

library(shiny)
library(tidyverse)
library(magick)
library(zoo)
library(patchwork)
library(shinyWidgets)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Paths - adjust if running from different location
base_path <- here::here()
data_path <- file.path(base_path, "TAM-SC-IsaacA")
output_path <- file.path(base_path, "output")

# Image enhancement settings
BRIGHTNESS_BOOST <- 80
CONTRAST_BOOST <- 25
GAMMA_CORRECTION <- 2.0
BLACK_POINT <- 2
WHITE_POINT <- 60
WINDOW_SIZE <- 5

# ==============================================================================
# SECTION CONFIGURATION (Stratigraphic Order)
# ==============================================================================

section_config <- tribble(
  ~section, ~optical_path, ~group, ~site, ~strat_order,
  # GROUP1 - TAM (oldest)

  "TAM-1-2-3B-A", "GROUP1/TAM-1-2-3B/TAM-1-2-3B-A", "GROUP1", "TAM", 1,
  "TAM-1-2-3B-B", "GROUP1/TAM-1-2-3B/TAM-1-2-3B-B", "GROUP1", "TAM", 2,
  "TAM-1-2-3B-C", "GROUP1/TAM-1-2-3B/TAM-1-2-3B-C", "GROUP1", "TAM", 3,
  # GROUP2 - TAM
  "TAM-3A-4-5CDE-A", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-A", "GROUP2", "TAM", 4,
  "TAM-3A-4-5CDE-B", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-B", "GROUP2", "TAM", 5,
  "TAM-3A-4-5CDE-RUN2-C", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-RUN2-C", "GROUP2", "TAM", 6,
  "TAM-3A-4-5CDE-RUN2-D", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-RUN2-D", "GROUP2", "TAM", 7,
  "TAM-3A-4-5CDE-RUN2-E", "GROUP2/TAM-3A-4-5CDE/TAM-3A-4-5CDE-RUN2-E", "GROUP2", "TAM", 8,
  # GROUP3 - TAM (youngest TAM)
  "TAM-5AB-6-7-A", "GROUP3/TAM-5AB-6-7/TAM-5AB-6-7-A", "GROUP3", "TAM", 9,
  "TAM-5AB-6-7-B-RUN2", "GROUP3/TAM-5AB-6-7/TAM-5AB-6-7-B-RUN2", "GROUP3", "TAM", 10,
  "TAM-5AB-6-7-C", "GROUP3/TAM-5AB-6-7/TAM-5AB-6-7-C", "GROUP3", "TAM", 11,
  # GROUP4 - SC (oldest SC)
  "SC-1ABC-2-3C-A-RUN1", "GROUP4/SC-1ABC-2-3C-RUN1/SC-1ABC-2-3C-A-RUN1", "GROUP4", "SC", 12,
  "SC-1ABC-2-3C-RUN2-B", "GROUP4/SC-1ABC-2-3C-RUN2/SC-1ABC-2-3C-RUN2-B", "GROUP4", "SC", 13,
  "SC-1ABC-2-3C-RUN2-C", "GROUP4/SC-1ABC-2-3C-RUN2/SC-1ABC-2-3C-RUN2-C", "GROUP4", "SC", 14,
  # GROUP5 - SC
  "SC-3AB-4ABCD-A", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-A", "GROUP5", "SC", 15,
  "SC-3AB-4ABCD-B", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-B", "GROUP5", "SC", 16,
  "SC-3AB-4ABCD-C", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-C", "GROUP5", "SC", 17,
  "SC-3AB-4ABCD-D", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-D", "GROUP5", "SC", 18,
  "SC-3AB-4ABCD-RUN2-D", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-RUN2-D", "GROUP5", "SC", 19,
  "SC-3AB-4ABCD-RUN2-E", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-RUN2-E", "GROUP5", "SC", 20,
  "SC-3AB-4ABCD-RUN2-F", "GROUP5/SC-3AB-4ABCD/SC-3AB-4ABCD-RUN2-F", "GROUP5", "SC", 21,
  # GROUP6 - SC
  "SC-5-6-7ABC-A", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-A", "GROUP6", "SC", 22,
  "SC-5-6-7ABC-B", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-B", "GROUP6", "SC", 23,
  "SC-5-6-7ABC-C", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-C", "GROUP6", "SC", 24,
  "SC-5-6-7ABC-D", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-D", "GROUP6", "SC", 25,
  "SC-5-6-7ABC-E", "GROUP6/SC-5-6-7ABC/SC-5-6-7ABC-E", "GROUP6", "SC", 26,
  # GROUP7 - SC (youngest)
  "SC8-A", "GROUP7/SC8-A/SC8-A", "GROUP7", "SC", 27
)

# Available proxies/elements - complete XRF suite
# Named vectors: display label = internal value
element_choices <- c(
  # Major elements
  "Al" = "Al", "Si" = "Si", "K" = "K", "Ca" = "Ca", "Ti" = "Ti", "Fe" = "Fe",
  # Minor elements
  "P" = "P", "S" = "S", "Mn" = "Mn", "Rb" = "Rb", "Sr" = "Sr", "Zr" = "Zr", "Ba" = "Ba",
  # Trace elements
  "V" = "V", "Cr" = "Cr", "Ni" = "Ni", "Cu" = "Cu", "Zn" = "Zn",
  "As" = "As", "Y" = "Y", "Nb" = "Nb", "Pb" = "Pb"
)

ratio_choices <- c(
  # Primary paleoenvironmental ratios
  "Ca/Ti (Carbonate input)" = "Ca_Ti",
  "Fe/Mn (Redox conditions)" = "Fe_Mn",
  "K/Ti (Chemical weathering)" = "K_Ti",
  "Zr/Rb (Grain size/energy)" = "Zr_Rb",
  # Additional ratios
  "Sr/Ca (Salinity/aragonite)" = "Sr_Ca",
  "Rb/Sr (Weathering intensity)" = "Rb_Sr",
  "Si/Al (Biogenic silica)" = "Si_Al",
  "Fe/Ti (Fe enrichment)" = "Fe_Ti",
  "Mn/Ti (Mn enrichment)" = "Mn_Ti",
  "K/Al (Clay mineralogy)" = "K_Al",
  "Sr/Ti (Marine influence)" = "Sr_Ti",
  "Ba/Sr (Productivity)" = "Ba_Sr",
  "Cu/Zn (Metal source)" = "Cu_Zn",
  "V/Cr (Paleoredox)" = "V_Cr"
)

# Default selections (Ca, Ti, Fe, Mn and their ratios)
default_elements <- c("Ca", "Ti", "Fe", "Mn")
default_ratios <- c("Ca_Ti", "Fe_Mn")

proxy_colors <- c(
  # Major elements
  "Al" = "#636363", "Si" = "#bdbdbd", "K" = "#9467bd", "Ca" = "#1f77b4",

  "Ti" = "#8c564b", "Fe" = "#d62728",
  # Minor elements
  "P" = "#98df8a", "S" = "#ffbb78", "Mn" = "#ff7f0e", "Rb" = "#7f7f7f",
  "Sr" = "#17becf", "Zr" = "#2ca02c", "Ba" = "#e377c2",
  # Trace elements
  "V" = "#bcbd22", "Cr" = "#aec7e8", "Ni" = "#c49c94", "Cu" = "#f7b6d2",
  "Zn" = "#c7c7c7", "As" = "#dbdb8d", "Y" = "#9edae5", "Nb" = "#393b79",
  "Pb" = "#843c39",
  # Ratios
  "Ca_Ti" = "#2166ac", "Fe_Mn" = "#7b3294", "K_Ti" = "#e08214",
  "Zr_Rb" = "#1a9850", "Sr_Ca" = "#d95f02", "Rb_Sr" = "#e7298a",
  "Si_Al" = "#5254a3", "Fe_Ti" = "#ad494a", "Mn_Ti" = "#637939",
  "K_Al" = "#8c6d31", "Sr_Ti" = "#31a354", "Ba_Sr" = "#756bb1",
  "Cu_Zn" = "#de9ed6", "V_Cr" = "#6baed6"
)

proxy_labels <- c(
  # Major elements
  "Al" = "Al (cps)", "Si" = "Si (cps)", "K" = "K (cps)", "Ca" = "Ca (cps)",
  "Ti" = "Ti (cps)", "Fe" = "Fe (cps)",
  # Minor elements
  "P" = "P (cps)", "S" = "S (cps)", "Mn" = "Mn (cps)", "Rb" = "Rb (cps)",
  "Sr" = "Sr (cps)", "Zr" = "Zr (cps)", "Ba" = "Ba (cps)",
  # Trace elements
  "V" = "V (cps)", "Cr" = "Cr (cps)", "Ni" = "Ni (cps)", "Cu" = "Cu (cps)",
  "Zn" = "Zn (cps)", "As" = "As (cps)", "Y" = "Y (cps)", "Nb" = "Nb (cps)",
  "Pb" = "Pb (cps)",
  # Ratios
  "Ca_Ti" = "Ca/Ti", "Fe_Mn" = "Fe/Mn", "K_Ti" = "K/Ti", "Zr_Rb" = "Zr/Rb",
  "Sr_Ca" = "Sr/Ca", "Rb_Sr" = "Rb/Sr", "Si_Al" = "Si/Al", "Fe_Ti" = "Fe/Ti",
  "Mn_Ti" = "Mn/Ti", "K_Al" = "K/Al", "Sr_Ti" = "Sr/Ti", "Ba_Sr" = "Ba/Sr",
  "Cu_Zn" = "Cu/Zn", "V_Cr" = "V/Cr"
)

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("Loading XRF data...")
xrf_data <- read_csv(file.path(output_path, "tables", "xrf_data_stacked.csv"),
                     show_col_types = FALSE) %>%
  mutate(
    # Additional ratios not in stacked data
    Sr_Ca = Sr / Ca,
    Si_Al = Si / Al,
    Fe_Ti = Fe / Ti,
    Mn_Ti = Mn / Ti,
    K_Al = K / Al,
    Sr_Ti = Sr / Ti,
    Ba_Sr = Ba / Sr,
    Cu_Zn = Cu / Zn,
    V_Cr = V / Cr
  )

xrf_valid <- xrf_data %>% filter(qc_pass, !excluded)

exclusion_zones <- read_csv(file.path(base_path, "data", "exclusion_zones.csv"),
                            show_col_types = FALSE) %>%
  filter(!is.na(exclude_start_mm))

message(sprintf("Loaded %d valid measurements from %d sections",
                nrow(xrf_valid), n_distinct(xrf_valid$section)))

# Get sections with valid data
valid_sections <- xrf_valid %>%
  group_by(section) %>%
  summarise(n_valid = n(), .groups = "drop") %>%
  filter(n_valid >= 10) %>%
  pull(section)

# Filter section_config to only include sections with valid data
section_config <- section_config %>%
  filter(section %in% valid_sections)

message(sprintf("Sections with sufficient data: %d", nrow(section_config)))

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
  list(xrf_start = xrf_start, xrf_end = xrf_end,
       optical_start = optical_start, optical_end = optical_end)
}

# Dynamic image enhancement function with adjustable parameters
brighten_image_dynamic <- function(img, brightness = 80, contrast = 25, gamma = 2.0) {
  img %>%
    image_level(black_point = 2, white_point = 60,
                mid_point = 1/gamma) %>%
    image_normalize() %>%
    image_modulate(brightness = 100 + brightness) %>%
    image_contrast(sharpen = contrast)
}

overlay_excluded_sections <- function(img, section_name, xrf_start, xrf_end, pixels_per_mm) {
  excl_zones <- exclusion_zones %>%
    filter(section == section_name) %>%
    arrange(exclude_start_mm)
  if (nrow(excl_zones) == 0) return(img)
  img_info <- image_info(img)
  for (i in seq_len(nrow(excl_zones))) {
    zone <- excl_zones[i, ]
    zone_start_px <- round((zone$exclude_start_mm - xrf_start) * pixels_per_mm)
    zone_end_px <- round((zone$exclude_end_mm - xrf_start) * pixels_per_mm)
    zone_start_px <- max(0, zone_start_px)
    zone_end_px <- min(img_info$width, zone_end_px)
    zone_width <- zone_end_px - zone_start_px
    if (zone_width > 0) {
      overlay <- image_blank(zone_width, img_info$height, color = "#E0E0E0")
      img <- image_composite(img, overlay, operator = "blend",
                              compose_args = "50", offset = sprintf("+%d+0", zone_start_px))
    }
  }
  img
}

process_core_image <- function(section_name, optical_path, data_start, data_end,
                               brightness = 80, contrast = 25, gamma = 2.0) {
  section_dir <- file.path(data_path, optical_path)
  img_path <- file.path(section_dir, "optical.tif")
  doc_path <- file.path(section_dir, "document.txt")

  if (!file.exists(img_path)) return(NULL)
  params <- read_document_params(doc_path)
  if (is.null(params)) return(NULL)

  img <- image_read(img_path)
  img_info <- image_info(img)
  optical_range <- params$optical_end - params$optical_start
  pixels_per_mm <- img_info$width / optical_range

  xrf_start_px <- max(0, round((data_start - params$optical_start) * pixels_per_mm))
  xrf_end_px <- min(img_info$width, round((data_end - params$optical_start) * pixels_per_mm))
  crop_width <- xrf_end_px - xrf_start_px
  if (crop_width <= 0) return(NULL)

  crop_geom <- sprintf("%dx%d+%d+0", crop_width, img_info$height, xrf_start_px)
  img_cropped <- image_crop(img, crop_geom)
  cropped_ppm <- crop_width / (data_end - data_start)

  img_bright <- brighten_image_dynamic(img_cropped, brightness, contrast, gamma)
  # Note: Exclusion zone shading is now handled by ggplot overlay (add_exclusion_overlay_core)
  # which responds to the "Shade Excluded Zones" toggle

  target_width <- max(crop_width * 4, 2000)
  img_hires <- image_resize(img_bright, sprintf("%dx", target_width))

  list(image = img_hires, section = section_name,
       xrf_start = data_start, xrf_end = data_end)
}

# ==============================================================================
# UI
# ==============================================================================

ui <- fluidPage(

  titlePanel("XRF Core Figures"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("1. Select Site"),
      radioButtons("site", NULL,
                   choices = c("Tamshiyacu (TAM)" = "TAM", "Santa Corina (SC)" = "SC"),
                   selected = "TAM"),

      hr(),

      h4("2. Select Sections"),
      helpText("Sections are listed in stratigraphic order (top = oldest)"),
      uiOutput("section_selector"),

      actionButton("select_all", "Select All", class = "btn-sm"),
      actionButton("select_none", "Clear All", class = "btn-sm"),

      hr(),

      h4("3. Select Elements"),
      pickerInput("elements", NULL,
                  choices = element_choices,
                  selected = default_elements,
                  multiple = TRUE,
                  options = pickerOptions(
                    actionsBox = TRUE,
                    liveSearch = TRUE,
                    size = 10,
                    title = "Choose elements..."
                  )),

      hr(),

      h4("4. Select Ratios"),
      pickerInput("ratios", NULL,
                  choices = ratio_choices,
                  selected = default_ratios,
                  multiple = TRUE,
                  options = pickerOptions(
                    actionsBox = TRUE,
                    liveSearch = TRUE,
                    size = 10,
                    title = "Choose ratios..."
                  )),

      # Custom ratio builder
      tags$details(
        tags$summary(tags$b("Custom Ratio")),
        fluidRow(
          column(5, selectInput("custom_num", "Numerator",
                                choices = c("", names(element_choices)), selected = "")),
          column(2, tags$div(style = "text-align: center; padding-top: 30px;", "/")),
          column(5, selectInput("custom_den", "Denominator",
                                choices = c("", names(element_choices)), selected = ""))
        ),
        actionButton("add_custom", "Add Custom Ratio", class = "btn-sm btn-info"),
        uiOutput("custom_ratio_display")
      ),

      hr(),

      h4("5. Display Options"),
      checkboxInput("show_core", "Show Core Image", value = TRUE),
      checkboxInput("sediment_filter", "Sediment Filter (collapse excluded zones)", value = FALSE),
      checkboxInput("show_exclusions", "Shade Excluded Zones", value = TRUE),
      checkboxInput("show_points", "Show Data Points", value = TRUE),
      checkboxInput("show_boundaries", "Show Section Boundaries", value = TRUE),
      sliderInput("smooth_window", "Smoothing Window", min = 1, max = 15, value = 5, step = 2),

      hr(),

      h4("6. Image Enhancement"),
      sliderInput("img_brightness", "Brightness", min = 0, max = 150, value = 80, step = 5),
      sliderInput("img_contrast", "Contrast", min = 0, max = 50, value = 25, step = 5),
      sliderInput("img_gamma", "Gamma", min = 0.5, max = 4.0, value = 2.0, step = 0.1),

      hr(),

      h4("7. Export"),
      numericInput("fig_width", "Width (inches)", value = 16, min = 8, max = 30),
      numericInput("fig_height", "Height (inches)", value = 10, min = 6, max = 30),
      numericInput("fig_dpi", "DPI", value = 300, min = 100, max = 600),
      downloadButton("download_plot", "Download Figure", class = "btn-primary"),

      hr(),

      h4("8. Download Data"),
      helpText("Download raw XRF data for selected sections"),
      checkboxInput("include_smoothed", "Include smoothed values", value = TRUE),
      checkboxInput("include_all_elements", "Include all elements (not just selected)", value = FALSE),
      downloadButton("download_data", "Download Data (CSV)", class = "btn-success")
    ),

    mainPanel(
      width = 9,

      tabsetPanel(
        tabPanel("Figure",
                 br(),
                 plotOutput("strat_plot", height = "1200px")
        ),
        tabPanel("Site Comparison",
                 br(),
                 fluidRow(
                   column(4,
                     h4("Distribution Plot Settings"),
                     selectInput("density_proxy", "Select Element or Ratio",
                                 choices = c("Elements" = "",
                                             setNames(names(element_choices), names(element_choices)),
                                             "Ratios" = "",
                                             ratio_choices),
                                 selected = "Fe_Mn"),
                     checkboxInput("density_log", "Log scale", value = TRUE),
                     numericInput("density_threshold", "Threshold line (optional)", value = NA),
                     textInput("density_threshold_label", "Threshold label", value = ""),
                     hr(),
                     h5("Custom Ratio"),
                     fluidRow(
                       column(5, selectInput("density_custom_num", "Numerator",
                                             choices = c("", names(element_choices)))),
                       column(2, tags$div(style = "text-align: center; padding-top: 30px;", "/")),
                       column(5, selectInput("density_custom_den", "Denominator",
                                             choices = c("", names(element_choices))))
                     ),
                     actionButton("density_use_custom", "Use Custom Ratio", class = "btn-sm btn-info"),
                     hr(),
                     numericInput("density_width", "Width (inches)", value = 7, min = 4, max = 14),
                     numericInput("density_height", "Height (inches)", value = 5, min = 3, max = 10),
                     downloadButton("download_density", "Download Figure", class = "btn-primary")
                   ),
                   column(8,
                     plotOutput("density_plot", height = "500px"),
                     br(),
                     verbatimTextOutput("density_stats")
                   )
                 )
        ),
        tabPanel("Data Summary",
                 br(),
                 verbatimTextOutput("data_summary"),
                 br(),
                 DT::dataTableOutput("section_table")
        ),
        tabPanel("Temporal Evolution",
                 br(),
                 fluidRow(
                   column(3,
                     h4("View Mode"),
                     radioButtons("temporal_view_mode", NULL,
                                  choices = c("Time Series" = "timeseries",
                                              "Core Comparison" = "cores"),
                                  selected = "timeseries"),

                     hr(),

                     h4("Site Selection"),
                     checkboxGroupInput("temporal_sites", NULL,
                                        choices = c("TAM (Tamshiyacu)" = "TAM",
                                                    "SC (Santa Corina)" = "SC"),
                                        selected = c("TAM", "SC")),

                     hr(),

                     h4("Proxy Selection"),
                     pickerInput("temporal_proxies", NULL,
                                 choices = c("Mn/Ti (Redox)" = "Mn_Ti",
                                             "Ca/Ti (Carbonate)" = "Ca_Ti",
                                             "Fe/Ti (Iron)" = "Fe_Ti",
                                             "Fe/Mn (Reducing)" = "Fe_Mn",
                                             "K/Ti (Weathering)" = "K_Ti",
                                             "Zr/Rb (Grain size)" = "Zr_Rb",
                                             "Sr/Ca (Salinity)" = "Sr_Ca"),
                                 selected = c("Mn_Ti"),
                                 multiple = TRUE,
                                 options = pickerOptions(
                                   actionsBox = TRUE,
                                   maxOptions = 4,
                                   maxOptionsText = "Max 4 proxies"
                                 )),

                     hr(),

                     conditionalPanel(
                       condition = "input.temporal_view_mode == 'timeseries'",
                       h4("Time Range (Ma)"),
                       fluidRow(
                         column(6, numericInput("temporal_age_min", "Youngest",
                                                value = 12.8, min = 12.5, max = 14.5, step = 0.1)),
                         column(6, numericInput("temporal_age_max", "Oldest",
                                                value = 14.5, min = 12.5, max = 14.5, step = 0.1))
                       ),
                       actionButton("temporal_reset_range", "Reset to Full Range", class = "btn-sm"),
                       actionButton("temporal_zoom_overlap", "Zoom to Overlap", class = "btn-sm btn-warning"),
                       hr()
                     ),

                     conditionalPanel(
                       condition = "input.temporal_view_mode == 'cores'",
                       h4("Section Selection"),
                       helpText("Select sections for side-by-side comparison"),
                       uiOutput("temporal_tam_section_selector"),
                       uiOutput("temporal_sc_section_selector"),
                       hr(),
                       h4("Image Enhancement"),
                       sliderInput("temporal_brightness", "Brightness", min = 0, max = 150, value = 80, step = 5),
                       sliderInput("temporal_contrast", "Contrast", min = 0, max = 50, value = 25, step = 5),
                       hr()
                     ),

                     h4("Display Options"),
                     checkboxInput("temporal_show_overlap", "Highlight Overlap Period", value = TRUE),
                     checkboxInput("temporal_show_points", "Show Data Points", value = TRUE),
                     checkboxInput("temporal_log_scale", "Log Scale Y-axis", value = TRUE),
                     sliderInput("temporal_smooth", "Smoothing Window", min = 1, max = 15, value = 5, step = 2),

                     hr(),

                     tags$details(
                       tags$summary(tags$b("Age Model Info")),
                       tags$div(
                         style = "font-size: 11px; padding: 5px;",
                         tags$b("TAM:"), " 12.935-13.446 Ma", tags$br(),
                         "SR: 778 cm/Ma (2.6 yr/step)", tags$br(), tags$br(),
                         tags$b("SC:"), " 13.275-14.298 Ma", tags$br(),
                         "SR: 414 cm/Ma (4.9 yr/step)", tags$br(), tags$br(),
                         tags$b("Overlap:"), " 13.275-13.446 Ma", tags$br(),
                         "Duration: 171 kyr"
                       )
                     ),

                     hr(),

                     h4("Export"),
                     numericInput("temporal_width", "Width (inches)", value = 12, min = 6, max = 20),
                     numericInput("temporal_height", "Height (inches)", value = 8, min = 4, max = 16),
                     downloadButton("download_temporal", "Download Figure", class = "btn-primary")
                   ),
                   column(9,
                     conditionalPanel(
                       condition = "input.temporal_view_mode == 'timeseries'",
                       h4("Temporal Evolution"),
                       plotOutput("temporal_plot", height = "550px"),
                       br(),
                       h4("Overlap Comparison (13.275-13.446 Ma)"),
                       plotOutput("overlap_plot", height = "350px"),
                       br(),
                       verbatimTextOutput("temporal_stats")
                     ),
                     conditionalPanel(
                       condition = "input.temporal_view_mode == 'cores'",
                       h4("Core Comparison"),
                       helpText("Side-by-side view of TAM and SC cores with vertical XRF profiles"),
                       plotOutput("temporal_core_plot", height = "900px"),
                       br(),
                       verbatimTextOutput("temporal_core_stats")
                     )
                   )
                 )
        ),
        tabPanel("Exclusion Editor",
                 br(),
                 fluidRow(
                   column(3,
                     h4("1. Select Section"),
                     selectInput("editor_section", "Section to edit:",
                                 choices = NULL),  # Populated dynamically
                     hr(),
                     h4("2. Mark Exclusion Zone"),
                     helpText("Click and drag on the plot to select a region, then click 'Add Zone'"),
                     verbatimTextOutput("brush_info"),
                     fluidRow(
                       column(6, actionButton("add_zone", "Add Zone", class = "btn-success btn-sm")),
                       column(6, actionButton("clear_brush", "Clear Selection", class = "btn-warning btn-sm"))
                     ),
                     hr(),
                     h4("3. Current Exclusions"),
                     helpText("Select rows and click 'Delete Selected' to remove"),
                     DT::dataTableOutput("exclusion_table"),
                     br(),
                     actionButton("delete_zone", "Delete Selected", class = "btn-danger btn-sm"),
                     hr(),
                     h4("4. Save Changes"),
                     helpText("Save changes to exclusion_zones.csv"),
                     actionButton("save_exclusions", "Save to CSV", class = "btn-primary"),
                     br(), br(),
                     verbatimTextOutput("save_status")
                   ),
                   column(9,
                     h4("Core Image with XRF Data"),
                     helpText("Gray shading = existing exclusion zones. Drag to select new zones."),
                     plotOutput("editor_plot", height = "800px",
                                brush = brushOpts(id = "editor_brush", direction = "y",
                                                  fill = "red", stroke = "darkred", opacity = 0.3))
                   )
                 )
        ),
        tabPanel("Help",
                 br(),
                 h4("How to Use"),
                 tags$ol(
                   tags$li("Select a site (TAM or SC)"),
                   tags$li("Choose sections to display (in stratigraphic order)"),
                   tags$li("Select elements and/or ratios to plot"),
                   tags$li("Adjust display options as needed"),
                   tags$li("Download the figure when satisfied")
                 ),
                 br(),
                 h4("Site Comparison Tab"),
                 p("Generate density distribution plots comparing TAM and SC for any element or ratio. Useful for visualizing differences between sites."),
                 br(),
                 h4("Gray Shading"),
                 p("Gray shaded areas indicate sections where XRF measurements were excluded due to uneven surface or foam fill.")
        )
      )
    )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {

  # ==============================================================================
  # EXCLUSION EDITOR: Reactive values for managing exclusion zones
  # ==============================================================================

  # Working copy of exclusion zones (allows edits before saving)
  working_exclusions <- reactiveVal({
    read_csv(file.path(base_path, "data", "exclusion_zones.csv"), show_col_types = FALSE)
  })

  # Track if there are unsaved changes
  unsaved_changes <- reactiveVal(FALSE)

  # Populate section selector for editor
  observe({
    choices <- setNames(section_config$section,
                        paste0(section_config$group, " - ", section_config$section))
    updateSelectInput(session, "editor_section", choices = choices)
  })

  # Dynamic section choices based on site
  output$section_selector <- renderUI({
    site_sections <- section_config %>%
      filter(site == input$site) %>%
      arrange(strat_order)

    choices <- setNames(site_sections$section,
                        paste0(site_sections$group, " - ", site_sections$section))

    checkboxGroupInput("sections", NULL,
                       choices = choices,
                       selected = choices[1:min(3, length(choices))])
  })

  # Select all/none buttons
  observeEvent(input$select_all, {
    site_sections <- section_config %>%
      filter(site == input$site) %>%
      pull(section)
    updateCheckboxGroupInput(session, "sections", selected = site_sections)
  })

  observeEvent(input$select_none, {
    updateCheckboxGroupInput(session, "sections", selected = character(0))
  })

  # Custom ratio management
  custom_ratios <- reactiveVal(list())

  observeEvent(input$add_custom, {
    req(input$custom_num, input$custom_den)
    req(input$custom_num != "", input$custom_den != "")
    req(input$custom_num != input$custom_den)

    num <- input$custom_num
    den <- input$custom_den
    ratio_name <- paste0(num, "_", den)
    ratio_label <- paste0(num, "/", den)

    # Add to custom ratios list
    current <- custom_ratios()
    current[[ratio_name]] <- list(num = num, den = den, label = ratio_label)
    custom_ratios(current)

    # Reset dropdowns
    updateSelectInput(session, "custom_num", selected = "")
    updateSelectInput(session, "custom_den", selected = "")
  })

  output$custom_ratio_display <- renderUI({
    ratios <- custom_ratios()
    if (length(ratios) == 0) return(NULL)

    tags$div(
      style = "margin-top: 10px;",
      tags$b("Active custom ratios:"),
      tags$ul(
        lapply(names(ratios), function(r) {
          tags$li(
            ratios[[r]]$label,
            actionLink(paste0("remove_", r), "(remove)", style = "color: red; font-size: 10px;")
          )
        })
      )
    )
  })

  # Remove custom ratio observers
  observe({
    ratios <- custom_ratios()
    lapply(names(ratios), function(r) {
      observeEvent(input[[paste0("remove_", r)]], {
        current <- custom_ratios()
        current[[r]] <- NULL
        custom_ratios(current)
      }, ignoreInit = TRUE, once = TRUE)
    })
  })

  # Get selected proxies
  selected_proxies <- reactive({
    custom <- names(custom_ratios())
    c(input$elements, input$ratios, custom)
  })

  # Process data for selected sections
  processed_data <- reactive({
    req(input$sections)
    req(length(selected_proxies()) > 0)

    withProgress(message = "Processing data...", {

      # Filter to selected sections (already excludes QC failures)
      data <- xrf_valid %>%
        filter(section %in% input$sections) %>%
        arrange(section, position_mm)

      if (nrow(data) == 0) return(NULL)

      # Get section order
      section_order <- section_config %>%
        filter(section %in% input$sections) %>%
        arrange(strat_order) %>%
        pull(section)

      # Get FULL data ranges (including excluded) for coordinate alignment
      # This ensures graph positions match image positions
      full_data_ranges <- xrf_data %>%
        filter(section %in% input$sections) %>%
        group_by(section) %>%
        summarise(data_start = min(position_mm), data_end = max(position_mm), .groups = "drop")

      # Calculate cumulative depth for stacking
      cumulative_offset <- 0
      result <- tibble()

      for (sect in section_order) {
        sect_data <- data %>% filter(section == sect)
        if (nrow(sect_data) == 0) next

        # Get full range for this section (for coordinate alignment with image)
        sect_full_range <- full_data_ranges %>% filter(section == sect)

        if (input$sediment_filter) {
          # SEDIMENT FILTER MODE: Only show valid sediment measurements
          # Use sequential positioning based on actual valid measurements
          # This ensures graph matches cropped image (both exclude same regions)

          sect_data <- sect_data %>%
            arrange(position_mm) %>%
            mutate(
              measurement_order = row_number(),
              # Each measurement at 3mm intervals (matching XRF scan resolution)
              cumulative_depth = (measurement_order - 1) * 3 + cumulative_offset,
              cumulative_depth_cm = cumulative_depth / 10
            )

          # Section span = number of valid measurements * 3mm spacing
          section_span <- nrow(sect_data) * 3
          cumulative_offset <- cumulative_offset + section_span + 30

        } else {
          # NORMAL MODE: Physical position preserved (gaps where excluded)
          # CRITICAL: Use FULL data range (same as image) for coordinate alignment
          # This ensures valid measurements align with their positions on the core image
          data_start <- sect_full_range$data_start
          data_end <- sect_full_range$data_end

          sect_data <- sect_data %>%
            mutate(
              # Position relative to full data start (not just valid data start)
              relative_depth = position_mm - data_start,
              cumulative_depth = relative_depth + cumulative_offset,
              cumulative_depth_cm = cumulative_depth / 10
            )
          # Use full data span for offset calculation
          cumulative_offset <- cumulative_offset + (data_end - data_start) + 30
        }

        # Calculate custom ratios
        custom <- custom_ratios()
        for (ratio_name in names(custom)) {
          num <- custom[[ratio_name]]$num
          den <- custom[[ratio_name]]$den
          if (num %in% names(sect_data) && den %in% names(sect_data)) {
            sect_data[[ratio_name]] <- sect_data[[num]] / sect_data[[den]]
          }
        }

        # Apply smoothing to all proxies
        for (proxy in selected_proxies()) {
          if (proxy %in% names(sect_data)) {
            smooth_col <- paste0(proxy, "_smooth")
            sect_data[[smooth_col]] <- zoo::rollmean(sect_data[[proxy]],
                                                      input$smooth_window,
                                                      fill = NA, align = "center")
          }
        }

        result <- bind_rows(result, sect_data)
      }

      result
    })
  })

  # Get section boundary positions for red lines
  section_boundaries <- reactive({
    req(input$sections)

    section_order <- section_config %>%
      filter(section %in% input$sections) %>%
      arrange(strat_order) %>%
      pull(section)

    cumulative_offset <- 0
    boundaries <- c()

    if (input$sediment_filter) {
      # SEDIMENT FILTER MODE: boundaries based on measurement counts
      # Must match processed_data logic: each measurement at 3mm intervals
      for (i in seq_along(section_order)) {
        sect <- section_order[i]
        sect_data <- xrf_valid %>% filter(section == sect)

        if (nrow(sect_data) == 0) next

        # Section length = number of valid measurements * 3mm
        section_length <- nrow(sect_data) * 3

        if (i < length(section_order)) {
          boundary_pos <- (cumulative_offset + section_length + 15) / 10
          boundaries <- c(boundaries, boundary_pos)
        }

        cumulative_offset <- cumulative_offset + section_length + 30
      }
    } else {
      # NORMAL MODE: boundaries based on FULL data range (same as image and graph)
      data_ranges <- xrf_data %>%
        filter(section %in% input$sections) %>%
        group_by(section) %>%
        summarise(data_start = min(position_mm), data_end = max(position_mm), .groups = "drop")

      for (i in seq_along(section_order)) {
        sect <- section_order[i]
        range_info <- data_ranges %>% filter(section == sect)
        if (nrow(range_info) == 0) next

        # Use full data span (matches image and graph positioning)
        section_length <- range_info$data_end - range_info$data_start

        if (i < length(section_order)) {
          boundary_pos <- (cumulative_offset + section_length + 15) / 10
          boundaries <- c(boundaries, boundary_pos)
        }

        cumulative_offset <- cumulative_offset + section_length + 30
      }
    }

    boundaries
  })

  # Get exclusion zones for selected sections
  section_exclusions <- reactive({
    req(input$sections)

    # No exclusion overlays in sediment filter mode (they're cropped out)
    if (input$sediment_filter) {
      return(tibble(start_cm = numeric(), end_cm = numeric()))
    }

    # Get data ranges per section
    data_ranges <- xrf_data %>%
      filter(section %in% input$sections) %>%
      group_by(section) %>%
      summarise(data_start = min(position_mm), data_end = max(position_mm), .groups = "drop")

    # Calculate cumulative offsets
    section_order <- section_config %>%
      filter(section %in% input$sections) %>%
      arrange(strat_order) %>%
      pull(section)

    cumulative_offset <- 0
    offsets <- tibble()

    for (sect in section_order) {
      range_info <- data_ranges %>% filter(section == sect)
      if (nrow(range_info) == 0) next

      offsets <- bind_rows(offsets, tibble(
        section = sect,
        offset = cumulative_offset,
        data_start = range_info$data_start,
        data_end = range_info$data_end
      ))

      cumulative_offset <- cumulative_offset + (range_info$data_end - range_info$data_start) + 30
    }

    # Join with exclusion zones and calculate positions
    # Clamp exclusion zones to actual data range (some zones extend beyond data)
    exclusion_zones %>%
      filter(section %in% input$sections) %>%
      left_join(offsets, by = "section") %>%
      mutate(
        # Clamp exclusion boundaries to data range
        excl_start_clamped = pmax(exclude_start_mm, data_start),
        excl_end_clamped = pmin(exclude_end_mm, data_end),
        # Calculate cumulative positions
        start_cm = (excl_start_clamped - data_start + offset) / 10,
        end_cm = (excl_end_clamped - data_start + offset) / 10
      ) %>%
      # Only keep zones that have positive width after clamping
      filter(!is.na(start_cm), !is.na(end_cm), end_cm > start_cm)
  })

  # Generate plot
  generate_plot <- reactive({
    req(processed_data())
    req(length(selected_proxies()) > 0)

    data <- processed_data()
    proxies <- selected_proxies()
    excl <- section_exclusions()

    depth_max <- max(data$cumulative_depth_cm, na.rm = TRUE)
    depth_min <- 0

    # Theme
    theme_strat <- function(base_size = 11) {
      theme_minimal(base_size = base_size) +
        theme(
          plot.title = element_text(size = base_size + 1, face = "bold"),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size - 1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
          legend.position = "none",
          plot.margin = margin(5, 10, 5, 5)
        )
    }

    # Function to add exclusion overlay (grey shading for styrofoam) - CORE IMAGE ONLY
    add_exclusion_overlay_core <- function(p) {
      if (input$show_exclusions && nrow(excl) > 0) {
        for (i in seq_len(nrow(excl))) {
          p <- p + annotate("rect",
                            xmin = -Inf, xmax = Inf,
                            ymin = excl$start_cm[i], ymax = excl$end_cm[i],
                            fill = "gray50", alpha = 0.6)
        }
      }
      p
    }

    # Function to add section boundary lines - CORE IMAGE ONLY
    add_section_boundaries_core <- function(p) {
      if (input$show_boundaries && length(boundaries) > 0) {
        p <- p + geom_hline(yintercept = boundaries, color = "red",
                            linewidth = 0.8, linetype = "solid")
      }
      p
    }

    # Get section boundaries for red lines
    boundaries <- section_boundaries()

    # Create plot list
    plot_list <- list()

    # Core image panel (if enabled)
    if (input$show_core) {
      # Process and stack core images
      withProgress(message = "Processing core images...", {

        section_order <- section_config %>%
          filter(section %in% input$sections) %>%
          arrange(strat_order) %>%
          pull(section)

        # Get data ranges (full range including excluded)
        data_ranges <- xrf_data %>%
          filter(section %in% input$sections) %>%
          group_by(section) %>%
          summarise(data_start = min(position_mm), data_end = max(position_mm), .groups = "drop")

        images <- list()
        for (sect in section_order) {
          sect_info <- section_config %>% filter(section == sect)
          range_info <- data_ranges %>% filter(section == sect)

          if (nrow(range_info) > 0) {
            img_result <- tryCatch(
              process_core_image(sect, sect_info$optical_path,
                                  range_info$data_start, range_info$data_end,
                                  brightness = input$img_brightness,
                                  contrast = input$img_contrast,
                                  gamma = input$img_gamma),
              error = function(e) NULL
            )
            if (!is.null(img_result)) {
              images[[sect]] <- img_result
            }
          }
        }

        if (length(images) > 0) {
          img_list <- list()

          if (input$sediment_filter) {
            # SEDIMENT FILTER MODE: Crop out excluded zones from images
            # CRITICAL: Base cropping on ACTUAL VALID MEASUREMENT POSITIONS
            # not on exclusion zone boundaries (which may not align exactly)
            #
            # The graph uses sequential positioning: each measurement at 3mm intervals
            # Image segments must be cropped to actual measurement ranges and scaled by count

            target_width <- 400  # Fixed width after rotation
            pixels_per_measurement <- 5  # Base scale factor for height

            for (i in seq_along(section_order)) {
              sect <- section_order[i]
              if (!(sect %in% names(images))) next

              img <- images[[sect]]$image
              img_info <- image_info(img)
              sect_range <- data_ranges %>% filter(section == sect)

              if (nrow(sect_range) == 0) next

              sect_start <- sect_range$data_start
              sect_end <- sect_range$data_end
              sect_length <- sect_end - sect_start
              px_per_mm <- img_info$width / sect_length

              # Get valid measurements for this section, sorted by position
              sect_valid <- xrf_valid %>%
                filter(section == sect) %>%
                arrange(position_mm)

              if (nrow(sect_valid) == 0) next

              # Find continuous groups of measurements (gaps > 10mm indicate separate segments)
              positions <- sect_valid$position_mm
              gaps <- diff(positions)
              gap_indices <- which(gaps > 10)  # Gaps larger than 10mm indicate foam

              # Build valid ranges from actual measurement positions
              valid_ranges <- list()
              start_idx <- 1

              for (gap_idx in gap_indices) {
                # Range from start_idx to gap_idx
                range_start <- positions[start_idx]
                range_end <- positions[gap_idx]
                n_meas <- gap_idx - start_idx + 1
                valid_ranges[[length(valid_ranges) + 1]] <- list(
                  start = range_start,
                  end = range_end,
                  n_measurements = n_meas
                )
                start_idx <- gap_idx + 1
              }

              # Add final range
              if (start_idx <= length(positions)) {
                range_start <- positions[start_idx]
                range_end <- positions[length(positions)]
                n_meas <- length(positions) - start_idx + 1
                valid_ranges[[length(valid_ranges) + 1]] <- list(
                  start = range_start,
                  end = range_end,
                  n_measurements = n_meas
                )
              }

              # Process each valid range: crop to actual measurement bounds
              for (vr in valid_ranges) {
                # Add small buffer (1.5mm = half measurement spacing) for image aesthetics
                crop_start_mm <- vr$start - 1.5
                crop_end_mm <- vr$end + 1.5

                crop_start_px <- max(0, round((crop_start_mm - sect_start) * px_per_mm))
                crop_end_px <- min(img_info$width, round((crop_end_mm - sect_start) * px_per_mm))
                crop_width <- crop_end_px - crop_start_px

                if (crop_width > 10 && vr$n_measurements > 0) {
                  crop_geom <- sprintf("%dx%d+%d+0", crop_width, img_info$height, crop_start_px)
                  cropped <- image_crop(img, crop_geom)

                  # Rotate (width becomes height)
                  rotated <- image_rotate(cropped, 90)

                  # Scale: height proportional to measurement count
                  # This ensures image segment matches graph height for same measurements
                  target_height <- vr$n_measurements * pixels_per_measurement
                  scaled <- image_resize(rotated, sprintf("%dx%d!", target_width, target_height))

                  img_list[[length(img_list) + 1]] <- scaled
                }
              }

              # Add gap between sections (not after last)
              if (i < length(section_order)) {
                if (length(img_list) > 0) {
                  # Gap is proportional: 30mm gap / 3mm per measurement = 10 measurements worth
                  gap_height <- 10 * pixels_per_measurement
                  gap_img <- image_blank(target_width, gap_height, color = "gray30")
                  img_list[[length(img_list) + 1]] <- gap_img
                }
              }
            }
          } else {
            # NORMAL MODE: Full images with gaps
            # CRITICAL: Scale each image height proportional to its depth span
            # This ensures image and graph alignment
            target_width <- 400
            pixels_per_mm <- 1.0  # Base scale factor

            for (i in seq_along(section_order)) {
              sect <- section_order[i]
              if (!(sect %in% names(images))) next

              img <- images[[sect]]$image
              sect_range <- data_ranges %>% filter(section == sect)
              if (nrow(sect_range) == 0) next

              sect_length <- sect_range$data_end - sect_range$data_start

              # Rotate first (width becomes height)
              rotated <- image_rotate(img, 90)

              # Scale: height proportional to section length in mm
              # This ensures image proportions match graph depth proportions
              target_height <- round(sect_length * pixels_per_mm)
              scaled <- image_resize(rotated, sprintf("%dx%d!", target_width, target_height))

              img_list[[length(img_list) + 1]] <- scaled

              # Add gap after each image except the last
              if (i < length(section_order)) {
                # Gap height = 30mm * pixels_per_mm
                gap_height <- round(30 * pixels_per_mm)
                gap_img <- image_blank(target_width, gap_height, color = "gray30")
                img_list[[length(img_list) + 1]] <- gap_img
              }
            }
          }

          if (length(img_list) > 0) {
            stacked <- image_append(do.call(c, img_list), stack = TRUE)
            core_raster <- as.raster(stacked)

            core_plot <- ggplot() +
              annotation_raster(core_raster,
                                xmin = 0, xmax = 1,
                                ymin = depth_min, ymax = depth_max) +
              scale_y_reverse(limits = c(depth_max, depth_min), expand = c(0, 0)) +
              scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
              labs(x = NULL, y = "Depth (cm)", title = "Core") +
              theme_strat() +
              theme(axis.text.x = element_blank(),
                    panel.background = element_rect(fill = "gray30", color = NA))

            # Add exclusion overlay and section boundaries to core panel ONLY
            core_plot <- add_exclusion_overlay_core(core_plot)
            plot_list$core <- add_section_boundaries_core(core_plot)
          }
        }
      })
    }

    # Custom ratio color palette for dynamic assignment
    custom_color_palette <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
                               "#ff7f00", "#a65628", "#f781bf", "#999999")

    # Create proxy panels
    for (i in seq_along(proxies)) {
      proxy <- proxies[i]
      smooth_col <- paste0(proxy, "_smooth")

      if (!proxy %in% names(data)) next

      # Handle colors and labels for custom ratios
      custom <- custom_ratios()
      if (proxy %in% names(custom)) {
        # Custom ratio - assign color from palette and use stored label
        custom_idx <- which(names(custom) == proxy)
        color <- custom_color_palette[((custom_idx - 1) %% length(custom_color_palette)) + 1]
        label <- custom[[proxy]]$label
        is_element <- FALSE
      } else {
        color <- proxy_colors[proxy]
        label <- proxy_labels[proxy]
        is_element <- proxy %in% element_choices
      }

      p <- ggplot(data, aes(y = cumulative_depth_cm))

      if (input$show_points) {
        p <- p + geom_point(aes_string(x = proxy), color = color, alpha = 0.4, size = 1)
      }

      if (smooth_col %in% names(data)) {
        p <- p + geom_path(aes_string(x = smooth_col), color = color, linewidth = 1.2, na.rm = TRUE)
      }


      if (is_element) {
        p <- p + scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()))
      } else {
        p <- p + scale_x_log10()
      }

      p <- p +
        scale_y_reverse(limits = c(depth_max, depth_min), expand = c(0.01, 0)) +
        labs(x = label, y = if (i == 1 && !input$show_core) "Depth (cm)" else NULL,
             title = proxy) +
        theme_strat()

      if (i > 1 || input$show_core) {
        p <- p + theme(axis.text.y = element_blank())
      }

      # Proxy graphs show only valid sediment data (no overlay, no boundaries)
      # Excluded data is already filtered out via xrf_valid
      plot_list[[proxy]] <- p
    }

    # Combine all panels
    if (length(plot_list) == 0) return(NULL)

    # Set widths
    n_panels <- length(plot_list)
    if (input$show_core && "core" %in% names(plot_list)) {
      widths <- c(0.6, rep(1, n_panels - 1))
    } else {
      widths <- rep(1, n_panels)
    }

    combined <- wrap_plots(plot_list, nrow = 1, widths = widths) +
      plot_annotation(
        title = sprintf("%s - Stratigraphic Section Comparison",
                        if (input$site == "TAM") "Tamshiyacu" else "Santa Corina"),
        subtitle = sprintf("%d sections | %d measurements | %s",
                           length(input$sections), nrow(data),
                           paste(proxies, collapse = ", ")),
        caption = "Gray shading = excluded zones (uneven surface/foam)"
      )

    combined
  })

  # Render plot
  output$strat_plot <- renderPlot({
    generate_plot()
  }, height = function() { input$fig_height * 50 })

  # Data summary
  output$data_summary <- renderPrint({
    req(processed_data())
    data <- processed_data()

    cat("=== DATA SUMMARY ===\n\n")
    cat(sprintf("Total measurements: %d\n", nrow(data)))
    cat(sprintf("Sections: %d\n", n_distinct(data$section)))
    cat(sprintf("Depth range: %.1f cm\n\n", max(data$cumulative_depth_cm, na.rm = TRUE)))

    cat("Selected proxies:\n")
    for (proxy in selected_proxies()) {
      if (proxy %in% names(data)) {
        vals <- data[[proxy]]
        cat(sprintf("  %s: mean=%.2f, sd=%.2f, range=[%.2f, %.2f]\n",
                    proxy, mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE),
                    min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)))
      }
    }
  })

  # Section table
  output$section_table <- DT::renderDataTable({
    req(input$sections)

    xrf_valid %>%
      filter(section %in% input$sections) %>%
      group_by(section, group) %>%
      summarise(
        n = n(),
        depth_mm = max(position_mm) - min(position_mm),
        mean_CaTi = round(mean(Ca_Ti, na.rm = TRUE), 2),
        mean_FeMn = round(mean(Fe_Mn, na.rm = TRUE), 1),
        pct_reducing = round(mean(Fe_Mn > 50, na.rm = TRUE) * 100, 1),
        .groups = "drop"
      ) %>%
      arrange(group, section)
  })

  # Download handler for figure
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("pebas_xrf_", input$site, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      p <- generate_plot()
      ggsave(file, p, width = input$fig_width, height = input$fig_height,
             dpi = input$fig_dpi, bg = "white")
    }
  )

  # Download handler for raw data
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("pebas_xrf_data_", input$site, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(processed_data())

      data <- processed_data()

      # Define base columns to always include
      base_cols <- c("section", "group", "position_mm", "cumulative_depth_cm")

      # Get element and ratio columns
      if (input$include_all_elements) {
        # Include all available elements and ratios
        element_cols <- names(element_choices)
        ratio_cols <- c("Ca_Ti", "Fe_Mn", "K_Ti", "Zr_Rb", "Sr_Ca", "Rb_Sr",
                        "Si_Al", "Fe_Ti", "Mn_Ti", "K_Al", "Sr_Ti", "Ba_Sr",
                        "Cu_Zn", "V_Cr")
        # Add any custom ratios
        custom <- custom_ratios()
        ratio_cols <- c(ratio_cols, names(custom))
      } else {
        # Only include selected proxies
        element_cols <- input$elements
        ratio_cols <- input$ratios
        # Add custom ratios
        custom <- custom_ratios()
        ratio_cols <- c(ratio_cols, names(custom))
      }

      # Build list of columns to export
      export_cols <- base_cols

      # Add elements that exist in data
      for (col in element_cols) {
        if (col %in% names(data)) {
          export_cols <- c(export_cols, col)
        }
      }

      # Add ratios that exist in data
      for (col in ratio_cols) {
        if (col %in% names(data)) {
          export_cols <- c(export_cols, col)
        }
      }

      # Add smoothed columns if requested
      if (input$include_smoothed) {
        proxies <- c(element_cols, ratio_cols)
        for (proxy in proxies) {
          smooth_col <- paste0(proxy, "_smooth")
          if (smooth_col %in% names(data)) {
            export_cols <- c(export_cols, smooth_col)
          }
        }
      }

      # Remove duplicates and select columns
      export_cols <- unique(export_cols)
      export_cols <- export_cols[export_cols %in% names(data)]

      export_data <- data %>%
        select(all_of(export_cols)) %>%
        arrange(section, position_mm)

      write_csv(export_data, file)
    }
  )

  # ==============================================================================
  # DENSITY PLOT (Site Comparison Tab)
  # ==============================================================================

  # Track custom ratio for density plot
  density_custom_ratio <- reactiveVal(NULL)
  use_custom_ratio <- reactiveVal(FALSE)

 observeEvent(input$density_use_custom, {
    req(input$density_custom_num, input$density_custom_den)
    req(input$density_custom_num != "", input$density_custom_den != "")
    req(input$density_custom_num != input$density_custom_den)

    ratio_name <- paste0(input$density_custom_num, "_", input$density_custom_den)
    density_custom_ratio(list(
      name = ratio_name,
      num = input$density_custom_num,
      den = input$density_custom_den,
      label = paste0(input$density_custom_num, "/", input$density_custom_den)
    ))
    use_custom_ratio(TRUE)
  })

  # Reset custom when selecting from dropdown
  observeEvent(input$density_proxy, {
    use_custom_ratio(FALSE)
  })

  # Generate density plot
  generate_density_plot <- reactive({
    # Prepare data with site labels
    plot_data <- xrf_valid %>%
      mutate(site = if_else(str_detect(section, "^TAM"), "TAM", "SC"))

    # Determine which proxy to use
    custom <- density_custom_ratio()
    if (use_custom_ratio() && !is.null(custom)) {
      # Use custom ratio
      plot_data <- plot_data %>%
        mutate(!!custom$name := .data[[custom$num]] / .data[[custom$den]])
      proxy_col <- custom$name
      proxy_label <- custom$label
    } else {
      # Use selected proxy from dropdown
      req(input$density_proxy)
      proxy_col <- input$density_proxy
      proxy_label <- proxy_labels[proxy_col]
      if (is.na(proxy_label)) proxy_label <- proxy_col
    }

    req(proxy_col %in% names(plot_data))

    # Remove NA and invalid values
    plot_data <- plot_data %>%
      filter(!is.na(.data[[proxy_col]]), .data[[proxy_col]] > 0)

    # Site colors
    site_colors <- c("TAM" = "#D55E00", "SC" = "#0072B2")

    # Build plot with rug marks
    p <- ggplot(plot_data, aes(x = .data[[proxy_col]], fill = site, color = site)) +
      geom_density(alpha = 0.4, linewidth = 0.8) +
      geom_rug(alpha = 0.3, length = unit(0.02, "npc"), sides = "b") +
      scale_fill_manual(values = site_colors, name = "Site") +
      scale_color_manual(values = site_colors, name = "Site") +
      labs(
        x = proxy_label,
        y = "Density",
        title = paste(proxy_label, "distribution by site"),
        subtitle = sprintf("TAM: n=%d, SC: n=%d",
                           sum(plot_data$site == "TAM"),
                           sum(plot_data$site == "SC"))
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray40")
      )

    # Log scale if requested
    if (input$density_log) {
      p <- p + scale_x_log10()
    }

    # Add threshold line if specified
    if (!is.na(input$density_threshold)) {
      p <- p + geom_vline(xintercept = input$density_threshold,
                          linetype = "dashed", color = "gray40", linewidth = 0.8)
      if (nchar(input$density_threshold_label) > 0) {
        p <- p + annotate("text", x = input$density_threshold, y = Inf,
                          label = input$density_threshold_label,
                          vjust = 2, hjust = -0.1, size = 3.5, color = "gray30")
      }
    }

    p
  })

  # Render density plot
  output$density_plot <- renderPlot({
    generate_density_plot()
  })

  # Density statistics
  output$density_stats <- renderPrint({
    plot_data <- xrf_valid %>%
      mutate(site = if_else(str_detect(section, "^TAM"), "TAM", "SC"))

    # Determine which proxy to use
    custom <- density_custom_ratio()
    if (use_custom_ratio() && !is.null(custom)) {
      plot_data <- plot_data %>%
        mutate(!!custom$name := .data[[custom$num]] / .data[[custom$den]])
      proxy_col <- custom$name
      proxy_label <- custom$label
    } else {
      req(input$density_proxy)
      proxy_col <- input$density_proxy
      proxy_label <- proxy_labels[proxy_col]
      if (is.na(proxy_label)) proxy_label <- proxy_col
    }

    req(proxy_col %in% names(plot_data))

    tam_vals <- plot_data %>% filter(site == "TAM") %>% pull(!!sym(proxy_col)) %>% na.omit()
    sc_vals <- plot_data %>% filter(site == "SC") %>% pull(!!sym(proxy_col)) %>% na.omit()

    cat("=== SITE COMPARISON ===\n\n")
    cat(sprintf("Proxy: %s\n\n", proxy_label))

    cat("TAM (Tamshiyacu):\n")
    cat(sprintf("  n = %d\n", length(tam_vals)))
    cat(sprintf("  median = %.2f\n", median(tam_vals)))
    cat(sprintf("  mean = %.2f\n", mean(tam_vals)))
    cat(sprintf("  range = [%.2f, %.2f]\n\n", min(tam_vals), max(tam_vals)))

    cat("SC (Santa Corina):\n")
    cat(sprintf("  n = %d\n", length(sc_vals)))
    cat(sprintf("  median = %.2f\n", median(sc_vals)))
    cat(sprintf("  mean = %.2f\n", mean(sc_vals)))
    cat(sprintf("  range = [%.2f, %.2f]\n", min(sc_vals), max(sc_vals)))
  })

  # Download density plot
  output$download_density <- downloadHandler(
    filename = function() {
      custom <- density_custom_ratio()
      if (use_custom_ratio() && !is.null(custom)) {
        proxy_name <- custom$name
      } else {
        proxy_name <- input$density_proxy
      }
      paste0("density_", proxy_name, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      p <- generate_density_plot()
      ggsave(file, p, width = input$density_width, height = input$density_height,
             dpi = 300, bg = "white")
    }
  )

  # ==============================================================================
  # TEMPORAL EVOLUTION: Server Logic
  # ==============================================================================

  # Age model parameters
  tam_age_top <- 12.935
  tam_age_bottom <- 13.446
  sc_age_top <- 13.275
  sc_age_bottom <- 14.298
  overlap_start <- 13.275
  overlap_end <- 13.446

  # Proxy labels and colors
  temporal_proxy_labels <- c(
    "Mn_Ti" = "Mn/Ti", "Ca_Ti" = "Ca/Ti", "Fe_Ti" = "Fe/Ti", "Fe_Mn" = "Fe/Mn",
    "K_Ti" = "K/Ti", "Zr_Rb" = "Zr/Rb", "Sr_Ca" = "Sr/Ca"
  )

  temporal_proxy_colors <- c(
    "Mn_Ti" = "#637939", "Ca_Ti" = "#2166ac", "Fe_Ti" = "#b2182b", "Fe_Mn" = "#7b3294",
    "K_Ti" = "#e08214", "Zr_Rb" = "#1a9850", "Sr_Ca" = "#d95f02"
  )

  site_colors <- c("TAM" = "#D55E00", "SC" = "#0072B2")

  # Reset time range button

  observeEvent(input$temporal_reset_range, {
    updateNumericInput(session, "temporal_age_min", value = 12.8)
    updateNumericInput(session, "temporal_age_max", value = 14.5)
  })

  # Zoom to overlap button
  observeEvent(input$temporal_zoom_overlap, {
    updateNumericInput(session, "temporal_age_min", value = 13.2)
    updateNumericInput(session, "temporal_age_max", value = 13.5)
  })

  # Calculate age for each measurement
  temporal_data <- reactive({
    req(input$temporal_proxies)

    data <- xrf_valid %>%
      mutate(site = if_else(str_detect(section, "^TAM"), "TAM", "SC"))

    # Calculate cumulative depth per site
    tam_data <- data %>%
      filter(site == "TAM") %>%
      arrange(strat_order, position_mm) %>%
      mutate(cumulative_depth_cm = cumulative_depth / 10)

    sc_data <- data %>%
      filter(site == "SC") %>%
      arrange(strat_order, position_mm) %>%
      mutate(cumulative_depth_cm = cumulative_depth / 10)

    # Calculate age from depth
    tam_data <- tam_data %>%
      mutate(
        depth_fraction = cumulative_depth_cm / max(cumulative_depth_cm, na.rm = TRUE),
        age_ma = tam_age_top + depth_fraction * (tam_age_bottom - tam_age_top),
        in_overlap = age_ma >= overlap_start & age_ma <= overlap_end
      )

    sc_data <- sc_data %>%
      mutate(
        depth_fraction = cumulative_depth_cm / max(cumulative_depth_cm, na.rm = TRUE),
        age_ma = sc_age_top + depth_fraction * (sc_age_bottom - sc_age_top),
        in_overlap = age_ma >= overlap_start & age_ma <= overlap_end
      )

    # Apply smoothing to all selected proxies
    smooth_window <- input$temporal_smooth
    for (proxy_col in input$temporal_proxies) {
      smooth_col <- paste0(proxy_col, "_smooth")
      if (proxy_col %in% names(tam_data)) {
        tam_data[[smooth_col]] <- zoo::rollmean(tam_data[[proxy_col]], smooth_window, fill = NA, align = "center")
        sc_data[[smooth_col]] <- zoo::rollmean(sc_data[[proxy_col]], smooth_window, fill = NA, align = "center")
      }
    }

    list(tam = tam_data, sc = sc_data)
  })

  # Generate temporal evolution plot
  generate_temporal_plot <- reactive({
    req(input$temporal_proxies)
    req(input$temporal_sites)
    req(length(input$temporal_sites) > 0)

    data <- temporal_data()

    # Filter by selected sites
    plot_data <- tibble()
    if ("TAM" %in% input$temporal_sites) plot_data <- bind_rows(plot_data, data$tam)
    if ("SC" %in% input$temporal_sites) plot_data <- bind_rows(plot_data, data$sc)

    if (nrow(plot_data) == 0) return(NULL)

    # Get time range
    age_min <- input$temporal_age_min
    age_max <- input$temporal_age_max

    # Filter to time range
    plot_data <- plot_data %>%
      filter(age_ma >= age_min, age_ma <= age_max)

    if (nrow(plot_data) == 0) return(NULL)

    # Create faceted plot for multiple proxies
    n_proxies <- length(input$temporal_proxies)

    if (n_proxies == 1) {
      # Single proxy - simple plot
      proxy_col <- input$temporal_proxies[1]
      smooth_col <- paste0(proxy_col, "_smooth")

      p <- ggplot(plot_data, aes(x = age_ma, color = site))

      if (input$temporal_show_overlap) {
        p <- p + annotate("rect", xmin = overlap_start, xmax = overlap_end,
                          ymin = -Inf, ymax = Inf, fill = "yellow", alpha = 0.2)
      }

      if (smooth_col %in% names(plot_data)) {
        p <- p + geom_line(aes(y = .data[[smooth_col]]), linewidth = 1.2, na.rm = TRUE)
      }

      if (input$temporal_show_points) {
        p <- p + geom_point(aes(y = .data[[proxy_col]]), alpha = 0.3, size = 0.8)
      }

      p <- p +
        scale_color_manual(values = site_colors, name = "Site") +
        scale_x_reverse(limits = c(age_max, age_min)) +
        labs(x = "Age (Ma)", y = temporal_proxy_labels[proxy_col],
             title = "Temporal Evolution") +
        theme_bw(base_size = 12) +
        theme(panel.grid.minor = element_blank(),
              legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
              legend.background = element_rect(fill = "white", color = "gray80"))

      if (input$temporal_log_scale) p <- p + scale_y_log10()

    } else {
      # Multiple proxies - faceted plot
      plot_long <- plot_data %>%
        select(age_ma, site, in_overlap, all_of(input$temporal_proxies),
               all_of(paste0(input$temporal_proxies, "_smooth"))) %>%
        pivot_longer(cols = all_of(input$temporal_proxies),
                     names_to = "proxy", values_to = "value") %>%
        mutate(proxy_label = temporal_proxy_labels[proxy])

      # Add smoothed values
      smooth_long <- plot_data %>%
        select(age_ma, site, all_of(paste0(input$temporal_proxies, "_smooth"))) %>%
        pivot_longer(cols = all_of(paste0(input$temporal_proxies, "_smooth")),
                     names_to = "proxy_smooth", values_to = "value_smooth") %>%
        mutate(proxy = str_remove(proxy_smooth, "_smooth"))

      plot_long <- plot_long %>%
        left_join(smooth_long %>% select(age_ma, site, proxy, value_smooth),
                  by = c("age_ma", "site", "proxy"))

      p <- ggplot(plot_long, aes(x = age_ma, color = site))

      if (input$temporal_show_overlap) {
        p <- p + annotate("rect", xmin = overlap_start, xmax = overlap_end,
                          ymin = -Inf, ymax = Inf, fill = "yellow", alpha = 0.2)
      }

      p <- p + geom_line(aes(y = value_smooth), linewidth = 1, na.rm = TRUE)

      if (input$temporal_show_points) {
        p <- p + geom_point(aes(y = value), alpha = 0.2, size = 0.5)
      }

      p <- p +
        facet_wrap(~proxy_label, scales = "free_y", ncol = 1) +
        scale_color_manual(values = site_colors, name = "Site") +
        scale_x_reverse(limits = c(age_max, age_min)) +
        labs(x = "Age (Ma)", y = NULL, title = "Temporal Evolution") +
        theme_bw(base_size = 11) +
        theme(panel.grid.minor = element_blank(),
              legend.position = "top",
              strip.background = element_rect(fill = "gray95"),
              strip.text = element_text(face = "bold"))

      if (input$temporal_log_scale) p <- p + scale_y_log10()
    }

    p
  })

  # Generate overlap comparison plot
  generate_overlap_plot <- reactive({
    req(input$temporal_proxies)
    req(input$temporal_sites)

    data <- temporal_data()

    # Use first proxy for overlap comparison
    proxy_col <- input$temporal_proxies[1]

    # Filter to overlap period only
    tam_overlap <- data$tam %>% filter(in_overlap)
    sc_overlap <- data$sc %>% filter(in_overlap)

    # Check which sites to show
    show_tam <- "TAM" %in% input$temporal_sites && nrow(tam_overlap) > 0
    show_sc <- "SC" %in% input$temporal_sites && nrow(sc_overlap) > 0

    if (!show_tam && !show_sc) return(NULL)

    plot_data <- tibble()
    if (show_tam) plot_data <- bind_rows(plot_data, tam_overlap)
    if (show_sc) plot_data <- bind_rows(plot_data, sc_overlap)

    if (nrow(plot_data) == 0) return(NULL)

    p <- ggplot(plot_data, aes(x = .data[[proxy_col]], fill = site, color = site)) +
      geom_density(alpha = 0.4, linewidth = 0.8) +
      geom_rug(alpha = 0.3, length = unit(0.02, "npc"), sides = "b") +
      scale_fill_manual(values = site_colors, name = "Site") +
      scale_color_manual(values = site_colors, name = "Site") +
      labs(x = temporal_proxy_labels[proxy_col], y = "Density",
           title = sprintf("Overlap Comparison: %s", temporal_proxy_labels[proxy_col]),
           subtitle = sprintf("TAM n=%d, SC n=%d", nrow(tam_overlap), nrow(sc_overlap))) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank(),
            legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
            legend.background = element_rect(fill = "white", color = "gray80"))

    if (input$temporal_log_scale) p <- p + scale_x_log10()

    p
  })

  # Render temporal plot
  output$temporal_plot <- renderPlot({
    p <- generate_temporal_plot()
    if (is.null(p)) {
      plot.new()
      text(0.5, 0.5, "Select at least one site and one proxy", cex = 1.5)
    } else {
      p
    }
  })

  # Render overlap plot
  output$overlap_plot <- renderPlot({
    p <- generate_overlap_plot()
    if (is.null(p)) {
      plot.new()
      text(0.5, 0.5, "No data in overlap period for selected sites", cex = 1.5)
    } else {
      p
    }
  })

  # Temporal statistics
  output$temporal_stats <- renderPrint({
    req(input$temporal_proxies)

    data <- temporal_data()
    proxy_col <- input$temporal_proxies[1]

    tam_overlap <- data$tam %>% filter(in_overlap)
    sc_overlap <- data$sc %>% filter(in_overlap)

    tam_vals <- tam_overlap[[proxy_col]] %>% na.omit()
    sc_vals <- sc_overlap[[proxy_col]] %>% na.omit()

    cat("=== OVERLAP STATISTICS ===\n")
    cat(sprintf("Proxy: %s | Period: 13.275-13.446 Ma (171 kyr)\n\n", temporal_proxy_labels[proxy_col]))

    if ("TAM" %in% input$temporal_sites && length(tam_vals) > 0) {
      cat(sprintf("TAM: n=%d, median=%.3f, mean=%.3f\n",
                  length(tam_vals), median(tam_vals), mean(tam_vals)))
    }

    if ("SC" %in% input$temporal_sites && length(sc_vals) > 0) {
      cat(sprintf("SC:  n=%d, median=%.3f, mean=%.3f\n",
                  length(sc_vals), median(sc_vals), mean(sc_vals)))
    }

    if (length(tam_vals) > 0 && length(sc_vals) > 0 &&
        "TAM" %in% input$temporal_sites && "SC" %in% input$temporal_sites) {
      test_result <- wilcox.test(tam_vals, sc_vals)
      cat(sprintf("\nWilcoxon test: p = %.2e %s\n",
                  test_result$p.value,
                  if (test_result$p.value < 0.05) "(significant)" else "(not significant)"))
    }
  })

  # Download temporal plot
  output$download_temporal <- downloadHandler(
    filename = function() {
      paste0("temporal_evolution_", Sys.Date(), ".png")
    },
    content = function(file) {
      if (input$temporal_view_mode == "cores") {
        p <- generate_temporal_core_plot()
        if (!is.null(p)) {
          ggsave(file, p, width = input$temporal_width, height = input$temporal_height,
                 dpi = 300, bg = "white")
        }
      } else {
        p1 <- generate_temporal_plot()
        p2 <- generate_overlap_plot()
        if (!is.null(p1) && !is.null(p2)) {
          combined <- p1 / p2 + plot_layout(heights = c(1.5, 1))
          ggsave(file, combined, width = input$temporal_width, height = input$temporal_height,
                 dpi = 300, bg = "white")
        } else if (!is.null(p1)) {
          ggsave(file, p1, width = input$temporal_width, height = input$temporal_height,
                 dpi = 300, bg = "white")
        }
      }
    }
  )

  # ==============================================================================
  # TEMPORAL EVOLUTION: Core Comparison View
  # ==============================================================================

  # Section selectors for core comparison
  output$temporal_tam_section_selector <- renderUI({
    tam_sections <- section_config %>%
      filter(site == "TAM") %>%
      arrange(strat_order)

    choices <- setNames(tam_sections$section,
                        paste0(tam_sections$group, " - ", tam_sections$section))

    selectInput("temporal_tam_section", "TAM Section",
                choices = choices,
                selected = choices[1])
  })

  output$temporal_sc_section_selector <- renderUI({
    sc_sections <- section_config %>%
      filter(site == "SC") %>%
      arrange(strat_order)

    choices <- setNames(sc_sections$section,
                        paste0(sc_sections$group, " - ", sc_sections$section))

    selectInput("temporal_sc_section", "SC Section",
                choices = choices,
                selected = choices[1])
  })

  # Generate vertical core plot for a single section
  generate_section_core_plot <- function(section_name, site_label, proxies, brightness = 80, contrast = 25) {

    # Get section info
    sect_info <- section_config %>% filter(section == section_name)
    if (nrow(sect_info) == 0) return(NULL)

    # Get section data
    sect_data <- xrf_valid %>%
      filter(section == section_name) %>%
      arrange(position_mm)

    if (nrow(sect_data) == 0) return(NULL)

    # Get data range
    data_start <- min(sect_data$position_mm)
    data_end <- max(sect_data$position_mm)
    depth_cm <- (data_end - data_start) / 10

    # Apply smoothing
    smooth_window <- input$temporal_smooth
    for (proxy_col in proxies) {
      smooth_col <- paste0(proxy_col, "_smooth")
      if (proxy_col %in% names(sect_data)) {
        sect_data[[smooth_col]] <- zoo::rollmean(sect_data[[proxy_col]], smooth_window, fill = NA, align = "center")
      }
    }

    # Theme for stratigraphic plots
    theme_strat <- function(base_size = 10) {
      theme_minimal(base_size = base_size) +
        theme(
          plot.title = element_text(size = base_size + 1, face = "bold"),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size - 1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
          legend.position = "none",
          plot.margin = margin(5, 5, 5, 5)
        )
    }

    plot_list <- list()

    # Core image panel
    img_result <- tryCatch({
      process_core_image(section_name, sect_info$optical_path, data_start, data_end,
                         brightness = brightness, contrast = contrast, gamma = 2.0)
    }, error = function(e) NULL)

    if (!is.null(img_result)) {
      # Rotate image 90 degrees for vertical display
      img_rotated <- image_rotate(img_result$image, 90)
      target_width <- 300
      target_height <- round(depth_cm * 8)  # Scale factor for display
      img_scaled <- image_resize(img_rotated, sprintf("%dx%d!", target_width, target_height))
      core_raster <- as.raster(img_scaled)

      core_plot <- ggplot() +
        annotation_raster(core_raster, xmin = 0, xmax = 1, ymin = 0, ymax = depth_cm) +
        scale_y_reverse(limits = c(depth_cm, 0), expand = c(0, 0)) +
        scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        labs(x = NULL, y = "Depth (cm)", title = site_label) +
        theme_strat() +
        theme(axis.text.x = element_blank(),
              panel.background = element_rect(fill = "gray30", color = NA))

      plot_list$core <- core_plot
    }

    # Proxy panels
    for (proxy in proxies) {
      if (!proxy %in% names(sect_data)) next

      smooth_col <- paste0(proxy, "_smooth")
      color <- temporal_proxy_colors[proxy]
      label <- temporal_proxy_labels[proxy]

      # Convert position to cm from start
      sect_data_plot <- sect_data %>%
        mutate(depth_cm = (position_mm - data_start) / 10)

      p <- ggplot(sect_data_plot, aes(y = depth_cm))

      if (input$temporal_show_points) {
        p <- p + geom_point(aes(x = .data[[proxy]]), color = color, alpha = 0.4, size = 1)
      }

      if (smooth_col %in% names(sect_data_plot)) {
        p <- p + geom_path(aes(x = .data[[smooth_col]]), color = color, linewidth = 1, na.rm = TRUE)
      }

      if (input$temporal_log_scale) {
        p <- p + scale_x_log10()
      }

      p <- p +
        scale_y_reverse(limits = c(depth_cm, 0), expand = c(0.01, 0)) +
        labs(x = label, y = NULL, title = proxy) +
        theme_strat() +
        theme(axis.text.y = element_blank())

      plot_list[[proxy]] <- p
    }

    if (length(plot_list) == 0) return(NULL)

    # Combine plots
    n_panels <- length(plot_list)
    if ("core" %in% names(plot_list)) {
      widths <- c(0.5, rep(1, n_panels - 1))
    } else {
      widths <- rep(1, n_panels)
    }

    wrap_plots(plot_list, nrow = 1, widths = widths)
  }

  # Generate side-by-side core comparison plot
  generate_temporal_core_plot <- reactive({
    req(input$temporal_proxies)

    proxies <- input$temporal_proxies

    plot_list <- list()

    # Generate TAM plot if selected
    if ("TAM" %in% input$temporal_sites && !is.null(input$temporal_tam_section)) {
      tam_plot <- generate_section_core_plot(
        input$temporal_tam_section,
        paste0("TAM: ", input$temporal_tam_section),
        proxies,
        brightness = input$temporal_brightness,
        contrast = input$temporal_contrast
      )
      if (!is.null(tam_plot)) {
        plot_list$tam <- tam_plot
      }
    }

    # Generate SC plot if selected
    if ("SC" %in% input$temporal_sites && !is.null(input$temporal_sc_section)) {
      sc_plot <- generate_section_core_plot(
        input$temporal_sc_section,
        paste0("SC: ", input$temporal_sc_section),
        proxies,
        brightness = input$temporal_brightness,
        contrast = input$temporal_contrast
      )
      if (!is.null(sc_plot)) {
        plot_list$sc <- sc_plot
      }
    }

    if (length(plot_list) == 0) return(NULL)

    # Stack TAM and SC vertically
    if (length(plot_list) == 2) {
      plot_list$tam / plot_list$sc +
        plot_annotation(title = "Core Comparison: TAM vs SC",
                        subtitle = paste("Proxies:", paste(temporal_proxy_labels[proxies], collapse = ", ")))
    } else if (length(plot_list) == 1) {
      plot_list[[1]] +
        plot_annotation(title = "Core Section",
                        subtitle = paste("Proxies:", paste(temporal_proxy_labels[proxies], collapse = ", ")))
    }
  })

  # Render core comparison plot
  output$temporal_core_plot <- renderPlot({
    p <- generate_temporal_core_plot()
    if (is.null(p)) {
      plot.new()
      text(0.5, 0.5, "Select at least one site and one proxy", cex = 1.5)
    } else {
      p
    }
  })

  # Core comparison statistics
  output$temporal_core_stats <- renderPrint({
    req(input$temporal_proxies)

    proxy_col <- input$temporal_proxies[1]
    cat("=== SECTION COMPARISON ===\n")
    cat(sprintf("Primary proxy: %s\n\n", temporal_proxy_labels[proxy_col]))

    if ("TAM" %in% input$temporal_sites && !is.null(input$temporal_tam_section)) {
      tam_data <- xrf_valid %>%
        filter(section == input$temporal_tam_section)

      if (nrow(tam_data) > 0 && proxy_col %in% names(tam_data)) {
        vals <- tam_data[[proxy_col]] %>% na.omit()
        cat(sprintf("TAM (%s):\n  n=%d, median=%.3f, mean=%.3f\n\n",
                    input$temporal_tam_section, length(vals), median(vals), mean(vals)))
      }
    }

    if ("SC" %in% input$temporal_sites && !is.null(input$temporal_sc_section)) {
      sc_data <- xrf_valid %>%
        filter(section == input$temporal_sc_section)

      if (nrow(sc_data) > 0 && proxy_col %in% names(sc_data)) {
        vals <- sc_data[[proxy_col]] %>% na.omit()
        cat(sprintf("SC (%s):\n  n=%d, median=%.3f, mean=%.3f\n",
                    input$temporal_sc_section, length(vals), median(vals), mean(vals)))
      }
    }
  })

  # ==============================================================================
  # EXCLUSION EDITOR: Server Logic
  # ==============================================================================

  # Get data for the selected section in editor
  editor_section_data <- reactive({
    req(input$editor_section)

    sect <- input$editor_section

    # Get XRF data for this section (all data, not just valid)
    sect_data <- xrf_data %>%
      filter(section == sect) %>%
      arrange(position_mm)

    if (nrow(sect_data) == 0) return(NULL)

    sect_data
  })

  # Get exclusion zones for the selected section
  editor_section_exclusions <- reactive({
    req(input$editor_section)

    working_exclusions() %>%
      filter(section == input$editor_section, !is.na(exclude_start_mm))
  })

  # Display brush selection info
  output$brush_info <- renderText({
    brush <- input$editor_brush
    if (is.null(brush)) {
      return("No selection. Drag vertically on the plot.")
    }

    # Convert brush y values (in mm) to selection range
    start_mm <- round(min(brush$ymin, brush$ymax), 1)
    end_mm <- round(max(brush$ymin, brush$ymax), 1)

    sprintf("Selected range:\n  Start: %.1f mm\n  End: %.1f mm\n  Length: %.1f mm",
            start_mm, end_mm, end_mm - start_mm)
  })

  # Clear brush selection
  observeEvent(input$clear_brush, {
    session$resetBrush("editor_brush")
  })

  # Add new exclusion zone from brush selection
  observeEvent(input$add_zone, {
    brush <- input$editor_brush
    req(brush)
    req(input$editor_section)

    start_mm <- round(min(brush$ymin, brush$ymax), 0)
    end_mm <- round(max(brush$ymin, brush$ymax), 0)

    # Don't add if selection is too small
    if (end_mm - start_mm < 5) {
      showNotification("Selection too small (< 5mm). Please select a larger region.",
                       type = "warning")
      return()
    }

    # Add new zone to working exclusions
    new_zone <- tibble(
      section = input$editor_section,
      exclude_start_mm = start_mm,
      exclude_end_mm = end_mm,
      notes = "foam"
    )

    current <- working_exclusions()
    updated <- bind_rows(current, new_zone) %>%
      arrange(section, exclude_start_mm)

    working_exclusions(updated)
    unsaved_changes(TRUE)

    # Clear the brush
    session$resetBrush("editor_brush")

    showNotification(sprintf("Added exclusion zone: %d - %d mm", start_mm, end_mm),
                     type = "message")
  })

  # Display exclusion table for current section
  output$exclusion_table <- DT::renderDataTable({
    excl <- editor_section_exclusions()

    if (is.null(excl) || nrow(excl) == 0) {
      return(DT::datatable(
        tibble(Message = "No exclusion zones for this section"),
        options = list(dom = 't'),
        rownames = FALSE
      ))
    }

    excl %>%
      select(Start = exclude_start_mm, End = exclude_end_mm, Notes = notes) %>%
      DT::datatable(
        options = list(dom = 't', pageLength = 20),
        rownames = FALSE,
        selection = 'multiple'
      )
  })

  # Delete selected exclusion zones
  observeEvent(input$delete_zone, {
    req(input$editor_section)
    selected_rows <- input$exclusion_table_rows_selected

    if (is.null(selected_rows) || length(selected_rows) == 0) {
      showNotification("No rows selected. Click on table rows to select them.",
                       type = "warning")
      return()
    }

    sect_excl <- editor_section_exclusions()
    if (nrow(sect_excl) == 0) return()

    # Get the zones to delete
    zones_to_delete <- sect_excl[selected_rows, ]

    # Remove from working exclusions
    current <- working_exclusions()
    updated <- current %>%
      anti_join(zones_to_delete, by = c("section", "exclude_start_mm", "exclude_end_mm"))

    working_exclusions(updated)
    unsaved_changes(TRUE)

    showNotification(sprintf("Deleted %d exclusion zone(s)", length(selected_rows)),
                     type = "message")
  })

  # Save exclusions to CSV
  observeEvent(input$save_exclusions, {
    csv_path <- file.path(base_path, "data", "exclusion_zones.csv")

    tryCatch({
      write_csv(working_exclusions(), csv_path)

      # Update the global exclusion_zones variable used by the main plots
      exclusion_zones <<- working_exclusions() %>%
        filter(!is.na(exclude_start_mm))

      unsaved_changes(FALSE)

      showNotification("Exclusion zones saved successfully!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error saving:", e$message), type = "error")
    })
  })

  # Save status display
  output$save_status <- renderText({
    if (unsaved_changes()) {
      "⚠ UNSAVED CHANGES"
    } else {
      "✓ All changes saved"
    }
  })

  # Generate the editor plot with core image and XRF data
  output$editor_plot <- renderPlot({
    req(input$editor_section)

    sect <- input$editor_section
    sect_data <- editor_section_data()
    req(sect_data)

    sect_excl <- editor_section_exclusions()

    # Get section info
    sect_info <- section_config %>% filter(section == sect)
    if (nrow(sect_info) == 0) return(NULL)

    # Get data range
    data_start <- min(sect_data$position_mm)
    data_end <- max(sect_data$position_mm)

    # Try to load and process core image
    img_result <- tryCatch({
      process_core_image(sect, sect_info$optical_path, data_start, data_end,
                         brightness = 80, contrast = 25, gamma = 2.0)
    }, error = function(e) NULL)

    # Create base plot
    p <- ggplot()

    # Add core image as background if available
    if (!is.null(img_result)) {
      # Rotate image 90 degrees clockwise for vertical display
      img_rotated <- image_rotate(img_result$image, 90)

      # Convert to raster array for ggplot
      img_array <- as.raster(img_rotated)

      # Add image as annotation raster (image spans x = 0 to 1, y = data range)
      p <- p + annotation_raster(img_array,
                                  xmin = 0, xmax = 1,
                                  ymin = data_end, ymax = data_start)
    }

    # Add exclusion zone shading (semi-transparent overlay on image)
    if (nrow(sect_excl) > 0) {
      for (i in seq_len(nrow(sect_excl))) {
        p <- p + annotate("rect",
                          xmin = -Inf, xmax = Inf,
                          ymin = sect_excl$exclude_start_mm[i],
                          ymax = sect_excl$exclude_end_mm[i],
                          fill = "red", alpha = 0.35)
      }
    }

    # Add XRF data points (Fe normalized) on right side of image
    if ("Fe" %in% names(sect_data)) {
      fe_max <- max(sect_data$Fe, na.rm = TRUE)
      sect_data <- sect_data %>%
        mutate(Fe_norm = 1.0 + (Fe / fe_max) * 0.4)  # Scale to 1.0-1.4 range (right of image)

      # Color points by QC status
      p <- p +
        geom_point(data = sect_data %>% filter(qc_pass, !excluded),
                   aes(x = Fe_norm, y = position_mm),
                   color = "darkgreen", alpha = 0.7, size = 1.5) +
        geom_point(data = sect_data %>% filter(!qc_pass | excluded),
                   aes(x = Fe_norm, y = position_mm),
                   color = "red", alpha = 0.5, size = 1.5)
    }

    # Add vertical line separating image from XRF data
    p <- p + geom_vline(xintercept = 1.0, color = "gray40", linewidth = 0.5)

    # Add labels and theme
    p <- p +
      scale_y_reverse(limits = c(data_end + 5, data_start - 5),
                      expand = c(0, 0)) +
      scale_x_continuous(limits = c(0, 1.5),
                         breaks = c(0.5, 1.2, 1.4),
                         labels = c("Core Image", "Fe", "")) +
      labs(
        title = sprintf("Section: %s", sect),
        subtitle = sprintf("Drag vertically on plot to select foam regions | %.0f - %.0f mm",
                           data_start, data_end),
        x = NULL,
        y = "Position (mm)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray80", linewidth = 0.3),
        plot.title = element_text(face = "bold"),
        axis.title = element_text(size = 11),
        panel.background = element_rect(fill = "gray95", color = NA)
      )

    # Add annotation for exclusion zones
    if (nrow(sect_excl) > 0) {
      p <- p +
        labs(caption = sprintf("Red shading = %d existing exclusion zone(s) | Green dots = valid XRF | Red dots = excluded",
                               nrow(sect_excl)))
    } else {
      p <- p +
        labs(caption = "Green dots = valid XRF data | Red dots = excluded/QC fail")
    }

    p
  }, res = 96)
}

# ==============================================================================
# RUN APP
# ==============================================================================

shinyApp(ui = ui, server = server)
