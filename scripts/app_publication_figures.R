# ==============================================================================
# Pebas-XRF: Publication Figure Generator
# ==============================================================================
# Interactive Shiny app for overlaying XRF data trends on core images
# Export publication-quality figures with customizable styling
# ==============================================================================

library(shiny)
library(tidyverse)
library(magick)
library(zoo)
library(grid)
library(gridExtra)
library(png)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

base_path <- here::here()
output_path <- file.path(base_path, "output")
fig_path <- file.path(output_path, "figures", "case_study")
export_path <- file.path(output_path, "figures", "publication")
dir.create(export_path, showWarnings = FALSE, recursive = TRUE)

# Load XRF data
xrf_data <- read_csv(file.path(output_path, "tables", "xrf_data_stacked.csv"),
                     show_col_types = FALSE)

# Available proxies
proxy_info <- list(
  Ca_Ti = list(name = "Ca/Ti (Carbonate)", color = "#1B7837", threshold = c(2, 5, 10)),
  Fe_Mn = list(name = "Fe/Mn (Redox)", color = "#762A83", threshold = 50),
  K_Ti = list(name = "K/Ti (Weathering)", color = "#E66101", threshold = NULL),
  Zr_Rb = list(name = "Zr/Rb (Grain Size)", color = "#B2182B", threshold = NULL),
  Fe = list(name = "Fe (Terrigenous)", color = "#8C510A", threshold = NULL),
  Ca = list(name = "Ca (Carbonate)", color = "#01665E", threshold = NULL)
)

# ==============================================================================
# UI
# ==============================================================================

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .well { background-color: #f8f9fa; }
      .sidebar { background-color: #f0f0f0; padding: 15px; border-radius: 5px; }
      h4 { color: #333; margin-top: 20px; }
      .export-btn { margin-top: 10px; width: 100%; }
    "))
  ),

  titlePanel("Pebas-XRF Publication Figure Generator"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      class = "sidebar",

      h4("Core Selection"),
      selectInput("group", "Select Core Group:",
                  choices = sort(unique(xrf_data$group)),
                  selected = "GROUP3"),

      hr(),
      h4("Proxy Display"),

      checkboxGroupInput("proxies", "Show Proxies:",
                         choices = setNames(names(proxy_info),
                                            sapply(proxy_info, function(x) x$name)),
                         selected = c("Ca_Ti", "Fe_Mn")),

      sliderInput("window_size", "Filter Window Size:",
                  min = 1, max = 15, value = 5, step = 2),

      hr(),
      h4("Image Settings"),

      sliderInput("brightness", "Image Brightness:",
                  min = 0, max = 100, value = 40),

      sliderInput("core_width", "Core Image Width:",
                  min = 50, max = 300, value = 150),

      checkboxInput("show_raw", "Show Raw Data Points", value = FALSE),

      hr(),
      h4("Overlay Settings"),

      sliderInput("overlay_alpha", "Curve Transparency:",
                  min = 0.3, max = 1.0, value = 0.8, step = 0.1),

      sliderInput("line_width", "Line Width:",
                  min = 0.5, max = 3.0, value = 1.2, step = 0.1),

      checkboxInput("show_thresholds", "Show Threshold Lines", value = TRUE),

      checkboxInput("show_facies", "Show Facies Column", value = TRUE),

      hr(),
      h4("Export Options"),

      numericInput("export_width", "Width (inches):", value = 10, min = 4, max = 20),
      numericInput("export_height", "Height (inches):", value = 12, min = 4, max = 24),
      numericInput("export_dpi", "Resolution (DPI):", value = 300, min = 150, max = 600),

      textInput("export_filename", "Filename:", value = "publication_figure"),

      actionButton("export_png", "Export PNG", class = "btn-primary export-btn"),
      actionButton("export_pdf", "Export PDF", class = "btn-success export-btn"),

      hr(),
      verbatimTextOutput("status")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Preview",
                 plotOutput("preview_plot", height = "800px")),
        tabPanel("Statistics",
                 verbatimTextOutput("stats_output"),
                 tableOutput("facies_table")),
        tabPanel("Help",
                 includeMarkdown_or_text())
      )
    )
  )
)

# Helper for help text
includeMarkdown_or_text <- function() {
  HTML("
    <h3>Publication Figure Generator - Help</h3>

    <h4>Core Selection</h4>
    <p>Select from 7 core groups: GROUP1-3 (Tamshiyacu) and GROUP4-7 (Santa Corina).</p>

    <h4>Proxy Display</h4>
    <ul>
      <li><strong>Ca/Ti</strong>: Carbonate vs terrigenous balance. Higher = more biogenic carbonate.</li>
      <li><strong>Fe/Mn</strong>: Redox indicator. Values >50 indicate reducing conditions.</li>
      <li><strong>K/Ti</strong>: Chemical weathering intensity proxy.</li>
      <li><strong>Zr/Rb</strong>: Grain size/depositional energy proxy.</li>
    </ul>

    <h4>Filter Window Size</h4>
    <p>Moving window filter to reduce noise. Larger = smoother curves. At 3mm step, window of 5 = 15mm effective resolution.</p>

    <h4>Image Settings</h4>
    <ul>
      <li><strong>Brightness</strong>: Adjust core image brightness (0-100%).</li>
      <li><strong>Core Width</strong>: Width of core image strip in pixels.</li>
      <li><strong>Raw Data</strong>: Show individual measurement points.</li>
    </ul>

    <h4>Thresholds</h4>
    <ul>
      <li><strong>Ca/Ti</strong>: 2, 5, 10 (Clastic/Mixed/Carbonate/Shell-rich)</li>
      <li><strong>Fe/Mn</strong>: 50 (Oxic/Reducing boundary)</li>
    </ul>

    <h4>Export</h4>
    <p>PNG for presentations/web, PDF for publication. Standard journal size: 7-10 inches wide, 300 DPI.</p>
  ")
}

# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {

  # Reactive: Get processed group data
  group_data <- reactive({
    req(input$group)

    data <- xrf_data %>%
      filter(group == input$group) %>%
      arrange(position_mm) %>%
      mutate(depth_cm = position_mm / 10)

    # Apply moving window filter
    if (input$window_size > 1) {
      for (proxy in names(proxy_info)) {
        if (proxy %in% names(data)) {
          data[[paste0(proxy, "_filt")]] <- zoo::rollmean(
            data[[proxy]], input$window_size, fill = NA, align = "center"
          )
        }
      }
    } else {
      for (proxy in names(proxy_info)) {
        if (proxy %in% names(data)) {
          data[[paste0(proxy, "_filt")]] <- data[[proxy]]
        }
      }
    }

    # Add facies classification
    data <- data %>%
      mutate(
        facies = case_when(
          Ca_Ti > 10 ~ "Shell-rich",
          Ca_Ti > 5 ~ "Carbonate",
          Ca_Ti > 2 ~ "Mixed",
          TRUE ~ "Clastic"
        ),
        facies = factor(facies, levels = c("Shell-rich", "Carbonate", "Mixed", "Clastic"))
      )

    data
  })

  # Reactive: Load and process core image
  core_image <- reactive({
    req(input$group)

    img_path <- file.path(fig_path, sprintf("core_optical_%s.png", input$group))

    if (!file.exists(img_path)) {
      return(NULL)
    }

    img <- image_read(img_path)

    # Apply brightness
    if (input$brightness != 40) {
      brightness_adj <- (input$brightness - 40) + 100
      img <- image_modulate(img, brightness = brightness_adj)
    }

    # Resize to target width
    img <- image_resize(img, sprintf("%dx", input$core_width))

    img
  })

  # Main plot generation
  generate_plot <- reactive({
    req(group_data())

    data <- group_data()
    depth_min <- min(data$position_mm)
    depth_max <- max(data$position_mm)

    facies_colors <- c("Shell-rich" = "#2166AC", "Carbonate" = "#67A9CF",
                       "Mixed" = "#D1E5F0", "Clastic" = "#B2182B")

    plot_list <- list()
    widths <- c()

    # Core image panel
    img <- core_image()
    if (!is.null(img)) {
      core_raster <- as.raster(img)

      plot_list$core <- ggplot() +
        annotation_raster(core_raster,
                          xmin = 0, xmax = 1,
                          ymin = depth_min/10, ymax = depth_max/10) +
        scale_y_reverse(limits = c(depth_max/10, depth_min/10)) +
        scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        labs(x = NULL, y = "Depth (cm)") +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid = element_blank(),
              plot.margin = margin(5, 2, 5, 5))

      widths <- c(widths, 0.5)
    }

    # Facies panel
    if (input$show_facies) {
      plot_list$facies <- ggplot(data, aes(y = depth_cm)) +
        geom_tile(aes(x = 0.5, fill = facies), width = 1, height = 0.3) +
        scale_fill_manual(values = facies_colors, name = "Facies") +
        scale_y_reverse() +
        labs(x = NULL, y = if (is.null(img)) "Depth (cm)" else NULL) +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid = element_blank(),
              legend.position = "none",
              axis.text.y = if (!is.null(img)) element_blank() else element_text(),
              plot.margin = margin(5, 2, 5, 2))

      widths <- c(widths, 0.25)
    }

    # Proxy panels
    for (proxy in input$proxies) {
      if (!(proxy %in% names(data))) next

      info <- proxy_info[[proxy]]
      filt_col <- paste0(proxy, "_filt")

      p <- ggplot(data, aes(y = depth_cm))

      # Raw data points
      if (input$show_raw) {
        p <- p + geom_point(aes_string(x = proxy),
                            alpha = 0.2, size = 0.8, color = info$color)
      }

      # Filtered curve
      p <- p + geom_path(aes_string(x = filt_col),
                         color = info$color,
                         linewidth = input$line_width,
                         alpha = input$overlay_alpha,
                         na.rm = TRUE)

      # Threshold lines
      if (input$show_thresholds && !is.null(info$threshold)) {
        for (thresh in info$threshold) {
          p <- p + geom_vline(xintercept = thresh,
                              linetype = "dashed",
                              color = "gray50",
                              linewidth = 0.4)
        }
      }

      # Special handling for Fe/Mn redox
      if (proxy == "Fe_Mn" && input$show_thresholds) {
        p <- p + geom_vline(xintercept = 50,
                            linetype = "dashed",
                            color = "red",
                            linewidth = 0.5)
      }

      p <- p +
        scale_y_reverse() +
        labs(x = info$name, y = NULL) +
        theme_minimal(base_size = 11) +
        theme(axis.text.y = element_blank(),
              plot.margin = margin(5, 5, 5, 2))

      plot_list[[proxy]] <- p
      widths <- c(widths, 1)
    }

    if (length(plot_list) == 0) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select at least one proxy") +
               theme_void())
    }

    # Combine plots
    combined <- wrap_plots(plot_list, nrow = 1, widths = widths)

    # Get site info
    site <- if (grepl("^GROUP[1-3]$", input$group)) "TAM" else "SC"
    site_name <- if (site == "TAM") "Tamshiyacu" else "Santa Corina"

    combined <- combined +
      plot_annotation(
        title = sprintf("%s %s: Geochemical Stratigraphy", site, input$group),
        subtitle = sprintf("%s, Pebas Formation | n=%d | %.1f cm depth",
                           site_name, nrow(data), (depth_max - depth_min)/10),
        theme = theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 11, color = "gray40")
        )
      )

    combined
  })

  # Preview output
  output$preview_plot <- renderPlot({
    generate_plot()
  }, res = 96)

  # Statistics output
  output$stats_output <- renderPrint({
    data <- group_data()

    cat("=== GROUP STATISTICS ===\n\n")
    cat(sprintf("Group: %s\n", input$group))
    cat(sprintf("Measurements: %d\n", nrow(data)))
    cat(sprintf("Depth range: %.1f - %.1f cm\n",
                min(data$depth_cm), max(data$depth_cm)))

    cat("\n--- Proxy Summary ---\n")
    for (proxy in input$proxies) {
      if (proxy %in% names(data)) {
        vals <- data[[proxy]]
        cat(sprintf("%s: mean=%.2f, sd=%.2f, range=[%.2f, %.2f]\n",
                    proxy,
                    mean(vals, na.rm = TRUE),
                    sd(vals, na.rm = TRUE),
                    min(vals, na.rm = TRUE),
                    max(vals, na.rm = TRUE)))
      }
    }

    cat(sprintf("\nReducing conditions (Fe/Mn > 50): %.1f%%\n",
                mean(data$Fe_Mn > 50, na.rm = TRUE) * 100))
  })

  # Facies table
  output$facies_table <- renderTable({
    data <- group_data()
    data %>%
      count(facies) %>%
      mutate(percent = n / sum(n) * 100) %>%
      arrange(desc(n))
  })

  # Export handlers
  observeEvent(input$export_png, {
    p <- generate_plot()
    filename <- paste0(input$export_filename, "_", input$group, ".png")
    filepath <- file.path(export_path, filename)

    ggsave(filepath, p,
           width = input$export_width,
           height = input$export_height,
           dpi = input$export_dpi,
           bg = "white")

    output$status <- renderText({
      sprintf("Exported: %s\nPath: %s", filename, filepath)
    })
  })

  observeEvent(input$export_pdf, {
    p <- generate_plot()
    filename <- paste0(input$export_filename, "_", input$group, ".pdf")
    filepath <- file.path(export_path, filename)

    ggsave(filepath, p,
           width = input$export_width,
           height = input$export_height,
           bg = "white")

    output$status <- renderText({
      sprintf("Exported: %s\nPath: %s", filename, filepath)
    })
  })

  # Status output
  output$status <- renderText({
    sprintf("Ready. Group: %s | Proxies: %s",
            input$group, paste(input$proxies, collapse = ", "))
  })
}

# ==============================================================================
# RUN APP
# ==============================================================================

shinyApp(ui = ui, server = server)
