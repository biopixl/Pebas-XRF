# ==============================================================================
# Pebas-XRF: Interactive QC Viewer (v2 - Aligned)
# ==============================================================================
# Shiny app for aligning XRF spectra with optical images
# Properly handles consecutive sections within the same core
# ==============================================================================

library(shiny)
library(tidyverse)
library(png)
library(grid)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

BASE_PATH <- here::here()
DATA_PATH <- file.path(BASE_PATH, "TAM-SC-IsaacA")
OUTPUT_PATH <- file.path(BASE_PATH, "output")
EXCLUSION_FILE <- file.path(BASE_PATH, "data", "exclusion_zones.csv")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Get optical image info for a section
get_optical_info <- function(section_name) {
  # Find the section folder
  section_dirs <- list.dirs(DATA_PATH, recursive = TRUE)
  section_dir <- section_dirs[basename(section_dirs) == section_name &
                               !grepl("/copied/", section_dirs)]

  if (length(section_dir) == 0) return(NULL)
  section_dir <- section_dir[1]

  # Read document.txt for positions
  doc_file <- file.path(section_dir, "document.txt")
  optical_file <- file.path(section_dir, "optical.tif")

  if (!file.exists(doc_file)) return(NULL)

  doc <- readLines(doc_file, warn = FALSE, encoding = "latin1")

  # Parse XRF scan coordinates
  xrf_line <- grep("^Start coordinate", doc, value = TRUE)
  if (length(xrf_line) > 0) {
    parts <- strsplit(xrf_line, "\t")[[1]]
    xrf_start <- as.numeric(parts[2])
    xrf_stop <- as.numeric(parts[4])
  } else {
    xrf_start <- NA
    xrf_stop <- NA
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
    section = section_name,
    optical_file = optical_file,
    optical_start = optical_start,
    optical_end = optical_end,
    xrf_start = xrf_start,
    xrf_stop = xrf_stop
  )
}

#' Get optical PNG for display (from pre-generated cache)
get_optical_png <- function(section_name) {
  info <- get_optical_info(section_name)
  if (is.null(info)) return(NULL)

  # Check for pre-generated PNG in optical_review directory
  png_file <- file.path(OUTPUT_PATH, "figures", "optical_review",
                        paste0("optical_", section_name, ".png"))

  if (file.exists(png_file)) {
    return(list(
      file = png_file,
      optical_start = info$optical_start,
      optical_end = info$optical_end,
      xrf_start = info$xrf_start,
      xrf_stop = info$xrf_stop
    ))
  }

  NULL
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

# Load XRF data
xrf_data <- tryCatch({
  read_csv(file.path(OUTPUT_PATH, "tables", "xrf_data_qc.csv"), show_col_types = FALSE)
}, error = function(e) tibble())

sections <- sort(unique(xrf_data$section))

# Load existing exclusions
load_exclusions <- function() {
  if (file.exists(EXCLUSION_FILE)) {
    read_csv(EXCLUSION_FILE, show_col_types = FALSE)
  } else {
    tibble(section = character(), exclude_start_mm = numeric(),
           exclude_end_mm = numeric(), notes = character())
  }
}

# ==============================================================================
# UI
# ==============================================================================

ui <- fluidPage(

  titlePanel("Pebas-XRF: QC Viewer - Optical/XRF Alignment"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      selectInput("section", "Select Section:", choices = sections),

      hr(),

      selectInput("element", "Element:",
                  choices = c("Fe", "Ti", "Ca", "K", "Si", "Mn", "cps"),
                  selected = "Fe"),

      checkboxInput("log_scale", "Log Scale", value = TRUE),

      hr(),

      h4("Mark Exclusion Zone"),
      p("Click on the aligned plot to set start/end positions"),

      fluidRow(
        column(6, numericInput("excl_start", "Start (mm):", value = NA, step = 1)),
        column(6, numericInput("excl_end", "End (mm):", value = NA, step = 1))
      ),

      textInput("excl_notes", "Notes:", placeholder = "foam, gap, etc."),

      fluidRow(
        column(6, actionButton("add_exclusion", "Add Zone", class = "btn-warning")),
        column(6, actionButton("clear_inputs", "Clear", class = "btn-secondary"))
      ),

      hr(),

      h4("Current Exclusions"),
      verbatimTextOutput("exclusion_list"),

      actionButton("remove_last", "Remove Last", class = "btn-danger btn-sm"),

      hr(),

      actionButton("save_all", "💾 Save All Exclusions",
                   class = "btn-success", style = "width:100%"),

      hr(),

      h5("Position Info"),
      verbatimTextOutput("position_info")
    ),

    mainPanel(
      width = 9,

      h4(textOutput("section_header")),

      # Combined aligned plot
      plotOutput("aligned_plot", height = "600px", click = "plot_click",
                 brush = brushOpts(id = "plot_brush", direction = "x")),

      hr(),

      fluidRow(
        column(6,
          h5("Data Summary"),
          verbatimTextOutput("data_summary")
        ),
        column(6,
          h5("Exclusion Zones for This Section"),
          tableOutput("exclusion_table")
        )
      )
    )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {

  # Reactive values
  exclusions <- reactiveVal(load_exclusions())

  # Get section data
  section_data <- reactive({
    req(input$section)
    xrf_data %>% filter(section == input$section)
  })

  # Get optical info
  optical_info <- reactive({
    req(input$section)
    get_optical_png(input$section)
  })

  # Section header
  output$section_header <- renderText({
    data <- section_data()
    info <- optical_info()

    if (nrow(data) > 0 && !is.null(info)) {
      sprintf("%s | XRF: %.0f-%.0f mm | Optical: %.0f-%.0f mm | N=%d",
              input$section, info$xrf_start, info$xrf_stop,
              info$optical_start, info$optical_end, nrow(data))
    } else {
      input$section
    }
  })

  # Position info
  output$position_info <- renderText({
    info <- optical_info()
    if (is.null(info)) return("No position info")

    sprintf("Optical range: %.0f - %.0f mm\nXRF scan range: %.0f - %.0f mm\nScale: %.2f mm/pixel",
            info$optical_start, info$optical_end,
            info$xrf_start, info$xrf_stop,
            (info$optical_end - info$optical_start) / 2000)  # Approximate
  })

  # Main aligned plot
  output$aligned_plot <- renderPlot({
    req(input$section, input$element)

    data <- section_data()
    info <- optical_info()
    sect_excl <- exclusions() %>% filter(section == input$section, !is.na(exclude_start_mm))

    if (nrow(data) == 0) {
      plot.new()
      text(0.5, 0.5, "No XRF data for this section", cex = 2)
      return()
    }

    # Set up plot layout: optical on top, XRF below, shared x-axis
    par(mfrow = c(2, 1), mar = c(1, 4, 2, 2), oma = c(3, 0, 0, 0))

    # Determine x-axis range (use XRF data range)
    x_min <- min(data$position_mm, na.rm = TRUE)
    x_max <- max(data$position_mm, na.rm = TRUE)
    x_range <- x_max - x_min
    x_min <- x_min - x_range * 0.02
    x_max <- x_max + x_range * 0.02

    # === PANEL 1: Optical Image ===
    if (!is.null(info) && file.exists(info$file)) {
      img <- tryCatch(readPNG(info$file), error = function(e) NULL)

      if (!is.null(img)) {
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

          # Plot cropped image
          plot(0, 0, type = "n", xlim = c(x_min, x_max), ylim = c(0, 1),
               xlab = "", ylab = "", xaxt = "n", yaxt = "n",
               main = "Optical Image (aligned)")

          rasterImage(img_crop, x_min, 0, x_max, 1)

          # Add exclusion zones as semi-transparent overlays
          if (nrow(sect_excl) > 0) {
            for (i in 1:nrow(sect_excl)) {
              rect(sect_excl$exclude_start_mm[i], 0,
                   sect_excl$exclude_end_mm[i], 1,
                   col = rgb(1, 0, 0, 0.4), border = "red", lwd = 2)
            }
          }

          # Add pending exclusion
          if (!is.na(input$excl_start) && !is.na(input$excl_end)) {
            rect(input$excl_start, 0, input$excl_end, 1,
                 col = rgb(1, 0.5, 0, 0.4), border = "orange", lwd = 2)
          }

          # Add position scale
          axis(1, at = pretty(c(x_min, x_max)), labels = FALSE, tck = -0.05)
        }
      }
    } else {
      plot.new()
      text(0.5, 0.5, "Optical image not available", cex = 1.5)
    }

    # === PANEL 2: XRF Profile ===
    y_vals <- data[[input$element]]

    if (input$log_scale && all(y_vals > 0, na.rm = TRUE)) {
      y_vals <- log10(y_vals)
      y_lab <- paste0("log10(", input$element, ")")
    } else {
      y_lab <- paste0(input$element, " (cps)")
    }

    plot(data$position_mm, y_vals, type = "l", col = "steelblue", lwd = 1.5,
         xlim = c(x_min, x_max), xlab = "", ylab = y_lab,
         main = paste(input$element, "Profile"))
    points(data$position_mm, y_vals, pch = 16, cex = 0.5, col = "steelblue")

    # Add exclusion zones
    if (nrow(sect_excl) > 0) {
      for (i in 1:nrow(sect_excl)) {
        rect(sect_excl$exclude_start_mm[i], par("usr")[3],
             sect_excl$exclude_end_mm[i], par("usr")[4],
             col = rgb(1, 0, 0, 0.3), border = NA)
        abline(v = c(sect_excl$exclude_start_mm[i], sect_excl$exclude_end_mm[i]),
               col = "red", lty = 2, lwd = 1.5)
        text(mean(c(sect_excl$exclude_start_mm[i], sect_excl$exclude_end_mm[i])),
             par("usr")[4], sect_excl$notes[i], pos = 1, col = "red", cex = 0.8)
      }
    }

    # Add pending exclusion
    if (!is.na(input$excl_start) && !is.na(input$excl_end)) {
      rect(input$excl_start, par("usr")[3], input$excl_end, par("usr")[4],
           col = rgb(1, 0.5, 0, 0.3), border = "orange", lwd = 2)
    }

    # Single click marker
    if (!is.na(input$excl_start) && is.na(input$excl_end)) {
      abline(v = input$excl_start, col = "orange", lwd = 2)
    }

    # Shared x-axis label
    mtext("Position (mm)", side = 1, outer = TRUE, line = 1.5)

    # Add grid
    grid(nx = NA, ny = NULL, col = "gray80", lty = 1)
  })

  # Handle plot clicks
  observeEvent(input$plot_click, {
    click_x <- round(input$plot_click$x)

    if (is.na(input$excl_start)) {
      updateNumericInput(session, "excl_start", value = click_x)
    } else if (is.na(input$excl_end)) {
      updateNumericInput(session, "excl_end", value = click_x)
    } else {
      updateNumericInput(session, "excl_start", value = click_x)
      updateNumericInput(session, "excl_end", value = NA)
    }
  })

  # Handle brush selection
  observeEvent(input$plot_brush, {
    brush <- input$plot_brush
    if (!is.null(brush)) {
      updateNumericInput(session, "excl_start", value = round(brush$xmin))
      updateNumericInput(session, "excl_end", value = round(brush$xmax))
    }
  })

  # Add exclusion
  observeEvent(input$add_exclusion, {
    req(input$excl_start, input$excl_end)

    start_val <- min(input$excl_start, input$excl_end)
    end_val <- max(input$excl_start, input$excl_end)

    new_excl <- tibble(
      section = input$section,
      exclude_start_mm = start_val,
      exclude_end_mm = end_val,
      notes = ifelse(input$excl_notes == "", "excluded", input$excl_notes)
    )

    exclusions(bind_rows(exclusions(), new_excl))

    updateNumericInput(session, "excl_start", value = NA)
    updateNumericInput(session, "excl_end", value = NA)
    updateTextInput(session, "excl_notes", value = "")

    showNotification(sprintf("Added: %.0f-%.0f mm", start_val, end_val), type = "message")
  })

  # Clear inputs
  observeEvent(input$clear_inputs, {
    updateNumericInput(session, "excl_start", value = NA)
    updateNumericInput(session, "excl_end", value = NA)
    updateTextInput(session, "excl_notes", value = "")
  })

  # Remove last exclusion for this section
  observeEvent(input$remove_last, {
    current <- exclusions()
    sect_rows <- which(current$section == input$section & !is.na(current$exclude_start_mm))

    if (length(sect_rows) > 0) {
      last_row <- max(sect_rows)
      exclusions(current[-last_row, ])
      showNotification("Removed last exclusion", type = "warning")
    }
  })

  # Save all exclusions
  observeEvent(input$save_all, {
    all_sections <- tibble(section = sections)

    final <- all_sections %>%
      left_join(exclusions() %>% filter(!is.na(exclude_start_mm)), by = "section")

    write_csv(final, EXCLUSION_FILE)

    n_zones <- sum(!is.na(final$exclude_start_mm))
    showNotification(sprintf("Saved %d exclusion zones", n_zones),
                     type = "message", duration = 5)
  })

  # Exclusion list display
  output$exclusion_list <- renderText({
    sect_excl <- exclusions() %>%
      filter(section == input$section, !is.na(exclude_start_mm))

    if (nrow(sect_excl) == 0) return("None defined")

    paste(sprintf("• %.0f-%.0f: %s",
                  sect_excl$exclude_start_mm,
                  sect_excl$exclude_end_mm,
                  sect_excl$notes), collapse = "\n")
  })

  # Exclusion table
  output$exclusion_table <- renderTable({
    exclusions() %>%
      filter(section == input$section, !is.na(exclude_start_mm)) %>%
      select(Start = exclude_start_mm, End = exclude_end_mm, Notes = notes)
  })

  # Data summary
  output$data_summary <- renderText({
    data <- section_data()
    if (nrow(data) == 0) return("No data")

    sprintf("Measurements: %d\nPosition: %.0f - %.0f mm\nQC pass: %d (%.1f%%)",
            nrow(data),
            min(data$position_mm, na.rm = TRUE),
            max(data$position_mm, na.rm = TRUE),
            sum(data$qc_pass, na.rm = TRUE),
            100 * mean(data$qc_pass, na.rm = TRUE))
  })
}

# ==============================================================================
# RUN
# ==============================================================================

shinyApp(ui = ui, server = server)
