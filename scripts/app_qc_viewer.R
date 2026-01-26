# ==============================================================================
# Pebas-XRF: Interactive QC Viewer
# ==============================================================================
# Shiny app for aligning XRF spectra with optical images
# and marking exclusion zones for data filtering
# ==============================================================================
#
# USAGE:
#   1. Open R/RStudio in the Pebas-XRF directory
#   2. Run: shiny::runApp("scripts/app_qc_viewer.R")
#   3. Select a section from the dropdown
#   4. Click on the plot to mark exclusion start/end points
#   5. Click "Save Exclusions" when done
#
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
OPTICAL_PATH <- file.path(BASE_PATH, "output", "figures", "optical_review")
OUTPUT_PATH <- file.path(BASE_PATH, "output")
EXCLUSION_FILE <- file.path(BASE_PATH, "data", "exclusion_zones.csv")

# ==============================================================================
# LOAD DATA
# ==============================================================================

# Load XRF data
xrf_data <- tryCatch({
  read_csv(file.path(OUTPUT_PATH, "tables", "xrf_data_qc.csv"), show_col_types = FALSE)
}, error = function(e) {
  # If processed data doesn't exist, return empty
  tibble()
})

# Get list of sections
sections <- sort(unique(xrf_data$section))

# Load existing exclusions
load_exclusions <- function() {
  if (file.exists(EXCLUSION_FILE)) {
    read_csv(EXCLUSION_FILE, show_col_types = FALSE)
  } else {
    tibble(
      section = character(),
      exclude_start_mm = numeric(),
      exclude_end_mm = numeric(),
      notes = character()
    )
  }
}

# ==============================================================================
# UI
# ==============================================================================

ui <- fluidPage(

  titlePanel("Pebas-XRF: QC Viewer & Exclusion Tool"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      # Section selector
      selectInput("section", "Select Section:", choices = sections),

      hr(),

      # Element selector for profile
      selectInput("element", "Element to Display:",
                  choices = c("Fe", "Ti", "Ca", "K", "Si", "Al", "cps"),
                  selected = "Fe"),

      # Log scale toggle
      checkboxInput("log_scale", "Log Scale", value = TRUE),

      hr(),

      h4("Exclusion Zones"),

      # Current exclusions display
      verbatimTextOutput("current_exclusions"),

      # Manual entry
      fluidRow(
        column(6, numericInput("excl_start", "Start (mm):", value = NA, step = 1)),
        column(6, numericInput("excl_end", "End (mm):", value = NA, step = 1))
      ),

      textInput("excl_notes", "Notes:", placeholder = "e.g., foam fill"),

      fluidRow(
        column(6, actionButton("add_exclusion", "Add Zone", class = "btn-warning")),
        column(6, actionButton("clear_exclusion", "Clear", class = "btn-secondary"))
      ),

      hr(),

      # Save button
      actionButton("save_exclusions", "Save All Exclusions", class = "btn-success",
                   style = "width: 100%;"),

      hr(),

      # Instructions
      h5("Instructions:"),
      tags$ol(
        tags$li("Select a section from dropdown"),
        tags$li("Review optical image and XRF profile"),
        tags$li("Enter start/end positions for exclusion zones"),
        tags$li("Add notes describing why (foam, gap, etc.)"),
        tags$li("Click 'Add Zone' to add exclusion"),
        tags$li("Click 'Save All Exclusions' when done")
      )
    ),

    mainPanel(
      width = 9,

      # Optical image
      h4(textOutput("section_title")),

      plotOutput("optical_plot", height = "200px"),

      hr(),

      # XRF profile with exclusions marked
      plotOutput("xrf_plot", height = "400px", click = "plot_click"),

      hr(),

      # Data summary
      fluidRow(
        column(4, verbatimTextOutput("data_summary")),
        column(8, tableOutput("exclusion_table"))
      )
    )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {

  # Reactive: current exclusions
  exclusions <- reactiveVal(load_exclusions())

  # Reactive: section data
  section_data <- reactive({
    req(input$section)
    xrf_data %>% filter(section == input$section)
  })

  # Reactive: optical image path
  optical_file <- reactive({
    req(input$section)
    file.path(OPTICAL_PATH, paste0("optical_", input$section, ".png"))
  })

  # Output: Section title
  output$section_title <- renderText({
    req(input$section)
    data <- section_data()
    if (nrow(data) > 0) {
      sprintf("%s | Position: %.0f - %.0f mm | %d measurements",
              input$section,
              min(data$position_mm, na.rm = TRUE),
              max(data$position_mm, na.rm = TRUE),
              nrow(data))
    } else {
      input$section
    }
  })

  # Output: Optical image plot
  output$optical_plot <- renderPlot({
    req(optical_file())

    if (file.exists(optical_file())) {
      img <- readPNG(optical_file())
      grid.raster(img)
    } else {
      plot.new()
      text(0.5, 0.5, "Optical image not found", cex = 1.5)
    }
  })

  # Output: XRF profile plot with exclusions
  output$xrf_plot <- renderPlot({
    req(input$section, input$element)

    data <- section_data()
    if (nrow(data) == 0) {
      plot.new()
      text(0.5, 0.5, "No data for this section", cex = 1.5)
      return()
    }

    # Get exclusions for this section
    sect_excl <- exclusions() %>%
      filter(section == input$section, !is.na(exclude_start_mm))

    # Create plot
    p <- ggplot(data, aes(x = position_mm, y = .data[[input$element]])) +
      geom_line(color = "steelblue", linewidth = 0.8) +
      geom_point(color = "steelblue", size = 1, alpha = 0.5) +
      labs(
        x = "Position (mm)",
        y = paste0(input$element, " (cps)"),
        title = paste0(input$element, " Profile - ", input$section)
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold")
      )

    # Add log scale if selected
    if (input$log_scale && input$element != "cps") {
      p <- p + scale_y_log10()
    }

    # Add exclusion zones as red shaded areas
    if (nrow(sect_excl) > 0) {
      for (i in 1:nrow(sect_excl)) {
        p <- p + annotate(
          "rect",
          xmin = sect_excl$exclude_start_mm[i],
          xmax = sect_excl$exclude_end_mm[i],
          ymin = -Inf, ymax = Inf,
          fill = "red", alpha = 0.3
        ) +
          annotate(
            "text",
            x = (sect_excl$exclude_start_mm[i] + sect_excl$exclude_end_mm[i]) / 2,
            y = Inf,
            label = sect_excl$notes[i],
            vjust = 1.5, color = "red", size = 3
          )
      }
    }

    # Add pending exclusion if values entered
    if (!is.na(input$excl_start) && !is.na(input$excl_end)) {
      p <- p + annotate(
        "rect",
        xmin = input$excl_start,
        xmax = input$excl_end,
        ymin = -Inf, ymax = Inf,
        fill = "orange", alpha = 0.3
      )
    }

    print(p)
  })

  # Output: Current exclusions for this section
  output$current_exclusions <- renderText({
    sect_excl <- exclusions() %>%
      filter(section == input$section, !is.na(exclude_start_mm))

    if (nrow(sect_excl) == 0) {
      "No exclusions defined"
    } else {
      paste(
        sprintf("%.0f-%.0f mm: %s",
                sect_excl$exclude_start_mm,
                sect_excl$exclude_end_mm,
                sect_excl$notes),
        collapse = "\n"
      )
    }
  })

  # Output: Data summary
  output$data_summary <- renderText({
    data <- section_data()
    if (nrow(data) == 0) return("No data")

    sprintf(
      "Section: %s\nMeasurements: %d\nPosition range: %.0f - %.0f mm\nQC pass: %d (%.1f%%)",
      input$section,
      nrow(data),
      min(data$position_mm, na.rm = TRUE),
      max(data$position_mm, na.rm = TRUE),
      sum(data$qc_pass, na.rm = TRUE),
      100 * mean(data$qc_pass, na.rm = TRUE)
    )
  })

  # Output: Exclusion table
  output$exclusion_table <- renderTable({
    exclusions() %>%
      filter(section == input$section, !is.na(exclude_start_mm)) %>%
      select(Start = exclude_start_mm, End = exclude_end_mm, Notes = notes)
  })

  # Handle plot clicks to set exclusion points
  observeEvent(input$plot_click, {
    click_x <- round(input$plot_click$x)

    if (is.na(input$excl_start)) {
      updateNumericInput(session, "excl_start", value = click_x)
    } else if (is.na(input$excl_end)) {
      updateNumericInput(session, "excl_end", value = click_x)
    } else {
      # Both filled, start new
      updateNumericInput(session, "excl_start", value = click_x)
      updateNumericInput(session, "excl_end", value = NA)
    }
  })

  # Add exclusion zone
  observeEvent(input$add_exclusion, {
    req(input$excl_start, input$excl_end)

    # Ensure start < end
    start_val <- min(input$excl_start, input$excl_end)
    end_val <- max(input$excl_start, input$excl_end)

    new_excl <- tibble(
      section = input$section,
      exclude_start_mm = start_val,
      exclude_end_mm = end_val,
      notes = ifelse(input$excl_notes == "", "manual exclusion", input$excl_notes)
    )

    # Add to exclusions (remove any existing for same section/range)
    current <- exclusions() %>%
      filter(!(section == input$section &
               exclude_start_mm == start_val &
               exclude_end_mm == end_val))

    exclusions(bind_rows(current, new_excl))

    # Clear inputs
    updateNumericInput(session, "excl_start", value = NA)
    updateNumericInput(session, "excl_end", value = NA)
    updateTextInput(session, "excl_notes", value = "")

    showNotification(
      sprintf("Added exclusion: %s %.0f-%.0f mm", input$section, start_val, end_val),
      type = "message"
    )
  })

  # Clear exclusion inputs
  observeEvent(input$clear_exclusion, {
    updateNumericInput(session, "excl_start", value = NA)
    updateNumericInput(session, "excl_end", value = NA)
    updateTextInput(session, "excl_notes", value = "")
  })

  # Save all exclusions
  observeEvent(input$save_exclusions, {
    # Ensure all sections are represented
    all_sections <- tibble(section = sections)

    final_excl <- all_sections %>%
      left_join(
        exclusions() %>% filter(!is.na(exclude_start_mm)),
        by = "section"
      )

    write_csv(final_excl, EXCLUSION_FILE)

    showNotification(
      sprintf("Saved %d exclusion zones to %s",
              sum(!is.na(final_excl$exclude_start_mm)),
              basename(EXCLUSION_FILE)),
      type = "message",
      duration = 5
    )
  })
}

# ==============================================================================
# RUN APP
# ==============================================================================

shinyApp(ui = ui, server = server)
