# ==============================================================================
# Pebas-XRF: Master Analysis Script
# ==============================================================================
# Run all analysis scripts in sequence
# TAM = Tamshiyacu, SC = Santa Corina (Pebas Formation, Peru)
# ==============================================================================

# Set working directory to project root
if (!require(here)) install.packages("here")
setwd(here::here())

message("=" |> rep(70) |> paste(collapse = ""))
message("PEBAS-XRF ANALYSIS PIPELINE")
message("=" |> rep(70) |> paste(collapse = ""))

# Check/install required packages
required_packages <- c(
  "tidyverse",   # Data manipulation and plotting
 "fs",          # File system operations
  "here",        # Project paths
  "patchwork",   # Combining plots
  "compositions" # Log-ratio transformations
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    message(sprintf("Installing package: %s", pkg))
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

message("\nAll required packages loaded.\n")

# Run scripts in order
scripts <- c(
  "01_data_import.R",
  "02_element_ratios.R",
  "03_stratigraphic_plots.R",
  "04_pca_analysis.R"
)

for (script in scripts) {
  script_path <- file.path("scripts", script)
  message("\n", "=" |> rep(70) |> paste(collapse = ""))
  message(sprintf("RUNNING: %s", script))
  message("=" |> rep(70) |> paste(collapse = ""), "\n")

  tryCatch({
    source(script_path)
    message(sprintf("\n[SUCCESS] %s completed", script))
  }, error = function(e) {
    message(sprintf("\n[ERROR] %s failed: %s", script, e$message))
  })
}

message("\n", "=" |> rep(70) |> paste(collapse = ""))
message("PIPELINE COMPLETE")
message("=" |> rep(70) |> paste(collapse = ""))
message("\nOutput files are in: output/")
message("  - tables/: CSV data files")
message("  - figures/: PNG plot files")
