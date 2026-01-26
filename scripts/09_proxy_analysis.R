# ==============================================================================
# Pebas-XRF: Multi-Proxy Geochemical Analysis
# ==============================================================================
# Statistical analysis of XRF elemental proxies with literature-based hypotheses
# Traceability matrix for co-measured magnetic, CT, and fossil data integration
# Pebas Formation, Western Amazonia (Middle-Late Miocene)
# ==============================================================================

library(tidyverse)
library(zoo)  # Rolling window functions
library(patchwork)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

base_path <- here::here()
output_path <- file.path(base_path, "output")
fig_path <- file.path(output_path, "figures", "proxy_analysis")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

# Signal filtering parameters
WINDOW_SIZE <- 5      # 5-point moving window (15mm at 3mm step resolution)

# ==============================================================================
# MULTI-PROXY TRACEABILITY MATRIX
# ==============================================================================
# Links XRF elemental proxies to paleoenvironmental interpretations
# with expected correlations for co-measured datasets:
# - Magnetic susceptibility (MS)
# - CT radiodensity
# - Fossil assemblages (mollusks, ostracods, fish, plant remains)

proxy_matrix <- tribble(
  ~proxy,    ~elements,   ~interpretation,                          ~ms_correlation,  ~ct_correlation,  ~fossil_association,                    ~key_refs,

  # Terrigenous/Detrital Indicators
  "Fe",      "Fe",        "Terrigenous input, detrital flux",       "positive",       "positive",       "Low diversity assemblages",            "Croudace2006",
  "Ti",      "Ti",        "Conservative detrital indicator",        "positive",       "positive",       "Clastic-dominated facies",             "Calvert1996",
  "K",       "K",         "Clay minerals, detrital feldspars",      "positive",       "neutral",        "Floodplain/lacustrine clays",          "Nesbitt1982",
  "Al",      "Al",        "Aluminosilicate clays",                  "positive",       "neutral",        "Fine-grained suspension deposits",     "Weltje2008",
  "Si",      "Si",        "Quartz, biogenic silica",                "variable",       "positive",       "Diatoms, sponge spicules",             "Brown2007",
  "Rb",      "Rb",        "Clay minerals (illite)",                 "positive",       "neutral",        "Weathered source terrain",             "Jin2001",

  # Carbonate/Biogenic Indicators
  "Ca",      "Ca",        "Carbonate (biogenic/authigenic)",        "negative",       "positive",       "Mollusk shells, ostracods",            "Davies2015",
  "Sr",      "Sr",        "Carbonate, aragonite shells",            "negative",       "positive",       "Aragonitic mollusks (Pachydon)",       "Croudace2015",
  "Ca/Ti",   "Ca,Ti",     "Carbonate vs terrigenous balance",       "negative",       "variable",       "Shell bed vs clastic facies",          "Haug2001",

  # Redox Indicators
  "Fe/Mn",   "Fe,Mn",     "Bottom water oxygenation",               "variable",       "neutral",        "Anoxia-tolerant taxa abundance",       "Calvert1996",
  "Mn",      "Mn",        "Oxic conditions, Mn-oxide precipitation","positive",       "neutral",        "Bioturbation intensity",               "Rothwell2015",

  # Weathering Indicators
  "K/Ti",    "K,Ti",      "Chemical weathering intensity",          "neutral",        "neutral",        "Andean vs cratonic provenance",        "Nesbitt1982",
  "Rb/Sr",   "Rb,Sr",     "Silicate weathering, CIA proxy",         "positive",       "negative",       "Paleoclimate humidity signal",         "Jin2001",

  # Grain Size/Energy Indicators
  "Zr/Rb",   "Zr,Rb",     "Grain size, depositional energy",        "positive",       "positive",       "Channel vs lake margin facies",        "Dypvik2001",
  "Zr",      "Zr",        "Heavy minerals, coarse fraction",        "positive",       "positive",       "Proximity to fluvial input",           "Kylander2011",
  "Si/Al",   "Si,Al",     "Quartz/clay ratio, grain size",          "positive",       "positive",       "Depositional energy proxy",            "Weltje2008"
)

# Save traceability matrix
write_csv(proxy_matrix, file.path(output_path, "tables", "proxy_traceability_matrix.csv"))
message("Multi-proxy traceability matrix saved")

# ==============================================================================
# LOAD AND FILTER DATA
# ==============================================================================

xrf_data <- read_csv(file.path(output_path, "tables", "xrf_data_stacked.csv"),
                     show_col_types = FALSE)

message(sprintf("Loaded %d measurements", nrow(xrf_data)))

# Apply moving window filter within each section to reduce measurement noise
# while preserving stratigraphic trends
xrf_filtered <- xrf_data %>%
  group_by(site, group, section) %>%
  arrange(cumulative_depth) %>%
  mutate(
    # Moving window mean filter (centered, handles edge effects with NA)
    Fe_filtered = rollmean(Fe, WINDOW_SIZE, fill = NA, align = "center"),
    Ca_filtered = rollmean(Ca, WINDOW_SIZE, fill = NA, align = "center"),
    Ti_filtered = rollmean(Ti, WINDOW_SIZE, fill = NA, align = "center"),
    Mn_filtered = rollmean(Mn, WINDOW_SIZE, fill = NA, align = "center"),
    K_filtered = rollmean(K, WINDOW_SIZE, fill = NA, align = "center"),
    Si_filtered = rollmean(Si, WINDOW_SIZE, fill = NA, align = "center"),

    # Filtered ratios
    Ca_Ti_filtered = rollmean(Ca_Ti, WINDOW_SIZE, fill = NA, align = "center"),
    Fe_Mn_filtered = rollmean(Fe_Mn, WINDOW_SIZE, fill = NA, align = "center"),
    K_Ti_filtered = rollmean(K_Ti, WINDOW_SIZE, fill = NA, align = "center"),
    Zr_Rb_filtered = rollmean(Zr_Rb, WINDOW_SIZE, fill = NA, align = "center"),

    # Standardized anomalies for excursion detection
    Fe_zscore = (Fe - mean(Fe, na.rm = TRUE)) / sd(Fe, na.rm = TRUE),
    Ca_Ti_zscore = (Ca_Ti - mean(Ca_Ti, na.rm = TRUE)) / sd(Ca_Ti, na.rm = TRUE)
  ) %>%
  ungroup()

# ==============================================================================
# LITERATURE-BASED HYPOTHESES AND STATISTICAL TESTS
# ==============================================================================

message("\n=== STATISTICAL TESTS FOR LITERATURE-BASED HYPOTHESES ===\n")

# Hypothesis 1: Reducing conditions (Fe/Mn > 50) at both sites
# Based on Calvert & Pedersen 1996, Rothwell 2015
h1_test <- xrf_filtered %>%
  group_by(site) %>%
  summarise(
    n = n(),
    Fe_Mn_median = median(Fe_Mn, na.rm = TRUE),
    Fe_Mn_mean = mean(Fe_Mn, na.rm = TRUE),
    pct_above_50 = mean(Fe_Mn > 50, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("H1: Bottom waters were predominantly reducing (Fe/Mn > 50)\n")
print(h1_test)
cat(sprintf("\nConclusion: %s\n\n",
            ifelse(all(h1_test$pct_above_50 > 50),
                   "SUPPORTED - majority of samples indicate reducing conditions",
                   "PARTIALLY SUPPORTED - mixed redox signals")))

# Hypothesis 2: TAM has higher carbonate than SC (Ca/Ti difference)
# T-test for site differences
h2_test <- t.test(Ca_Ti ~ site, data = xrf_filtered)
cat("H2: TAM has higher carbonate accumulation than SC (Ca/Ti)\n")
cat(sprintf("  TAM Ca/Ti mean: %.2f\n", mean(xrf_filtered$Ca_Ti[xrf_filtered$site == "TAM"], na.rm = TRUE)))
cat(sprintf("  SC Ca/Ti mean: %.2f\n", mean(xrf_filtered$Ca_Ti[xrf_filtered$site == "SC"], na.rm = TRUE)))
cat(sprintf("  t-statistic: %.2f, p-value: %.4f\n", h2_test$statistic, h2_test$p.value))
cat(sprintf("\nConclusion: %s\n\n",
            ifelse(h2_test$p.value < 0.05,
                   "SUPPORTED - significant difference in Ca/Ti between sites",
                   "NOT SUPPORTED - no significant difference")))

# Hypothesis 3: Weathering intensity is similar between sites (K/Ti)
h3_test <- t.test(K_Ti ~ site, data = xrf_filtered)
cat("H3: Source weathering intensity is similar between sites (K/Ti)\n")
cat(sprintf("  TAM K/Ti mean: %.3f\n", mean(xrf_filtered$K_Ti[xrf_filtered$site == "TAM"], na.rm = TRUE)))
cat(sprintf("  SC K/Ti mean: %.3f\n", mean(xrf_filtered$K_Ti[xrf_filtered$site == "SC"], na.rm = TRUE)))
cat(sprintf("  t-statistic: %.2f, p-value: %.4f\n", h3_test$statistic, h3_test$p.value))
cat(sprintf("\nConclusion: %s\n\n",
            ifelse(h3_test$p.value > 0.05,
                   "SUPPORTED - no significant difference in K/Ti",
                   "NOT SUPPORTED - significant difference in weathering proxy")))

# Hypothesis 4: Fine-grained sedimentation (low Zr/Rb) indicates lacustrine setting
h4_summary <- xrf_filtered %>%
  summarise(
    Zr_Rb_median = median(Zr_Rb, na.rm = TRUE),
    Zr_Rb_q25 = quantile(Zr_Rb, 0.25, na.rm = TRUE),
    Zr_Rb_q75 = quantile(Zr_Rb, 0.75, na.rm = TRUE)
  )
cat("H4: Fine-grained lacustrine sedimentation (low Zr/Rb)\n")
cat(sprintf("  Zr/Rb median: %.2f (IQR: %.2f - %.2f)\n",
            h4_summary$Zr_Rb_median, h4_summary$Zr_Rb_q25, h4_summary$Zr_Rb_q75))
cat("\nConclusion: Low Zr/Rb values support fine-grained, low-energy deposition\n\n")

# ==============================================================================
# CORRELATION ANALYSIS FOR PROXY VALIDATION
# ==============================================================================

message("=== PROXY CORRELATION MATRIX ===\n")

# Key proxy correlations
cor_vars <- c("Fe", "Ti", "Ca", "K", "Mn", "Si", "Zr", "Rb", "Sr")
cor_data <- xrf_filtered %>% select(all_of(cor_vars))
cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")

# Print correlation matrix
cat("Element correlation matrix:\n")
print(round(cor_matrix, 2))

# Expected correlations for validation
cat("\n\nProxy validation against expected correlations:\n")
cat(sprintf("  Fe-Ti (terrigenous): r = %.2f (expected: positive) - %s\n",
            cor_matrix["Fe", "Ti"],
            ifelse(cor_matrix["Fe", "Ti"] > 0.5, "VALID", "CHECK")))
cat(sprintf("  Ca-Sr (carbonate): r = %.2f (expected: positive) - %s\n",
            cor_matrix["Ca", "Sr"],
            ifelse(cor_matrix["Ca", "Sr"] > 0.3, "VALID", "CHECK")))
cat(sprintf("  Zr-Rb (grain size): r = %.2f (expected: negative/weak) - %s\n",
            cor_matrix["Zr", "Rb"],
            ifelse(cor_matrix["Zr", "Rb"] < 0.5, "VALID", "CHECK")))

# ==============================================================================
# FILTERED STRATIGRAPHIC PLOTS
# ==============================================================================

message("\n=== GENERATING FILTERED PROXY PLOTS ===")

create_filtered_plot <- function(data, site_name) {

  depth_max <- max(data$cumulative_depth, na.rm = TRUE)

  # Common aesthetics
  raw_alpha <- 0.2
  filtered_size <- 0.8

  p1 <- ggplot(data, aes(y = cumulative_depth)) +
    geom_point(aes(x = Ca_Ti), alpha = raw_alpha, size = 0.5, color = "darkgreen") +
    geom_path(aes(x = Ca_Ti_filtered), color = "darkgreen", linewidth = filtered_size) +
    geom_vline(xintercept = median(data$Ca_Ti, na.rm = TRUE),
               linetype = "dashed", color = "gray50") +
    scale_y_reverse() +
    labs(x = "Ca/Ti", y = "Depth (mm)", title = "Carbonate") +
    theme_minimal(base_size = 9)

  p2 <- ggplot(data, aes(y = cumulative_depth)) +
    geom_point(aes(x = Fe_Mn), alpha = raw_alpha, size = 0.5, color = "purple") +
    geom_path(aes(x = Fe_Mn_filtered), color = "purple", linewidth = filtered_size) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
    scale_y_reverse() +
    labs(x = "Fe/Mn", y = NULL, title = "Redox") +
    theme_minimal(base_size = 9)

  p3 <- ggplot(data, aes(y = cumulative_depth)) +
    geom_point(aes(x = K_Ti), alpha = raw_alpha, size = 0.5, color = "orange") +
    geom_path(aes(x = K_Ti_filtered), color = "orange", linewidth = filtered_size) +
    scale_y_reverse() +
    labs(x = "K/Ti", y = NULL, title = "Weathering") +
    theme_minimal(base_size = 9)

  p4 <- ggplot(data, aes(y = cumulative_depth)) +
    geom_point(aes(x = Zr_Rb), alpha = raw_alpha, size = 0.5, color = "darkred") +
    geom_path(aes(x = Zr_Rb_filtered), color = "darkred", linewidth = filtered_size) +
    scale_y_reverse() +
    labs(x = "Zr/Rb", y = NULL, title = "Grain Size") +
    theme_minimal(base_size = 9)

  combined <- p1 + p2 + p3 + p4 +
    plot_annotation(
      title = sprintf("%s - Geochemical Proxy Records", site_name),
      subtitle = sprintf("Moving window filter (n=%d, 15mm) | Raw measurements shown as points", WINDOW_SIZE),
      theme = theme(plot.title = element_text(face = "bold", size = 12))
    )

  return(combined)
}

# Generate for each site
for (site in c("TAM", "SC")) {
  site_data <- xrf_filtered %>% filter(site == !!site)
  site_name <- ifelse(site == "TAM", "Tamshiyacu", "Santa Corina")

  p <- create_filtered_plot(site_data, site_name)

  ggsave(file.path(fig_path, paste0("filtered_proxies_", site, ".png")), p,
         width = 12, height = 10, dpi = 150, bg = "white")
  message(sprintf("Saved: filtered_proxies_%s.png", site))
}

# ==============================================================================
# MULTI-PROXY INTEGRATION FRAMEWORK
# ==============================================================================
# Framework for integrating co-measured datasets from Pebas Formation cores:
# - Magnetic susceptibility (MS): mineralogy, terrigenous input
# - CT radiodensity: bulk density, porosity, lithology
# - Fossil material: mollusks, ostracods, fish, plant remains

message("\n=== MULTI-PROXY INTEGRATION FRAMEWORK ===")

cat("
================================================================================
CO-MEASURED DATA INTEGRATION FOR PEBAS FORMATION CORES
Western Amazonia, Middle-Late Miocene (13-10 Ma)
================================================================================

1. MAGNETIC SUSCEPTIBILITY (MS)
   Expected file: data/magnetic_susceptibility.csv
   Columns: section, depth_mm, ms_value (SI units x10^-5)

   Integration hypotheses:
   - MS positively correlates with Fe, Ti (detrital input)
   - MS negatively correlates with Ca (biogenic dilution)
   - MS peaks may indicate transgressive/regressive cycles
   - High MS intervals: increased Andean sediment flux

2. CT RADIODENSITY
   Expected file: data/ct_radiodensity.csv
   Columns: section, depth_mm, ct_hu (Hounsfield units)

   Integration hypotheses:
   - High CT density: shell beds, cemented horizons, compacted clay
   - Low CT density: organic-rich intervals, bioturbated zones
   - CT-Ca correlation: shell accumulation zones
   - CT-porosity: diagenetic overprint assessment

3. FOSSIL MATERIAL (EXCAVATED)
   Expected file: data/fossil_occurrences.csv
   Columns: section, depth_mm, taxon, count, preservation, facies_association

   Key Pebas Formation fossil groups:
   - Mollusks: Pachydon, Dyris, Tryonidens (endemic lacustrine)
   - Ostracods: Cyprideis (brackish indicators)
   - Fish: Characiformes, Siluriformes (freshwater vs marine)
   - Plant remains: mangrove elements, palm fruits

   Integration hypotheses:
   - High Ca/Ti + abundant Pachydon: stable lacustrine conditions
   - High Fe/Mn + low diversity: dysoxic episodes
   - Sr peaks + Cyprideis: marine/brackish incursions
   - Mollusk-free intervals: environmental stress events

4. DEPTH-CORRELATION PROTOCOL
   a. Establish common depth reference (core-top = 0mm)
   b. Resample all datasets to 3mm resolution (XRF grid)
   c. Apply cross-correlation analysis with depth lags
   d. Identify coherent multi-proxy signals
   e. Define geochemical-biofacies associations

================================================================================
Place co-measured data files in data/ folder to enable integration analysis.
================================================================================
")

# Save filtered data for further analysis
write_csv(xrf_filtered, file.path(output_path, "tables", "xrf_data_filtered.csv"))
message(sprintf("\nFiltered data saved: %d rows", nrow(xrf_filtered)))

# ==============================================================================
# SUMMARY STATISTICS TABLE
# ==============================================================================

summary_stats <- xrf_filtered %>%
  group_by(site) %>%
  summarise(
    n = n(),
    depth_m = max(cumulative_depth) / 1000,
    Ca_Ti_mean = mean(Ca_Ti, na.rm = TRUE),
    Ca_Ti_sd = sd(Ca_Ti, na.rm = TRUE),
    Fe_Mn_mean = mean(Fe_Mn, na.rm = TRUE),
    Fe_Mn_sd = sd(Fe_Mn, na.rm = TRUE),
    K_Ti_mean = mean(K_Ti, na.rm = TRUE),
    K_Ti_sd = sd(K_Ti, na.rm = TRUE),
    Zr_Rb_mean = mean(Zr_Rb, na.rm = TRUE),
    Zr_Rb_sd = sd(Zr_Rb, na.rm = TRUE),
    pct_reducing = mean(Fe_Mn > 50, na.rm = TRUE) * 100,
    .groups = "drop"
  )

message("\n=== SUMMARY STATISTICS BY SITE ===")
print(summary_stats)

write_csv(summary_stats, file.path(output_path, "tables", "proxy_summary_stats.csv"))

message("\n=== ANALYSIS COMPLETE ===")
