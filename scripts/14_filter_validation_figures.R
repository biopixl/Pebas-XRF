# ==============================================================================
# Pebas-XRF: Filter Validation Figures
# ==============================================================================
# Generate spectral evidence supporting data filtering decisions:
# 1. Foam zone identification (anomalous elemental signatures)
# 2. Facies threshold validation (Ca/Ti distribution analysis)
# 3. Redox threshold validation (Fe/Mn systematics)
# 4. Window size optimization (signal-to-noise analysis)
# 5. QC filter effectiveness (MSE, CPS thresholds)
# ==============================================================================

library(tidyverse)
library(patchwork)
library(zoo)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

base_path <- here::here()
output_path <- file.path(base_path, "output")
fig_path <- file.path(output_path, "figures", "filter_validation")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

# Load data
xrf_data <- read_csv(file.path(output_path, "tables", "xrf_data_stacked.csv"),
                     show_col_types = FALSE)

# Load exclusion zones
exclusion_zones <- read_csv(file.path(base_path, "data", "exclusion_zones.csv"),
                            show_col_types = FALSE) %>%
  filter(!is.na(exclude_start_mm))

message(sprintf("Loaded %d measurements, %d exclusion zones",
                nrow(xrf_data), nrow(exclusion_zones)))

# ==============================================================================
# 1. FOAM ZONE SPECTRAL SIGNATURES
# ==============================================================================

message("\n=== Generating foam zone validation figures ===")

# Compare elemental signatures in foam vs sediment
foam_comparison <- xrf_data %>%
  mutate(
    zone_type = ifelse(excluded, "Foam/Gap", "Sediment"),
    zone_type = factor(zone_type, levels = c("Sediment", "Foam/Gap"))
  ) %>%
  filter(!is.na(zone_type))

# Key elements that distinguish foam from sediment
diagnostic_elements <- c("Ca", "Ti", "Fe", "Mn", "K", "Rb", "Sr", "Zr", "Mo_inc", "Mo_coh")

# Create comparison boxplots
foam_boxplots <- foam_comparison %>%
  select(zone_type, all_of(diagnostic_elements)) %>%
  pivot_longer(cols = -zone_type, names_to = "element", values_to = "counts") %>%
  mutate(element = factor(element, levels = diagnostic_elements)) %>%
  ggplot(aes(x = zone_type, y = counts, fill = zone_type)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  facet_wrap(~element, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c("Sediment" = "#2166AC", "Foam/Gap" = "#B2182B")) +
  scale_y_log10() +
  labs(
    title = "Spectral Signatures: Sediment vs Excluded Zones",
    subtitle = "Log-scale element counts distinguish foam/gaps from valid sediment measurements",
    x = NULL, y = "Counts (log scale)",
    fill = "Zone Type"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(fig_path, "fig_V1_foam_spectral_signatures.png"), foam_boxplots,
       width = 14, height = 8, dpi = 300, bg = "white")
ggsave(file.path(fig_path, "fig_V1_foam_spectral_signatures.pdf"), foam_boxplots,
       width = 14, height = 8, bg = "white")

# Statistical summary
foam_stats <- foam_comparison %>%
  group_by(zone_type) %>%
  summarise(
    n = n(),
    across(all_of(diagnostic_elements),
           list(median = ~median(., na.rm = TRUE),
                iqr = ~IQR(., na.rm = TRUE)),
           .names = "{.col}_{.fn}")
  )

write_csv(foam_stats, file.path(fig_path, "foam_zone_statistics.csv"))
message("  Saved: fig_V1_foam_spectral_signatures.png/pdf")

# ==============================================================================
# 2. INC/COH RATIO FOR FOAM DETECTION
# ==============================================================================

# Inc/Coh ratio is diagnostic for organic matter and matrix composition
inc_coh_comparison <- foam_comparison %>%
  mutate(inc_coh = Mo_inc / Mo_coh) %>%
  filter(is.finite(inc_coh))

p_inc_coh <- ggplot(inc_coh_comparison, aes(x = inc_coh, fill = zone_type)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = c(2, 4), linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Sediment" = "#2166AC", "Foam/Gap" = "#B2182B")) +
  labs(
    title = "Incoherent/Coherent Scattering Ratio: Foam Detection Criterion",
    subtitle = "Inc/Coh ratio reflects matrix composition; foam shows anomalous values",
    x = "Mo Inc/Coh Ratio", y = "Density",
    fill = "Zone Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_path, "fig_V2_inc_coh_foam_detection.png"), p_inc_coh,
       width = 10, height = 6, dpi = 300, bg = "white")
message("  Saved: fig_V2_inc_coh_foam_detection.png")

# ==============================================================================
# 3. FACIES THRESHOLD VALIDATION (Ca/Ti)
# ==============================================================================

message("\n=== Generating facies threshold validation ===")

# Filter to valid sediment only
sediment_data <- xrf_data %>% filter(!excluded)

# Ca/Ti distribution with proposed thresholds
p_cati_dist <- ggplot(sediment_data, aes(x = Ca_Ti)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100,
                 fill = "#2166AC", alpha = 0.7, color = "white") +
  geom_density(color = "darkblue", linewidth = 1) +
  geom_vline(xintercept = c(2, 5, 10), linetype = "dashed",
             color = c("#B2182B", "#F4A582", "#2166AC"), linewidth = 1) +
  annotate("text", x = 1, y = 0.15, label = "Clastic\n(<2)", color = "#B2182B", size = 3.5) +
  annotate("text", x = 3.5, y = 0.15, label = "Mixed\n(2-5)", color = "#F4A582", size = 3.5) +
  annotate("text", x = 7.5, y = 0.15, label = "Carbonate\n(5-10)", color = "#4393C3", size = 3.5) +
  annotate("text", x = 15, y = 0.15, label = "Shell-rich\n(>10)", color = "#2166AC", size = 3.5) +
  scale_x_log10(limits = c(0.5, 50)) +
  labs(
    title = "Ca/Ti Ratio Distribution with Facies Thresholds",
    subtitle = sprintf("n = %d measurements | Thresholds based on natural breaks in distribution",
                       nrow(sediment_data)),
    x = "Ca/Ti (log scale)", y = "Density"
  ) +
  theme_minimal(base_size = 12)

# Ca/Ti by site comparison
p_cati_site <- ggplot(sediment_data, aes(x = Ca_Ti, fill = site)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = c(2, 5, 10), linetype = "dashed", color = "gray40") +
  scale_x_log10(limits = c(0.5, 50)) +
  scale_fill_manual(values = c("TAM" = "#2166AC", "SC" = "#B2182B"),
                    labels = c("TAM" = "Tamshiyacu", "SC" = "Santa Corina")) +
  labs(
    title = "Ca/Ti Distribution by Site",
    subtitle = "Thresholds applicable to both localities despite different mean values",
    x = "Ca/Ti (log scale)", y = "Density",
    fill = "Site"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

p_facies_combined <- p_cati_dist / p_cati_site +
  plot_annotation(tag_levels = "A")

ggsave(file.path(fig_path, "fig_V3_facies_threshold_validation.png"), p_facies_combined,
       width = 10, height = 10, dpi = 300, bg = "white")
ggsave(file.path(fig_path, "fig_V3_facies_threshold_validation.pdf"), p_facies_combined,
       width = 10, height = 10, bg = "white")
message("  Saved: fig_V3_facies_threshold_validation.png/pdf")

# Facies proportions by site
facies_props <- sediment_data %>%
  mutate(facies = case_when(
    Ca_Ti > 10 ~ "Shell-rich",
    Ca_Ti > 5 ~ "Carbonate",
    Ca_Ti > 2 ~ "Mixed",
    TRUE ~ "Clastic"
  )) %>%
  count(site, facies) %>%
  group_by(site) %>%
  mutate(prop = n / sum(n) * 100) %>%
  ungroup()

write_csv(facies_props, file.path(fig_path, "facies_proportions_by_site.csv"))

# ==============================================================================
# 4. REDOX THRESHOLD VALIDATION (Fe/Mn)
# ==============================================================================

message("\n=== Generating redox threshold validation ===")

# Fe/Mn distribution with threshold
p_femn_dist <- ggplot(sediment_data, aes(x = Fe_Mn)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100,
                 fill = "#762A83", alpha = 0.7, color = "white") +
  geom_density(color = "purple4", linewidth = 1) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red", linewidth = 1.2) +
  annotate("text", x = 25, y = 0.015, label = "Oxic\n(<50)", color = "#F1A340", size = 4) +
  annotate("text", x = 100, y = 0.015, label = "Reducing\n(>50)", color = "#762A83", size = 4) +
  scale_x_log10(limits = c(5, 500)) +
  labs(
    title = "Fe/Mn Ratio Distribution with Redox Threshold",
    subtitle = "Threshold at Fe/Mn = 50 separates oxic and reducing facies (Calvert & Pedersen, 1996)",
    x = "Fe/Mn (log scale)", y = "Density"
  ) +
  theme_minimal(base_size = 12)

# Fe/Mn by site
p_femn_site <- ggplot(sediment_data, aes(x = Fe_Mn, fill = site)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red", linewidth = 1) +
  scale_x_log10(limits = c(5, 500)) +
  scale_fill_manual(values = c("TAM" = "#2166AC", "SC" = "#B2182B"),
                    labels = c("TAM" = "Tamshiyacu", "SC" = "Santa Corina")) +
  labs(
    title = "Fe/Mn Distribution by Site",
    subtitle = "TAM shows more reducing conditions (higher Fe/Mn) than SC",
    x = "Fe/Mn (log scale)", y = "Density",
    fill = "Site"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Fe vs Mn crossplot colored by site
p_fe_mn_cross <- ggplot(sediment_data, aes(x = Mn, y = Fe, color = site)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(intercept = 0, slope = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 8000, y = 200000, label = "Fe/Mn = 50", color = "red", angle = 30) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = c("TAM" = "#2166AC", "SC" = "#B2182B"),
                     labels = c("TAM" = "Tamshiyacu", "SC" = "Santa Corina")) +
  labs(
    title = "Fe-Mn Systematics",
    subtitle = "Parallel trends confirm redox control; TAM offset to higher Fe/Mn",
    x = "Mn (cps)", y = "Fe (cps)",
    color = "Site"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

p_redox_combined <- (p_femn_dist | p_femn_site) / p_fe_mn_cross +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(fig_path, "fig_V4_redox_threshold_validation.png"), p_redox_combined,
       width = 12, height = 10, dpi = 300, bg = "white")
ggsave(file.path(fig_path, "fig_V4_redox_threshold_validation.pdf"), p_redox_combined,
       width = 12, height = 10, bg = "white")
message("  Saved: fig_V4_redox_threshold_validation.png/pdf")

# Redox statistics by site
redox_stats <- sediment_data %>%
  group_by(site) %>%
  summarise(
    n = n(),
    Fe_Mn_median = median(Fe_Mn, na.rm = TRUE),
    Fe_Mn_q25 = quantile(Fe_Mn, 0.25, na.rm = TRUE),
    Fe_Mn_q75 = quantile(Fe_Mn, 0.75, na.rm = TRUE),
    pct_reducing = mean(Fe_Mn > 50, na.rm = TRUE) * 100
  )

write_csv(redox_stats, file.path(fig_path, "redox_statistics_by_site.csv"))

# ==============================================================================
# 5. WINDOW SIZE OPTIMIZATION
# ==============================================================================

message("\n=== Generating window size optimization ===")

# Test different window sizes on a representative section
test_section <- sediment_data %>%
  filter(section == "TAM-5AB-6-7-B-RUN2") %>%
  arrange(position_mm)

window_sizes <- c(1, 3, 5, 7, 11)

window_comparison <- map_dfr(window_sizes, function(ws) {
  test_section %>%
    mutate(
      Ca_Ti_filt = if(ws == 1) Ca_Ti else zoo::rollmean(Ca_Ti, ws, fill = NA, align = "center"),
      window_size = paste0("Window = ", ws)
    ) %>%
    select(position_mm, Ca_Ti, Ca_Ti_filt, window_size)
}) %>%
  mutate(window_size = factor(window_size,
                               levels = paste0("Window = ", window_sizes)))

p_window <- ggplot(window_comparison, aes(y = position_mm / 10)) +
  geom_point(aes(x = Ca_Ti), alpha = 0.2, size = 0.5, color = "gray50") +
  geom_path(aes(x = Ca_Ti_filt), color = "#2166AC", linewidth = 0.8, na.rm = TRUE) +
  facet_wrap(~window_size, nrow = 1) +
  scale_y_reverse() +
  labs(
    title = "Window Size Comparison for Moving Average Filter",
    subtitle = "TAM-5AB-6-7-B-RUN2 | Window = 5 (15mm) balances noise reduction with feature preservation",
    x = "Ca/Ti", y = "Depth (cm)"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(fig_path, "fig_V5_window_size_optimization.png"), p_window,
       width = 14, height = 8, dpi = 300, bg = "white")
ggsave(file.path(fig_path, "fig_V5_window_size_optimization.pdf"), p_window,
       width = 14, height = 8, bg = "white")
message("  Saved: fig_V5_window_size_optimization.png/pdf")

# ==============================================================================
# 6. QC FILTER EFFECTIVENESS
# ==============================================================================

message("\n=== Generating QC filter effectiveness ===")

# Load raw QC data
xrf_qc <- read_csv(file.path(output_path, "tables", "xrf_data_qc.csv"),
                   show_col_types = FALSE)

# QC filter summary
qc_summary <- xrf_qc %>%
  summarise(
    n_total = n(),
    n_pass_mse = sum(qc_mse, na.rm = TRUE),
    n_pass_cps = sum(qc_cps, na.rm = TRUE),
    n_pass_surface = sum(qc_surface, na.rm = TRUE),
    n_pass_all = sum(qc_pass, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "filter", values_to = "n") %>%
  mutate(
    pct = n / first(n) * 100,
    filter = factor(filter,
                    levels = c("n_total", "n_pass_mse", "n_pass_cps", "n_pass_surface", "n_pass_all"),
                    labels = c("Total", "MSE ≤ 10", "CPS ≥ 20k", "Surface < 8", "All Filters"))
  )

# MSE distribution
p_mse <- ggplot(xrf_qc, aes(x = MSE, fill = qc_mse)) +
  geom_histogram(bins = 100, alpha = 0.7) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                    labels = c("TRUE" = "Pass", "FALSE" = "Fail")) +
  labs(
    title = "MSE Threshold (≤ 10)",
    subtitle = "Mean squared error of spectral fit",
    x = "MSE", y = "Count", fill = "QC Status"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

# CPS distribution
p_cps <- ggplot(xrf_qc, aes(x = cps, fill = qc_cps)) +
  geom_histogram(bins = 100, alpha = 0.7) +
  geom_vline(xintercept = 20000, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                    labels = c("TRUE" = "Pass", "FALSE" = "Fail")) +
  scale_x_log10() +
  labs(
    title = "CPS Threshold (≥ 20,000)",
    subtitle = "Total counts per second (signal strength)",
    x = "CPS (log scale)", y = "Count", fill = "QC Status"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

# Sample surface distribution
p_surface <- ggplot(xrf_qc, aes(x = sample_surface, fill = qc_surface)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  geom_vline(xintercept = 8, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                    labels = c("TRUE" = "Pass", "FALSE" = "Fail")) +
  labs(
    title = "Sample Surface Threshold (< 8)",
    subtitle = "Distance from detector (mm)",
    x = "Sample Surface (mm)", y = "Count", fill = "QC Status"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

# QC filter cascade
p_qc_cascade <- ggplot(qc_summary, aes(x = filter, y = pct, fill = filter)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct, n)),
            vjust = -0.3, size = 3.5) +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.8) +
  scale_y_continuous(limits = c(0, 110)) +
  labs(
    title = "QC Filter Cascade",
    subtitle = "Progressive data retention through quality filters",
    x = NULL, y = "Retention (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

p_qc_combined <- (p_mse | p_cps | p_surface) / p_qc_cascade +
  plot_layout(heights = c(1, 0.8)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(fig_path, "fig_V6_qc_filter_effectiveness.png"), p_qc_combined,
       width = 14, height = 10, dpi = 300, bg = "white")
ggsave(file.path(fig_path, "fig_V6_qc_filter_effectiveness.pdf"), p_qc_combined,
       width = 14, height = 10, bg = "white")
message("  Saved: fig_V6_qc_filter_effectiveness.png/pdf")

write_csv(qc_summary, file.path(fig_path, "qc_filter_summary.csv"))

# ==============================================================================
# 7. ALIGNMENT VERIFICATION - STRATIGRAPHIC CORRELATION
# ==============================================================================

message("\n=== Generating alignment verification ===")

# Show alignment of key features across panels for GROUP3
g3_data <- sediment_data %>%
  filter(group == "GROUP3") %>%
  arrange(position_mm) %>%
  mutate(
    depth_cm = position_mm / 10,
    Ca_Ti_filt = zoo::rollmean(Ca_Ti, 5, fill = NA, align = "center"),
    Fe_Mn_filt = zoo::rollmean(Fe_Mn, 5, fill = NA, align = "center")
  )

# Mark key stratigraphic features
key_features <- tibble(
  depth_cm = c(5, 35, 55, 75, 95),
  label = c("High carbonate", "Transition", "Carbonate peak", "Reducing peak", "Low carbonate")
)

p_alignment <- ggplot(g3_data, aes(y = depth_cm)) +
  # Ca/Ti panel
  geom_path(aes(x = Ca_Ti_filt / 20), color = "#2166AC", linewidth = 0.8, na.rm = TRUE) +
  # Fe/Mn panel (scaled)
  geom_path(aes(x = Fe_Mn_filt / 100 + 1), color = "#762A83", linewidth = 0.8, na.rm = TRUE) +
  # Feature markers
  geom_hline(yintercept = key_features$depth_cm, linetype = "dotted", color = "gray50") +
  scale_y_reverse() +
  labs(
    title = "Stratigraphic Feature Alignment: GROUP3",
    subtitle = "Blue = Ca/Ti (scaled), Purple = Fe/Mn (scaled) | Dotted lines = key stratigraphic markers",
    x = "Scaled proxy values", y = "Depth (cm)"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(fig_path, "fig_V7_alignment_verification.png"), p_alignment,
       width = 8, height = 10, dpi = 300, bg = "white")
message("  Saved: fig_V7_alignment_verification.png")

# ==============================================================================
# SUMMARY
# ==============================================================================

message("\n=== FILTER VALIDATION COMPLETE ===")
message(sprintf("Figures saved to: %s", fig_path))
message("\nGenerated figures:")
message("  V1: Foam spectral signatures")
message("  V2: Inc/Coh foam detection")
message("  V3: Facies threshold validation")
message("  V4: Redox threshold validation")
message("  V5: Window size optimization")
message("  V6: QC filter effectiveness")
message("  V7: Alignment verification")
