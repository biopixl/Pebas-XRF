# ==============================================================================
# Pebas-XRF: QC/QA Report for Retained XRF Data
# ==============================================================================
# Quality control and quality assurance figures for the analyzed samples
# ==============================================================================

library(tidyverse)
library(patchwork)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

base_path <- here::here()
output_path <- file.path(base_path, "output")
qc_path <- file.path(output_path, "figures", "qc_qa")
dir.create(qc_path, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# LOAD DATA
# ==============================================================================

xrf_data <- read_csv(file.path(output_path, "tables", "xrf_data_qc.csv"),
                     show_col_types = FALSE)

# Filter to only QC-passed data (the samples of interest)
xrf_clean <- xrf_data %>% filter(qc_pass == TRUE)

message(sprintf("QC/QA analysis for %d measurements from %d sections",
                nrow(xrf_clean), n_distinct(xrf_clean$section)))

# Key elements for analysis
elements <- c("Fe", "Ti", "Ca", "K", "Si", "Mn", "Sr", "Rb", "Zr")

# ==============================================================================
# 1. INSTRUMENT PERFORMANCE: Total counts (cps) distribution
# ==============================================================================

p1 <- ggplot(xrf_clean, aes(x = cps, fill = core_series)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "A. Total X-ray Counts Distribution",
    subtitle = "Higher counts = better signal quality",
    x = "Total counts per second (cps)",
    y = "Number of measurements",
    fill = "Core"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = c(0.85, 0.85))

# ==============================================================================
# 2. FIT QUALITY: MSE distribution
# ==============================================================================

p2 <- ggplot(xrf_clean, aes(x = MSE, fill = core_series)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  annotate("text", x = 10.5, y = Inf, label = "QC threshold",
           hjust = 0, vjust = 2, color = "red", size = 3) +
  labs(
    title = "B. Spectral Fit Quality (MSE)",
    subtitle = "Lower MSE = better spectral fit",
    x = "Mean Squared Error",
    y = "Number of measurements",
    fill = "Core"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = c(0.85, 0.85))

# ==============================================================================
# 3. ELEMENT DETECTION: Counts above detection limit
# ==============================================================================

detection_stats <- xrf_clean %>%
  summarise(across(all_of(elements), ~sum(. > 0, na.rm = TRUE) / n() * 100)) %>%
  pivot_longer(everything(), names_to = "Element", values_to = "pct_detected") %>%
  mutate(Element = factor(Element, levels = elements))

p3 <- ggplot(detection_stats, aes(x = Element, y = pct_detected)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "darkgreen") +
  annotate("text", x = 0.5, y = 96, label = "95% threshold",
           hjust = 0, color = "darkgreen", size = 3) +
  labs(
    title = "C. Element Detection Rate",
    subtitle = "% of measurements with signal above zero",
    x = "Element",
    y = "Detection rate (%)"
  ) +
  theme_minimal(base_size = 11) +
  ylim(0, 105)

# ==============================================================================
# 4. DATA COMPLETENESS: Measurements per section
# ==============================================================================

section_counts <- xrf_clean %>%
  group_by(section, core_series) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(core_series, section)

p4 <- ggplot(section_counts, aes(x = reorder(section, -n), y = n, fill = core_series)) +
  geom_col(alpha = 0.8) +
  labs(
    title = "D. Measurements per Section (after QC)",
    subtitle = sprintf("Total: %d measurements across %d sections",
                       nrow(xrf_clean), n_distinct(xrf_clean$section)),
    x = NULL,
    y = "Number of measurements",
    fill = "Core"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = c(0.9, 0.85)
  )

# ==============================================================================
# 5. ELEMENT CORRELATIONS: Key diagnostic ratios
# ==============================================================================

p5 <- ggplot(xrf_clean, aes(x = Ti, y = Fe, color = core_series)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_x_log10() + scale_y_log10() +
  labs(
    title = "E. Fe-Ti Correlation (Detrital Proxy)",
    subtitle = "Strong correlation indicates consistent detrital signal",
    x = "Ti (cps)", y = "Fe (cps)",
    color = "Core"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = c(0.15, 0.85))

# Calculate R² for annotation
r2_tam <- cor(log10(xrf_clean$Ti[xrf_clean$core_series == "TAM"]),
              log10(xrf_clean$Fe[xrf_clean$core_series == "TAM"]),
              use = "complete.obs")^2
r2_sc <- cor(log10(xrf_clean$Ti[xrf_clean$core_series == "SC"]),
             log10(xrf_clean$Fe[xrf_clean$core_series == "SC"]),
             use = "complete.obs")^2

p5 <- p5 + annotate("text", x = max(xrf_clean$Ti, na.rm = TRUE) * 0.5,
                    y = min(xrf_clean$Fe, na.rm = TRUE) * 2,
                    label = sprintf("R² TAM: %.2f\nR² SC: %.2f", r2_tam, r2_sc),
                    hjust = 1, size = 3)

# ==============================================================================
# 6. Ca-Ti: Carbonate vs Detrital
# ==============================================================================

p6 <- ggplot(xrf_clean, aes(x = Ti, y = Ca, color = core_series)) +
  geom_point(alpha = 0.3, size = 1) +
  scale_x_log10() + scale_y_log10() +
  labs(
    title = "F. Ca-Ti Relationship",
    subtitle = "Carbonate (Ca) vs detrital (Ti) end members",
    x = "Ti (cps)", y = "Ca (cps)",
    color = "Core"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = c(0.15, 0.85))

# ==============================================================================
# COMBINE AND SAVE
# ==============================================================================

combined <- (p1 + p2) / (p3 + p4) / (p5 + p6) +
  plot_annotation(
    title = "XRF Data Quality Control / Quality Assurance",
    subtitle = sprintf("Pebas Formation cores: TAM (Tamshiyacu) and SC (Santa Corina) | N = %d QC-passed measurements",
                       nrow(xrf_clean)),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )
  )

ggsave(file.path(qc_path, "qc_qa_summary.png"), combined,
       width = 12, height = 14, dpi = 150, bg = "white")

message(sprintf("\nQC/QA summary saved to: %s", file.path(qc_path, "qc_qa_summary.png")))

# ==============================================================================
# SUMMARY STATISTICS TABLE
# ==============================================================================

qc_summary <- xrf_clean %>%
  group_by(core_series) %>%
  summarise(
    n_sections = n_distinct(section),
    n_measurements = n(),
    cps_median = median(cps, na.rm = TRUE),
    cps_iqr = IQR(cps, na.rm = TRUE),
    MSE_median = median(MSE, na.rm = TRUE),
    MSE_max = max(MSE, na.rm = TRUE),
    Fe_median = median(Fe, na.rm = TRUE),
    Ti_median = median(Ti, na.rm = TRUE),
    Ca_median = median(Ca, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(qc_summary, file.path(output_path, "tables", "qc_qa_summary.csv"))

message("\n=== QC/QA SUMMARY ===")
print(qc_summary)

message("\n=== ELEMENT DETECTION RATES ===")
print(detection_stats)
