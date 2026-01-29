# ==============================================================================
# Pebas-XRF: Proxy Evaluation Figures for Manuscript
# ==============================================================================
# Generates diagnostic figures supporting proxy selection for the
# Miocene Pebas Formation lacustrine-wetland setting
# ==============================================================================

library(tidyverse)
library(patchwork)
library(ggrepel)
library(scales)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

base_path <- here::here()
output_path <- file.path(base_path, "output")
fig_path <- file.path(output_path, "figures", "proxy_evaluation")

# Create output directory
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

# Load ratio data
xrf <- read_csv(file.path(output_path, "tables", "xrf_data_ratios.csv"),
                show_col_types = FALSE) %>%
  filter(qc_pass, !excluded)

message(sprintf("Loaded %d measurements for proxy evaluation", nrow(xrf)))

# Define consistent theme for publication
theme_pub <- theme_minimal(base_size = 10) +

theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey70"),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    legend.position = "bottom"
  )

# Color palette for core series
core_colors <- c("TAM" = "#E69F00", "SC" = "#0072B2")

# ==============================================================================
# 2. FIGURE 1: REDUNDANCY ANALYSIS - Ca/Ti vs Ca
# ==============================================================================

message("Generating Figure 1: Redundancy analysis...")

# Panel A: Ca vs Ca/Ti scatter
p1a <- ggplot(xrf, aes(x = Ca, y = Ca_Ti, color = core_series)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10() +
  scale_color_manual(values = core_colors) +
  labs(
    title = "A) Ca/Ti is redundant with Ca",
    subtitle = sprintf("r = %.2f (ratio provides no additional information)",
                       cor(xrf$Ca, xrf$Ca_Ti, use = "complete.obs")),
    x = "Ca (cps)",
    y = "Ca/Ti",
    color = "Core"
  ) +
  theme_pub +
  theme(legend.position = "none")

# Panel B: Fe vs Fe/Mn scatter (useful ratio)
p1b <- ggplot(xrf, aes(x = Fe, y = Fe_Mn, color = core_series)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10() +
  scale_color_manual(values = core_colors) +
  labs(
    title = "B) Fe/Mn adds information beyond Fe",
    subtitle = sprintf("r = %.2f (Mn normalization differentiates signal)",
                       cor(xrf$Fe, xrf$Fe_Mn, use = "complete.obs")),
    x = "Fe (cps)",
    y = "Fe/Mn",
    color = "Core"
  ) +
  theme_pub +
  theme(legend.position = "none")

# Panel C: Zr vs Zr/Rb scatter
p1c <- ggplot(xrf, aes(x = Zr, y = Zr_Rb, color = core_series)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10() +
  scale_color_manual(values = core_colors) +
  labs(
    title = "C) Zr/Rb adds information beyond Zr",
    subtitle = sprintf("r = %.2f (Rb normalization improves grain size signal)",
                       cor(xrf$Zr, xrf$Zr_Rb, use = "complete.obs")),
    x = "Zr (cps)",
    y = "Zr/Rb",
    color = "Core"
  ) +
  theme_pub +
  theme(legend.position = "none")

# Panel D: Summary bar chart of correlations
redundancy_data <- tibble(
  proxy = c("Ca/Ti", "Ba/Ti", "Zr/Rb", "Fe/Ti", "K/Ti", "Fe/Mn", "Rb/Sr"),
  r = c(
    cor(xrf$Ca, xrf$Ca_Ti, use = "complete.obs"),
    cor(xrf$Ba, xrf$Ba_Ti, use = "complete.obs"),
    cor(xrf$Zr, xrf$Zr_Rb, use = "complete.obs"),
    cor(xrf$Fe, xrf$Fe_Ti, use = "complete.obs"),
    cor(xrf$K, xrf$K_Ti, use = "complete.obs"),
    cor(xrf$Fe, xrf$Fe_Mn, use = "complete.obs"),
    cor(xrf$Rb, xrf$Rb_Sr, use = "complete.obs")
  )
) %>%
  mutate(
    status = ifelse(abs(r) > 0.85, "Redundant", "Useful"),
    proxy = factor(proxy, levels = proxy[order(abs(r), decreasing = TRUE)])
  )

p1d <- ggplot(redundancy_data, aes(x = proxy, y = abs(r), fill = status)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "red") +
  annotate("text", x = 6.5, y = 0.88, label = "Redundancy threshold",
           size = 3, color = "red") +
  scale_fill_manual(values = c("Redundant" = "#D55E00", "Useful" = "#009E73")) +
  labs(
    title = "D) Proxy redundancy summary",
    subtitle = "Ratios with |r| > 0.85 vs numerator are redundant",
    x = "Elemental Ratio",
    y = "|Correlation with numerator|",
    fill = ""
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine panels
fig1 <- (p1a | p1b) / (p1c | p1d) +
  plot_annotation(
    title = "Figure S1: Proxy Redundancy Analysis",
    subtitle = "Empirical evaluation of ratio utility in Pebas Formation XRF data",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )

ggsave(file.path(fig_path, "fig_S1_redundancy_analysis.pdf"), fig1,
       width = 10, height = 8)
ggsave(file.path(fig_path, "fig_S1_redundancy_analysis.png"), fig1,
       width = 10, height = 8, dpi = 300)

# ==============================================================================
# 3. FIGURE 2: ELEMENT VARIABILITY (CV Analysis)
# ==============================================================================

message("Generating Figure 2: Element variability...")

# Calculate CV for each element
cv_data <- xrf %>%
  summarise(across(
    c(Ca, Ti, Fe, Mn, K, Rb, Sr, Zr, Al, Si, Ba),
    list(
      mean = ~mean(., na.rm = TRUE),
      sd = ~sd(., na.rm = TRUE),
      cv = ~100 * sd(., na.rm = TRUE) / mean(., na.rm = TRUE)
    )
  )) %>%
  pivot_longer(everything(), names_to = c("element", "stat"), names_sep = "_") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    detection = case_when(
      element %in% c("Fe", "Ti", "Ca", "Mn", "Rb", "Sr", "Zr") ~ "Excellent",
      element %in% c("K", "Ba") ~ "Good",
      element %in% c("Al", "Si") ~ "Marginal"
    ),
    element = factor(element, levels = element[order(cv, decreasing = TRUE)])
  )

p2 <- ggplot(cv_data, aes(x = element, y = cv, fill = detection)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", cv)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c(
    "Excellent" = "#009E73",
    "Good" = "#E69F00",
    "Marginal" = "#D55E00"
  )) +
  labs(
    title = "Figure S2: Element Variability and Detection Quality",
    subtitle = "Coefficient of Variation (CV) indicates signal dynamic range; color indicates Mo tube detection",
    x = "Element",
    y = "Coefficient of Variation (%)",
    fill = "Mo Tube Detection"
  ) +
  theme_pub +
  theme(legend.position = "right")

ggsave(file.path(fig_path, "fig_S2_element_variability.pdf"), p2,
       width = 8, height = 5)
ggsave(file.path(fig_path, "fig_S2_element_variability.png"), p2,
       width = 8, height = 5, dpi = 300)

# ==============================================================================
# 4. FIGURE 3: ELEMENT CORRELATION MATRIX
# ==============================================================================

message("Generating Figure 3: Correlation matrix...")

# Calculate correlation matrix
cor_matrix <- xrf %>%
  select(Ca, Ti, Fe, Mn, K, Rb, Sr, Zr) %>%
  cor(use = "pairwise.complete.obs")

# Convert to long format for ggplot
cor_long <- as.data.frame(cor_matrix) %>%
  rownames_to_column("element1") %>%
  pivot_longer(-element1, names_to = "element2", values_to = "r") %>%
  mutate(
    element1 = factor(element1, levels = c("Ca", "Sr", "Ti", "Fe", "Mn", "K", "Rb", "Zr")),
    element2 = factor(element2, levels = c("Ca", "Sr", "Ti", "Fe", "Mn", "K", "Rb", "Zr"))
  )

p3 <- ggplot(cor_long, aes(x = element1, y = element2, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
  scale_fill_gradient2(
    low = "#0072B2", mid = "white", high = "#D55E00",
    midpoint = 0, limits = c(-1, 1)
  ) +
  labs(
    title = "Figure S3: Element Correlation Matrix",
    subtitle = "Pairwise correlations identify geochemical associations",
    x = "", y = "",
    fill = "r"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave(file.path(fig_path, "fig_S3_correlation_matrix.pdf"), p3,
       width = 7, height = 6)
ggsave(file.path(fig_path, "fig_S3_correlation_matrix.png"), p3,
       width = 7, height = 6, dpi = 300)

# ==============================================================================
# 5. FIGURE 4: PROXY CROSS-PLOTS WITH INTERPRETIVE ANNOTATIONS
# ==============================================================================

message("Generating Figure 4: Interpretive cross-plots...")

# Panel A: Ca vs Ti (carbonate-detrital mixing)
p4a <- ggplot(xrf, aes(x = Ti, y = Ca, color = core_series)) +
  geom_point(alpha = 0.3, size = 0.8) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  scale_color_manual(values = core_colors) +
  annotate("text", x = 25000, y = 5000, label = "Terrigenous\ndominance",
           size = 3, fontface = "italic") +
  annotate("text", x = 500, y = 300000, label = "Carbonate\ndominance",
           size = 3, fontface = "italic") +
  labs(
    title = "A) Ca vs Ti: Carbonate-Detrital Mixing",
    x = "Ti (cps) - Detrital flux",
    y = "Ca (cps) - Carbonate signal",
    color = "Core"
  ) +
  theme_pub +
  theme(legend.position = "none")

# Panel B: Fe vs Mn (redox systematics)
p4b <- ggplot(xrf, aes(x = Mn, y = Fe, color = core_series)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_abline(intercept = 0, slope = 50, linetype = "dashed", color = "grey50") +
  geom_abline(intercept = 0, slope = 100, linetype = "dashed", color = "grey50") +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  scale_color_manual(values = core_colors) +
  annotate("text", x = 8000, y = 300000, label = "Fe/Mn = 100\n(reducing)",
           size = 2.5, fontface = "italic") +
  annotate("text", x = 8000, y = 150000, label = "Fe/Mn = 50\n(oxic)",
           size = 2.5, fontface = "italic") +
  labs(
    title = "B) Fe vs Mn: Redox Systematics",
    subtitle = "TAM shows more reducing conditions than SC",
    x = "Mn (cps)",
    y = "Fe (cps)",
    color = "Core"
  ) +
  theme_pub +
  theme(legend.position = "none")

# Panel C: Zr vs Rb (grain size)
p4c <- ggplot(xrf, aes(x = Rb, y = Zr, color = core_series)) +
  geom_point(alpha = 0.3, size = 0.8) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  scale_color_manual(values = core_colors) +
  annotate("text", x = 150, y = 3000, label = "Coarse\n(high energy)",
           size = 3, fontface = "italic") +
  annotate("text", x = 800, y = 200, label = "Fine\n(low energy)",
           size = 3, fontface = "italic") +
  labs(
    title = "C) Zr vs Rb: Grain Size Systematics",
    x = "Rb (cps) - Clay fraction",
    y = "Zr (cps) - Heavy minerals",
    color = "Core"
  ) +
  theme_pub +
  theme(legend.position = "none")

# Panel D: Sr vs Ca (carbonate mineralogy)
p4d <- ggplot(xrf, aes(x = Ca, y = Sr, color = core_series)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  scale_color_manual(values = core_colors) +
  annotate("text", x = 300000, y = 15000, label = "High Sr/Ca\n(aragonite)",
           size = 3, fontface = "italic") +
  annotate("text", x = 300000, y = 2000, label = "Low Sr/Ca\n(calcite)",
           size = 3, fontface = "italic") +
  labs(
    title = "D) Sr vs Ca: Carbonate Mineralogy",
    subtitle = sprintf("r = %.2f", cor(xrf$Sr, xrf$Ca, use = "complete.obs")),
    x = "Ca (cps)",
    y = "Sr (cps)",
    color = "Core"
  ) +
  theme_pub

# Combine panels
fig4 <- (p4a | p4b) / (p4c | p4d) +
  plot_annotation(
    title = "Figure S4: Proxy Cross-Plots with Interpretive Framework",
    subtitle = "Geochemical relationships supporting proxy interpretation in Pebas Formation",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )

ggsave(file.path(fig_path, "fig_S4_crossplots.pdf"), fig4,
       width = 10, height = 8)
ggsave(file.path(fig_path, "fig_S4_crossplots.png"), fig4,
       width = 10, height = 8, dpi = 300)

# ==============================================================================
# 6. FIGURE 5: PROXY DISTRIBUTIONS BY CORE
# ==============================================================================

message("Generating Figure 5: Proxy distributions...")

# Prepare data in long format
proxy_long <- xrf %>%
  select(core_series, Ca, Ti, Fe_Mn, Zr_Rb, Sr, Fe) %>%
  pivot_longer(-core_series, names_to = "proxy", values_to = "value") %>%
  mutate(
    proxy = factor(proxy, levels = c("Ca", "Ti", "Fe_Mn", "Zr_Rb", "Sr", "Fe")),
    proxy_label = recode(proxy,
      "Ca" = "Ca (carbonate)",
      "Ti" = "Ti (detrital)",
      "Fe_Mn" = "Fe/Mn (redox)",
      "Zr_Rb" = "Zr/Rb (grain size)",
      "Sr" = "Sr (mineralogy)",
      "Fe" = "Fe (lateritic)"
    )
  )

p5 <- ggplot(proxy_long, aes(x = value, fill = core_series)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~proxy_label, scales = "free", ncol = 3) +
  scale_x_log10(labels = label_comma()) +
  scale_fill_manual(values = core_colors) +
  labs(
    title = "Figure S5: Proxy Distributions by Core Series",
    subtitle = "Density plots showing geochemical differences between TAM and SC cores",
    x = "Value (cps or ratio)",
    y = "Density",
    fill = "Core Series"
  ) +
  theme_pub

ggsave(file.path(fig_path, "fig_S5_distributions.pdf"), p5,
       width = 10, height = 6)
ggsave(file.path(fig_path, "fig_S5_distributions.png"), p5,
       width = 10, height = 6, dpi = 300)

# ==============================================================================
# 7. FIGURE 6: STATISTICAL SUMMARY TABLE
# ==============================================================================

message("Generating statistical summary...")

# Summary statistics by core series
summary_stats <- xrf %>%
  group_by(core_series) %>%
  summarise(
    n = n(),
    across(
      c(Ca, Ti, Fe, Mn, Sr, Zr, Rb, Fe_Mn, Zr_Rb),
      list(
        median = ~median(., na.rm = TRUE),
        q25 = ~quantile(., 0.25, na.rm = TRUE),
        q75 = ~quantile(., 0.75, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  )

write_csv(summary_stats, file.path(fig_path, "proxy_summary_statistics.csv"))

# ==============================================================================
# 8. COMBINED MAIN FIGURE: RECOMMENDED PROXY PANEL
# ==============================================================================

message("Generating combined main figure...")

# Select one representative section for stratigraphic example
example_section <- xrf %>%
  filter(section == "TAM-1-2-3B-A") %>%
  arrange(position_mm)

# Create 4-panel stratigraphic profile
p_ca <- ggplot(example_section, aes(x = Ca, y = position_mm)) +
  geom_path(color = "#0072B2", linewidth = 0.5) +
  scale_y_reverse() +
  scale_x_log10(labels = label_comma()) +
  labs(x = "Ca (cps)", y = "Position (mm)") +
  theme_pub +
  theme(axis.title.y = element_blank())

p_ti <- ggplot(example_section, aes(x = Ti, y = position_mm)) +
  geom_path(color = "#E69F00", linewidth = 0.5) +
  scale_y_reverse() +
  scale_x_log10(labels = label_comma()) +
  labs(x = "Ti (cps)", y = "") +
  theme_pub +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())

p_femn <- ggplot(example_section, aes(x = Fe_Mn, y = position_mm)) +
  geom_path(color = "#D55E00", linewidth = 0.5) +
  scale_y_reverse() +
  scale_x_log10() +
  labs(x = "Fe/Mn", y = "") +
  theme_pub +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())

p_zrrb <- ggplot(example_section, aes(x = Zr_Rb, y = position_mm)) +
  geom_path(color = "#009E73", linewidth = 0.5) +
  scale_y_reverse() +
  scale_x_log10() +
  labs(x = "Zr/Rb", y = "") +
  theme_pub +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())

fig_strat <- p_ca | p_ti | p_femn | p_zrrb

fig_main <- fig_strat +
  plot_annotation(
    title = "Figure 2: Recommended Proxy Suite - Stratigraphic Example",
    subtitle = "Section TAM-1-2-3B-A showing independent signals from each proxy",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )

ggsave(file.path(fig_path, "fig_main_stratigraphic_example.pdf"), fig_main,
       width = 10, height = 6)
ggsave(file.path(fig_path, "fig_main_stratigraphic_example.png"), fig_main,
       width = 10, height = 6, dpi = 300)

# ==============================================================================
# 9. COMPLETE
# ==============================================================================

message("\n=== PROXY EVALUATION FIGURES COMPLETE ===")
message(sprintf("Figures saved to: %s", fig_path))
message("\nGenerated figures:")
message("  - fig_S1_redundancy_analysis.pdf/png")
message("  - fig_S2_element_variability.pdf/png")
message("  - fig_S3_correlation_matrix.pdf/png")
message("  - fig_S4_crossplots.pdf/png")
message("  - fig_S5_distributions.pdf/png")
message("  - fig_main_stratigraphic_example.pdf/png")
message("  - proxy_summary_statistics.csv")
