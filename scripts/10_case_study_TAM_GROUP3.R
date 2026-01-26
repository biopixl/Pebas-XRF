# ==============================================================================
# Pebas-XRF: Case Study - TAM GROUP3 High-Resolution Geochemical Record
# ==============================================================================
# Detailed analysis of the highest-quality continuous XRF record
# Tamshiyacu core, GROUP3 (sections TAM-5AB-6-7-A/B/C)
# 303 measurements at 3.7mm resolution over 1.12m depth
# ==============================================================================

library(tidyverse)
library(zoo)
library(patchwork)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

base_path <- here::here()
output_path <- file.path(base_path, "output")
fig_path <- file.path(output_path, "figures", "case_study")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

WINDOW_SIZE <- 5  # 5-point moving window (18.5mm at 3.7mm step)

# ==============================================================================
# LOAD AND PREPARE DATA
# ==============================================================================

xrf_data <- read_csv(file.path(output_path, "tables", "xrf_data_stacked.csv"),
                     show_col_types = FALSE)

# Extract GROUP3 - highest quality continuous record
g3 <- xrf_data %>%
  filter(group == "GROUP3") %>%
  arrange(position_mm)

message(sprintf("GROUP3: %d measurements, %.0f-%.0f mm depth",
                nrow(g3), min(g3$position_mm), max(g3$position_mm)))

# ==============================================================================
# SIGNAL PROCESSING
# ==============================================================================

g3_processed <- g3 %>%
  mutate(
    # Depth in cm for clearer presentation
    depth_cm = position_mm / 10,

    # Moving window filter
    Fe_filt = rollmean(Fe, WINDOW_SIZE, fill = NA, align = "center"),
    Ca_filt = rollmean(Ca, WINDOW_SIZE, fill = NA, align = "center"),
    Ti_filt = rollmean(Ti, WINDOW_SIZE, fill = NA, align = "center"),
    Mn_filt = rollmean(Mn, WINDOW_SIZE, fill = NA, align = "center"),
    K_filt = rollmean(K, WINDOW_SIZE, fill = NA, align = "center"),

    Ca_Ti_filt = rollmean(Ca_Ti, WINDOW_SIZE, fill = NA, align = "center"),
    Fe_Mn_filt = rollmean(Fe_Mn, WINDOW_SIZE, fill = NA, align = "center"),
    K_Ti_filt = rollmean(K_Ti, WINDOW_SIZE, fill = NA, align = "center"),
    Zr_Rb_filt = rollmean(Zr_Rb, WINDOW_SIZE, fill = NA, align = "center"),

    # Geochemical zonation
    facies = case_when(
      Ca_Ti > 10 ~ "Shell-rich",
      Ca_Ti > 5 ~ "Carbonate",
      Ca_Ti > 2 ~ "Mixed",
      TRUE ~ "Clastic"
    ),
    facies = factor(facies, levels = c("Shell-rich", "Carbonate", "Mixed", "Clastic")),

    # Redox classification
    redox = ifelse(Fe_Mn > 50, "Reducing", "Oxic")
  )

# ==============================================================================
# FIGURE 1: HIGH-RESOLUTION STRATIGRAPHIC COLUMN
# ==============================================================================

message("Generating Figure 1: Stratigraphic column...")

# Color palette for facies
facies_colors <- c("Shell-rich" = "#2166AC", "Carbonate" = "#67A9CF",
                   "Mixed" = "#D1E5F0", "Clastic" = "#B2182B")

# Create multi-panel stratigraphic figure
p_facies <- ggplot(g3_processed, aes(y = depth_cm)) +
  geom_tile(aes(x = 0.5, fill = facies), width = 1, height = 0.4) +
  scale_fill_manual(values = facies_colors, name = "Facies") +
  scale_y_reverse() +
  labs(x = NULL, y = "Depth (cm)") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")

p_cati <- ggplot(g3_processed, aes(y = depth_cm)) +
  geom_point(aes(x = Ca_Ti), alpha = 0.3, size = 0.8, color = "darkgreen") +
  geom_path(aes(x = Ca_Ti_filt), color = "darkgreen", linewidth = 0.6) +
  geom_vline(xintercept = c(2, 5, 10), linetype = "dashed", color = "gray60", linewidth = 0.3) +
  scale_y_reverse() +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20)) +
  labs(x = "Ca/Ti", y = NULL, title = "Carbonate") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_blank())

p_femn <- ggplot(g3_processed, aes(y = depth_cm)) +
  geom_ribbon(aes(xmin = 0, xmax = Fe_Mn_filt,
                  fill = ifelse(Fe_Mn_filt > 50, "Reducing", "Oxic")),
              alpha = 0.4) +
  geom_path(aes(x = Fe_Mn_filt), color = "purple", linewidth = 0.6) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_fill_manual(values = c("Reducing" = "purple", "Oxic" = "orange"), guide = "none") +
  scale_y_reverse() +
  labs(x = "Fe/Mn", y = NULL, title = "Redox") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_blank())

p_fe <- ggplot(g3_processed, aes(y = depth_cm)) +
  geom_point(aes(x = Fe), alpha = 0.3, size = 0.8, color = "brown") +
  geom_path(aes(x = Fe_filt), color = "brown", linewidth = 0.6) +
  scale_y_reverse() +
  scale_x_log10() +
  labs(x = "Fe (cps)", y = NULL, title = "Terrigenous") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_blank())

p_kti <- ggplot(g3_processed, aes(y = depth_cm)) +
  geom_point(aes(x = K_Ti), alpha = 0.3, size = 0.8, color = "orange") +
  geom_path(aes(x = K_Ti_filt), color = "darkorange", linewidth = 0.6) +
  scale_y_reverse() +
  labs(x = "K/Ti", y = NULL, title = "Weathering") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_blank())

p_zrrb <- ggplot(g3_processed, aes(y = depth_cm)) +
  geom_point(aes(x = Zr_Rb), alpha = 0.3, size = 0.8, color = "darkred") +
  geom_path(aes(x = Zr_Rb_filt), color = "darkred", linewidth = 0.6) +
  scale_y_reverse() +
  labs(x = "Zr/Rb", y = NULL, title = "Grain Size") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_blank())

fig1 <- p_facies + p_cati + p_femn + p_fe + p_kti + p_zrrb +
  plot_layout(widths = c(0.6, 1, 1, 1, 1, 1), guides = "collect") +
  plot_annotation(
    title = "TAM GROUP3: High-Resolution Geochemical Stratigraphy",
    subtitle = "Tamshiyacu, Pebas Formation | 303 measurements at 3.7mm resolution | 1.12m continuous record",
    caption = "Dashed lines: Ca/Ti facies thresholds (2, 5, 10); Fe/Mn=50 oxic/reducing boundary",
    theme = theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      plot.caption = element_text(size = 8, color = "gray50"),
      legend.position = "bottom"
    )
  )

ggsave(file.path(fig_path, "fig1_stratigraphy_GROUP3.png"), fig1,
       width = 14, height = 10, dpi = 300, bg = "white")
ggsave(file.path(fig_path, "fig1_stratigraphy_GROUP3.pdf"), fig1,
       width = 14, height = 10, bg = "white")
message("  Saved: fig1_stratigraphy_GROUP3.png/pdf")

# ==============================================================================
# FIGURE 2: ELEMENT CORRELATION MATRIX
# ==============================================================================

message("Generating Figure 2: Correlation matrix...")

cor_vars <- c("Fe", "Ca", "Ti", "K", "Mn", "Si", "Zr", "Rb", "Sr")
cor_data <- g3 %>% select(all_of(cor_vars))
cor_mat <- cor(cor_data, use = "pairwise.complete.obs")

# Convert to long format for ggplot
cor_long <- cor_mat %>%
  as.data.frame() %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "r") %>%
  mutate(
    var1 = factor(var1, levels = cor_vars),
    var2 = factor(var2, levels = rev(cor_vars)),
    label = sprintf("%.2f", r)
  )

fig2 <- ggplot(cor_long, aes(x = var1, y = var2, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                       midpoint = 0, limits = c(-1, 1), name = "Correlation") +
  labs(
    title = "Element Correlation Matrix - TAM GROUP3",
    subtitle = "Pearson correlations (n=303) | Strong Fe-Ti-Rb terrigenous association",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave(file.path(fig_path, "fig2_correlation_GROUP3.png"), fig2,
       width = 8, height = 7, dpi = 300, bg = "white")
message("  Saved: fig2_correlation_GROUP3.png")

# ==============================================================================
# FIGURE 3: BIVARIATE PROXY RELATIONSHIPS
# ==============================================================================

message("Generating Figure 3: Proxy relationships...")

p3a <- ggplot(g3_processed, aes(x = Ca_Ti, y = Fe_Mn, color = depth_cm)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 5, linetype = "dashed", color = "gray50") +
  scale_color_viridis_c(name = "Depth (cm)", direction = -1) +
  labs(x = "Ca/Ti (carbonate)", y = "Fe/Mn (redox)",
       title = "Carbonate vs Redox") +
  theme_minimal(base_size = 10) +
  annotate("text", x = 15, y = 30, label = "Oxic\nCarbonate", size = 3, color = "gray40") +
  annotate("text", x = 15, y = 80, label = "Reducing\nCarbonate", size = 3, color = "gray40") +
  annotate("text", x = 1, y = 30, label = "Oxic\nClastic", size = 3, color = "gray40") +
  annotate("text", x = 1, y = 80, label = "Reducing\nClastic", size = 3, color = "gray40")

p3b <- ggplot(g3_processed, aes(x = K_Ti, y = Zr_Rb, color = facies)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = facies_colors, name = "Facies") +
  labs(x = "K/Ti (weathering)", y = "Zr/Rb (grain size)",
       title = "Weathering vs Grain Size") +
  theme_minimal(base_size = 10)

p3c <- ggplot(g3_processed, aes(x = Fe, y = Ti, color = facies)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
  scale_color_manual(values = facies_colors, name = "Facies") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Fe (cps)", y = "Ti (cps)",
       title = "Fe-Ti Terrigenous Association") +
  theme_minimal(base_size = 10)

p3d <- ggplot(g3_processed, aes(x = Ca, y = Sr, color = facies)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
  scale_color_manual(values = facies_colors, name = "Facies") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Ca (cps)", y = "Sr (cps)",
       title = "Ca-Sr Carbonate Association") +
  theme_minimal(base_size = 10)

fig3 <- (p3a + p3b) / (p3c + p3d) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Bivariate Proxy Relationships - TAM GROUP3",
    subtitle = "Geochemical facies discrimination and elemental associations",
    theme = theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      legend.position = "right"
    )
  )

ggsave(file.path(fig_path, "fig3_bivariate_GROUP3.png"), fig3,
       width = 12, height = 10, dpi = 300, bg = "white")
message("  Saved: fig3_bivariate_GROUP3.png")

# ==============================================================================
# FIGURE 4: DEPTH PROFILE WITH INTERPRETED ZONES
# ==============================================================================

message("Generating Figure 4: Interpreted zonation...")

# Identify distinct intervals based on sustained proxy values
g3_zones <- g3_processed %>%
  mutate(
    zone_id = cumsum(c(1, abs(diff(as.numeric(facies))) > 0)),
    interval = cut(depth_cm, breaks = seq(0, 120, by = 20), labels = paste0(seq(0, 100, by = 20), "-", seq(20, 120, by = 20), " cm"))
  )

# Summarize by interval
interval_summary <- g3_zones %>%
  group_by(interval) %>%
  summarise(
    Ca_Ti_mean = mean(Ca_Ti, na.rm = TRUE),
    Fe_Mn_mean = mean(Fe_Mn, na.rm = TRUE),
    pct_reducing = mean(Fe_Mn > 50, na.rm = TRUE) * 100,
    dominant_facies = names(which.max(table(facies))),
    n = n(),
    .groups = "drop"
  )

cat("\n=== INTERVAL SUMMARY ===\n")
print(interval_summary)

# Create interpreted zonation figure
fig4 <- ggplot(g3_processed, aes(y = depth_cm)) +
  # Background shading for facies
  geom_tile(aes(x = 3, fill = facies), width = 6, alpha = 0.3) +
  # Proxy curves
  geom_path(aes(x = Ca_Ti_filt / 5), color = "darkgreen", linewidth = 0.8) +
  geom_path(aes(x = Fe_Mn_filt / 25), color = "purple", linewidth = 0.8) +
  # Reference lines
  geom_vline(xintercept = 50/25, linetype = "dashed", color = "red", linewidth = 0.4) +
  # Scales

  scale_y_reverse() +
  scale_fill_manual(values = facies_colors, name = "Facies") +
  scale_x_continuous(
    name = "Ca/Ti (green)",
    sec.axis = sec_axis(~.*25, name = "Fe/Mn (purple)")
  ) +
  labs(
    y = "Depth (cm)",
    title = "Integrated Geochemical Zonation - TAM GROUP3",
    subtitle = "Carbonate (Ca/Ti) and redox (Fe/Mn) proxies with facies interpretation",
    caption = "Red dashed line: Fe/Mn=50 oxic/reducing threshold"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(fig_path, "fig4_zonation_GROUP3.png"), fig4,
       width = 8, height = 10, dpi = 300, bg = "white")
message("  Saved: fig4_zonation_GROUP3.png")

# ==============================================================================
# STATISTICAL SUMMARY TABLE
# ==============================================================================

message("\nGenerating summary statistics...")

stats_table <- g3_processed %>%
  group_by(facies) %>%
  summarise(
    n = n(),
    pct = n() / nrow(g3_processed) * 100,
    Ca_Ti_mean = mean(Ca_Ti, na.rm = TRUE),
    Ca_Ti_sd = sd(Ca_Ti, na.rm = TRUE),
    Fe_Mn_mean = mean(Fe_Mn, na.rm = TRUE),
    Fe_Mn_sd = sd(Fe_Mn, na.rm = TRUE),
    K_Ti_mean = mean(K_Ti, na.rm = TRUE),
    Zr_Rb_mean = mean(Zr_Rb, na.rm = TRUE),
    pct_reducing = mean(Fe_Mn > 50, na.rm = TRUE) * 100,
    .groups = "drop"
  )

write_csv(stats_table, file.path(output_path, "tables", "case_study_GROUP3_stats.csv"))

cat("\n=== FACIES STATISTICS ===\n")
print(stats_table)

# ==============================================================================
# EXPORT DATA FOR CASE STUDY
# ==============================================================================

write_csv(g3_processed, file.path(output_path, "tables", "case_study_GROUP3_data.csv"))
message(sprintf("\nData exported: %d rows", nrow(g3_processed)))

message("\n=== CASE STUDY FIGURES COMPLETE ===")
message(sprintf("Output directory: %s", fig_path))
