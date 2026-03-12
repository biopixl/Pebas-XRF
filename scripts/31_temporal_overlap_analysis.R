# Temporal Overlap Analysis
# Generates figures comparing TAM and SC for the 171 kyr overlap period
# TAM: 12.935-13.446 Ma, SR=778 cm/Ma, overlap = bottom 133 cm
# SC: 13.275-14.298 Ma, SR=414 cm/Ma, overlap = top 71 cm

library(tidyverse)
library(patchwork)

# Load data
df <- read_csv("output/tables/xrf_data_stacked.csv", show_col_types=FALSE)

# Calculate Mn/Ti ratio if not present
if(!"Mn_Ti" %in% names(df)) {
  df <- df %>% mutate(Mn_Ti = Mn / Ti)
}

# Age model parameters
tam_age_top <- 12.935  # Ma
tam_age_bottom <- 13.446  # Ma
tam_thickness <- 3974  # mm (397.4 cm)
tam_sed_rate <- 778  # cm/Ma = 7.78 mm/kyr

sc_age_top <- 13.275  # Ma
sc_age_bottom <- 14.298  # Ma
sc_thickness <- 4235  # mm (423.5 cm)
sc_sed_rate <- 414  # cm/Ma = 4.14 mm/kyr

# Overlap period
overlap_start <- 13.275  # Ma (younger end)
overlap_end <- 13.446  # Ma (older end)
overlap_duration <- overlap_end - overlap_start  # 0.171 Myr

# Overlap depths
tam_overlap_depth <- overlap_duration * 778 * 10  # mm from bottom = 1330 mm
sc_overlap_depth <- overlap_duration * 414 * 10   # mm from top = 708 mm

cat("Temporal Overlap Analysis\n")
cat("========================\n")
cat(sprintf("Overlap period: %.3f - %.3f Ma (%.0f kyr)\n", overlap_start, overlap_end, overlap_duration*1000))
cat(sprintf("TAM overlap: bottom %.0f mm of core\n", tam_overlap_depth))
cat(sprintf("SC overlap: top %.0f mm of core\n", sc_overlap_depth))

# Get cumulative depths for each site
# For TAM: need to identify which sections are at the bottom
# For SC: need to identify which sections are at the top

# Get section depth ranges
section_info <- df %>%
  group_by(core_series, section) %>%
  summarise(
    min_depth = min(depth_mm, na.rm=TRUE),
    max_depth = max(depth_mm, na.rm=TRUE),
    n = n(),
    mean_Mn_Ti = mean(Mn_Ti, na.rm=TRUE),
    mean_Ca_Ti = mean(Ca_Ti, na.rm=TRUE),
    mean_Fe_Mn = mean(Fe_Mn, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  arrange(core_series, section)

print(section_info, n=30)

# Since we don't have explicit stratigraphic ordering, we'll use section naming conventions
# and focus on the clearest comparison: representative sections from each site

# For the overlap comparison, let's select:
# - TAM sections that might be near the base (look at naming)
# - SC sections that might be near the top (look at naming)

# Based on typical core naming, lower numbers may be near top/younger
# SC-1ABC... likely near top (youngest SC)
# TAM-5AB... or TAM-1... structure varies

# Let's analyze all sections and identify best candidates for comparison
cat("\n\nSection Statistics:\n")
section_info %>%
  arrange(core_series, desc(mean_Mn_Ti)) %>%
  print(n=30)

# Create comparison figures for publication

# Figure 1: Site-level Mn/Ti distribution comparison
fig1 <- df %>%
  filter(!is.na(Mn_Ti), Mn_Ti < 5) %>%  # Remove extreme outliers
  ggplot(aes(x = Mn_Ti, fill = core_series, color = core_series)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_rug(alpha = 0.2, sides = "b") +
  scale_fill_manual(values = c("TAM" = "#D55E00", "SC" = "#0072B2"), name = "Site") +
  scale_color_manual(values = c("TAM" = "#D55E00", "SC" = "#0072B2"), name = "Site") +
  labs(
    x = "Mn/Ti ratio",
    y = "Density",
    title = "Mn/Ti Distribution by Site",
    subtitle = "TAM: 12.9-13.4 Ma | SC: 13.3-14.3 Ma | Overlap: 13.3-13.4 Ma (171 kyr)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80"),
    panel.grid.minor = element_blank()
  ) +
  annotate("text", x = 0.18, y = Inf, label = "TAM median\n0.18",
           vjust = 2, hjust = 0.5, size = 3, color = "#D55E00") +
  annotate("text", x = 0.37, y = Inf, label = "SC median\n0.37",
           vjust = 2, hjust = 0.5, size = 3, color = "#0072B2")

ggsave("output/figures/temporal_overlap/fig1_MnTi_density.png", fig1,
       width = 8, height = 6, dpi = 300)
ggsave("output/figures/temporal_overlap/fig1_MnTi_density.pdf", fig1,
       width = 8, height = 6)

# Figure 2: Temporal framework diagram
# Create age-depth model visualization
age_depth <- tibble(
  site = c("TAM", "TAM", "SC", "SC"),
  age = c(tam_age_top, tam_age_bottom, sc_age_top, sc_age_bottom),
  depth = c(0, tam_thickness, 0, sc_thickness)
)

fig2 <- ggplot() +
  # TAM line
  geom_segment(aes(x = tam_age_top, xend = tam_age_bottom, y = 0, yend = tam_thickness),
               color = "#D55E00", linewidth = 2) +
  # SC line
  geom_segment(aes(x = sc_age_top, xend = sc_age_bottom, y = 0, yend = sc_thickness),
               color = "#0072B2", linewidth = 2) +
  # Overlap region
  annotate("rect", xmin = overlap_start, xmax = overlap_end,
           ymin = -200, ymax = 4500, fill = "gray80", alpha = 0.5) +
  # Labels
  annotate("text", x = tam_age_top + 0.1, y = tam_thickness/2,
           label = "TAM\n778 cm/Ma", color = "#D55E00", hjust = 0, size = 4) +
  annotate("text", x = sc_age_bottom - 0.1, y = sc_thickness/2,
           label = "SC\n414 cm/Ma", color = "#0072B2", hjust = 1, size = 4) +
  annotate("text", x = (overlap_start + overlap_end)/2, y = 4300,
           label = "OVERLAP\n171 kyr", hjust = 0.5, size = 3.5, fontface = "bold") +
  scale_x_reverse(limits = c(14.5, 12.7), breaks = seq(13, 14.5, 0.25)) +
  scale_y_continuous(limits = c(-200, 4500)) +
  labs(
    x = "Age (Ma)",
    y = "Depth (mm)",
    title = "Core Chronology: TAM vs SC",
    subtitle = "Shaded region = temporal overlap (13.275-13.446 Ma)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("output/figures/temporal_overlap/fig2_age_depth_model.png", fig2,
       width = 8, height = 6, dpi = 300)
ggsave("output/figures/temporal_overlap/fig2_age_depth_model.pdf", fig2,
       width = 8, height = 6)

# Figure 3: Case study - Dramatic oxygenation event (SC-1ABC-2-3C-A-RUN1)
case1_data <- df %>%
  filter(section == "SC-1ABC-2-3C-A-RUN1", !is.na(Mn_Ti))

fig3 <- case1_data %>%
  ggplot(aes(x = depth_mm, y = Mn_Ti)) +
  geom_line(color = "#0072B2", linewidth = 0.8) +
  geom_point(color = "#0072B2", size = 1, alpha = 0.5) +
  geom_hline(yintercept = 0.37, linetype = "dashed", color = "gray50") +
  annotate("text", x = max(case1_data$depth_mm), y = 0.37,
           label = "SC median", hjust = 1.1, vjust = -0.5, size = 3) +
  # Highlight the dramatic transition
  annotate("rect", xmin = 350, xmax = 380, ymin = 0, ymax = max(case1_data$Mn_Ti, na.rm=TRUE),
           fill = "yellow", alpha = 0.3) +
  annotate("text", x = 365, y = max(case1_data$Mn_Ti, na.rm=TRUE) * 0.9,
           label = "+833%\ntransition", hjust = 0.5, size = 3, fontface = "bold") +
  labs(
    x = "Depth (mm)",
    y = "Mn/Ti ratio",
    title = "Case Study 1: Dramatic Oxygenation Event",
    subtitle = "SC-1ABC-2-3C-A-RUN1 | Most extreme Mn/Ti transition in dataset"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  coord_flip()

ggsave("output/figures/temporal_overlap/fig3_case_oxygenation.png", fig3,
       width = 8, height = 6, dpi = 300)
ggsave("output/figures/temporal_overlap/fig3_case_oxygenation.pdf", fig3,
       width = 8, height = 6)

# Figure 4: Case study - Dysoxic baseline (TAM-1-2-3B-A)
case2_data <- df %>%
  filter(section == "TAM-1-2-3B-A", !is.na(Mn_Ti))

fig4 <- case2_data %>%
  ggplot(aes(x = depth_mm)) +
  geom_line(aes(y = Mn_Ti), color = "#D55E00", linewidth = 0.8) +
  geom_hline(yintercept = 0.18, linetype = "dashed", color = "gray50") +
  annotate("text", x = max(case2_data$depth_mm), y = 0.18,
           label = "TAM median", hjust = 1.1, vjust = -0.5, size = 3) +
  # Highlight dysoxic interval
  annotate("rect", xmin = 243, xmax = 444, ymin = 0, ymax = 0.15,
           fill = "#D55E00", alpha = 0.2) +
  annotate("text", x = 343, y = 0.05,
           label = "Extreme dysoxia\nMn/Ti = 0.078", hjust = 0.5, size = 3) +
  labs(
    x = "Depth (mm)",
    y = "Mn/Ti ratio",
    title = "Case Study 2: Persistent Dysoxic Baseline",
    subtitle = "TAM-1-2-3B-A | Pachydon-dominated interval with extreme reducing conditions"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  coord_flip()

ggsave("output/figures/temporal_overlap/fig4_case_dysoxia.png", fig4,
       width = 8, height = 6, dpi = 300)
ggsave("output/figures/temporal_overlap/fig4_case_dysoxia.pdf", fig4,
       width = 8, height = 6)

# Figure 5: Concordant proxy event (dual proxy)
case3_data <- df %>%
  filter(section == "SC-3AB-4ABCD-RUN2-F", !is.na(Mn_Ti), !is.na(Ca_Ti))

fig5a <- case3_data %>%
  ggplot(aes(x = depth_mm, y = Mn_Ti)) +
  geom_line(color = "#0072B2", linewidth = 0.8) +
  geom_point(color = "#0072B2", size = 1, alpha = 0.5) +
  annotate("rect", xmin = 1240, xmax = 1280, ymin = 0, ymax = max(case3_data$Mn_Ti, na.rm=TRUE),
           fill = "yellow", alpha = 0.3) +
  labs(x = "", y = "Mn/Ti", title = "Mn/Ti (redox proxy)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank()) +
  coord_flip()

fig5b <- case3_data %>%
  ggplot(aes(x = depth_mm, y = Ca_Ti)) +
  geom_line(color = "#009E73", linewidth = 0.8) +
  geom_point(color = "#009E73", size = 1, alpha = 0.5) +
  annotate("rect", xmin = 1240, xmax = 1280, ymin = 0, ymax = max(case3_data$Ca_Ti, na.rm=TRUE),
           fill = "yellow", alpha = 0.3) +
  labs(x = "Depth (mm)", y = "Ca/Ti", title = "Ca/Ti (carbonate proxy)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank()) +
  coord_flip()

fig5 <- fig5a / fig5b +
  plot_annotation(
    title = "Case Study 3: Concordant Oxic-Carbonate Event",
    subtitle = "SC-3AB-4ABCD-RUN2-F @ 1258mm | Independent proxies move together (+5.8σ Mn/Ti, +2.0σ Ca/Ti)"
  )

ggsave("output/figures/temporal_overlap/fig5_case_concordant.png", fig5,
       width = 8, height = 8, dpi = 300)
ggsave("output/figures/temporal_overlap/fig5_case_concordant.pdf", fig5,
       width = 8, height = 8)

# Figure 6: Box plot comparison by site
fig6 <- df %>%
  filter(!is.na(Mn_Ti), Mn_Ti < 5) %>%
  ggplot(aes(x = core_series, y = Mn_Ti, fill = core_series)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(alpha = 0.1, width = 0.2, size = 0.5) +
  scale_fill_manual(values = c("TAM" = "#D55E00", "SC" = "#0072B2")) +
  labs(
    x = "Site",
    y = "Mn/Ti ratio",
    title = "Mn/Ti Distribution by Site",
    subtitle = "TAM (younger, 12.9-13.4 Ma) vs SC (older, 13.3-14.3 Ma)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  stat_summary(fun = median, geom = "text",
               aes(label = sprintf("median = %.2f", after_stat(y))),
               vjust = -1.5, size = 3.5)

ggsave("output/figures/temporal_overlap/fig6_boxplot_comparison.png", fig6,
       width = 6, height = 6, dpi = 300)
ggsave("output/figures/temporal_overlap/fig6_boxplot_comparison.pdf", fig6,
       width = 6, height = 6)

# Figure 7: Environmental modes scatter
fig7 <- df %>%
  filter(!is.na(Mn_Ti), !is.na(Ca_Ti), Mn_Ti < 5, Ca_Ti < 60) %>%
  ggplot(aes(x = Mn_Ti, y = Ca_Ti, color = core_series)) +
  geom_point(alpha = 0.3, size = 1) +
  scale_color_manual(values = c("TAM" = "#D55E00", "SC" = "#0072B2"), name = "Site") +
  # Add mode labels
  annotate("text", x = 1.5, y = 40, label = "Mode 1\nOxic + Carbonate",
           hjust = 0.5, size = 3, fontface = "bold") +
  annotate("text", x = 0.1, y = 5, label = "Mode 2\nDysoxic Baseline",
           hjust = 0, size = 3, fontface = "bold") +
  annotate("text", x = 0.1, y = 40, label = "Mode 3\nShells + Reducing",
           hjust = 0, size = 3, fontface = "bold") +
  labs(
    x = "Mn/Ti (redox proxy)",
    y = "Ca/Ti (carbonate proxy)",
    title = "Environmental Modes: Mn/Ti vs Ca/Ti",
    subtitle = "Three distinct paleoenvironmental states in the Pebas system"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = c(0.9, 0.15),
    legend.background = element_rect(fill = "white", color = "gray80"),
    panel.grid.minor = element_blank()
  )

ggsave("output/figures/temporal_overlap/fig7_environmental_modes.png", fig7,
       width = 8, height = 6, dpi = 300)
ggsave("output/figures/temporal_overlap/fig7_environmental_modes.pdf", fig7,
       width = 8, height = 6)

cat("\n\nFigures saved to output/figures/temporal_overlap/\n")
cat("Files generated:\n")
list.files("output/figures/temporal_overlap/", pattern = "fig")
