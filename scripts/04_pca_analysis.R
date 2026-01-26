# ==============================================================================
# Pebas-XRF: PCA Analysis with CLR Transformation
# ==============================================================================
# Compositional data analysis using centered log-ratio transformation
# TAM = Tamshiyacu, SC = Santa Corina
# ==============================================================================

library(tidyverse)
library(compositions)  # For log-ratio transformations
library(patchwork)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

base_path <- here::here()
output_path <- file.path(base_path, "output")

xrf <- read_csv(file.path(output_path, "tables", "xrf_data_ratios.csv"),
                show_col_types = FALSE)

message(sprintf("Loaded %d measurements", nrow(xrf)))

# ==============================================================================
# 2. PREPARE COMPOSITIONAL DATA
# ==============================================================================

# Elements for PCA (well-detected with Mo tube, geochemically meaningful)
pca_elements <- c("Al", "Si", "K", "Ca", "Ti", "Fe", "Mn", "Rb", "Sr", "Zr")

# Check which elements are available
available_elements <- intersect(pca_elements, names(xrf))
message(sprintf("Using %d elements for PCA: %s",
                length(available_elements),
                paste(available_elements, collapse = ", ")))

# Extract element matrix (remove zeros and NAs for log-ratio)
xrf_elements <- xrf %>%
  select(all_of(available_elements)) %>%
  # Replace zeros with small value (detection limit proxy)
  mutate(across(everything(), ~ifelse(. <= 0, 0.1, .))) %>%
  # Remove rows with any NA
  drop_na()

# Get metadata for non-NA rows
xrf_meta <- xrf %>%
  select(group, core_series, section, position_mm, depth_mm) %>%
  slice(which(complete.cases(xrf[, available_elements])))

message(sprintf("PCA dataset: %d complete observations", nrow(xrf_elements)))

# ==============================================================================
# 3. CENTERED LOG-RATIO (CLR) TRANSFORMATION
# ==============================================================================

message("\nApplying CLR transformation...")

# Convert to compositional data object
xrf_comp <- acomp(xrf_elements)

# Apply CLR transformation
xrf_clr <- clr(xrf_comp)

# Convert back to data frame
xrf_clr_df <- as.data.frame(xrf_clr)
names(xrf_clr_df) <- paste0(names(xrf_clr_df), "_clr")

# Combine with metadata
xrf_pca_data <- bind_cols(xrf_meta, xrf_clr_df)

message("CLR transformation complete")

# ==============================================================================
# 4. PRINCIPAL COMPONENT ANALYSIS
# ==============================================================================

message("\nRunning PCA...")

# PCA on CLR-transformed data
pca_result <- prcomp(xrf_clr, scale. = TRUE, center = TRUE)

# Summary
pca_summary <- summary(pca_result)
print(pca_summary)

# Variance explained
var_explained <- pca_summary$importance[2, ] * 100
cum_var <- pca_summary$importance[3, ] * 100

var_df <- tibble(
  PC = paste0("PC", 1:length(var_explained)),
  Variance = var_explained,
  Cumulative = cum_var
)

print(var_df)

# ==============================================================================
# 5. EXTRACT PCA SCORES AND LOADINGS
# ==============================================================================

# Scores (sample positions in PC space)
scores <- as.data.frame(pca_result$x)
names(scores) <- paste0("PC", 1:ncol(scores))

# Combine scores with metadata
xrf_scores <- bind_cols(xrf_meta, scores)

# Loadings (variable contributions to PCs)
loadings <- as.data.frame(pca_result$rotation)
loadings$element <- rownames(loadings)

# Save results
write_csv(xrf_scores, file.path(output_path, "tables", "pca_scores.csv"))
write_csv(loadings, file.path(output_path, "tables", "pca_loadings.csv"))
write_csv(var_df, file.path(output_path, "tables", "pca_variance.csv"))

# ==============================================================================
# 6. PCA VISUALIZATION
# ==============================================================================

message("\nGenerating PCA plots...")

# --- Scree plot ---
p_scree <- ggplot(var_df[1:8, ], aes(x = fct_inorder(PC), y = Variance)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_line(aes(group = 1), color = "darkred", linewidth = 1) +
  geom_point(color = "darkred", size = 3) +
  geom_text(aes(label = sprintf("%.1f%%", Variance)),
            vjust = -0.5, size = 3) +
  labs(
    title = "PCA Scree Plot",
    subtitle = "Variance explained by each principal component",
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

# --- Biplot (PC1 vs PC2) ---
# Scale loadings for visualization
loading_scale <- 4

loadings_plot <- loadings %>%
  mutate(
    PC1_scaled = PC1 * loading_scale,
    PC2_scaled = PC2 * loading_scale
  )

p_biplot <- ggplot() +
  # Sample scores
  geom_point(data = xrf_scores,
             aes(x = PC1, y = PC2, color = core_series),
             alpha = 0.5, size = 1.5) +
  # Loading arrows
  geom_segment(data = loadings_plot,
               aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.8) +
  # Loading labels
  geom_text(data = loadings_plot,
            aes(x = PC1_scaled * 1.15, y = PC2_scaled * 1.15, label = element),
            fontface = "bold", size = 3.5) +
  scale_color_manual(
    values = c("TAM" = "#2166ac", "SC" = "#b2182b"),
    labels = c("TAM" = "Tamshiyacu", "SC" = "Santa Corina")
  ) +
  labs(
    title = "PCA Biplot: Element Composition",
    subtitle = sprintf("PC1 (%.1f%%) vs PC2 (%.1f%%)",
                       var_explained[1], var_explained[2]),
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    color = "Core Series"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# --- PC scores by section ---
p_pc1_sections <- xrf_scores %>%
  ggplot(aes(x = PC1, y = position_mm, color = core_series)) +
  geom_path(linewidth = 0.5) +
  scale_y_reverse() +
  facet_wrap(~section, scales = "free", ncol = 5) +
  scale_color_manual(
    values = c("TAM" = "#2166ac", "SC" = "#b2182b"),
    labels = c("TAM" = "Tamshiyacu", "SC" = "Santa Corina")
  ) +
  labs(
    title = "PC1 Stratigraphic Profiles by Section",
    x = "PC1 Score",
    y = "Position (mm)",
    color = "Core Series"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 7),
    legend.position = "bottom"
  )

# --- Combined figure ---
p_combined <- (p_scree | p_biplot) / p_pc1_sections +
  plot_layout(heights = c(1, 1.5))

ggsave(file.path(output_path, "figures", "pca_analysis.png"),
       p_combined, width = 14, height = 12, dpi = 150, bg = "white")

# Save individual plots
ggsave(file.path(output_path, "figures", "pca_biplot.png"),
       p_biplot, width = 10, height = 8, dpi = 150, bg = "white")

ggsave(file.path(output_path, "figures", "pca_scree.png"),
       p_scree, width = 8, height = 5, dpi = 150, bg = "white")

# ==============================================================================
# 7. LOADING INTERPRETATION
# ==============================================================================

message("\n" , paste(rep("=", 60), collapse = ""))
message("PCA LOADING INTERPRETATION")
message(paste(rep("=", 60), collapse = ""))

# PC1 interpretation
pc1_loadings <- loadings %>%
  select(element, PC1) %>%
  arrange(desc(abs(PC1)))

message("\nPC1 Loadings (strongest to weakest):")
print(pc1_loadings)

# PC2 interpretation
pc2_loadings <- loadings %>%
  select(element, PC2) %>%
  arrange(desc(abs(PC2)))

message("\nPC2 Loadings (strongest to weakest):")
print(pc2_loadings)

# Interpretation guide
message("\n--- INTERPRETATION GUIDE ---")
message("PC1 likely represents: [Interpret based on loading pattern]")
message("  - If Fe, Ti, K, Al load positively: Terrigenous/detrital input")
message("  - If Ca, Sr load oppositely: Carbonate vs siliciclastic")
message("")
message("PC2 likely represents: [Interpret based on loading pattern]")
message("  - If Mn loads strongly: Redox conditions")
message("  - If Zr/Rb separate: Grain size/energy")

# ==============================================================================
# 8. CLUSTER ANALYSIS (OPTIONAL)
# ==============================================================================

message("\nRunning hierarchical clustering on PC scores...")

# Use first 3 PCs for clustering
pc_for_cluster <- xrf_scores %>%
  select(PC1, PC2, PC3) %>%
  as.matrix()

# Hierarchical clustering
hc <- hclust(dist(pc_for_cluster), method = "ward.D2")

# Cut tree into groups
xrf_scores$cluster <- cutree(hc, k = 4)

# Cluster summary
cluster_summary <- xrf_scores %>%
  group_by(cluster, core_series) %>%
  summarise(
    n = n(),
    PC1_mean = mean(PC1),
    PC2_mean = mean(PC2),
    .groups = "drop"
  )

print(cluster_summary)

# Save clustered data
write_csv(xrf_scores, file.path(output_path, "tables", "pca_scores_clustered.csv"))

message("\n", paste(rep("=", 60), collapse = ""))
message("PCA ANALYSIS COMPLETE")
message(paste(rep("=", 60), collapse = ""))
message(sprintf("\nOutput files saved to: %s", output_path))
message("  - pca_scores.csv: Sample PC scores")
message("  - pca_loadings.csv: Variable loadings")
message("  - pca_variance.csv: Variance explained")
message("  - pca_analysis.png: Combined figure")
