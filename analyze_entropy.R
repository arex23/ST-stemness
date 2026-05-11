library(Seurat)
library(ggplot2)
library(patchwork)

source("shannon_entropy.R")

spatial_obj <- Load10X_Spatial(
  data.dir = "data/raw_data",
  filename = "GSM9372563_MNG_Sample_01_raw_feature_bc_matrix.h5",
  assay = "Spatial",
  filter.matrix = TRUE
)

spatial_obj <- calculate_shannon_entropy(spatial_obj, assay = "Spatial", layer = "counts")
spatial_obj$normalized_entropy <- spatial_obj$shannon_entropy / log2(spatial_obj$nFeature_Spatial)

entropy_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10)
)

p_norm <- SpatialFeaturePlot(spatial_obj, features = "normalized_entropy") +
  scale_fill_viridis_c(option = "magma", name = "Normalized\nEntropy") +
  ggtitle("Spatial Distribution of Normalized Shannon Entropy") +
  entropy_theme

p_raw <- SpatialFeaturePlot(spatial_obj, features = "shannon_entropy") +
  scale_fill_viridis_c(option = "magma", name = "Shannon\nEntropy") +
  ggtitle("Spatial Distribution of Shannon Entropy") +
  entropy_theme

dir.create("results", showWarnings = FALSE)
ggsave("results/spatial_normalized_entropy_plot.png", plot = p_norm, width = 8, height = 7, dpi = 300)
ggsave("results/spatial_entropy_plot.png", plot = p_raw, width = 8, height = 7, dpi = 300)