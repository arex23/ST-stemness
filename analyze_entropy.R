library(Seurat)
library(ggplot2)
library(patchwork)

source("shannon_entropy.R")
source("entropy_correlation.R")

samples <- list.dirs("data", full.names = FALSE, recursive = FALSE)
samples <- samples[samples != ""]

analyze_sample <- function(sample_name) {
  cat("Analyzing sample:", sample_name, "\n")
  
  data_dir <- file.path("data", sample_name, "raw_data")
  if (!dir.exists(data_dir)) {
    data_dir <- file.path("data", sample_name)
  }
  
  h5_file <- list.files(path = data_dir, pattern = "\\.h5$", full.names = FALSE)[1]
  
  spatial_obj <- Load10X_Spatial(
    data.dir = data_dir,
    filename = h5_file,
    assay = "Spatial",
    filter.matrix = TRUE
  )
  
  spatial_obj <- calculate_shannon_entropy(spatial_obj, assay = "Spatial", layer = "counts")
  
  p_raw <- SpatialFeaturePlot(spatial_obj, features = "shannon_entropy") +
    scale_fill_viridis_c(option = "magma", name = "Shannon\nEntropy") +
    ggtitle(paste("Spatial Distribution of Shannon Entropy -", sample_name)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  if (!dir.exists("results")) dir.create("results", recursive = TRUE)
  ggsave(file.path("results", paste0(sample_name, "_spatial_entropy_plot.png")), plot = p_raw, width = 8, height = 7, dpi = 300)
  
  corr_res <- calculate_entropy_correlations(spatial_obj, sample_name = sample_name, output_dir = "results")
  invisible(corr_res)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  arg <- args[1]
  run_samples <- if (tolower(arg) == "all") samples else arg[arg %in% samples]
} else {
  run_samples <- samples
}

for (s in run_samples) {
  analyze_sample(s)
}
