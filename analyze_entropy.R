library(Seurat)
library(ggplot2)
library(patchwork)

source("shannon_entropy.R")

samples <- list.dirs("data", full.names = FALSE, recursive = FALSE)
samples <- samples[samples != ""]

if (length(samples) == 0) {
  stop("No samples found in data folder")
}

analyze_sample <- function(sample_name) {
  print(paste("Analyzing sample:", sample_name))
  
  path1 <- file.path("data", sample_name, "raw_data")
  path2 <- file.path("data", sample_name)
  
  data_dir <- NULL
  if (dir.exists(path1) && dir.exists(file.path(path1, "spatial"))) {
    data_dir <- path1
  } else if (dir.exists(path2) && dir.exists(file.path(path2, "spatial"))) {
    data_dir <- path2
  }
  
  if (is.null(data_dir)) {
    warning(paste("Data not found for sample:", sample_name))
    return(NULL)
  }
  
  h5_file <- list.files(path = data_dir, pattern = "\\.h5$", full.names = FALSE)[1]
  
  spatial_obj <- Load10X_Spatial(
    data.dir = data_dir,
    filename = h5_file,
    assay = "Spatial",
    filter.matrix = TRUE
  )
  
  spatial_obj <- calculate_shannon_entropy(spatial_obj, assay = "Spatial", layer = "counts")
  spatial_obj$normalized_entropy <- spatial_obj$shannon_entropy / log2(spatial_obj$nFeature_Spatial)
  
  theme_entropy <- theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
  
  p_norm <- SpatialFeaturePlot(spatial_obj, features = "normalized_entropy") +
    scale_fill_viridis_c(option = "magma", name = "Normalized\nEntropy") +
    ggtitle(paste("Spatial Distribution of Normalized Shannon Entropy -", sample_name)) +
    theme_entropy
  
  p_raw <- SpatialFeaturePlot(spatial_obj, features = "shannon_entropy") +
    scale_fill_viridis_c(option = "magma", name = "Shannon\nEntropy") +
    ggtitle(paste("Spatial Distribution of Shannon Entropy -", sample_name)) +
    theme_entropy
  
  dir.create("results", showWarnings = FALSE)
  ggsave(file.path("results", paste0(sample_name, "_spatial_normalized_entropy_plot.png")), plot = p_norm, width = 8, height = 7, dpi = 300)
  ggsave(file.path("results", paste0(sample_name, "_spatial_entropy_plot.png")), plot = p_raw, width = 8, height = 7, dpi = 300)
}

args <- commandArgs(trailingOnly = TRUE)
run_samples <- c()

if (length(args) > 0) {
  arg <- args[1]
  if (tolower(arg) == "all") {
    run_samples <- samples
  } else if (arg %in% samples) {
    run_samples <- arg
  } else {
    stop("Invalid sample name")
  }
} else if (interactive()) {
  print("Available samples:")
  print("0: All samples")
  for (i in 1:length(samples)) {
    print(paste0(i, ": ", samples[i]))
  }
  
  choice <- readline("Select sample (number/name/all): ")
  choice <- trimws(choice)
  
  if (choice == "0" || tolower(choice) == "all" || choice == "") {
    run_samples <- samples
  } else {
    idx <- suppressWarnings(as.numeric(choice))
    if (!is.na(idx) && idx >= 1 && idx <= length(samples)) {
      run_samples <- samples[idx]
    } else if (choice %in% samples) {
      run_samples <- choice
    } else {
      stop("Invalid selection")
    }
  }
} else {
  run_samples <- samples
}

for (s in run_samples) {
  analyze_sample(s)
}


