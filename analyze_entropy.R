library(Seurat)
library(ggplot2)
library(patchwork)
if (requireNamespace("SpaNorm", quietly = TRUE)) {
  library(SpaNorm)
}

source("shannon_entropy.R")
source("entropy_correlation.R")

samples <- list.dirs("data", full.names = FALSE, recursive = FALSE)
samples <- samples[samples != ""]

apply_spanorm_normalization <- function(seurat_obj, sample_p = 0.1) {
  if (!requireNamespace("SpaNorm", quietly = TRUE)) {
    cat("Warning: SpaNorm is not installed in the active environment. Skipping SpaNorm step.\n")
    return(seurat_obj)
  }
  if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
    cat("Warning: SpatialExperiment is not installed in the active environment. Skipping SpaNorm step.\n")
    return(seurat_obj)
  }
  
  cat("Applying SpaNorm normalization via SpatialExperiment...\n")
  
  counts_mat <- Seurat::GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")
  coords_df <- Seurat::GetTissueCoordinates(seurat_obj)
  
  coord_cols <- intersect(c("x", "y", "imagecol", "imagerow", "row", "col"), colnames(coords_df))
  if (length(coord_cols) < 2) {
    num_cols <- names(which(sapply(coords_df, is.numeric)))
    coord_cols <- num_cols[1:2]
  }
  coords_mat <- as.matrix(coords_df[, coord_cols[1:2]])
  
  common_spots <- intersect(colnames(counts_mat), rownames(coords_mat))
  if (length(common_spots) == 0) {
    if (ncol(counts_mat) == nrow(coords_mat)) {
      rownames(coords_mat) <- colnames(counts_mat)
      common_spots <- colnames(counts_mat)
    } else {
      stop("Could not align spatial coordinates with Seurat object spot barcodes.")
    }
  }
  
  counts_mat <- counts_mat[, common_spots, drop = FALSE]
  coords_mat <- coords_mat[common_spots, , drop = FALSE]
  
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts_mat),
    spatialCoords = coords_mat
  )
  
  # Run SpaNorm with sample.p to reduce memory footprint
  spe <- SpaNorm::SpaNorm(spe, sample.p = sample_p)
  
  avail_assays <- SummarizedExperiment::assayNames(spe)
  target_assay <- if ("logcounts" %in% avail_assays) {
    "logcounts"
  } else if ("normcounts" %in% avail_assays) {
    "normcounts"
  } else {
    avail_assays[length(avail_assays)]
  }
  
  norm_mat <- SummarizedExperiment::assay(spe, target_assay)
  
  tryCatch({
    seurat_obj[["Spatial"]]$data <- norm_mat
  }, error = function(e) {
    seurat_obj <- Seurat::SetAssayData(seurat_obj, assay = "Spatial", slot = "data", new.data = norm_mat)
  })
  
  rm(spe, counts_mat, coords_mat, norm_mat)
  gc()
  
  return(seurat_obj)
}

analyze_sample <- function(sample_name) {
  cat("\n==========================================\n")
  cat("Analyzing sample:", sample_name, "\n")
  cat("==========================================\n")
  
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
  
  cat(sprintf("Original data loaded: %d spots, %d genes\n", ncol(spatial_obj), nrow(spatial_obj)))
  
  # Aggressive QC Filtering
  cat("Applying aggressive QC filtering...\n")
  min_counts <- 500
  min_features <- 250
  valid_spots <- colnames(spatial_obj)[spatial_obj$nCount_Spatial >= min_counts & spatial_obj$nFeature_Spatial >= min_features]
  
  raw_counts <- Seurat::GetAssayData(spatial_obj, assay = "Spatial", layer = "counts")
  # Require a gene to be present in at least 2% of all spots
  min_spots_per_gene <- max(20, ceiling(0.02 * ncol(spatial_obj)))
  expressed_genes <- rownames(raw_counts)[Matrix::rowSums(raw_counts > 0) >= min_spots_per_gene]
  
  spatial_obj <- subset(spatial_obj, cells = valid_spots, features = expressed_genes)
  cat(sprintf("QC complete: retained %d spots and %d genes\n", ncol(spatial_obj), nrow(spatial_obj)))
  
  rm(raw_counts, valid_spots, expressed_genes)
  gc()
  
  # Step 1: Normalize spatial data using SpaNorm
  spatial_obj <- apply_spanorm_normalization(spatial_obj)
  
  # Step 2: Calculate Shannon entropy on normalized data
  spatial_obj <- calculate_shannon_entropy(spatial_obj, assay = "Spatial", layer = "data")
  
  p_raw <- SpatialFeaturePlot(spatial_obj, features = "shannon_entropy") +
    scale_fill_viridis_c(option = "magma", name = "Shannon\nEntropy") +
    ggtitle(paste("Spatial Distribution of Shannon Entropy -", sample_name)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  he_dir <- file.path("results", "he_graphs")
  stat_dir <- file.path("results", "statistical_tests")
  
  if (!dir.exists(he_dir)) dir.create(he_dir, recursive = TRUE)
  if (!dir.exists(stat_dir)) dir.create(stat_dir, recursive = TRUE)
  
  ggsave(file.path(he_dir, paste0(sample_name, "_spatial_entropy_plot.png")), plot = p_raw, width = 8, height = 7, dpi = 300)
  
  corr_res <- calculate_entropy_correlations(spatial_obj, sample_name = sample_name, output_dir = stat_dir)
  
  rm(spatial_obj, p_raw)
  gc()
  
  cat("Analysis complete for sample:", sample_name, "\n")
  invisible(corr_res)
}

# Process one sample at a time via command-line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("Usage: Rscript analyze_entropy.R <sample_name>\n")
  cat("Available samples:", paste(samples, collapse = ", "), "\n")
  cat("No sample specified. Defaulting to first sample:", samples[1], "\n")
  target_sample <- samples[1]
} else {
  arg_sample <- args[1]
  if (arg_sample %in% samples) {
    target_sample <- arg_sample
  } else {
    stop(sprintf("Sample '%s' not found. Available samples: %s", arg_sample, paste(samples, collapse = ", ")))
  }
}

analyze_sample(target_sample)
