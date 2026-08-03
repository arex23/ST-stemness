library(Seurat)
library(ggplot2)
library(patchwork)

source("shannon_entropy.R")
source("entropy_correlation.R")

# checks data folder to find the sample folders to run the analysis on multiple
samples <- list.dirs("data", full.names = FALSE, recursive = FALSE)
samples <- samples[samples != ""]

if (length(samples) == 0) {
  stop("No samples found in data folder")
}

# function for the analysis of one sample, to be applied to all the folders
analyze_sample <- function(sample_name) {
  print(paste("Analyzing sample:", sample_name))
  
  # save path of the raw_data folder and then check if "spatial" folder exists
  # if it exists it becomes the data directory
  
  path1 <- file.path("data", sample_name, "raw_data")
  path2 <- file.path("data", sample_name)
  
  data_dir <- NULL
  if (dir.exists(path1) && dir.exists(file.path(path1, "spatial"))) {
    data_dir <- path1
    
    # questo elif è un po' inutile perchè controlla se c'è spatial fuori da raw
    # data che non dovrebbe succedere infatti valuterò di cavarlo
  } else if (dir.exists(path2) && dir.exists(file.path(path2, "spatial"))) {
    data_dir <- path2
  }
  
  if (is.null(data_dir)) {
    warning(paste("Data not found for sample:", sample_name))
    return(NULL)
  }
  
  # find the h5 file which contains all the reads
  h5_file <- list.files(path = data_dir, pattern = "\\.h5$", full.names = FALSE)[1]
  
  # load it with Seurat
  spatial_obj <- Load10X_Spatial(
    data.dir = data_dir,
    filename = h5_file,
    assay = "Spatial",
    filter.matrix = TRUE
  )
  
  # calculate entropy per spot (RAW)
  spatial_obj <- calculate_shannon_entropy(spatial_obj, assay = "Spatial", layer = "counts")
  
  # various plot settings
  theme_entropy <- theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
  
  p_raw <- SpatialFeaturePlot(spatial_obj, features = "shannon_entropy") +
    scale_fill_viridis_c(option = "magma", name = "Shannon\nEntropy") +
    ggtitle(paste("Spatial Distribution of Shannon Entropy -", sample_name)) +
    theme_entropy
  ggsave(file.path("results", paste0(sample_name, "_spatial_entropy_plot.png")), plot = p_raw, width = 8, height = 7, dpi = 300)
  
  # entropy correlation function to see if there is correlation between counts 
  # or features and entropy
  
  corr_res <- calculate_entropy_correlations(spatial_obj, sample_name = sample_name, output_dir = "results")
  return(invisible(corr_res))
}

# no idea di tutto sto blocco ma penso sia per applicare l'analisi a più sample
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


