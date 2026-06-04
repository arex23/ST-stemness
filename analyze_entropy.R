library(Seurat)
library(ggplot2)
library(patchwork)

source("shannon_entropy.R")

# 1. Check data directory
if (!dir.exists("data")) {
  stop("The 'data' directory does not exist.")
}

# 2. Get list of available samples
available_samples <- list.dirs("data", full.names = FALSE, recursive = FALSE)
available_samples <- available_samples[available_samples != ""]

if (length(available_samples) == 0) {
  stop("No sample directories found in the 'data' folder.")
}

# Function to run the entropy analysis for a specific sample
analyze_sample <- function(sample_name) {
  message(paste("\n=== Analyzing sample:", sample_name, "==="))
  
  # Locate data directory and H5 file
  paths_to_check <- c(
    file.path("data", sample_name, "raw_data"),
    file.path("data", sample_name)
  )
  
  data_dir <- NULL
  h5_file <- NULL
  
  for (path in paths_to_check) {
    if (dir.exists(path)) {
      spatial_exists <- dir.exists(file.path(path, "spatial"))
      h5_files <- list.files(path = path, pattern = "\\.h5$", full.names = FALSE)
      
      if (spatial_exists && length(h5_files) > 0) {
        data_dir <- path
        h5_file <- h5_files[1]
        if (length(h5_files) > 1) {
          message(paste("Multiple .h5 files found in", path, "- using the first one:", h5_file))
        }
        break
      }
    }
  }
  
  if (is.null(data_dir)) {
    warning(paste("Could not find valid 10X spatial data for sample:", sample_name, 
                  "\nMake sure the directory contains both a 'spatial' folder and a '.h5' file."))
    return(NULL)
  }
  
  message(paste("Found data in:", data_dir))
  message(paste("Using matrix file:", h5_file))
  
  # Load the spatial data
  spatial_obj <- tryCatch({
    Load10X_Spatial(
      data.dir = data_dir,
      filename = h5_file,
      assay = "Spatial",
      filter.matrix = TRUE
    )
  }, error = function(e) {
    warning(paste("Error loading 10X spatial data for", sample_name, ":", e$message))
    return(NULL)
  })
  
  if (is.null(spatial_obj)) return(NULL)
  
  # Calculate Shannon entropy
  spatial_obj <- calculate_shannon_entropy(spatial_obj, assay = "Spatial", layer = "counts")
  spatial_obj$normalized_entropy <- spatial_obj$shannon_entropy / log2(spatial_obj$nFeature_Spatial)
  
  entropy_theme <- theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
  
  p_norm <- SpatialFeaturePlot(spatial_obj, features = "normalized_entropy") +
    scale_fill_viridis_c(option = "magma", name = "Normalized\nEntropy") +
    ggtitle(paste("Spatial Distribution of Normalized Shannon Entropy -", sample_name)) +
    entropy_theme
  
  p_raw <- SpatialFeaturePlot(spatial_obj, features = "shannon_entropy") +
    scale_fill_viridis_c(option = "magma", name = "Shannon\nEntropy") +
    ggtitle(paste("Spatial Distribution of Shannon Entropy -", sample_name)) +
    entropy_theme
  
  # Create results directory if needed
  dir.create("results", showWarnings = FALSE)
  
  norm_plot_path <- file.path("results", paste0(sample_name, "_spatial_normalized_entropy_plot.png"))
  raw_plot_path <- file.path("results", paste0(sample_name, "_spatial_entropy_plot.png"))
  
  ggsave(norm_plot_path, plot = p_norm, width = 8, height = 7, dpi = 300)
  ggsave(raw_plot_path, plot = p_raw, width = 8, height = 7, dpi = 300)
  
  message(paste("Completed analysis for", sample_name))
  message(paste("Saved plots to:\n  -", norm_plot_path, "\n  -", raw_plot_path))
}

# Determine which sample(s) to analyze
selected_samples <- c()

# Read command line arguments if present
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  arg <- args[1]
  if (tolower(arg) == "all") {
    selected_samples <- available_samples
  } else if (arg %in% available_samples) {
    selected_samples <- arg
  } else {
    stop(paste("Invalid sample specified in command line arguments:", arg, 
               "\nAvailable options: all,", paste(available_samples, collapse = ", ")))
  }
} else if (interactive()) {
  # Interactive mode selection
  cat("Available samples for analysis:\n")
  cat("0: All samples\n")
  for (i in seq_along(available_samples)) {
    cat(paste0(i, ": ", available_samples[i], "\n"))
  }
  
  # Prompt for input
  choice <- readline(prompt = "Select a sample (number or name, or 'all'): ")
  choice <- trimws(choice)
  
  if (tolower(choice) %in% c("all", "0", "")) {
    selected_samples <- available_samples
  } else {
    # Check if number
    num_choice <- suppressWarnings(as.numeric(choice))
    if (!is.na(num_choice) && num_choice >= 1 && num_choice <= length(available_samples)) {
      selected_samples <- available_samples[num_choice]
    } else if (choice %in% available_samples) {
      selected_samples <- choice
    } else {
      stop("Invalid selection. Aborting.")
    }
  }
} else {
  # Non-interactive and no arguments: default to all samples
  message("Non-interactive mode and no command-line arguments provided. Analyzing all samples.")
  selected_samples <- available_samples
}

# Run analysis on selected samples
for (sample in selected_samples) {
  analyze_sample(sample)
}

