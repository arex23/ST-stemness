calculate_entropy_correlations <- function(seurat_obj,
                                            entropy_cols = NULL,
                                            count_col = NULL,
                                            feature_col = NULL,
                                            sample_name = NULL,
                                            output_dir = "results",
                                            save_outputs = TRUE) {
  
  meta <- seurat_obj@meta.data
  
  # Auto-detect count column if not specified
  if (is.null(count_col)) {
    possible_counts <- grep("^nCount", colnames(meta), value = TRUE)
    if (length(possible_counts) > 0) {
      count_col <- possible_counts[1]
    } else {
      stop("Could not automatically find an nCounts column (e.g. nCount_Spatial). Please specify 'count_col'.")
    }
  }
  
  # Auto-detect feature column if not specified
  if (is.null(feature_col)) {
    possible_features <- grep("^nFeature", colnames(meta), value = TRUE)
    if (length(possible_features) > 0) {
      feature_col <- possible_features[1]
    } else {
      stop("Could not automatically find an nFeatures column (e.g. nFeature_Spatial). Please specify 'feature_col'.")
    }
  }
  
  # Auto-detect entropy columns if not specified
  if (is.null(entropy_cols)) {
    candidates <- c("shannon_entropy")
    entropy_cols <- candidates[candidates %in% colnames(meta)]
    if (length(entropy_cols) == 0) {
      stop("No entropy metadata columns found in Seurat object. Run calculate_shannon_entropy first.")
    }
  } else {
    entropy_cols <- entropy_cols[entropy_cols %in% colnames(meta)]
    if (length(entropy_cols) == 0) {
      stop("Specified entropy_cols not found in Seurat object metadata.")
    }
  }
  
  sample_prefix <- if (!is.null(sample_name) && nchar(sample_name) > 0) sample_name else "Sample"
  
  summary_rows <- list()
  plot_list <- list()
  
  targets <- list(
    list(col = count_col, label = "nCounts"),
    list(col = feature_col, label = "nFeatures")
  )
  
  for (e_col in entropy_cols) {
    e_label <- if (e_col == "shannon_entropy") {
      "Shannon Entropy"
    } else {
      e_col
    }
    
    for (t in targets) {
      t_col <- t$col
      t_label <- t$label
      
      x <- meta[[t_col]]
      y <- meta[[e_col]]
      
      valid_idx <- !is.na(x) & !is.na(y)
      x_clean <- x[valid_idx]
      y_clean <- y[valid_idx]
      
      # Pearson correlation
      p_test <- cor.test(x_clean, y_clean, method = "pearson")
      p_r <- unname(p_test$estimate)
      p_pval <- p_test$p.value
      
      # Spearman correlation
      s_test <- suppressWarnings(cor.test(x_clean, y_clean, method = "spearman", exact = FALSE))
      s_rho <- unname(s_test$estimate)
      s_pval <- s_test$p.value
      
      format_p <- function(p) {
        if (p < 0.0001) {
          sprintf("%.2e", p)
        } else {
          sprintf("%.4f", p)
        }
      }
      
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Sample = sample_prefix,
        Entropy_Metric = e_label,
        Target_Variable = t_label,
        Pearson_r = round(p_r, 4),
        Pearson_p = format_p(p_pval),
        Spearman_rho = round(s_rho, 4),
        Spearman_p = format_p(s_pval),
        stringsAsFactors = FALSE
      )
      
      anno_text <- sprintf(
        "Pearson r: %.4f (p = %s)\nSpearman \u03c1: %.4f (p = %s)",
        p_r, format_p(p_pval), s_rho, format_p(s_pval)
      )
      
      p_scatter <- ggplot2::ggplot(meta, ggplot2::aes(x = .data[[t_col]], y = .data[[e_col]])) +
        ggplot2::geom_point(alpha = 0.4, color = "#2c3e50", size = 1.2) +
        ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "#e74c3c", fill = "#e74c3c", alpha = 0.2, se = TRUE) +
        ggplot2::labs(
          title = paste(e_label, "vs", t_label),
          subtitle = paste("Sample:", sample_prefix),
          x = paste(t_label, paste0("(", t_col, ")")),
          y = e_label
        ) +
        ggplot2::annotate(
          "label",
          x = Inf, y = -Inf,
          label = anno_text,
          hjust = 1.05, vjust = -0.2,
          size = 3.5,
          fill = ggplot2::alpha("white", 0.85),
          label.padding = grid::unit(0.4, "lines")
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
          plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5, color = "gray30"),
          axis.title = ggplot2::element_text(face = "bold", size = 10),
          panel.grid.minor = ggplot2::element_blank()
        )
      
      plot_list[[paste(e_col, t_label, sep = "_")]] <- p_scatter
    }
  }
  
  summary_table <- do.call(rbind, summary_rows)
  
  cat("\n======================================================================\n")
  cat(sprintf("  Entropy Correlation Summary Table - Sample: %s\n", sample_prefix))
  cat("======================================================================\n")
  if (requireNamespace("knitr", quietly = TRUE)) {
    print(knitr::kable(summary_table, format = "simple"))
  } else {
    print(summary_table)
  }
  cat("======================================================================\n\n")
  
  combined_plot <- NULL
  if (length(plot_list) > 0) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      if (requireNamespace("gridExtra", quietly = TRUE)) {
        table_grob <- gridExtra::tableGrob(
          summary_table,
          rows = NULL,
          theme = gridExtra::ttheme_minimal(
            core = list(fg_params = list(cex = 0.8)),
            colhead = list(fg_params = list(cex = 0.85, fontface = "bold"))
          )
        )
        table_plot <- patchwork::wrap_elements(table_grob)
        
        scatter_grid <- patchwork::wrap_plots(plot_list, ncol = 2)
        combined_plot <- scatter_grid / table_plot + patchwork::plot_layout(heights = c(3, 1))
      } else {
        combined_plot <- patchwork::wrap_plots(plot_list, ncol = 2)
      }
    }
  }
  
  if (save_outputs) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    csv_file <- file.path(output_dir, paste0(sample_prefix, "_entropy_correlations.csv"))
    write.csv(summary_table, csv_file, row.names = FALSE)
    cat(sprintf("Saved correlation summary table to: %s\n", csv_file))
    
    if (!is.null(combined_plot)) {
      plot_file <- file.path(output_dir, paste0(sample_prefix, "_entropy_correlations_plot.png"))
      n_plots <- length(plot_list)
      h <- if (n_plots > 2) 10 else 6
      ggplot2::ggsave(plot_file, plot = combined_plot, width = 10, height = h, dpi = 300)
      cat(sprintf("Saved correlation plots to: %s\n", plot_file))
    }
  }
  
  return(list(
    summary_table = summary_table,
    plots = plot_list,
    combined_plot = combined_plot
  ))
}
