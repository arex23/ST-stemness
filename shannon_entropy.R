calculate_shannon_entropy <- function(seurat_obj, assay = NULL, layer = "counts") {
  
  if (is.null(assay)) assay <- Seurat::DefaultAssay(seurat_obj)
  
  expr_mat <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  
  if (is.null(expr_mat) || nrow(expr_mat) == 0 || ncol(expr_mat) == 0)
    stop(sprintf("No data found in assay '%s', layer '%s'.", assay, layer))
  
  col_sums <- Matrix::colSums(expr_mat)
  
  if (inherits(expr_mat, "dgCMatrix")) {
    counts_per_cell <- diff(expr_mat@p)
    p <- expr_mat@x / rep(col_sums, counts_per_cell)
    entropy_components <- -(p * log2(p))
    col_indices <- rep(seq_len(ncol(expr_mat)), counts_per_cell)
    cell_entropies <- tapply(entropy_components, col_indices, sum)
    
    full_entropy <- numeric(ncol(expr_mat))
    names(full_entropy) <- colnames(expr_mat)
    full_entropy[as.numeric(names(cell_entropies))] <- cell_entropies
  } else {
    full_entropy <- apply(expr_mat, 2, function(x) {
      x <- x[x > 0]
      if (length(x) == 0) return(0)
      p <- x / sum(x)
      -sum(p * log2(p))
    })
  }
  
  Seurat::AddMetaData(seurat_obj, metadata = full_entropy, col.name = "shannon_entropy")
}
