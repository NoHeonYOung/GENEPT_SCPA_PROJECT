# Shared Phase 5/7 gene-masking algebra.
#
# Scientific interpretation is deliberately outside this helper.  These
# functions only encode the frozen matrix operations used by Phase 5 and the
# Phase 7 synthetic benchmark.

masking_score <- function(raw_p, raw_p_clip = 1e-300) {
  -log10(max(as.numeric(raw_p), as.numeric(raw_p_clip)))
}

masking_q_score <- function(adjusted_p) {
  value <- as.numeric(adjusted_p)
  if (!is.finite(value) || value < 0 || value > 1) stop("Invalid adjusted p-value")
  if (value == 0) Inf else sqrt(-log10(value))
}

vanilla_zero_mask_pair <- function(xa, xb, gene_index) {
  if (gene_index < 1L || gene_index > ncol(xa) || ncol(xa) != ncol(xb)) {
    stop("Invalid Vanilla masking gene index or matrix alignment")
  }
  xa_masked <- xa
  xb_masked <- xb
  xa_masked[, gene_index] <- 0
  xb_masked[, gene_index] <- 0
  list(a = xa_masked, b = xb_masked)
}

genept_non_l2_project_pair <- function(xa, xb, embeddings) {
  if (ncol(xa) != nrow(embeddings) || ncol(xb) != nrow(embeddings)) {
    stop("Expression columns must align to GenePT embedding rows")
  }
  list(a = xa %*% embeddings, b = xb %*% embeddings)
}

genept_non_l2_subtraction_mask_pair <- function(
    za, zb, xa, xb, embeddings, gene_index) {
  if (gene_index < 1L || gene_index > ncol(xa) ||
      ncol(xa) != nrow(embeddings) || ncol(xb) != nrow(embeddings)) {
    stop("Invalid GenePT masking gene index or matrix alignment")
  }
  list(
    a = za - tcrossprod(xa[, gene_index], embeddings[gene_index, ]),
    b = zb - tcrossprod(xb[, gene_index], embeddings[gene_index, ])
  )
}
