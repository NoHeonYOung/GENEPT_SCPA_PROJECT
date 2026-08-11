#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve script location")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
source(file.path(dirname(script_path), "scpa_core_adapter.R"))

argument <- function(name, required = TRUE, default = NULL) {
  index <- which(args == name)
  if (length(index) == 1L && index < length(args)) return(args[[index + 1L]])
  prefix <- paste0(name, "=")
  inline <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(inline) == 1L) return(sub(paste0("^", prefix), "", inline[[1]]))
  if (required) stop("Missing argument: ", name)
  default
}

input_h5 <- normalizePath(argument("--input-h5"))
manifest_path <- normalizePath(argument("--manifest"))
output_csv <- argument("--output-csv")
output_json <- argument("--output-json")
run_l2 <- identical(argument("--run-l2-sensitivity", required = FALSE, default = "false"), "true")

required <- c("rhdf5", "jsonlite", "Matrix", "SCPA", "multicross")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) stop("Missing Phase 4 R packages: ", paste(missing, collapse = ", "))

manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
pathways <- manifest$pathways
if (length(pathways) == 0L) stop("No eligible pathways in manifest")
gene_names <- as.character(rhdf5::h5read(input_h5, "gene_names"))
# h5py and R use opposite in-memory dimension order for 2D datasets.
cd4 <- t(rhdf5::h5read(input_h5, "expression/cd4"))
cd8 <- t(rhdf5::h5read(input_h5, "expression/cd8"))
embeddings <- t(rhdf5::h5read(input_h5, "embeddings"))
if (ncol(cd4) != length(gene_names) || ncol(cd8) != length(gene_names) ||
    nrow(embeddings) != length(gene_names)) stop("HDF5 gene axes are not aligned")

effective_rank <- function(singular_values, rows, columns) {
  if (length(singular_values) == 0L || max(singular_values) == 0) return(0L)
  tolerance <- max(rows, columns) * max(singular_values) * .Machine$double.eps
  as.integer(sum(singular_values > tolerance))
}

summarize_singular <- function(values) {
  list(
    max = unname(max(values)),
    median = unname(stats::median(values)),
    min = unname(min(values)),
    min_nonzero = if (any(values > 0)) unname(min(values[values > 0])) else 0
  )
}

rows <- vector("list", length(pathways))
rank_qc <- vector("list", length(pathways))
warnings <- character()
l2_rows <- if (run_l2) vector("list", length(pathways)) else NULL
core_started <- proc.time()[["elapsed"]]
cat(
  sprintf(
    "[Phase4 core] Starting %d pathways; branches/pathway=%d\n",
    length(pathways), if (run_l2) 3L else 2L
  )
)
flush.console()

for (i in seq_along(pathways)) {
  pathway <- pathways[[i]]
  indices <- as.integer(unlist(pathway$global_gene_indices)) + 1L
  expected_genes <- as.character(unlist(pathway$paired_genes))
  if (!identical(gene_names[indices], expected_genes)) {
    stop("Gene-order assertion failed for pathway: ", pathway$pathway)
  }
  x4 <- cd4[, indices, drop = FALSE]
  x8 <- cd8[, indices, drop = FALSE]
  ep <- embeddings[indices, , drop = FALSE]
  z4 <- x4 %*% ep
  z8 <- x8 %*% ep
  prefix <- sprintf("[Phase4 core %03d/%03d] %s", i, length(pathways), pathway$pathway)
  cat(prefix, " | Vanilla MCM...\n", sep = "")
  flush.console()
  vanilla <- run_mcm_raw(x4, x8)
  cat(prefix, " | GenePT non-L2 MCM...\n", sep = "")
  flush.console()
  genept <- run_mcm_raw(z4, z8)
  if (length(vanilla$warnings) > 0L) {
    warnings <- c(warnings, paste0(pathway$pathway, ":vanilla:", vanilla$warnings))
  }
  if (length(genept$warnings) > 0L) {
    warnings <- c(warnings, paste0(pathway$pathway, ":genept:", genept$warnings))
  }

  embedding_svd <- svd(ep, nu = min(nrow(ep), ncol(ep)), nv = 0L)
  embedding_sv <- embedding_svd$d
  embedding_rank <- effective_rank(embedding_sv, nrow(ep), ncol(ep))
  # E = U D V'; rank(XE) equals rank(X U D) over non-zero singular directions.
  nonzero <- seq_len(embedding_rank)
  reduced_projected <- rbind(x4, x8) %*%
    embedding_svd$u[, nonzero, drop = FALSE] %*%
    diag(embedding_svd$d[nonzero], nrow = embedding_rank)
  projected_sv <- svd(reduced_projected, nu = 0L, nv = 0L)$d
  projected_rank <- effective_rank(projected_sv, nrow(reduced_projected), ncol(reduced_projected))

  rows[[i]] <- data.frame(
    pathway = pathway$pathway,
    vanilla_raw_p = vanilla$raw_p,
    genept_raw_p = genept$raw_p,
    embedding_rank = embedding_rank,
    projected_rank = projected_rank,
    stringsAsFactors = FALSE
  )
  rank_qc[[i]] <- list(
    pathway = pathway$pathway,
    pathway_gene_count = length(indices),
    embedding_rank = embedding_rank,
    projected_rank = projected_rank,
    rank_deficiency_vs_1536 = 1536L - projected_rank,
    embedding_singular_values = summarize_singular(embedding_sv),
    projected_singular_values = summarize_singular(projected_sv)
  )
  if (run_l2) {
    row_l2 <- function(x) {
      norms <- sqrt(rowSums(x * x))
      nonzero_rows <- norms > 0
      x[nonzero_rows, ] <- x[nonzero_rows, , drop = FALSE] / norms[nonzero_rows]
      x
    }
    cat(prefix, " | GenePT L2 sensitivity MCM...\n", sep = "")
    flush.console()
    l2_result <- run_mcm_raw(row_l2(z4), row_l2(z8))
    l2_rows[[i]] <- data.frame(pathway = pathway$pathway, l2_raw_p = l2_result$raw_p)
    if (length(l2_result$warnings) > 0L) {
      warnings <- c(warnings, paste0(pathway$pathway, ":genept_l2:", l2_result$warnings))
    }
  }
  elapsed <- proc.time()[["elapsed"]] - core_started
  eta <- if (i > 0L) elapsed / i * (length(pathways) - i) else NA_real_
  cat(
    sprintf(
      "%s | done | elapsed %.1f min | ETA %.1f min\n",
      prefix, elapsed / 60, eta / 60
    )
  )
  flush.console()
}

results <- do.call(rbind, rows)
vanilla_corrected <- scpa_pathway_qvalues(results$vanilla_raw_p)
genept_corrected <- scpa_pathway_qvalues(results$genept_raw_p)
results$vanilla_adjusted_p <- vanilla_corrected$adjusted_p
results$vanilla_qval <- vanilla_corrected$qval
results$genept_adjusted_p <- genept_corrected$adjusted_p
results$genept_qval <- genept_corrected$qval
rank_tie_aware <- function(qval) {
  rank(-qval, ties.method = "average")
}
results$vanilla_rank <- rank_tie_aware(results$vanilla_qval)
results$genept_rank <- rank_tie_aware(results$genept_qval)
results$rank_delta <- results$genept_rank - results$vanilla_rank
if (run_l2) {
  l2 <- do.call(rbind, l2_rows)
  corrected <- scpa_pathway_qvalues(l2$l2_raw_p)
  l2$l2_adjusted_p <- corrected$adjusted_p
  l2$l2_qval <- corrected$qval
  l2$l2_rank <- rank_tie_aware(l2$l2_qval)
  results <- merge(results, l2, by = "pathway", all.x = TRUE, sort = FALSE)
}

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, output_csv, row.names = FALSE, quote = TRUE)
payload <- list(
  scpa_version = as.character(utils::packageVersion("SCPA")),
  multicross_version = as.character(utils::packageVersion("multicross")),
  implementation = "multicross::mcm via Phase 3 SCPA-core adapter",
  raw_p_source = "multicross::mcm result[[1]]",
  multiple_testing = "stats::p.adjust(method='bonferroni', n=eligible_pathway_count)",
  qval_formula = "sqrt(-log10(adjusted_p))",
  log_base = 10,
  eligible_pathway_count = nrow(results),
  l2_sensitivity_run = run_l2,
  warnings = unname(unique(warnings[nzchar(warnings)])),
  effective_rank = rank_qc
)
dir.create(dirname(output_json), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(payload, output_json, auto_unbox = TRUE, pretty = TRUE, null = "null")
cat("PHASE4_CORE status=PASS pathways=", nrow(results), " output=", output_csv, "\n", sep = "")
