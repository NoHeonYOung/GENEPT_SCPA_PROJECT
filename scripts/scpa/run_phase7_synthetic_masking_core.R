#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve script location")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]))
source(file.path(dirname(script_path), "scpa_core_adapter.R"))

argument <- function(name, required = TRUE, default = NULL) {
  index <- which(args == name)
  if (length(index) == 1L && index < length(args)) return(args[[index + 1L]])
  if (required) stop("Missing argument: ", name)
  default
}

input_h5 <- normalizePath(argument("--input-h5"))
manifest_path <- normalizePath(argument("--manifest"))
checkpoint_dir <- argument("--checkpoint-dir")
output_json <- argument("--output-json")
max_experiments <- as.integer(argument("--max-experiments", FALSE, "0"))

required <- c("rhdf5", "jsonlite", "SCPA", "multicross")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) stop("Missing Phase 7 R packages: ", paste(missing, collapse = ", "))

manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
if (!isTRUE(manifest$execution_gate$production_scpa_allowed)) {
  stop("Phase 7 production SCPA remains locked in the synthetic manifest")
}
experiments <- manifest$experiments
if (max_experiments > 0L) experiments <- utils::head(experiments, max_experiments)
if (length(experiments) == 0L) stop("No Phase 7 experiments")

raw_clip <- 1e-300
score <- function(p) -log10(max(as.numeric(p), raw_clip))
run_checked <- function(xa, xb, label) {
  result <- tryCatch(run_mcm_raw(xa, xb), error = function(condition) {
    stop(label, ": ", conditionMessage(condition))
  })
  if (!is.finite(result$raw_p) || length(result$warnings) > 0L) {
    stop(label, ": invalid result or warning")
  }
  result$raw_p
}
atomic_csv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".partial.", Sys.getpid())
  utils::write.csv(data, temporary, row.names = FALSE, quote = TRUE)
  if (!file.rename(temporary, path)) stop("Could not commit checkpoint")
}

dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
completed <- 0L
reused <- 0L
total_mcm <- 0L
started <- proc.time()[["elapsed"]]
for (experiment_index in seq_along(experiments)) {
  experiment <- experiments[[experiment_index]]
  destination <- file.path(checkpoint_dir, paste0(experiment$experiment_id, "_masking.csv"))
  if (file.exists(destination)) {
    cached <- utils::read.csv(destination, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(cached) != as.integer(experiment$analysis_gene_count)) {
      stop("Incompatible Phase 7 checkpoint: ", destination)
    }
    completed <- completed + 1L
    reused <- reused + 1L
    total_mcm <- total_mcm + 2L * (nrow(cached) + 1L)
    next
  }
  root <- sub("^/", "", experiment$expression_h5_group)
  genes <- as.character(rhdf5::h5read(input_h5, paste0(root, "/gene_names")))
  xa <- t(rhdf5::h5read(input_h5, paste0(root, "/condition_A/expression")))
  xb <- t(rhdf5::h5read(input_h5, paste0(root, "/condition_B/expression")))
  ep <- t(rhdf5::h5read(input_h5, paste0(root, "/embeddings")))
  expected <- as.character(unlist(experiment$analysis_genes))
  if (!identical(genes, expected) || ncol(xa) != length(genes) || nrow(ep) != length(genes)) {
    stop("Phase 7 HDF5 gene order mismatch: ", experiment$experiment_id)
  }
  za <- xa %*% ep
  zb <- xb %*% ep
  vanilla_full <- run_checked(xa, xb, "Phase 7 Vanilla full")
  genept_full <- run_checked(za, zb, "Phase 7 GenePT full")
  total_mcm <- total_mcm + 2L
  rows <- vector("list", length(genes))
  for (gene_index in seq_along(genes)) {
    xa_masked <- xa; xb_masked <- xb
    xa_masked[, gene_index] <- 0
    xb_masked[, gene_index] <- 0
    vanilla_masked <- run_checked(xa_masked, xb_masked, "Phase 7 Vanilla mask")
    genept_masked <- run_checked(
      za - tcrossprod(xa[, gene_index], ep[gene_index, ]),
      zb - tcrossprod(xb[, gene_index], ep[gene_index, ]),
      "Phase 7 GenePT mask"
    )
    total_mcm <- total_mcm + 2L
    rows[[gene_index]] <- data.frame(
      experiment_id = experiment$experiment_id,
      pathway = experiment$pathway,
      gene = genes[[gene_index]], gene_index = gene_index - 1L,
      vanilla_raw_p_full = vanilla_full, vanilla_raw_p_masked = vanilla_masked,
      vanilla_delta_score = score(vanilla_full) - score(vanilla_masked),
      genept_raw_p_full = genept_full, genept_raw_p_masked = genept_masked,
      genept_delta_score = score(genept_full) - score(genept_masked),
      stringsAsFactors = FALSE
    )
  }
  result <- do.call(rbind, rows)
  result$vanilla_signed_rank <- rank(-result$vanilla_delta_score, ties.method = "average")
  result$genept_signed_rank <- rank(-result$genept_delta_score, ties.method = "average")
  result$vanilla_absolute_rank <- rank(-abs(result$vanilla_delta_score), ties.method = "average")
  result$genept_absolute_rank <- rank(-abs(result$genept_delta_score), ties.method = "average")
  atomic_csv(result, destination)
  completed <- completed + 1L
  elapsed <- proc.time()[["elapsed"]] - started
  cat(sprintf(
    "[Phase7 SCPA %03d/%03d] experiment=%s genes=%d elapsed=%.1fmin\n",
    experiment_index, length(experiments), experiment$experiment_id, length(genes), elapsed / 60
  ))
  flush.console()
}

payload <- list(
  status = "PASS", backend = "multicross::mcm",
  completed_experiments = completed, reused_checkpoints = reused,
  mcm_count = total_mcm, warnings = character(), failed_mcm_calls = 0L,
  ground_truth_read = FALSE, non_l2_genept = TRUE
)
dir.create(dirname(output_json), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(payload, output_json, auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(sprintf("PHASE7_SCPA_CORE status=PASS experiments=%d mcm=%d reused=%d\n", completed, total_mcm, reused))
