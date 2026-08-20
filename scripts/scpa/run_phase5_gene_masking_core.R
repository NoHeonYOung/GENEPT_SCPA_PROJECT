#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve script location")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
source(file.path(dirname(script_path), "scpa_core_adapter.R"))
source(file.path(dirname(script_path), "gene_masking_lib.R"))

argument <- function(name, required = TRUE, default = NULL) {
  index <- which(args == name)
  if (length(index) == 1L && index < length(args)) return(args[[index + 1L]])
  inline <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(inline) == 1L) return(sub(paste0("^", name, "="), "", inline[[1]]))
  if (required) stop("Missing argument: ", name)
  default
}

has_flag <- function(name) any(args == name)

input_h5 <- normalizePath(argument("--input-h5"))
manifest_path <- normalizePath(argument("--manifest"))
checkpoint_dir <- argument("--checkpoint-dir")
output_json <- argument("--output-json")
max_targets <- as.integer(argument("--max-targets", FALSE, "0"))
max_genes <- as.integer(argument("--max-genes", FALSE, "0"))
preflight_only <- has_flag("--preflight-only")

required <- c("rhdf5", "jsonlite", "SCPA", "multicross")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) stop("Missing Phase 5 R packages: ", paste(missing, collapse = ", "))

manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
pathways <- manifest$pathways
targets <- manifest$phase5$targets
if (max_targets > 0L) targets <- utils::head(targets, max_targets)
if (length(targets) == 0L) stop("Empty Phase 5 target universe")

gene_names <- as.character(rhdf5::h5read(input_h5, "gene_names"))
embeddings <- t(rhdf5::h5read(input_h5, "embeddings"))
group_names <- names(manifest$groups)
groups <- setNames(lapply(group_names, function(name) {
  t(rhdf5::h5read(input_h5, paste0("expression/", name)))
}), group_names)
if (nrow(embeddings) != length(gene_names) ||
    any(vapply(groups, ncol, integer(1)) != length(gene_names))) {
  stop("Phase 5 HDF5 gene axes are not aligned")
}

pathway_lookup <- setNames(pathways, vapply(pathways, function(x) x$pathway, character(1)))
safe_p <- function(p) max(as.numeric(p), as.numeric(manifest$phase5$raw_p_clip))
score <- function(p) -log10(safe_p(p))
adjusted <- function(p) min(safe_p(p) * as.integer(manifest$phase5$eligible_pathway_count), 1)
qval <- function(p_adj) masking_q_score(p_adj)

run_checked <- function(xa, xb, label) {
  value <- tryCatch(run_mcm_raw(xa, xb), error = function(condition) {
    stop(label, ": ", conditionMessage(condition))
  })
  if (!is.finite(value$raw_p)) stop(label, ": non-finite raw p")
  if (length(value$warnings) > 0L) {
    stop(label, ": runtime warning: ", paste(value$warnings, collapse = " | "))
  }
  value
}

masking_toy_tests <- function() {
  set.seed(20260810L)
  xa <- matrix(stats::rnorm(40L * 5L), nrow = 40L)
  xb <- matrix(stats::rnorm(40L * 5L, mean = 0.2), nrow = 40L)
  xa_zero <- xa; xb_zero <- xb
  xa_zero[, 3L] <- 0; xb_zero[, 3L] <- 0
  vanilla_distance_difference <- max(abs(
    as.numeric(stats::dist(rbind(xa_zero, xb_zero))) -
      as.numeric(stats::dist(rbind(xa[, -3L], xb[, -3L])))
  ))
  vanilla_zero <- run_checked(xa_zero, xb_zero, "toy vanilla zero mask")
  vanilla_remove <- run_checked(xa[, -3L], xb[, -3L], "toy vanilla physical removal")
  embedding <- matrix(stats::rnorm(5L * 7L), nrow = 5L)
  za <- xa %*% embedding
  zb <- xb %*% embedding
  za_subtract <- za - tcrossprod(xa[, 3L], embedding[3L, ])
  zb_subtract <- zb - tcrossprod(xb[, 3L], embedding[3L, ])
  za_direct <- xa[, -3L, drop = FALSE] %*% embedding[-3L, , drop = FALSE]
  zb_direct <- xb[, -3L, drop = FALSE] %*% embedding[-3L, , drop = FALSE]
  genept_matrix_difference <- max(abs(c(za_subtract - za_direct, zb_subtract - zb_direct)))
  list(
    vanilla_zero_vs_removal_max_distance_difference = vanilla_distance_difference,
    vanilla_zero_vs_removal_raw_p_difference = abs(vanilla_zero$raw_p - vanilla_remove$raw_p),
    vanilla_zero_vs_removal_pass = vanilla_distance_difference <= 1e-12 &&
      abs(vanilla_zero$raw_p - vanilla_remove$raw_p) <= 1e-12,
    genept_subtraction_vs_direct_max_difference = genept_matrix_difference,
    genept_subtraction_vs_direct_pass = genept_matrix_difference <= 1e-10
  )
}

baseline_check <- function(target) {
  pathway <- pathway_lookup[[target$pathway]]
  indices <- as.integer(unlist(pathway$global_gene_indices)) + 1L
  expected_genes <- as.character(unlist(pathway$paired_genes))
  if (!identical(gene_names[indices], expected_genes)) {
    stop("Baseline gene order mismatch: ", target$comparison, ":", target$pathway)
  }
  xa <- groups[[target$group_a]][, indices, drop = FALSE]
  xb <- groups[[target$group_b]][, indices, drop = FALSE]
  ep <- embeddings[indices, , drop = FALSE]
  vanilla <- run_checked(xa, xb, "baseline vanilla")
  genept <- run_checked(xa %*% ep, xb %*% ep, "baseline genept")
  list(
    comparison = target$comparison,
    pathway = target$pathway,
    vanilla_expected = as.numeric(target$vanilla_raw_p),
    vanilla_observed = vanilla$raw_p,
    vanilla_absolute_difference = abs(vanilla$raw_p - as.numeric(target$vanilla_raw_p)),
    genept_expected = as.numeric(target$genept_raw_p),
    genept_observed = genept$raw_p,
    genept_absolute_difference = abs(genept$raw_p - as.numeric(target$genept_raw_p)),
    passed = abs(vanilla$raw_p - as.numeric(target$vanilla_raw_p)) <= 1e-12 &&
      abs(genept$raw_p - as.numeric(target$genept_raw_p)) <= 1e-12,
    warnings = unique(c(vanilla$warnings, genept$warnings))
  )
}

toy <- masking_toy_tests()
representative_targets <- lapply(unique(vapply(targets, function(x) x$comparison, character(1))), function(id) {
  Filter(function(x) identical(x$comparison, id), targets)[[1L]]
})
baseline_preflight <- lapply(representative_targets, baseline_check)
preflight_pass <- isTRUE(toy$vanilla_zero_vs_removal_pass) &&
  isTRUE(toy$genept_subtraction_vs_direct_pass) &&
  all(vapply(baseline_preflight, function(x) isTRUE(x$passed), logical(1))) &&
  !any(lengths(lapply(baseline_preflight, function(x) x$warnings)) > 0L)

if (preflight_only) {
  payload <- list(
    status = if (preflight_pass) "PASS" else "FAIL",
    scpa_version = as.character(utils::packageVersion("SCPA")),
    multicross_version = as.character(utils::packageVersion("multicross")),
    raw_p_source = "multicross::mcm result[[1]]",
    toy = toy,
    baseline_reproduction = baseline_preflight,
    same_cells_across_branches = isTRUE(manifest$same_cells_across_branches),
    same_pathways_and_gene_order = isTRUE(manifest$pathway_universe$same_across_comparisons) &&
      isTRUE(manifest$pathway_universe$gene_order_identical_across_comparisons),
    warnings = unname(unique(unlist(lapply(baseline_preflight, function(x) x$warnings))))
  )
  dir.create(dirname(output_json), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(payload, output_json, auto_unbox = TRUE, pretty = TRUE, null = "null")
  if (!preflight_pass) stop("Phase 5 preflight failed")
  cat("PHASE5_PREFLIGHT status=PASS baseline_comparisons=3\n")
  quit(save = "no", status = 0L)
}

if (!preflight_pass) stop("Phase 5 preflight failed before masking")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
run_started <- proc.time()[["elapsed"]]
target_total <- length(targets)
total_gene_branches <- sum(vapply(targets, function(target) {
  n <- as.integer(target$n_paired_genes)
  if (max_genes > 0L) n <- min(n, max_genes)
  n * 2L
}, integer(1)))
completed_branches <- 0L
completed_targets <- 0L
reused_targets <- 0L
runtime_warnings <- character()
baseline_all <- vector("list", target_total)

cat(sprintf(
  "[Phase5 core] Starting %d target pathway-comparisons, %d gene-mask branch evaluations\n",
  target_total, total_gene_branches
))
flush.console()

for (target_index in seq_along(targets)) {
  target <- targets[[target_index]]
  pathway <- pathway_lookup[[target$pathway]]
  indices <- as.integer(unlist(pathway$global_gene_indices)) + 1L
  genes <- as.character(unlist(pathway$paired_genes))
  if (max_genes > 0L) {
    indices <- utils::head(indices, max_genes)
    genes <- utils::head(genes, max_genes)
  }
  if (!identical(gene_names[indices], genes)) {
    stop("Gene-order assertion failed: ", target$comparison, ":", target$pathway)
  }
  checkpoint_path <- file.path(
    checkpoint_dir, target$comparison, paste0(target$pathway, "_gene_masking.csv")
  )
  if (file.exists(checkpoint_path)) {
    cached <- utils::read.csv(checkpoint_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(cached) != length(genes) || !identical(as.character(cached$gene), genes)) {
      stop("Incompatible Phase 5 checkpoint: ", checkpoint_path)
    }
    completed_targets <- completed_targets + 1L
    reused_targets <- reused_targets + 1L
    completed_branches <- completed_branches + length(genes) * 2L
    cat(sprintf(
      "[Target %02d/%02d] %s | %s | RESUME checkpoint reused\n",
      target_index, target_total, target$comparison, target$pathway
    ))
    flush.console()
    next
  }

  cat(sprintf(
    "[Target %02d/%02d] comparison=%s pathway=%s genes=%d START\n",
    target_index, target_total, target$comparison, target$pathway, length(genes)
  ))
  flush.console()
  xa <- groups[[target$group_a]][, indices, drop = FALSE]
  xb <- groups[[target$group_b]][, indices, drop = FALSE]
  ep <- embeddings[indices, , drop = FALSE]
  projected <- genept_non_l2_project_pair(xa, xb, ep)
  za <- projected$a
  zb <- projected$b
  baseline <- baseline_check(target)
  baseline_all[[target_index]] <- baseline
  if (!isTRUE(baseline$passed)) stop("Phase 4B baseline reproduction failed")
  runtime_warnings <- c(runtime_warnings, baseline$warnings)
  pathway_mean_embedding <- colMeans(ep)
  pathway_mean_norm <- sqrt(sum(pathway_mean_embedding * pathway_mean_embedding))
  target_rows <- vector("list", length(genes))

  for (gene_index in seq_along(genes)) {
    gene <- genes[[gene_index]]
    vanilla_mask <- vanilla_zero_mask_pair(xa, xb, gene_index)
    xa_masked <- vanilla_mask$a
    xb_masked <- vanilla_mask$b
    genept_mask <- genept_non_l2_subtraction_mask_pair(
      za, zb, xa, xb, ep, gene_index
    )
    za_masked <- genept_mask$a
    zb_masked <- genept_mask$b

    cat(sprintf(
      "[Target %02d/%02d gene %03d/%03d] %s | %s | branch=vanilla START\n",
      target_index, target_total, gene_index, length(genes), target$comparison, gene
    ))
    flush.console()
    vanilla <- run_checked(xa_masked, xb_masked, "vanilla gene mask")
    completed_branches <- completed_branches + 1L
    elapsed <- proc.time()[["elapsed"]] - run_started
    eta <- elapsed / completed_branches * (total_gene_branches - completed_branches)
    cat(sprintf("[Progress %04d/%04d] branch=vanilla DONE | elapsed %.1f min | ETA %.1f min\n",
                completed_branches, total_gene_branches, elapsed / 60, eta / 60))
    flush.console()

    cat(sprintf(
      "[Target %02d/%02d gene %03d/%03d] %s | %s | branch=genept_non_l2 START\n",
      target_index, target_total, gene_index, length(genes), target$comparison, gene
    ))
    flush.console()
    genept <- run_checked(za_masked, zb_masked, "GenePT gene mask")
    completed_branches <- completed_branches + 1L
    elapsed <- proc.time()[["elapsed"]] - run_started
    eta <- elapsed / completed_branches * (total_gene_branches - completed_branches)
    cat(sprintf("[Progress %04d/%04d] branch=genept_non_l2 DONE | elapsed %.1f min | ETA %.1f min\n",
                completed_branches, total_gene_branches, elapsed / 60, eta / 60))
    flush.console()
    runtime_warnings <- c(runtime_warnings, vanilla$warnings, genept$warnings)

    vanilla_adj_masked <- adjusted(vanilla$raw_p)
    genept_adj_masked <- adjusted(genept$raw_p)
    vanilla_full_sig <- as.numeric(target$vanilla_adjusted_p) < 0.05
    genept_full_sig <- as.numeric(target$genept_adjusted_p) < 0.05
    vanilla_masked_sig <- vanilla_adj_masked < 0.05
    genept_masked_sig <- genept_adj_masked < 0.05
    embedding_norm <- sqrt(sum(ep[gene_index, ] * ep[gene_index, ]))
    cosine <- if (embedding_norm > 0 && pathway_mean_norm > 0) {
      sum(ep[gene_index, ] * pathway_mean_embedding) / (embedding_norm * pathway_mean_norm)
    } else {
      0
    }
    target_rows[[gene_index]] <- data.frame(
      comparison = target$comparison,
      group_a = target$group_a,
      group_b = target$group_b,
      pathway = target$pathway,
      detection_state = target$detection_state,
      n_paired_genes = as.integer(target$n_paired_genes),
      gene = gene,
      gene_index = gene_index,
      vanilla_raw_p_full = as.numeric(target$vanilla_raw_p),
      vanilla_raw_p_masked = vanilla$raw_p,
      vanilla_score_full = score(target$vanilla_raw_p),
      vanilla_score_masked = score(vanilla$raw_p),
      vanilla_delta_score = score(target$vanilla_raw_p) - score(vanilla$raw_p),
      genept_raw_p_full = as.numeric(target$genept_raw_p),
      genept_raw_p_masked = genept$raw_p,
      genept_score_full = score(target$genept_raw_p),
      genept_score_masked = score(genept$raw_p),
      genept_delta_score = score(target$genept_raw_p) - score(genept$raw_p),
      vanilla_adjusted_p_full = as.numeric(target$vanilla_adjusted_p),
      vanilla_adjusted_p_masked = vanilla_adj_masked,
      genept_adjusted_p_full = as.numeric(target$genept_adjusted_p),
      genept_adjusted_p_masked = genept_adj_masked,
      vanilla_qval_full = as.numeric(target$vanilla_qval),
      vanilla_qval_masked = qval(vanilla_adj_masked),
      genept_qval_full = as.numeric(target$genept_qval),
      genept_qval_masked = qval(genept_adj_masked),
      vanilla_significant_full = vanilla_full_sig,
      vanilla_significant_masked = vanilla_masked_sig,
      genept_significant_full = genept_full_sig,
      genept_significant_masked = genept_masked_sig,
      vanilla_detection_flip = vanilla_full_sig != vanilla_masked_sig,
      genept_detection_flip = genept_full_sig != genept_masked_sig,
      mean_expression_group_A = mean(xa[, gene_index]),
      mean_expression_group_B = mean(xb[, gene_index]),
      expression_difference = mean(xb[, gene_index]) - mean(xa[, gene_index]),
      detection_fraction_group_A = mean(xa[, gene_index] > 0),
      detection_fraction_group_B = mean(xb[, gene_index] > 0),
      embedding_norm = embedding_norm,
      embedding_to_pathway_mean_cosine = cosine,
      stringsAsFactors = FALSE
    )
  }
  result <- do.call(rbind, target_rows)
  dir.create(dirname(checkpoint_path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(checkpoint_path, ".partial")
  utils::write.csv(result, temporary, row.names = FALSE, quote = TRUE)
  if (!file.rename(temporary, checkpoint_path)) stop("Could not commit checkpoint")
  completed_targets <- completed_targets + 1L
  cat(sprintf("[Target %02d/%02d] %s | %s | CHECKPOINT SAVED\n",
              target_index, target_total, target$comparison, target$pathway))
  flush.console()
}

payload <- list(
  status = "PASS",
  scpa_version = as.character(utils::packageVersion("SCPA")),
  multicross_version = as.character(utils::packageVersion("multicross")),
  raw_p_source = "multicross::mcm result[[1]]",
  raw_p_clip = as.numeric(manifest$phase5$raw_p_clip),
  eligible_pathway_count = as.integer(manifest$phase5$eligible_pathway_count),
  target_count = target_total,
  completed_target_count = completed_targets,
  reused_checkpoint_count = reused_targets,
  gene_mask_evaluation_count = total_gene_branches,
  failed_mcm_calls = 0L,
  warnings = unname(unique(runtime_warnings)),
  toy = toy,
  baseline_preflight = baseline_preflight,
  baseline_reproduction = Filter(Negate(is.null), baseline_all),
  checkpoint_resume = TRUE
)
dir.create(dirname(output_json), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(payload, output_json, auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(sprintf(
  "PHASE5_CORE status=PASS targets=%d gene_mask_branch_evaluations=%d reused=%d\n",
  completed_targets, total_gene_branches, reused_targets
))
