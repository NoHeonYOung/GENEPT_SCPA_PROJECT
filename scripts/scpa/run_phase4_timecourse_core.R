#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve script location")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
source(file.path(dirname(script_path), "scpa_core_adapter.R"))

argument <- function(name, required = TRUE, default = NULL) {
  index <- which(args == name)
  if (length(index) == 1L && index < length(args)) return(args[[index + 1L]])
  inline <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(inline) == 1L) return(sub(paste0("^", name, "="), "", inline[[1]]))
  if (required) stop("Missing argument: ", name)
  default
}

input_h5 <- normalizePath(argument("--input-h5"))
manifest_path <- normalizePath(argument("--manifest"))
output_csv <- argument("--output-csv")
output_json <- argument("--output-json")
max_pathways <- as.integer(argument("--max-pathways", FALSE, "0"))
max_comparisons <- as.integer(argument("--max-comparisons", FALSE, "0"))

required <- c("rhdf5", "jsonlite", "SCPA", "multicross")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) stop("Missing Phase 4B/C R packages: ", paste(missing, collapse = ", "))

manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
audit_pathways <- manifest$pathways
pathways <- audit_pathways
comparisons <- manifest$comparisons
if (max_pathways > 0L) pathways <- utils::head(pathways, max_pathways)
if (max_comparisons > 0L) comparisons <- utils::head(comparisons, max_comparisons)
if (length(pathways) == 0L || length(comparisons) == 0L) stop("Empty pathway/comparison universe")

gene_names <- as.character(rhdf5::h5read(input_h5, "gene_names"))
embeddings <- t(rhdf5::h5read(input_h5, "embeddings"))
group_names <- names(manifest$groups)
groups <- setNames(lapply(group_names, function(name) {
  t(rhdf5::h5read(input_h5, paste0("expression/", name)))
}), group_names)
if (nrow(embeddings) != length(gene_names) ||
    any(vapply(groups, ncol, integer(1)) != length(gene_names))) {
  stop("HDF5 time-course gene axes are not aligned")
}

rank_tie_aware <- function(qval) rank(-qval, ties.method = "average")
row_l2 <- function(x) {
  norms <- sqrt(rowSums(x * x))
  nonzero <- norms > 0
  x[nonzero, ] <- x[nonzero, , drop = FALSE] / norms[nonzero]
  x
}

official_raw_p <- function(xa, xb, genes, pathway_name, seed) {
  sample_a <- t(xa)
  sample_b <- t(xb)
  rownames(sample_a) <- genes
  rownames(sample_b) <- genes
  colnames(sample_a) <- paste0("A", seq_len(ncol(sample_a)))
  colnames(sample_b) <- paste0("B", seq_len(ncol(sample_b)))
  pathway <- data.frame(
    Pathway = rep(pathway_name, length(genes)),
    Genes = genes,
    stringsAsFactors = FALSE
  )
  set.seed(seed)
  runtime_warnings <- character()
  result <- withCallingHandlers(
    suppressMessages(SCPA::compare_pathways(
      samples = list(sample_a, sample_b),
      pathways = list(pathway),
      downsample = min(ncol(sample_a), ncol(sample_b)),
      min_genes = 1,
      max_genes = 500,
      parallel = FALSE
    )),
    warning = function(condition) {
      runtime_warnings <<- c(runtime_warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(raw_p = as.numeric(result$Pval[[1]]), warnings = unique(runtime_warnings))
}

# Cross-check three representative pathways on the CD4 0h-vs-24h positive control.
positive <- Filter(function(x) identical(x$id, "cd4_0h_vs_24h"), manifest$comparisons)
if (length(positive) != 1L) stop("Missing frozen CD4 0h-vs-24h positive control")
check_indices <- unique(as.integer(round(seq(1, length(audit_pathways), length.out = min(3L, length(audit_pathways))))))
crosscheck <- vector("list", length(check_indices))
for (j in seq_along(check_indices)) {
  pathway <- audit_pathways[[check_indices[[j]]]]
  indices <- as.integer(unlist(pathway$global_gene_indices)) + 1L
  genes <- as.character(unlist(pathway$paired_genes))
  xa <- groups[[positive[[1]]$group_a]][, indices, drop = FALSE]
  xb <- groups[[positive[[1]]$group_b]][, indices, drop = FALSE]
  adapter <- run_mcm_raw(xa, xb)
  official <- official_raw_p(xa, xb, genes, pathway$pathway, as.integer(manifest$seed))
  difference <- abs(adapter$raw_p - official$raw_p)
  crosscheck[[j]] <- list(
    pathway = pathway$pathway,
    adapter_raw_p = adapter$raw_p,
    official_raw_p = official$raw_p,
    absolute_difference = difference,
    tolerance = 1e-12,
    passed = difference <= 1e-12,
    warnings = unique(c(adapter$warnings, official$warnings))
  )
}
if (!all(vapply(crosscheck, function(x) isTRUE(x$passed), logical(1)))) {
  stop("Official SCPA raw-p cross-check failed")
}

all_results <- vector("list", length(comparisons))
runtime_warnings <- character()
run_started <- proc.time()[["elapsed"]]
total_units <- length(comparisons) * length(pathways)
completed <- 0L
total_branch_units <- total_units * 3L
branch_completed <- 0L

# The CSV is written after every completed comparison. On restart, validated
# completed comparisons are reused so a later failure does not restart the
# entire production run.
if (file.exists(output_csv)) {
  checkpoint <- utils::read.csv(output_csv, stringsAsFactors = FALSE, check.names = FALSE)
  allowed <- vapply(comparisons, function(x) x$id, character(1))
  if (!all(checkpoint$comparison %in% allowed)) {
    stop("Checkpoint contains comparisons outside the requested comparison set")
  }
  expected_pathways <- vapply(pathways, function(x) x$pathway, character(1))
  for (comparison_index in seq_along(comparisons)) {
    cached <- checkpoint[checkpoint$comparison == comparisons[[comparison_index]]$id, , drop = FALSE]
    if (nrow(cached) == 0L) next
    if (nrow(cached) != length(pathways) || !setequal(cached$pathway, expected_pathways)) {
      stop("Incomplete or incompatible checkpoint for ", comparisons[[comparison_index]]$id)
    }
    cached <- cached[match(expected_pathways, cached$pathway), , drop = FALSE]
    all_results[[comparison_index]] <- cached
    completed <- completed + length(pathways)
    branch_completed <- branch_completed + length(pathways) * 3L
  }
}
session_branch_completed <- 0L
session_branch_total <- total_branch_units - branch_completed

run_branch <- function(branch, expression, comparison_id, pathway_name) {
  cat(sprintf(
    "[Branch %04d/%04d] comparison=%s pathway=%s branch=%s START\n",
    branch_completed + 1L, total_branch_units, comparison_id, pathway_name, branch
  ))
  flush.console()
  value <- force(expression)
  branch_completed <<- branch_completed + 1L
  session_branch_completed <<- session_branch_completed + 1L
  elapsed <- proc.time()[["elapsed"]] - run_started
  eta <- if (session_branch_completed > 0L) {
    elapsed / session_branch_completed * (session_branch_total - session_branch_completed)
  } else {
    NA_real_
  }
  cat(sprintf(
    "[Branch %04d/%04d] comparison=%s pathway=%s branch=%s DONE | elapsed %.1f min | ETA %.1f min\n",
    branch_completed, total_branch_units, comparison_id, pathway_name, branch,
    elapsed / 60, eta / 60
  ))
  flush.console()
  value
}
cat(sprintf(
  "[Phase4B core] Starting %d comparisons x %d pathways x 3 branches\n",
  length(comparisons), length(pathways)
))
flush.console()

for (comparison_index in seq_along(comparisons)) {
  comparison <- comparisons[[comparison_index]]
  if (!is.null(all_results[[comparison_index]])) {
    cat(sprintf(
      "[Comparison %d/%d] %s RESUME: reusing completed checkpoint\n",
      comparison_index, length(comparisons), comparison$id
    ))
    flush.console()
    next
  }
  xa_all <- groups[[comparison$group_a]]
  xb_all <- groups[[comparison$group_b]]
  comparison_rows <- vector("list", length(pathways))
  cat(sprintf(
    "[Comparison %d/%d] %s (%s vs %s)\n",
    comparison_index, length(comparisons), comparison$id,
    comparison$group_a, comparison$group_b
  ))
  flush.console()
  for (pathway_index in seq_along(pathways)) {
    pathway <- pathways[[pathway_index]]
    indices <- as.integer(unlist(pathway$global_gene_indices)) + 1L
    expected <- as.character(unlist(pathway$paired_genes))
    if (!identical(gene_names[indices], expected)) {
      stop("Gene-order assertion failed for ", comparison$id, ":", pathway$pathway)
    }
    xa <- xa_all[, indices, drop = FALSE]
    xb <- xb_all[, indices, drop = FALSE]
    ep <- embeddings[indices, , drop = FALSE]
    za <- xa %*% ep
    zb <- xb %*% ep
    vanilla <- run_branch("vanilla", run_mcm_raw(xa, xb), comparison$id, pathway$pathway)
    genept <- run_branch("genept_non_l2", run_mcm_raw(za, zb), comparison$id, pathway$pathway)
    l2 <- run_branch("genept_row_l2", run_mcm_raw(row_l2(za), row_l2(zb)), comparison$id, pathway$pathway)
    if (length(vanilla$warnings)) runtime_warnings <- c(runtime_warnings, paste(comparison$id, pathway$pathway, "vanilla", vanilla$warnings, sep = ":"))
    if (length(genept$warnings)) runtime_warnings <- c(runtime_warnings, paste(comparison$id, pathway$pathway, "genept", genept$warnings, sep = ":"))
    if (length(l2$warnings)) runtime_warnings <- c(runtime_warnings, paste(comparison$id, pathway$pathway, "l2", l2$warnings, sep = ":"))
    comparison_rows[[pathway_index]] <- data.frame(
      comparison = comparison$id,
      comparison_class = comparison$comparison_class,
      group_a = comparison$group_a,
      group_b = comparison$group_b,
      pathway = pathway$pathway,
      source_database = pathway$source_database,
      n_primary_paired_genes = length(indices),
      vanilla_raw_p = vanilla$raw_p,
      genept_raw_p = genept$raw_p,
      l2_raw_p = l2$raw_p,
      stringsAsFactors = FALSE
    )
    completed <- completed + 1L
    elapsed <- proc.time()[["elapsed"]] - run_started
    eta <- elapsed / completed * (total_units - completed)
    cat(sprintf(
      "[Phase4B %04d/%04d] %s | %s | done | elapsed %.1f min | ETA %.1f min\n",
      completed, total_units, comparison$id, pathway$pathway, elapsed / 60, eta / 60
    ))
    flush.console()
  }
  result <- do.call(rbind, comparison_rows)
  for (method in c("vanilla", "genept", "l2")) {
    corrected <- scpa_pathway_qvalues(result[[paste0(method, "_raw_p")]])
    result[[paste0(method, "_adjusted_p")]] <- corrected$adjusted_p
    result[[paste0(method, "_qval")]] <- corrected$qval
    result[[paste0(method, "_rank")]] <- rank_tie_aware(corrected$qval)
  }
  result$rank_delta <- result$genept_rank - result$vanilla_rank
  all_results[[comparison_index]] <- result
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  completed_results <- Filter(Negate(is.null), all_results)
  utils::write.csv(do.call(rbind, completed_results), output_csv, row.names = FALSE, quote = TRUE)
}

results <- do.call(rbind, all_results)
payload <- list(
  status = "PASS",
  scpa_version = as.character(utils::packageVersion("SCPA")),
  multicross_version = as.character(utils::packageVersion("multicross")),
  raw_p_source = "multicross::mcm result[[1]]",
  multiple_testing = "Bonferroni within each comparison/method over frozen pathway universe",
  qval_formula = "sqrt(-log10(adjusted_p))",
  log_base = 10,
  tie_method = "average",
  comparison_count = length(comparisons),
  pathway_count = length(pathways),
  official_scpa_crosscheck = crosscheck,
  warnings = unname(unique(runtime_warnings))
)
dir.create(dirname(output_json), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(payload, output_json, auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(sprintf(
  "PHASE4_TIMECOURSE_CORE status=PASS comparisons=%d pathways=%d rows=%d\n",
  length(comparisons), length(pathways), nrow(results)
))
