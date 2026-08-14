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
has_flag <- function(name) any(args == name)

input_h5 <- normalizePath(argument("--input-h5"))
manifest_path <- normalizePath(argument("--manifest"))
checkpoint_dir <- argument("--checkpoint-dir")
output_json <- argument("--output-json")
mode <- argument("--mode", FALSE, "preflight")
cores <- max(1L, as.integer(argument("--cores", FALSE, "1")))
max_targets <- as.integer(argument("--max-targets", FALSE, "0"))
max_reps <- as.integer(argument("--max-reps", FALSE, "0"))
max_genes <- as.integer(argument("--max-genes", FALSE, "0"))

required <- c("rhdf5", "jsonlite", "SCPA", "multicross")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) stop("Missing Phase 6 R packages: ", paste(missing, collapse = ", "))

manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
gene_names <- as.character(rhdf5::h5read(input_h5, "gene_names"))
embeddings <- t(rhdf5::h5read(input_h5, "embeddings"))
group_names <- names(manifest$groups)
groups <- setNames(lapply(group_names, function(name) {
  t(rhdf5::h5read(input_h5, paste0("expression/", name)))
}), group_names)
pathways <- manifest$pathways
pathway_lookup <- setNames(pathways, vapply(pathways, function(x) x$pathway, character(1)))
targets <- if (mode == "gene") manifest$phase6$gene_targets else manifest$phase6$pathway_targets
if (max_targets > 0L) targets <- utils::head(targets, max_targets)
if (nrow(embeddings) != length(gene_names) ||
    any(vapply(groups, ncol, integer(1)) != length(gene_names))) {
  stop("Phase 6 HDF5 gene axes are not aligned")
}

raw_clip <- as.numeric(manifest$phase6$raw_p_clip)
eligible_count <- as.integer(manifest$phase6$eligible_pathway_count)
safe_p <- function(p) max(as.numeric(p), raw_clip)
score <- function(p) -log10(safe_p(p))
adjusted <- function(p) min(safe_p(p) * eligible_count, 1)
qval <- function(p_adj) if (p_adj == 0) Inf else sqrt(-log10(p_adj))

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

seed_for <- function(control, pathway_index, replicate_id, robustness = FALSE) {
  base <- if (control == "permuted") {
    as.integer(manifest$phase6$seeds$permuted)
  } else {
    as.integer(manifest$phase6$seeds$random)
  }
  if (robustness) base <- base + 200000L
  as.integer(base + pathway_index * 1000L + replicate_id)
}

make_control <- function(ep, control, seed) {
  set.seed(seed)
  p <- nrow(ep)
  true_norms <- sqrt(rowSums(ep * ep))
  if (control == "permuted") {
    permutation <- sample.int(p)
    attempts <- 1L
    while (mean(permutation != seq_len(p)) <= 0.9 && attempts < 100L) {
      permutation <- sample.int(p)
      attempts <- attempts + 1L
    }
    if (identical(permutation, seq_len(p)) || mean(permutation != seq_len(p)) <= 0.9) {
      stop("Could not generate a valid correspondence-destroying permutation")
    }
    controlled <- ep[permutation, , drop = FALSE]
    return(list(
      matrix = controlled, mapping_changed_fraction = mean(permutation != seq_len(p)),
      permutation = permutation, norm_max_difference = max(abs(sort(sqrt(rowSums(controlled * controlled))) - sort(true_norms))),
      vector_multiset_preserved = identical(sort(permutation), seq_len(p))
    ))
  }
  if (control == "random") {
    controlled <- matrix(stats::rnorm(p * ncol(ep)), nrow = p, ncol = ncol(ep))
    generated_norms <- sqrt(rowSums(controlled * controlled))
    if (any(generated_norms == 0)) stop("Zero-norm random direction")
    controlled <- controlled / generated_norms * true_norms
    return(list(
      matrix = controlled, mapping_changed_fraction = 1,
      permutation = integer(), norm_max_difference = max(abs(sqrt(rowSums(controlled * controlled)) - true_norms)),
      vector_multiset_preserved = FALSE
    ))
  }
  stop("Unknown control: ", control)
}

target_data <- function(target, source_groups = groups) {
  pathway <- pathway_lookup[[target$pathway]]
  indices <- as.integer(unlist(pathway$global_gene_indices)) + 1L
  genes <- as.character(unlist(pathway$paired_genes))
  if (!identical(gene_names[indices], genes)) stop("Gene order mismatch: ", target$pathway)
  list(
    indices = indices, genes = genes,
    xa = source_groups[[target$group_a]][, indices, drop = FALSE],
    xb = source_groups[[target$group_b]][, indices, drop = FALSE],
    ep = embeddings[indices, , drop = FALSE]
  )
}

atomic_csv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".partial.", Sys.getpid())
  utils::write.csv(data, temporary, row.names = FALSE, quote = TRUE)
  if (!file.rename(temporary, path)) stop("Could not atomically commit checkpoint: ", path)
}

atomic_progress <- function(payload, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".partial.", Sys.getpid())
  jsonlite::write_json(payload, temporary, auto_unbox = TRUE, pretty = FALSE, null = "null")
  if (!file.rename(temporary, path)) stop("Could not atomically commit progress marker: ", path)
}

preflight <- function() {
  target <- targets[[1L]]
  d <- target_data(target)
  true_result <- run_checked(d$xa %*% d$ep, d$xb %*% d$ep, "True reconstruction")
  true_difference <- abs(true_result$raw_p - as.numeric(target$genept_raw_p))
  perm_seed <- seed_for("permuted", 1L, 1L)
  random_seed <- seed_for("random", 1L, 1L)
  perm1 <- make_control(d$ep, "permuted", perm_seed)
  perm1_again <- make_control(d$ep, "permuted", perm_seed)
  perm2 <- make_control(d$ep, "permuted", perm_seed + 1L)
  random1 <- make_control(d$ep, "random", random_seed)
  random1_again <- make_control(d$ep, "random", random_seed)
  random2 <- make_control(d$ep, "random", random_seed + 1L)
  direct <- d$xa %*% perm1$matrix
  accumulated <- matrix(0, nrow(d$xa), ncol(perm1$matrix))
  for (g in seq_len(ncol(d$xa))) accumulated <- accumulated + tcrossprod(d$xa[, g], perm1$matrix[g, ])
  checks <- list(
    true_baseline_reproduced = true_difference <= 1e-12,
    permuted_vector_multiset_preserved = isTRUE(perm1$vector_multiset_preserved),
    permuted_mapping_changed_gt_0_9 = perm1$mapping_changed_fraction > 0.9,
    permuted_not_identity = !identical(perm1$permutation, seq_len(nrow(d$ep))),
    random_dimension_1536 = ncol(random1$matrix) == 1536L,
    row_norms_preserved = perm1$norm_max_difference <= 1e-12 && random1$norm_max_difference <= 1e-10,
    same_seed_deterministic = identical(perm1$matrix, perm1_again$matrix) && identical(random1$matrix, random1_again$matrix),
    different_seed_changes_control = !identical(perm1$matrix, perm2$matrix) && !identical(random1$matrix, random2$matrix),
    same_cells_and_gene_order = isTRUE(manifest$same_cells_across_branches) && identical(gene_names[d$indices], d$genes),
    no_l2_normalization = isFALSE(manifest$phase6$l2_normalization),
    direct_projection_sanity = max(abs(direct - accumulated)) <= 1e-10
  )
  list(
    status = if (all(unlist(checks))) "PASS" else "FAIL", checks = checks,
    true_raw_p_expected = as.numeric(target$genept_raw_p),
    true_raw_p_observed = true_result$raw_p,
    true_raw_p_absolute_difference = true_difference,
    permutation_changed_fraction = perm1$mapping_changed_fraction,
    permuted_norm_max_difference = perm1$norm_max_difference,
    random_norm_max_difference = random1$norm_max_difference,
    direct_projection_max_difference = max(abs(direct - accumulated)),
    warnings = character()
  )
}

pre <- preflight()
if (mode == "preflight") {
  payload <- c(pre, list(
    scpa_version = as.character(utils::packageVersion("SCPA")),
    multicross_version = as.character(utils::packageVersion("multicross")),
    raw_p_source = "multicross::mcm result[[1]]"
  ))
  dir.create(dirname(output_json), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(payload, output_json, auto_unbox = TRUE, pretty = TRUE, null = "null")
  if (payload$status != "PASS") stop("Phase 6 preflight failed")
  cat("PHASE6_PREFLIGHT status=PASS controls=permuted,random dimension=1536 no_l2=true\n")
  quit(save = "no", status = 0L)
}
if (pre$status != "PASS") stop("Phase 6 preflight failed before production")

control_types <- c("permuted", "random")
rep_count <- if (mode == "pathway") as.integer(manifest$phase6$replicates$pathway_per_control) else as.integer(manifest$phase6$replicates$gene_per_control)
if (max_reps > 0L) rep_count <- min(rep_count, max_reps)
tasks <- expand.grid(
  target_index = seq_along(targets), control = control_types,
  replicate_id = seq_len(rep_count), stringsAsFactors = FALSE
)

checkpoint_path <- function(task) file.path(
  checkpoint_dir, sprintf("target_%02d", task$target_index),
  sprintf("%s_rep_%03d.csv", task$control, task$replicate_id)
)

pathway_worker <- function(task_index) {
  task <- tasks[task_index, ]
  destination <- checkpoint_path(task)
  if (file.exists(destination)) return(list(reused = TRUE, path = destination))
  target <- targets[[task$target_index]]
  d <- target_data(target)
  seed <- seed_for(task$control, task$target_index, task$replicate_id)
  control <- make_control(d$ep, task$control, seed)
  result <- run_checked(d$xa %*% control$matrix, d$xb %*% control$matrix, "Phase 6A control")
  adj <- adjusted(result$raw_p)
  row <- data.frame(
    target_index = task$target_index, comparison = target$comparison,
    group_a = target$group_a, group_b = target$group_b, pathway = target$pathway,
    detection_state = target$detection_state, n_paired_genes = length(d$genes),
    control = task$control, replicate_id = task$replicate_id, seed = seed,
    mapping_changed_fraction = control$mapping_changed_fraction,
    vector_multiset_preserved = control$vector_multiset_preserved,
    norm_max_difference = control$norm_max_difference,
    raw_p = result$raw_p, adjusted_p = adj, qval = qval(adj), score = score(result$raw_p),
    significant = adj < 0.05, stringsAsFactors = FALSE
  )
  atomic_csv(row, destination)
  if (task_index == 1L || task_index %% 50L == 0L || task_index == nrow(tasks)) {
    cat(sprintf("[Phase6A %04d/%04d] target=%02d control=%s rep=%03d DONE pid=%d\n",
                task_index, nrow(tasks), task$target_index, task$control, task$replicate_id, Sys.getpid()))
    flush.console()
  }
  list(reused = FALSE, path = destination)
}

gene_worker <- function(task_index) {
  task <- tasks[task_index, ]
  destination <- checkpoint_path(task)
  if (file.exists(destination)) return(list(reused = TRUE, path = destination))
  progress_path <- paste0(destination, ".progress.json")
  target <- targets[[task$target_index]]
  d <- target_data(target)
  if (max_genes > 0L) {
    keep <- seq_len(min(max_genes, length(d$genes)))
    d$indices <- d$indices[keep]; d$genes <- d$genes[keep]
    d$xa <- d$xa[, keep, drop = FALSE]; d$xb <- d$xb[, keep, drop = FALSE]
    d$ep <- d$ep[keep, , drop = FALSE]
  }
  seed <- seed_for(task$control, task$target_index, task$replicate_id)
  control <- make_control(d$ep, task$control, seed)
  za <- d$xa %*% control$matrix
  zb <- d$xb %*% control$matrix
  full <- run_checked(za, zb, "Phase 6B control full")
  full_adj <- adjusted(full$raw_p)
  atomic_progress(list(
    mode = "gene", target_index = task$target_index, control = task$control,
    replicate_id = task$replicate_id, current_gene = "FULL_BASELINE_DONE",
    completed_mcm = 1L, total_mcm = length(d$genes) + 1L, pid = Sys.getpid()
  ), progress_path)
  rows <- vector("list", length(d$genes))
  started <- proc.time()[["elapsed"]]
  for (gene_index in seq_along(d$genes)) {
    masked <- run_checked(
      za - tcrossprod(d$xa[, gene_index], control$matrix[gene_index, ]),
      zb - tcrossprod(d$xb[, gene_index], control$matrix[gene_index, ]),
      "Phase 6B control gene mask"
    )
    masked_adj <- adjusted(masked$raw_p)
    rows[[gene_index]] <- data.frame(
      target_index = task$target_index, comparison = target$comparison,
      pathway = target$pathway, detection_state = target$detection_state,
      control = task$control, replicate_id = task$replicate_id, seed = seed,
      n_paired_genes = length(d$genes), gene = d$genes[[gene_index]], gene_index = gene_index,
      control_raw_p_full = full$raw_p, control_adjusted_p_full = full_adj,
      control_qval_full = qval(full_adj), control_score_full = score(full$raw_p),
      control_raw_p_masked = masked$raw_p, control_adjusted_p_masked = masked_adj,
      control_qval_masked = qval(masked_adj), control_score_masked = score(masked$raw_p),
      control_delta_score = score(full$raw_p) - score(masked$raw_p),
      control_significant_full = full_adj < 0.05,
      control_significant_masked = masked_adj < 0.05,
      control_detection_flip = (full_adj < 0.05) != (masked_adj < 0.05),
      mapping_changed_fraction = control$mapping_changed_fraction,
      norm_max_difference = control$norm_max_difference, stringsAsFactors = FALSE
    )
    if (gene_index %% 5L == 0L || gene_index == length(d$genes)) {
      atomic_progress(list(
        mode = "gene", target_index = task$target_index, control = task$control,
        replicate_id = task$replicate_id, current_gene = d$genes[[gene_index]],
        completed_mcm = gene_index + 1L, total_mcm = length(d$genes) + 1L,
        pid = Sys.getpid()
      ), progress_path)
    }
    if (gene_index == length(d$genes)) {
      elapsed <- proc.time()[["elapsed"]] - started
      cat(sprintf("[Phase6B task %03d/%03d DONE] target=%02d control=%s rep=%02d genes=%d elapsed=%.1fmin pid=%d\n",
                  task_index, nrow(tasks), task$target_index, task$control,
                  task$replicate_id, length(d$genes), elapsed / 60, Sys.getpid()))
      flush.console()
    }
  }
  atomic_csv(do.call(rbind, rows), destination)
  if (file.exists(progress_path)) unlink(progress_path)
  list(reused = FALSE, path = destination)
}

robustness_worker <- function(task_index) {
  # Robustness tasks are explicitly rebuilt as replicate x target x representation.
  robust_reps <- if (max_reps > 0L) min(as.integer(manifest$phase6$replicates$resampling), max_reps) else as.integer(manifest$phase6$replicates$resampling)
  representations <- c("true", "permuted", "random")
  robust_tasks <- expand.grid(target_index = seq_along(targets), representation = representations,
                              resample_id = seq_len(robust_reps), stringsAsFactors = FALSE)
  task <- robust_tasks[task_index, ]
  destination <- file.path(checkpoint_dir, sprintf("resample_%02d", task$resample_id),
                           sprintf("target_%02d_%s.csv", task$target_index, task$representation))
  if (file.exists(destination)) return(list(reused = TRUE, path = destination))
  target <- targets[[task$target_index]]
  source_groups <- setNames(lapply(group_names, function(name) {
    t(rhdf5::h5read(input_h5, sprintf("resampling/rep_%02d/expression/%s", task$resample_id, name)))
  }), group_names)
  d <- target_data(target, source_groups)
  if (task$representation == "true") {
    projection <- d$ep; seed <- NA_integer_; changed <- 0; norm_diff <- 0
  } else {
    seed <- seed_for(task$representation, task$target_index, task$resample_id, robustness = TRUE)
    control <- make_control(d$ep, task$representation, seed)
    projection <- control$matrix; changed <- control$mapping_changed_fraction; norm_diff <- control$norm_max_difference
  }
  result <- run_checked(d$xa %*% projection, d$xb %*% projection, "Phase 6 robustness")
  adj <- adjusted(result$raw_p)
  row <- data.frame(
    resample_id = task$resample_id, target_index = task$target_index,
    comparison = target$comparison, pathway = target$pathway,
    representation = task$representation, seed = seed,
    mapping_changed_fraction = changed, norm_max_difference = norm_diff,
    raw_p = result$raw_p, adjusted_p = adj, qval = qval(adj), score = score(result$raw_p),
    significant = adj < 0.05, stringsAsFactors = FALSE
  )
  atomic_csv(row, destination)
  if (task_index == 1L || task_index %% 10L == 0L || task_index == nrow(robust_tasks)) {
    cat(sprintf("[Phase6R %03d/%03d] resample=%02d target=%02d representation=%s DONE pid=%d\n",
                task_index, nrow(robust_tasks), task$resample_id, task$target_index, task$representation, Sys.getpid()))
    flush.console()
  }
  list(reused = FALSE, path = destination)
}

if (mode == "robustness") {
  robust_reps <- if (max_reps > 0L) min(as.integer(manifest$phase6$replicates$resampling), max_reps) else as.integer(manifest$phase6$replicates$resampling)
  task_count <- length(targets) * 3L * robust_reps
  worker <- robustness_worker
} else if (mode == "pathway") {
  task_count <- nrow(tasks); worker <- pathway_worker
} else if (mode == "gene") {
  task_count <- nrow(tasks); worker <- gene_worker
} else stop("Unknown Phase 6 mode: ", mode)

dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
# Progress markers are ephemeral and may describe interrupted workers. Completed
# CSV checkpoints are the sole resume authority, so stale markers are removed.
stale_progress <- list.files(checkpoint_dir, pattern = "\\.progress\\.json$", recursive = TRUE, full.names = TRUE)
if (length(stale_progress) > 0L) unlink(stale_progress)
cat(sprintf("[Phase6 core] mode=%s tasks=%d cores=%d START\n", mode, task_count, cores)); flush.console()
started <- proc.time()[["elapsed"]]
indices <- seq_len(task_count)
if (cores > 1L && .Platform$OS.type != "windows") {
  # Schedule one task at a time so that workers which finish short pathways can
  # immediately take another task. This prevents the 170-gene target from being
  # stranded on only two workers near the end of Phase 6B.
  outcomes <- parallel::mclapply(indices, worker, mc.cores = cores, mc.preschedule = FALSE)
} else {
  outcomes <- lapply(indices, worker)
}
bad <- vapply(outcomes, inherits, logical(1), "try-error")
if (any(bad)) stop("Phase 6 worker failure: ", paste(as.character(outcomes[bad]), collapse = " | "))
payload <- list(
  status = "PASS", mode = mode, task_count = task_count,
  completed_task_count = length(outcomes),
  reused_checkpoint_count = sum(vapply(outcomes, function(x) isTRUE(x$reused), logical(1))),
  failed_mcm_calls = 0L, warnings = character(), cores = cores,
  elapsed_seconds = proc.time()[["elapsed"]] - started,
  checkpoint_resume = TRUE,
  scpa_version = as.character(utils::packageVersion("SCPA")),
  multicross_version = as.character(utils::packageVersion("multicross"))
)
dir.create(dirname(output_json), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(payload, output_json, auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(sprintf("PHASE6_CORE status=PASS mode=%s tasks=%d reused=%d elapsed=%.1fmin\n",
            mode, task_count, payload$reused_checkpoint_count, payload$elapsed_seconds / 60))
