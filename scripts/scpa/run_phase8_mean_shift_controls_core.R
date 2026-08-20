#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve script location")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]))
source(file.path(dirname(script_path), "scpa_core_adapter.R"))
source(file.path(dirname(script_path), "gene_masking_lib.R"))

argument <- function(name, required = TRUE, default = NULL) {
  index <- which(args == name)
  if (length(index) == 1L && index < length(args)) return(args[[index + 1L]])
  if (required) stop("Missing argument: ", name)
  default
}

expression_h5 <- normalizePath(argument("--expression-h5"))
controls_h5 <- normalizePath(argument("--controls-h5"))
manifest_path <- normalizePath(argument("--manifest"))
checkpoint_dir <- argument("--checkpoint-dir")
output_json <- argument("--output-json")
cores <- max(1L, as.integer(argument("--cores", FALSE, "1")))
max_experiments <- max(0L, as.integer(argument("--max-experiments", FALSE, "0")))
progress_every_genes <- max(1L, as.integer(argument("--progress-every-genes", FALSE, "10")))

required <- c("rhdf5", "jsonlite", "SCPA", "multicross")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) stop("Missing Phase 8 R packages: ", paste(missing, collapse = ", "))

manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
if (!identical(manifest$scope, "mean_shift_mechanism_decomposition") ||
    !identical(manifest$control_protocol, "phase8_deranged_and_norm_matched_random_v1")) {
  stop("Manifest is not the frozen Phase 8 mechanism decomposition")
}
if (!isFALSE(manifest$anti_leakage$ground_truth_file_read) ||
    !isFALSE(manifest$anti_leakage$ground_truth_labels_parsed) ||
    !isFALSE(manifest$anti_leakage$truth_fields_in_execution_manifest)) {
  stop("Phase 8 anti-leakage manifest gate failed")
}
all_experiments <- manifest$experiments
manifest_experiment_count <- length(all_experiments)
experiments <- all_experiments
if (max_experiments > 0L) experiments <- utils::head(experiments, max_experiments)
if (length(experiments) == 0L) stop("No Phase 8 experiments")
raw_clip <- as.numeric(manifest$masking$raw_p_clip)

run_checked <- function(xa, xb, label) {
  value <- tryCatch(run_mcm_raw(xa, xb), error = function(condition) {
    stop(label, ": ", conditionMessage(condition))
  })
  if (!is.finite(value$raw_p)) stop(label, ": non-finite raw p")
  if (length(value$warnings) > 0L) {
    stop(label, ": runtime warning: ", paste(value$warnings, collapse = " | "))
  }
  as.numeric(value$raw_p)
}

atomic_csv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".partial.", Sys.getpid())
  utils::write.csv(data, temporary, row.names = FALSE, quote = TRUE)
  if (!file.rename(temporary, path)) stop("Could not commit checkpoint: ", path)
}

checkpoint_path <- function(experiment) {
  file.path(checkpoint_dir, paste0(experiment$experiment_id, "_controls.csv"))
}

checkpoint_valid <- function(experiment) {
  path <- checkpoint_path(experiment)
  if (!file.exists(path)) return(FALSE)
  cached <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(condition) NULL
  )
  if (is.null(cached)) return(FALSE)
  nrow(cached) == as.integer(experiment$analysis_gene_count) &&
    length(unique(cached$experiment_id)) == 1L &&
    identical(as.character(unique(cached$experiment_id)), experiment$experiment_id) &&
    identical(as.character(unique(cached$control_protocol)), manifest$control_protocol) &&
    identical(as.character(unique(cached$phase8_experiment_input_sha256)),
              experiment$phase8_experiment_input_sha256) &&
    identical(as.character(unique(cached$permuted_sha256)), experiment$permuted_sha256) &&
    identical(as.character(unique(cached$random_sha256)), experiment$random_sha256) &&
    all(c("permuted_delta_score", "random_delta_score") %in% names(cached)) &&
    all(is.finite(cached$permuted_delta_score)) && all(is.finite(cached$random_delta_score))
}

read_control <- function(path) {
  value <- t(rhdf5::h5read(controls_h5, sub("^/", "", path)))
  storage.mode(value) <- "double"
  value
}

run_experiment <- function(experiment) {
  id <- experiment$experiment_id
  destination <- checkpoint_path(experiment)
  mcm_per_experiment <- 2L * (as.integer(experiment$analysis_gene_count) + 1L)
  if (checkpoint_valid(experiment)) {
    return(list(id = id, reused = TRUE, mcm = mcm_per_experiment))
  }
  genes <- as.character(rhdf5::h5read(expression_h5, sub("^/", "", experiment$gene_names_h5)))
  xa <- t(rhdf5::h5read(expression_h5, sub("^/", "", experiment$condition_a_h5)))
  xb <- t(rhdf5::h5read(expression_h5, sub("^/", "", experiment$condition_b_h5)))
  controls <- list(
    permuted = read_control(experiment$permuted_h5),
    random = read_control(experiment$random_h5)
  )
  if (length(genes) != as.integer(experiment$analysis_gene_count) ||
      ncol(xa) != length(genes) || ncol(xb) != length(genes) ||
      any(vapply(controls, nrow, integer(1)) != length(genes)) ||
      any(vapply(controls, ncol, integer(1)) != 1536L)) {
    stop("HDF5 axes are not aligned for ", id)
  }
  result <- data.frame(
    experiment_id = rep(id, length(genes)),
    draw_id = rep(as.integer(experiment$draw_id), length(genes)),
    pathway = rep(experiment$pathway, length(genes)),
    perturbation_type = rep(experiment$perturbation_type, length(genes)),
    perturbation_strength = rep(as.numeric(experiment$perturbation_strength), length(genes)),
    perturbation_seed = rep(as.integer(experiment$perturbation_seed), length(genes)),
    gene = genes,
    gene_index = seq_along(genes) - 1L,
    control_protocol = rep(manifest$control_protocol, length(genes)),
    phase8_experiment_input_sha256 = rep(experiment$phase8_experiment_input_sha256, length(genes)),
    permuted_seed = rep(as.integer(experiment$permuted_seed), length(genes)),
    random_seed = rep(as.integer(experiment$random_seed), length(genes)),
    permuted_sha256 = rep(experiment$permuted_sha256, length(genes)),
    random_sha256 = rep(experiment$random_sha256, length(genes)),
    stringsAsFactors = FALSE
  )
  cat(sprintf("[Phase8 worker pid=%d] %s START genes=%d controls=2\n",
              Sys.getpid(), id, length(genes)))
  flush.console()
  for (control_name in names(controls)) {
    embeddings <- controls[[control_name]]
    projected <- genept_non_l2_project_pair(xa, xb, embeddings)
    full <- run_checked(projected$a, projected$b, paste(id, control_name, "full"))
    masked_p <- numeric(length(genes))
    for (gene_index in seq_along(genes)) {
      masked <- genept_non_l2_subtraction_mask_pair(
        projected$a, projected$b, xa, xb, embeddings, gene_index
      )
      masked_p[[gene_index]] <- run_checked(
        masked$a, masked$b, paste(id, control_name, "mask", genes[[gene_index]])
      )
      if (gene_index %% progress_every_genes == 0L || gene_index == length(genes)) {
        cat(sprintf("[Phase8 worker pid=%d] %s control=%s genes=%d/%d\n",
                    Sys.getpid(), id, control_name, gene_index, length(genes)))
        flush.console()
      }
    }
    result[[paste0(control_name, "_raw_p_full")]] <- full
    result[[paste0(control_name, "_raw_p_masked")]] <- masked_p
    delta <- masking_score(full, raw_clip) - vapply(masked_p, masking_score, numeric(1), raw_p_clip = raw_clip)
    result[[paste0(control_name, "_delta_score")]] <- delta
    result[[paste0(control_name, "_signed_rank_average")]] <- rank(-delta, ties.method = "average")
    order_index <- order(-delta, genes)
    strict <- integer(length(order_index)); strict[order_index] <- seq_along(order_index)
    result[[paste0(control_name, "_signed_rank")]] <- strict
  }
  atomic_csv(result, destination)
  cat(sprintf("[Phase8 worker pid=%d] %s DONE checkpoint=%s\n", Sys.getpid(), id, destination))
  flush.console()
  list(id = id, reused = FALSE, mcm = mcm_per_experiment)
}

dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
total <- length(experiments)
total_mcm <- sum(vapply(experiments, function(x) {
  2L * (as.integer(x$analysis_gene_count) + 1L)
}, integer(1)))
valid_checkpoints <- vapply(experiments, checkpoint_valid, logical(1))
pending <- experiments[!valid_checkpoints]
reused <- sum(valid_checkpoints)
completed <- reused
completed_mcm <- sum(vapply(experiments[valid_checkpoints], function(x) {
  2L * (as.integer(x$analysis_gene_count) + 1L)
}, integer(1)))
new_mcm <- 0L
started <- proc.time()[["elapsed"]]
cat(sprintf(paste0(
  "[Phase8 core] experiments=%d pending=%d reused=%d mcm=%d ",
  "permuted=%d random=%d cores=%d checkpoint_resume=TRUE\n"
), total, length(pending), reused, total_mcm, total_mcm %/% 2L, total_mcm %/% 2L, cores))
flush.console()

if (length(pending) > 0L) {
  indices <- seq_along(pending)
  for (batch_start in seq(1L, length(pending), by = cores)) {
    batch_indices <- indices[batch_start:min(length(pending), batch_start + cores - 1L)]
    batch <- pending[batch_indices]
    results <- if (cores > 1L && .Platform$OS.type != "windows") {
      parallel::mclapply(batch, run_experiment, mc.cores = min(cores, length(batch)),
                         mc.preschedule = FALSE)
    } else {
      lapply(batch, run_experiment)
    }
    failed <- vapply(results, inherits, logical(1), what = "try-error")
    if (any(failed)) stop("Phase 8 worker failed: ", paste(results[failed], collapse = " | "))
    completed <- completed + length(results)
    batch_mcm <- sum(vapply(results, function(x) as.integer(x$mcm), integer(1)))
    completed_mcm <- completed_mcm + batch_mcm
    new_mcm <- new_mcm + sum(vapply(results, function(x) {
      if (isTRUE(x$reused)) 0L else as.integer(x$mcm)
    }, integer(1)))
    reused <- reused + sum(vapply(results, function(x) isTRUE(x$reused), logical(1)))
    elapsed <- proc.time()[["elapsed"]] - started
    eta <- if (new_mcm > 0L) elapsed / new_mcm * (total_mcm - completed_mcm) else 0
    cat(sprintf(paste0(
      "[Phase8 progress] experiments=%d/%d mcm=%d/%d remaining=%d ",
      "elapsed=%.1fmin ETA=%.1fmin reused=%d\n"
    ), completed, total, completed_mcm, total_mcm, total_mcm - completed_mcm,
    elapsed / 60, eta / 60, reused))
    flush.console()
  }
}

control_units <- manifest$control_units
fixed_points <- vapply(control_units, function(x) as.integer(x$permutation_fixed_point_count), integer(1))
perm_norm_diff <- vapply(control_units, function(x) as.numeric(x$permuted_sorted_row_norm_max_abs_difference), numeric(1))
random_norm_diff <- vapply(control_units, function(x) as.numeric(x$random_corresponding_row_norm_max_abs_difference), numeric(1))
nonfinite <- vapply(control_units, function(x) as.integer(x$nonfinite_count), integer(1))
payload <- list(
  status = if (total < manifest_experiment_count) "SMOKE_PASS" else "PASS",
  scope = manifest$scope,
  backend = "multicross::mcm",
  scpa_version = as.character(utils::packageVersion("SCPA")),
  multicross_version = as.character(utils::packageVersion("multicross")),
  manifest_experiment_count = manifest_experiment_count,
  selected_experiment_count = total,
  partial_run = total < manifest_experiment_count,
  completed_experiments = completed,
  reused_checkpoints = reused,
  mcm_count = completed_mcm,
  mcm_count_permuted = completed_mcm %/% 2L,
  mcm_count_random = completed_mcm %/% 2L,
  new_mcm_this_invocation = new_mcm,
  failed_mcm_calls = 0L,
  new_vanilla_mcm = 0L,
  new_true_genept_mcm = 0L,
  masking_library = "scripts/scpa/gene_masking_lib.R",
  projection_l2 = FALSE,
  checkpoint_resume = TRUE,
  cores = cores,
  control_qc = list(
    control_unit_count = length(control_units),
    permuted_fixed_point_max = max(fixed_points),
    permuted_row_multiset_preserved_all = all(vapply(
      control_units, function(x) isTRUE(x$permuted_vector_multiset_preserved), logical(1)
    )),
    permuted_sorted_row_norm_max_abs_difference = max(perm_norm_diff),
    random_dimension = 1536L,
    random_corresponding_row_norm_max_abs_difference = max(random_norm_diff),
    nonfinite_count = sum(nonfinite),
    same_control_paths_for_three_states = TRUE
  ),
  ground_truth_read_during_control_generation = FALSE,
  ground_truth_read_during_masking = FALSE
)
dir.create(dirname(output_json), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(payload, output_json, auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(sprintf("PHASE8_CORE status=%s experiments=%d/%d mcm=%d reused=%d new_this_run=%d\n",
            payload$status, completed, manifest_experiment_count, completed_mcm, reused, new_mcm))
