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

input_h5 <- normalizePath(argument("--input-h5"))
manifest_path <- normalizePath(argument("--manifest"))
checkpoint_dir <- argument("--checkpoint-dir")
output_json <- argument("--output-json")
cores <- max(1L, as.integer(argument("--cores", FALSE, "1")))
max_experiments <- max(0L, as.integer(argument("--max-experiments", FALSE, "0")))
progress_every_genes <- max(1L, as.integer(argument("--progress-every-genes", FALSE, "10")))

required <- c("rhdf5", "jsonlite", "SCPA", "multicross")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) stop("Missing Phase 7 R packages: ", paste(missing, collapse = ", "))

manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
if (!identical(manifest$scope, "llm_free_vanilla_vs_genept_ground_truth_recovery") ||
    isTRUE(manifest$llm_backend_present)) {
  stop("Manifest is not the frozen LLM-free Phase 7 benchmark")
}
experiments <- manifest$experiments
manifest_experiment_count <- length(experiments)
if (max_experiments > 0L) experiments <- utils::head(experiments, max_experiments)
if (length(experiments) == 0L) stop("No Phase 7 experiments")
raw_clip <- as.numeric(manifest$masking_implementation$raw_p_clip)

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
  file.path(checkpoint_dir, paste0(experiment$experiment_id, "_masking.csv"))
}

checkpoint_valid <- function(experiment) {
  path <- checkpoint_path(experiment)
  if (!file.exists(path)) return(FALSE)
  cached <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  nrow(cached) == as.integer(experiment$analysis_gene_count) &&
    length(unique(cached$experiment_id)) == 1L &&
    identical(as.character(unique(cached$experiment_id)), experiment$experiment_id) &&
    "experiment_input_sha256" %in% names(cached) &&
    identical(as.character(unique(cached$experiment_input_sha256)),
              experiment$experiment_input_sha256) &&
    "masking_protocol_version" %in% names(cached) &&
    identical(as.character(unique(cached$masking_protocol_version)),
              manifest$masking_implementation$protocol_version)
}

run_experiment <- function(experiment) {
  id <- experiment$experiment_id
  destination <- checkpoint_path(experiment)
  if (checkpoint_valid(experiment)) {
    return(list(id = id, reused = TRUE,
                mcm = 2L * (as.integer(experiment$analysis_gene_count) + 1L)))
  }
  genes <- as.character(rhdf5::h5read(input_h5, sub("^/", "", experiment$gene_names_h5)))
  xa <- t(rhdf5::h5read(input_h5, sub("^/", "", experiment$condition_a_h5)))
  xb <- t(rhdf5::h5read(input_h5, sub("^/", "", experiment$condition_b_h5)))
  ep <- t(rhdf5::h5read(input_h5, sub("^/", "", experiment$embeddings_h5)))
  if (length(genes) != as.integer(experiment$analysis_gene_count) ||
      ncol(xa) != length(genes) || ncol(xb) != length(genes) || nrow(ep) != length(genes)) {
    stop("HDF5 axes are not aligned for ", id)
  }
  projected <- genept_non_l2_project_pair(xa, xb, ep)
  vanilla_full <- run_checked(xa, xb, paste(id, "Vanilla full"))
  genept_full <- run_checked(projected$a, projected$b, paste(id, "GenePT full"))
  rows <- vector("list", length(genes))
  cat(sprintf("[Phase7 worker pid=%d] %s START genes=%d\n", Sys.getpid(), id, length(genes)))
  flush.console()
  for (gene_index in seq_along(genes)) {
    vanilla_mask <- vanilla_zero_mask_pair(xa, xb, gene_index)
    genept_mask <- genept_non_l2_subtraction_mask_pair(
      projected$a, projected$b, xa, xb, ep, gene_index
    )
    vanilla_masked <- run_checked(vanilla_mask$a, vanilla_mask$b, paste(id, "Vanilla mask"))
    genept_masked <- run_checked(genept_mask$a, genept_mask$b, paste(id, "GenePT mask"))
    rows[[gene_index]] <- data.frame(
      experiment_id = id, draw_id = as.integer(experiment$draw_id),
      experiment_input_sha256 = experiment$experiment_input_sha256,
      masking_protocol_version = manifest$masking_implementation$protocol_version,
      pathway = experiment$pathway, perturbation_type = experiment$perturbation_type,
      perturbation_strength = as.numeric(experiment$perturbation_strength),
      perturbation_seed = as.integer(experiment$perturbation_seed),
      gene = genes[[gene_index]], gene_index = gene_index - 1L,
      vanilla_raw_p_full = vanilla_full, vanilla_raw_p_masked = vanilla_masked,
      vanilla_delta_score = masking_score(vanilla_full, raw_clip) -
        masking_score(vanilla_masked, raw_clip),
      genept_raw_p_full = genept_full, genept_raw_p_masked = genept_masked,
      genept_delta_score = masking_score(genept_full, raw_clip) -
        masking_score(genept_masked, raw_clip),
      stringsAsFactors = FALSE
    )
    if (gene_index %% progress_every_genes == 0L || gene_index == length(genes)) {
      cat(sprintf("[Phase7 worker pid=%d] %s genes=%d/%d\n",
                  Sys.getpid(), id, gene_index, length(genes)))
      flush.console()
    }
  }
  result <- do.call(rbind, rows)
  for (method in c("vanilla", "genept")) {
    delta <- result[[paste0(method, "_delta_score")]]
    result[[paste0(method, "_signed_rank_average")]] <- rank(-delta, ties.method = "average")
    order_index <- order(-delta, result$gene)
    strict <- integer(length(order_index)); strict[order_index] <- seq_along(order_index)
    result[[paste0(method, "_signed_rank")]] <- strict
  }
  atomic_csv(result, destination)
  cat(sprintf("[Phase7 worker pid=%d] %s DONE checkpoint=%s\n", Sys.getpid(), id, destination))
  flush.console()
  list(id = id, reused = FALSE, mcm = 2L * (length(genes) + 1L))
}

dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
total <- length(experiments)
total_mcm <- sum(vapply(experiments, function(x) {
  2L * (as.integer(x$analysis_gene_count) + 1L)
}, integer(1)))
valid_checkpoints <- vapply(experiments, checkpoint_valid, logical(1))
pending <- experiments[!valid_checkpoints]
reused <- sum(valid_checkpoints)
reused_mcm <- sum(vapply(experiments[valid_checkpoints], function(x) {
  2L * (as.integer(x$analysis_gene_count) + 1L)
}, integer(1)))
completed <- reused
completed_mcm <- reused_mcm
new_mcm <- 0L
started <- proc.time()[["elapsed"]]
cat(sprintf(paste0(
  "[Phase7 LLM-free core] experiments=%d pending=%d reused=%d ",
  "mcm=%d cores=%d checkpoint_resume=TRUE\n"
), total, length(pending), reused, total_mcm, cores))
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
    if (any(failed)) stop("Phase 7 worker failed: ", paste(results[failed], collapse = " | "))
    completed <- completed + length(results)
    batch_mcm <- sum(vapply(results, function(x) as.integer(x$mcm), integer(1)))
    completed_mcm <- completed_mcm + batch_mcm
    new_mcm <- new_mcm + sum(vapply(results, function(x) {
      if (isTRUE(x$reused)) 0L else as.integer(x$mcm)
    }, integer(1)))
    reused <- reused + sum(vapply(results, function(x) isTRUE(x$reused), logical(1)))
    elapsed <- proc.time()[["elapsed"]] - started
    eta <- if (new_mcm > 0L) elapsed / new_mcm * (total_mcm - completed_mcm) else NA_real_
    cat(sprintf("[Phase7 progress] experiments=%d/%d mcm=%d/%d elapsed=%.1fmin ETA=%.1fmin reused=%d\n",
                completed, total, completed_mcm, total_mcm, elapsed / 60, eta / 60, reused))
    flush.console()
  }
}

payload <- list(
  status = if (total < manifest_experiment_count) "SMOKE_PASS" else "PASS",
  scope = manifest$scope, backend = "multicross::mcm",
  scpa_version = as.character(utils::packageVersion("SCPA")),
  multicross_version = as.character(utils::packageVersion("multicross")),
  manifest_experiment_count = manifest_experiment_count,
  selected_experiment_count = total,
  partial_run = total < manifest_experiment_count,
  completed_experiments = completed, reused_checkpoints = reused,
  mcm_count = completed_mcm, failed_mcm_calls = 0L,
  masking_library = "scripts/scpa/gene_masking_lib.R",
  vanilla_mask = "zero expression column",
  genept_mask = "non-L2 exact outer-product subtraction",
  ground_truth_read = FALSE, checkpoint_resume = TRUE, cores = cores
)
dir.create(dirname(output_json), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(payload, output_json, auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(sprintf("PHASE7_LLMFREE_CORE status=%s experiments=%d/%d mcm=%d reused=%d\n",
            payload$status, completed, manifest_experiment_count, completed_mcm, reused))
