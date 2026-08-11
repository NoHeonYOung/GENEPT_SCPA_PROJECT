script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve reproduce_scpa_phase1b_reference.R location")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))

source(file.path(project_root, "scripts", "data", "phase1_validation_lib.R"))
source(file.path(dirname(script_path), "phase1b_scpa_lib.R"))
source(file.path(dirname(script_path), "phase1b_reference_lib.R"))
source(file.path(dirname(script_path), "phase1b_plot_lib.R"))
require_phase1b_packages()

config <- yaml::read_yaml(file.path(project_root, "config", "genept_scpa.yaml"))
if (config$project$active_phase != 1L || config$phase1$status != "in_progress") {
  stop("Phase 1 must be active and in_progress")
}
if (config$phase1$dataset_gate_status != "passed") stop("Phase 1A dataset gate is not passed")

processed_dir <- file.path(project_root, "data", "processed", "genept_scpa", "phase1")
interim_dir <- file.path(project_root, "data", "interim", "genept_scpa")
sampling_dir <- file.path(interim_dir, "phase1b_sampling")
figure_dir <- file.path(processed_dir, "figures")
qc_path <- file.path(interim_dir, "phase1b_scpa_qc.json")
summary_path <- file.path(interim_dir, "phase1b_reproduction_summary.md")
if (!file.exists(qc_path)) stop("Missing completed Hour-only Phase 1B QC JSON: ", qc_path)

qc <- jsonlite::fromJSON(qc_path, simplifyVector = FALSE)
results <- load_phase1b_hour_results(processed_dir)
hour_analysis_ids <- vapply(phase1b_analysis_specs(), function(x) x$id, character(1))
hour_analyses_pass <- all(vapply(
  hour_analysis_ids,
  function(id) identical(qc$analyses[[id]]$status, "PASS"),
  logical(1)
))
if (!hour_analyses_pass) stop("All four retained Hour-only analyses must have PASS QC")

spec <- phase1b_reference_spec()
scpa_config <- config$phase1$scpa
seed <- as.integer(scpa_config$seed)
pathway_path <- file.path(project_root, config$phase1$phase1b$pathways$file)
if (!file.exists(pathway_path)) stop("Missing official pathway file: ", pathway_path)
observed_pathway_sha <- sha256_file(pathway_path)
if (!identical(observed_pathway_sha, config$phase1$phase1b$pathways$sha256)) {
  stop("Official pathway SHA-256 mismatch")
}

qc$statistics <- list(
  primary_statistic = "SCPA qval",
  qval_convention = "larger qval = stronger multivariate pathway difference",
  pathway_ranking = "descending_qval",
  rank_1 = "largest_qval",
  qval_zero_interpretation = "weakest_end"
)
qc$analyses[[spec$id]] <- list(
  status = "RUNNING",
  population_1_definition = "Cell_Type=Resting AND Hour=0",
  population_2_definition = "Cell_Type=Activated AND Hour=24",
  parameter_tuning = FALSE
)

reference_error <- NULL
reference_outcome <- tryCatch({
  rds_path <- file.path(
    project_root, "data", "raw", "genept_scpa", config$phase1$dataset$extracted_filename
  )
  loaded <- load_phase1b_object(rds_path)
  object <- loaded$object
  metadata <- object[[]]
  cell_ids <- colnames(object)
  labels <- resolve_reference_population_labels(metadata)
  selection <- select_reference_cells(
    metadata = metadata,
    cell_ids = cell_ids,
    labels = labels,
    downsample = as.integer(scpa_config$downsample),
    seed = seed
  )
  sampling_files <- write_reference_sampling_files(selection, sampling_dir)
  matrices <- extract_reference_matrices(
    object = object,
    selection = selection,
    labels = labels,
    assay = config$phase1$phase1b$expression$assay,
    pseudocount = config$phase1$phase1b$expression$pseudocount
  )
  if (ncol(matrices$resting_0h) != length(selection$population_1) ||
      ncol(matrices$activated_24h) != length(selection$population_2)) {
    stop("Reference extraction cell counts do not match the saved sampling IDs")
  }
  rm(object, metadata)
  invisible(gc())

  outcome <- run_one_scpa_analysis(
    matrices = matrices,
    pathway_file = pathway_path,
    scpa_config = scpa_config,
    seed = seed,
    pairwise = TRUE
  )
  failures <- outcome$failures
  expected_rows <- as.integer(qc$pathways$analyzed_pathway_count)
  if (nrow(outcome$result) != expected_rows) failures <- c(failures, "result_row_count")
  if (length(failures) > 0L) {
    stop("Reference result validation failed: ", paste(unique(failures), collapse = ", "))
  }

  output_path <- file.path(processed_dir, spec$output)
  write_csv_atomic(outcome$result, output_path)
  figure_path <- render_phase1b_reference_figure(outcome$result, figure_dir)
  target_summary <- summarize_reference_targets(outcome$result)
  list(
    result = outcome$result,
    elapsed_seconds = outcome$elapsed_seconds,
    warnings = outcome$warnings,
    labels = labels,
    selection = selection,
    sampling_files = sampling_files,
    matrices = lapply(matrices, function(matrix) list(genes = nrow(matrix), cells = ncol(matrix))),
    output_path = output_path,
    figure_path = figure_path,
    target_summary = target_summary,
    compatibility = loaded$compatibility
  )
}, error = function(error_condition) {
  reference_error <<- conditionMessage(error_condition)
  NULL
})

if (is.null(reference_outcome)) {
  qc$analyses[[spec$id]] <- list(
    status = "FAIL",
    population_1_definition = "Cell_Type=Resting AND Hour=0",
    population_2_definition = "Cell_Type=Activated AND Hour=24",
    parameter_tuning = FALSE,
    error = reference_error
  )
  qc$gate <- list(
    status = "NEEDS_REVIEW",
    failed_checks = c("reference_resting0_vs_activated24:execution_or_population_error"),
    warnings = list(),
    criteria = list(
      hour_only_analyses_pass = hour_analyses_pass,
      qval_interpretation_corrected = TRUE,
      reference_execution_pass = FALSE,
      parameter_tuning_absent = TRUE
    )
  )
  qc$generated_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  write_qc_json(qc, qc_path)
  write_phase1b_summary(results, qc, summary_path)
  cat(
    "PHASE1B_REFERENCE_SUMMARY status=NEEDS_REVIEW error=",
    shQuote(reference_error),
    " qc_json=", normalizePath(qc_path), "\n",
    sep = ""
  )
  quit(save = "no", status = 1L)
}

target_summary <- reference_outcome$target_summary
targets_present <- all(vapply(target_summary$targets, function(x) isTRUE(x$present), logical(1)))
no_runtime_warnings <- length(reference_outcome$warnings) == 0L
finite_qval_count <- sum(is.finite(reference_outcome$result$qval))
finite_fc_count <- sum(is.finite(reference_outcome$result$FC))
all_finite <- finite_qval_count == nrow(reference_outcome$result) &&
  finite_fc_count == nrow(reference_outcome$result)

qc$expression$matrix_dimensions$reference_resting_0h <- reference_outcome$matrices$resting_0h
qc$expression$matrix_dimensions$reference_activated_24h <- reference_outcome$matrices$activated_24h
qc$analyses[[spec$id]] <- list(
  status = if (all_finite) "PASS" else "FAIL",
  population_1_definition = paste0(
    "Cell_Type=", reference_outcome$labels$resting_label, " AND Hour=0"
  ),
  population_2_definition = paste0(
    "Cell_Type=", reference_outcome$labels$activated_label, " AND Hour=24"
  ),
  metadata_columns = list(
    cell_type = reference_outcome$labels$cell_type_column,
    hour = reference_outcome$labels$hour_column
  ),
  resolved_metadata_labels = list(
    population_1_cell_type = reference_outcome$labels$resting_label,
    population_2_cell_type = reference_outcome$labels$activated_label
  ),
  full_cell_count_population_1 = reference_outcome$selection$full_counts$population_1,
  full_cell_count_population_2 = reference_outcome$selection$full_counts$population_2,
  actual_sampled_cells_population_1 = length(reference_outcome$selection$population_1),
  actual_sampled_cells_population_2 = length(reference_outcome$selection$population_2),
  seed = seed,
  assay = config$phase1$phase1b$expression$assay,
  layer = config$phase1$phase1b$expression$layer_or_slot,
  extraction_function = config$phase1$phase1b$expression$extraction_function,
  pseudocount = config$phase1$phase1b$expression$pseudocount,
  pathway_count = nrow(reference_outcome$result),
  finite_qval_count = finite_qval_count,
  finite_fc_count = finite_fc_count,
  elapsed_seconds = reference_outcome$elapsed_seconds,
  output_file = reference_outcome$output_path,
  sampling_files = reference_outcome$sampling_files,
  figure_file = reference_outcome$figure_path,
  warnings = reference_outcome$warnings,
  target_summary = target_summary,
  qualitative_agreement = target_summary$qualitative_agreement,
  parameter_tuning = FALSE
)
qc$compatibility$reference_run <- reference_outcome$compatibility
existing_figure_files <- unlist(qc$visualizations$files %||% list(), use.names = FALSE)
qc$visualizations$status <- "PASS"
qc$visualizations$files <- as.list(unique(c(existing_figure_files, reference_outcome$figure_path)))
qc$visualizations$reference_scope <- "official_two_population_qualitative_reproduction"

criteria <- list(
  hour_only_analyses_pass = hour_analyses_pass,
  qval_interpretation_corrected = TRUE,
  reference_execution_pass = identical(qc$analyses[[spec$id]]$status, "PASS"),
  reference_qval_and_fc_finite = all_finite,
  official_workflow_preserved = TRUE,
  arachidonic_targets_present = targets_present,
  qualitative_comparison_recorded = isTRUE(target_summary$comparison_recorded),
  no_critical_runtime_warning_or_error = no_runtime_warnings,
  parameter_tuning_absent = TRUE
)
gate_pass <- all(vapply(criteria, isTRUE, logical(1)))
failed_checks <- names(criteria)[!vapply(criteria, isTRUE, logical(1))]
qc$gate <- list(
  status = if (gate_pass) "PASS" else "NEEDS_REVIEW",
  failed_checks = unname(failed_checks),
  warnings = unname(reference_outcome$warnings),
  criteria = criteria,
  pass_basis = "successful pipeline execution plus qualitative reference reproduction, not numerical identity"
)
qc$generated_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
results[[spec$id]] <- reference_outcome$result
write_qc_json(qc, qc_path)
write_phase1b_summary(results, qc, summary_path)

cat(
  "PHASE1B_REFERENCE_SUMMARY",
  paste0("status=", qc$gate$status),
  paste0("population_1_full=", reference_outcome$selection$full_counts$population_1),
  paste0("population_2_full=", reference_outcome$selection$full_counts$population_2),
  paste0("sampled_each=", length(reference_outcome$selection$population_1), "/", length(reference_outcome$selection$population_2)),
  paste0("targets_present=", targets_present),
  paste0("qualitative_agreement=", target_summary$qualitative_agreement),
  paste0("qc_json=", normalizePath(qc_path)),
  sep = " "
)
cat("\n")
quit(save = "no", status = if (gate_pass) 0L else 1L)
