args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve reproduce_scpa_phase1b.R location")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))

source(file.path(project_root, "scripts", "data", "phase1_validation_lib.R"))
source(file.path(dirname(script_path), "phase1b_scpa_lib.R"))
source(file.path(dirname(script_path), "phase1b_plot_lib.R"))
require_phase1b_packages()
selection <- parse_phase1b_args(args)
specs <- selected_analysis_specs(selection)

config <- yaml::read_yaml(file.path(project_root, "config", "genept_scpa.yaml"))
if (config$project$active_phase != 1L || config$phase1$status != "in_progress") {
  stop("Phase 1 must be active and in_progress")
}
if (config$phase1$dataset_gate_status != "passed") stop("Phase 1A dataset gate is not passed")
if (config$phase1$stage != "phase1b_ready_for_user_run") {
  stop("Phase 1B config stage is not ready_for_user_run")
}

phase1a_qc_path <- file.path(project_root, "data", "interim", "genept_scpa", "phase1_dataset_qc.json")
phase1a_qc <- jsonlite::fromJSON(phase1a_qc_path, simplifyVector = FALSE)
if (!identical(phase1a_qc$gate$status, "PASS")) stop("Phase 1A QC JSON does not report PASS")

rds_path <- file.path(project_root, "data", "raw", "genept_scpa", config$phase1$dataset$extracted_filename)
pathway_relative <- config$phase1$phase1b$pathways$file
pathway_path <- file.path(project_root, pathway_relative)
if (!file.exists(pathway_path)) stop("Missing official pathway file: ", pathway_path)
observed_pathway_sha <- sha256_file(pathway_path)
if (!identical(observed_pathway_sha, config$phase1$phase1b$pathways$sha256)) {
  stop("Official pathway SHA-256 mismatch")
}

loaded <- load_phase1b_object(rds_path)
object <- loaded$object
metadata <- object[[]]
cell_ids <- colnames(object)
hours <- as.integer(unlist(config$phase1$phase1b$time_values))
scpa_config <- config$phase1$scpa
seed <- as.integer(scpa_config$seed)

sampling <- select_timepoint_cells(
  metadata = metadata,
  cell_ids = cell_ids,
  hours = hours,
  downsample = as.integer(scpa_config$downsample),
  seed = seed
)
hour_names <- paste0(hours, "h")
expected_counts <- vapply(
  hour_names,
  function(hour_name) as.integer(phase1a_qc$metadata$cells_per_timepoint[[hour_name]]),
  integer(1)
)
observed_counts <- vapply(sampling$input_counts[hour_names], as.integer, integer(1))
if (!identical(unname(observed_counts), unname(expected_counts))) {
  stop("Current Hour cell counts do not match the passed Phase 1A QC report")
}
selected_hours <- sort(unique(unlist(lapply(specs, function(x) x$hours))))
sampling_dir <- file.path(project_root, "data", "interim", "genept_scpa", "phase1b_sampling")
sampling_files <- write_sampling_files(sampling$selected, specs, sampling_dir)

matrices <- extract_phase1b_matrices(
  object = object,
  selected_cells = sampling$selected,
  hours = selected_hours,
  assay = config$phase1$phase1b$expression$assay,
  pseudocount = config$phase1$phase1b$expression$pseudocount
)
rm(object, metadata)
invisible(gc())

pathways <- read_pathways_official(pathway_path)
pathway_qc <- inspect_pathways(
  pathways,
  dataset_genes = rownames(matrices[[1]]),
  min_genes = as.integer(scpa_config$min_genes),
  max_genes = as.integer(scpa_config$max_genes)
)

processed_dir <- file.path(project_root, "data", "processed", "genept_scpa", "phase1")
interim_dir <- file.path(project_root, "data", "interim", "genept_scpa")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(interim_dir, recursive = TRUE, showWarnings = FALSE)

results <- list()
analysis_qc <- list()
failed_checks <- character()
warnings <- character()

for (spec in specs) {
  matrix_names <- paste0(spec$hours, "h")
  analysis_matrices <- matrices[matrix_names]
  outcome <- tryCatch(
    run_one_scpa_analysis(
      analysis_matrices,
      pathway_file = pathway_path,
      scpa_config = scpa_config,
      seed = seed,
      pairwise = spec$pairwise
    ),
    error = function(error_condition) error_condition
  )
  if (inherits(outcome, "error")) {
    failed_checks <- c(failed_checks, paste0(spec$id, ":execution_error"))
    analysis_qc[[spec$id]] <- list(status = "FAIL", error = conditionMessage(outcome))
    next
  }
  analysis_failures <- outcome$failures
  if (nrow(outcome$result) != pathway_qc$analyzed_pathway_count) {
    analysis_failures <- c(analysis_failures, "result_row_count")
  }
  if (length(analysis_failures) > 0L) {
    failed_checks <- c(failed_checks, paste0(spec$id, ":", analysis_failures))
  }
  if (length(outcome$warnings) > 0L) {
    warnings <- c(warnings, paste0(spec$id, ":", outcome$warnings))
  }
  output_path <- file.path(processed_dir, spec$output)
  write_csv_atomic(outcome$result, output_path)
  results[[spec$id]] <- outcome$result
  analysis_qc[[spec$id]] <- list(
    status = if (length(analysis_failures) == 0L) "PASS" else "FAIL",
    input_cell_counts = sampling$input_counts[matrix_names],
    actual_cell_counts = as.list(vapply(analysis_matrices, ncol, integer(1))),
    result_rows = nrow(outcome$result),
    finite_qval_count = sum(is.finite(outcome$result$qval)),
    finite_fc_count = if (spec$pairwise) sum(is.finite(outcome$result$FC)) else NULL,
    warnings = outcome$warnings,
    elapsed_seconds = outcome$elapsed_seconds,
    output_file = output_path,
    sampling_files = unname(sampling_files[grepl(paste0("^", spec$id, "_"), names(sampling_files))])
  )
}

all_requested <- selection == "all"
visualization_qc <- list(
  status = if (all_requested) "SKIPPED" else "NOT_REQUESTED",
  files = list(),
  error = NULL
)
if (all_requested && length(results) == 4L && length(failed_checks) == 0L) {
  figure_dir <- file.path(processed_dir, "figures")
  figure_outcome <- tryCatch(
    render_phase1b_figures(
      results,
      figure_dir,
      top_n = as.integer(config$phase1$phase1b$visualization$heatmap_top_n)
    ),
    error = function(error_condition) error_condition
  )
  if (inherits(figure_outcome, "error")) {
    failed_checks <- c(failed_checks, "visualization:execution_error")
    visualization_qc <- list(
      status = "FAIL",
      files = list(),
      error = conditionMessage(figure_outcome)
    )
  } else {
    visualization_qc <- list(
      status = "PASS",
      files = as.list(unname(figure_outcome)),
      error = NULL,
      comparison_scope = "qualitative_only_hour_groups_differ_from_paper_populations"
    )
  }
}

all_succeeded <- all_requested && length(results) == 4L && length(failed_checks) == 0L
if (!all_requested) warnings <- c(warnings, "partial_analysis_selection_not_phase1b_complete")
gate_status <- if (all_succeeded) {
  "READY_FOR_GPT_REVIEW"
} else if (length(failed_checks) > 0L) {
  "FAIL"
} else {
  "INCOMPLETE"
}

matrix_dimensions <- lapply(matrices, function(matrix) list(genes = nrow(matrix), cells = ncol(matrix)))
qc <- list(
  phase = "Phase 1B - Vanilla SCPA reproduction",
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  dataset = list(
    accession = config$phase1$dataset$accession,
    total_cells = phase1a_qc$object$cells,
    total_features = phase1a_qc$object$features,
    population_definition = config$phase1$phase1b$population_definition,
    full_timepoint_cell_counts = sampling$input_counts
  ),
  expression = list(
    assay = config$phase1$phase1b$expression$assay,
    layer_or_slot = config$phase1$phase1b$expression$layer_or_slot,
    extraction_function = config$phase1$phase1b$expression$extraction_function,
    pseudocount = config$phase1$phase1b$expression$pseudocount,
    matrix_dimensions = matrix_dimensions
  ),
  pathways = c(
    list(
      input_file = pathway_path,
      source_ref = config$phase1$phase1b$pathways$source_ref,
      sha256 = observed_pathway_sha,
      min_genes = scpa_config$min_genes,
      max_genes = scpa_config$max_genes
    ),
    pathway_qc
  ),
  analyses = analysis_qc,
  visualizations = visualization_qc,
  sampling = list(
    mechanism = config$phase1$phase1b$sampling$mechanism,
    same_timepoint_cells_reused = TRUE,
    canonical_cell_id_files = unname(sampling_files[grepl("^canonical_", names(sampling_files))])
  ),
  compatibility = loaded$compatibility,
  reproducibility = list(
    seed = seed,
    scpa_version = as.character(utils::packageVersion("SCPA")),
    R_version = R.version.string,
    Seurat_version = as.character(utils::packageVersion("Seurat")),
    SeuratObject_version = as.character(utils::packageVersion("SeuratObject")),
    parallel = isTRUE(scpa_config$parallel),
    cores = scpa_config$cores %||% NULL
  ),
  gate = list(
    status = gate_status,
    failed_checks = unname(unique(failed_checks)),
    warnings = unname(unique(warnings))
  )
)

qc_path <- file.path(interim_dir, "phase1b_scpa_qc.json")
summary_path <- file.path(interim_dir, "phase1b_reproduction_summary.md")
write_qc_json(qc, qc_path)
write_phase1b_summary(results, qc, summary_path)

cat(
  "PHASE1B_SUMMARY",
  paste0("status=", gate_status),
  paste0("analyses_completed=", length(results)),
  paste0("failed_checks=", length(unique(failed_checks))),
  paste0("qc_json=", normalizePath(qc_path)),
  paste0("summary_md=", normalizePath(summary_path)),
  sep = " "
)
cat("\n")
quit(save = "no", status = if (length(failed_checks) == 0L) 0L else 1L)
