script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve plot_phase1b_results.R location")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))

source(file.path(project_root, "scripts", "data", "phase1_validation_lib.R"))
source(file.path(dirname(script_path), "phase1b_scpa_lib.R"))
source(file.path(dirname(script_path), "phase1b_reference_lib.R"))
source(file.path(dirname(script_path), "phase1b_plot_lib.R"))
require_phase1b_packages()
config <- yaml::read_yaml(file.path(project_root, "config", "genept_scpa.yaml"))

processed_dir <- file.path(project_root, "data", "processed", "genept_scpa", "phase1")
results <- load_phase1b_hour_results(processed_dir)
reference_spec <- phase1b_reference_spec()
reference_path <- file.path(processed_dir, reference_spec$output)
if (file.exists(reference_path)) {
  reference_result <- utils::read.csv(reference_path, check.names = FALSE, stringsAsFactors = FALSE)
  failures <- validate_scpa_result(reference_result, pairwise = TRUE)
  if (length(failures) > 0L) {
    stop("Invalid reference result CSV ", reference_spec$output, ": ", paste(failures, collapse = ", "))
  }
  results[[reference_spec$id]] <- reference_result
}

figure_dir <- file.path(processed_dir, "figures")
files <- render_phase1b_figures(
  results,
  figure_dir,
  top_n = as.integer(config$phase1$phase1b$visualization$heatmap_top_n)
)
if (!is.null(results[[reference_spec$id]])) {
  files <- c(
    files,
    reference_png = render_phase1b_reference_figure(results[[reference_spec$id]], figure_dir)
  )
}

qc_path <- file.path(project_root, "data", "interim", "genept_scpa", "phase1b_scpa_qc.json")
summary_path <- file.path(project_root, "data", "interim", "genept_scpa", "phase1b_reproduction_summary.md")
if (!file.exists(qc_path)) stop("Missing Phase 1B QC JSON: ", qc_path)
qc <- jsonlite::fromJSON(qc_path, simplifyVector = FALSE)
qc$statistics <- list(
  primary_statistic = "SCPA qval",
  qval_convention = "larger qval = stronger multivariate pathway difference",
  pathway_ranking = "descending_qval",
  rank_1 = "largest_qval",
  qval_zero_interpretation = "weakest_end"
)
if (is.null(results[[reference_spec$id]])) {
  qc$gate$status <- "REFERENCE_READY_FOR_USER_RUN"
  qc$gate$failed_checks <- c("reference_resting0_vs_activated24:pending_user_run")
  qc$gate$pass_basis <- "reference run pending; numerical identity is not a PASS requirement"
}
qc$generated_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
write_qc_json(qc, qc_path)
write_phase1b_summary(results, qc, summary_path)

cat(
  "PHASE1B_FIGURES status=PASS files=", length(files),
  " reference_present=", !is.null(results[[reference_spec$id]]),
  " output_dir=", normalizePath(figure_dir), "\n",
  sep = ""
)
