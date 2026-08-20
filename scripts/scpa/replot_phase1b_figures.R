script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve replot script location")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))

source(file.path(project_root, "scripts", "data", "phase1_validation_lib.R"))
source(file.path(dirname(script_path), "phase1b_scpa_lib.R"))
source(file.path(dirname(script_path), "phase1b_reference_lib.R"))
source(file.path(dirname(script_path), "phase1b_plot_lib.R"))
require_phase1b_packages()

config <- yaml::read_yaml(file.path(project_root, "config", "genept_scpa.yaml"))
processed_dir <- file.path(project_root, "data", "processed", "genept_scpa", "phase1")
figure_dir <- file.path(processed_dir, "figures")
results <- load_phase1b_hour_results(processed_dir)
reference_spec <- phase1b_reference_spec()
reference_path <- file.path(processed_dir, reference_spec$output)

if (file.exists(reference_path)) {
  reference_result <- utils::read.csv(reference_path, check.names = FALSE, stringsAsFactors = FALSE)
  failures <- validate_scpa_result(reference_result, pairwise = TRUE)
  if (length(failures) > 0L) stop("Invalid reference result: ", paste(failures, collapse = ", "))
  results[[reference_spec$id]] <- reference_result
}

files <- render_phase1b_figures(
  results,
  figure_dir,
  top_n = as.integer(config$phase1$phase1b$visualization$heatmap_top_n)
)
if (!is.null(results[[reference_spec$id]])) {
  files <- c(files, reference_png = render_phase1b_reference_figure(results[[reference_spec$id]], figure_dir))
}

cat(
  "PHASE1B_REPLOT status=PASS figures=", length(files),
  " analysis_rerun=false output_dir=", normalizePath(figure_dir), "\n",
  sep = ""
)
