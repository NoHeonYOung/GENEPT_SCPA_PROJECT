args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1L) args[[1]] else ".")
source(file.path(project_root, "scripts", "data", "phase1_validation_lib.R"))
source(file.path(project_root, "scripts", "scpa", "phase1b_scpa_lib.R"))
source(file.path(project_root, "scripts", "scpa", "phase1b_reference_lib.R"))
source(file.path(project_root, "scripts", "scpa", "phase1b_plot_lib.R"))

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}
expect_error <- function(expression, pattern) {
  message <- tryCatch({ force(expression); NULL }, error = function(e) conditionMessage(e))
  assert_true(!is.null(message) && grepl(pattern, message), paste0("Expected error matching: ", pattern))
}

require_phase1b_packages()
config <- yaml::read_yaml(file.path(project_root, "config", "genept_scpa.yaml"))
assert_true(config$phase1$scpa$seed == 20260810L, "Phase 1B seed is not frozen")
assert_true(config$phase1$dataset_gate_status == "passed", "Phase 1A gate must remain passed")
assert_true(config$phase1$status == "passed", "Phase 1B must remain passed")
assert_true(config$phase1$stage == "phase1b_passed", "Phase 1B stage is incorrect")
official_pathway_file <- file.path(project_root, config$phase1$phase1b$pathways$file)
assert_true(file.exists(official_pathway_file), "Official Phase 1B pathway file is missing")
assert_true(
  sha256_file(official_pathway_file) == config$phase1$phase1b$pathways$sha256,
  "Official Phase 1B pathway SHA-256 mismatch"
)
official_pathways <- read_pathways_official(official_pathway_file)
assert_true(length(official_pathways) == 243L, "Official pathway collection must contain 243 pathways")

assert_true(parse_phase1b_args(character()) == "all", "Default analysis must be all")
assert_true(parse_phase1b_args(c("--analysis", "global")) == "global", "Plain CLI analysis parsing failed")
assert_true(parse_phase1b_args("--analysis=0_vs_12") == "0_vs_12", "Equals CLI analysis parsing failed")
expect_error(parse_phase1b_args("--analysis=invalid"), "Invalid --analysis")

ranking_fixture <- data.frame(
  Pathway = c("ZERO", "MAX", "MID"),
  Pval = c(1, 0.01, 0.2),
  adjPval = c(1, 0.02, 0.3),
  qval = c(0, 9, 3),
  FC = c(0, 1, -1),
  stringsAsFactors = FALSE
)
ranked_fixture <- rank_scpa_result(ranking_fixture)
assert_true(ranked_fixture$Pathway[[1]] == "MAX", "Rank 1 must have the maximum SCPA qval")
assert_true(ranked_fixture$Rank[[1]] == 1L, "Descending qval rank numbering failed")
assert_true(tail(ranked_fixture$Pathway, 1L) == "ZERO", "qval=0 must be at the weakest end")
ranking_markdown <- markdown_result_table(ranking_fixture)
assert_true(
  regexpr("MAX", ranking_markdown)[[1]] < regexpr("ZERO", ranking_markdown)[[1]],
  "Markdown result table is not sorted by descending qval"
)

metadata <- data.frame(
  Hour = rep(c(0L, 12L, 24L), each = 7L),
  row.names = paste0("cell", seq_len(21L))
)
sampling_a <- select_timepoint_cells(metadata, rownames(metadata), c(0L, 12L, 24L), 5L, 11L)
sampling_b <- select_timepoint_cells(metadata, rownames(metadata), c(0L, 12L, 24L), 5L, 11L)
assert_true(all(vapply(sampling_a$selected, length, integer(1)) == 5L), "Hour filtering/downsampling failed")
assert_true(identical(sampling_a$selected, sampling_b$selected), "Sampling is not seed deterministic")
expect_error(
  select_timepoint_cells(metadata[metadata$Hour != 24L, , drop = FALSE], rownames(metadata)[metadata$Hour != 24L], c(0L, 12L, 24L), 5L, 11L),
  "No cells found for Hour == 24"
)

reference_metadata <- metadata
reference_metadata$Cell_Type <- "Intermediate"
reference_metadata$Cell_Type[reference_metadata$Hour == 0L] <- "Resting"
reference_metadata$Cell_Type[reference_metadata$Hour == 24L] <- "Activated"
reference_labels <- resolve_reference_population_labels(reference_metadata)
assert_true(reference_labels$resting_label == "Resting", "Resting metadata label was not resolved")
assert_true(reference_labels$activated_label == "Activated", "Activated metadata label was not resolved")
reference_selection_a <- select_reference_cells(
  reference_metadata, rownames(reference_metadata), reference_labels, downsample = 4L, seed = 29L
)
reference_selection_b <- select_reference_cells(
  reference_metadata, rownames(reference_metadata), reference_labels, downsample = 4L, seed = 29L
)
assert_true(identical(reference_selection_a, reference_selection_b), "Reference sampling is not deterministic")
assert_true(length(reference_selection_a$population_1) == 4L, "Resting 0 h filtering failed")
assert_true(length(reference_selection_a$population_2) == 4L, "Activated 24 h filtering failed")
expect_error(resolve_reference_population_labels(metadata), "Missing required metadata column.*Cell_Type")

mock_counts <- matrix(
  seq_len(8L * 21L) %% 5L,
  nrow = 8L,
  dimnames = list(paste0("G", 1:8), rownames(metadata))
)
mock_object <- SeuratObject::CreateSeuratObject(mock_counts, project = "phase1b_mock")
mock_object <- SeuratObject::SetAssayData(
  mock_object,
  assay = "RNA",
  layer = "data",
  new.data = log1p(mock_counts)
)
mock_object$Hour <- metadata$Hour
mock_object$Cell_Type <- reference_metadata$Cell_Type
extracted <- extract_phase1b_matrices(
  mock_object,
  selected_cells = sampling_a$selected,
  hours = c(0L, 12L, 24L),
  assay = "RNA",
  pseudocount = 0.001
)
assert_true(all(vapply(extracted, nrow, integer(1)) == 8L), "RNA data feature extraction failed")
assert_true(all(vapply(extracted, ncol, integer(1)) == 5L), "Hour-filtered cell extraction failed")

reference_extracted <- extract_reference_matrices(
  mock_object,
  selection = reference_selection_a,
  labels = reference_labels,
  assay = "RNA",
  pseudocount = 0.001
)
assert_true(
  identical(vapply(reference_extracted, ncol, integer(1)), c(resting_0h = 4L, activated_24h = 4L)),
  "Dual metadata-filtered reference extraction failed"
)

test_root <- tempfile("phase1b_tests_")
dir.create(test_root)
on.exit(unlink(test_root, recursive = TRUE), add = TRUE)
missing_path <- file.path(test_root, "missing.csv")
expect_error(read_pathways_official(missing_path), ".")

pathway_file <- file.path(test_root, "mock_pathways.csv")
writeLines(
  c(
    "PATH_A,G1,G2,G3,G4",
    "PATH_B,G3,G4,G5,G6",
    "PATH_C,G5,G6,G7,G8"
  ),
  pathway_file
)
pathways <- read_pathways_official(pathway_file)
pathway_qc <- inspect_pathways(pathways, paste0("G", 1:8), min_genes = 2L, max_genes = 10L)
assert_true(pathway_qc$input_pathway_count == 3L, "Input pathway count is wrong")
assert_true(pathway_qc$analyzed_pathway_count == 3L, "Analyzed pathway count is wrong")
assert_true(pathway_qc$unmatched_unique_genes == 0L, "Unexpected unmatched genes")

set.seed(9)
make_matrix <- function(name, shift) {
  matrix(
    stats::rnorm(8L * 20L, mean = shift),
    nrow = 8L,
    dimnames = list(paste0("G", 1:8), paste0(name, "_cell", 1:20))
  )
}
matrices <- list(
  `0h` = make_matrix("h0", 0),
  `12h` = make_matrix("h12", 0.2),
  `24h` = make_matrix("h24", 0.5)
)
mock_scpa_config <- list(
  downsample = 10L,
  min_genes = 2L,
  max_genes = 10L,
  parallel = FALSE,
  cores = NULL
)

global <- run_one_scpa_analysis(matrices, pathway_file, mock_scpa_config, seed = 17L, pairwise = FALSE)
assert_true(length(global$failures) == 0L, "Three-sample SCPA schema/finite checks failed")
assert_true(is.character(global$warnings), "Global runtime warnings were not captured")
assert_true(identical(colnames(global$result), c("Pathway", "Pval", "adjPval", "qval")), "Global schema is wrong")

pairwise <- run_one_scpa_analysis(matrices[c("0h", "12h")], pathway_file, mock_scpa_config, seed = 17L, pairwise = TRUE)
assert_true(length(pairwise$failures) == 0L, "Two-sample SCPA schema/finite checks failed")
assert_true(is.character(pairwise$warnings), "Pairwise runtime warnings were not captured")
assert_true(identical(colnames(pairwise$result), c("Pathway", "Pval", "adjPval", "qval", "FC")), "Pairwise schema is wrong")

reference_fixture <- data.frame(
  Pathway = c(
    "REACTOME_ARACHIDONIC_ACID_METABOLISM",
    "KEGG_ARACHIDONIC_ACID_METABOLISM",
    "MOCK_OTHER"
  ),
  Pval = c(0.001, 0.01, 0.5),
  adjPval = c(0.003, 0.02, 0.5),
  qval = c(8, 5, 0),
  FC = c(0.4, -1.2, 0),
  stringsAsFactors = FALSE
)
assert_true(length(validate_scpa_result(reference_fixture, pairwise = TRUE)) == 0L, "Reference schema failed")
reference_targets <- summarize_reference_targets(reference_fixture)
assert_true(
  all(vapply(reference_targets$targets, function(x) isTRUE(x$present), logical(1))),
  "Arachidonic reference targets were not summarized"
)
assert_true(
  reference_targets$targets$REACTOME_ARACHIDONIC_ACID_METABOLISM$rank == 1L,
  "Reference target qval rank is wrong"
)

figure_results <- list(
  global_0_12_24 = global$result,
  pairwise_0_vs_12 = pairwise$result,
  pairwise_12_vs_24 = pairwise$result,
  pairwise_0_vs_24 = pairwise$result
)
figure_dir <- file.path(test_root, "figures")
figure_files <- render_phase1b_figures(figure_results, figure_dir, top_n = 3L)
assert_true(length(figure_files) == 6L, "Phase 1B figure output count is wrong")
assert_true(all(file.exists(figure_files)), "One or more Phase 1B figures were not generated")
assert_true(all(file.info(figure_files)$size > 0L), "One or more Phase 1B figures are empty")
reference_figure <- render_phase1b_reference_figure(reference_fixture, figure_dir)
assert_true(file.exists(reference_figure) && file.info(reference_figure)$size > 0L, "Reference figure was not generated")

qc <- list(
  dataset = list(accession = "MOCK", total_cells = 60L, total_features = 8L),
  expression = list(
    assay = "RNA",
    layer_or_slot = "data",
    matrix_dimensions = list(`0h` = list(genes = 8L, cells = 20L))
  ),
  pathways = pathway_qc,
  analyses = list(
    global_0_12_24 = list(result_rows = nrow(global$result), finite_qval_count = sum(is.finite(global$result$qval))),
    pairwise_0_vs_12 = list(result_rows = nrow(pairwise$result), finite_qval_count = sum(is.finite(pairwise$result$qval)), finite_fc_count = sum(is.finite(pairwise$result$FC))),
    reference_resting0_vs_activated24 = list(
      status = "PASS",
      population_1_definition = "Cell_Type=Resting AND Hour=0",
      population_2_definition = "Cell_Type=Activated AND Hour=24",
      full_cell_count_population_1 = 7L,
      full_cell_count_population_2 = 7L,
      actual_sampled_cells_population_1 = 4L,
      actual_sampled_cells_population_2 = 4L,
      qualitative_agreement = "QUALITATIVELY_CONSISTENT",
      parameter_tuning = FALSE
    )
  ),
  compatibility = list(serialized_seurat_version = "mock", current_object_version = "mock"),
  reproducibility = list(seed = 17L, scpa_version = as.character(packageVersion("SCPA"))),
  gate = list(status = "INCOMPLETE", failed_checks = character(), warnings = "mock")
)
qc_path <- file.path(test_root, "phase1b_scpa_qc.json")
write_qc_json(qc, qc_path)
roundtrip <- jsonlite::fromJSON(qc_path, simplifyVector = FALSE)
required_sections <- c("dataset", "expression", "pathways", "analyses", "reproducibility", "gate")
assert_true(all(required_sections %in% names(roundtrip)), "Phase 1B QC JSON schema is incomplete")
reference_qc_fields <- c(
  "population_1_definition", "population_2_definition",
  "full_cell_count_population_1", "full_cell_count_population_2",
  "actual_sampled_cells_population_1", "actual_sampled_cells_population_2",
  "qualitative_agreement", "parameter_tuning"
)
assert_true(
  all(reference_qc_fields %in% names(roundtrip$analyses$reference_resting0_vs_activated24)),
  "Reference QC JSON schema is incomplete"
)

summary_path <- file.path(test_root, "phase1b_reproduction_summary.md")
write_phase1b_summary(
  list(
    global_0_12_24 = global$result,
    pairwise_0_vs_12 = pairwise$result,
    reference_resting0_vs_activated24 = reference_fixture
  ),
  qc,
  summary_path
)
assert_true(file.exists(summary_path) && file.info(summary_path)$size > 0L, "Reproduction summary was not generated")

hour_output_names <- vapply(phase1b_analysis_specs(), function(x) x$output, character(1))
assert_true(
  all(file.exists(file.path(project_root, "data", "processed", "genept_scpa", "phase1", hour_output_names))),
  "One or more retained Hour-only result files are missing"
)
assert_true(
  !phase1b_reference_spec()$output %in% hour_output_names,
  "Reference output must not overwrite an Hour-only result"
)

cat("Phase 1B mock tests: PASS\n")
