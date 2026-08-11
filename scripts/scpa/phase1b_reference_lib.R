phase1b_reference_spec <- function() {
  list(
    id = "reference_resting0_vs_activated24",
    output = "05_reference_resting0_vs_activated24.csv",
    population_1 = list(cell_type = "Resting", hour = 0L),
    population_2 = list(cell_type = "Activated", hour = 24L)
  )
}

match_metadata_value <- function(values, requested, column) {
  observed <- unique(as.character(values[!is.na(values)]))
  matches <- observed[tolower(trimws(observed)) == tolower(trimws(requested))]
  if (length(matches) != 1L) {
    stop(
      "Could not uniquely resolve metadata value '", requested, "' in ", column,
      ". Observed values: ", paste(sort(observed), collapse = ", ")
    )
  }
  matches[[1]]
}

resolve_reference_population_labels <- function(metadata) {
  required <- c("Hour", "Cell_Type")
  missing <- setdiff(required, colnames(metadata))
  if (length(missing) > 0L) {
    stop("Missing required metadata column(s): ", paste(missing, collapse = ", "))
  }
  list(
    hour_column = "Hour",
    cell_type_column = "Cell_Type",
    resting_label = match_metadata_value(metadata$Cell_Type, "Resting", "Cell_Type"),
    activated_label = match_metadata_value(metadata$Cell_Type, "Activated", "Cell_Type")
  )
}

select_reference_cells <- function(metadata, cell_ids, labels, downsample, seed) {
  if (nrow(metadata) != length(cell_ids)) stop("Metadata and cell IDs are not aligned")
  if (!is.null(rownames(metadata)) && !identical(rownames(metadata), cell_ids)) {
    stop("Metadata row names and cell IDs are not in the same order")
  }
  spec <- phase1b_reference_spec()
  population_1_candidates <- cell_ids[
    as.character(metadata[[labels$hour_column]]) == as.character(spec$population_1$hour) &
      as.character(metadata[[labels$cell_type_column]]) == labels$resting_label
  ]
  population_2_candidates <- cell_ids[
    as.character(metadata[[labels$hour_column]]) == as.character(spec$population_2$hour) &
      as.character(metadata[[labels$cell_type_column]]) == labels$activated_label
  ]
  if (length(population_1_candidates) == 0L) {
    stop("No cells found for Cell_Type=", labels$resting_label, " AND Hour=0")
  }
  if (length(population_2_candidates) == 0L) {
    stop("No cells found for Cell_Type=", labels$activated_label, " AND Hour=24")
  }
  set.seed(seed)
  population_1 <- sample(
    population_1_candidates,
    min(length(population_1_candidates), as.integer(downsample)),
    replace = FALSE
  )
  population_2 <- sample(
    population_2_candidates,
    min(length(population_2_candidates), as.integer(downsample)),
    replace = FALSE
  )
  list(
    population_1 = population_1,
    population_2 = population_2,
    full_counts = list(
      population_1 = length(population_1_candidates),
      population_2 = length(population_2_candidates)
    )
  )
}

write_lines_atomic <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  writeLines(lines, temporary)
  if (!file.rename(temporary, path)) stop("Could not atomically write text file: ", path)
  path
}

write_reference_sampling_files <- function(selection, directory) {
  list(
    population_1 = write_lines_atomic(
      selection$population_1,
      file.path(directory, "reference_resting0_cells.txt")
    ),
    population_2 = write_lines_atomic(
      selection$population_2,
      file.path(directory, "reference_activated24_cells.txt")
    )
  )
}

extract_reference_matrices <- function(object, selection, labels, assay, pseudocount) {
  selected <- unique(c(selection$population_1, selection$population_2))
  subset_object <- object[, selected]
  spec <- phase1b_reference_spec()
  population_1 <- SCPA::seurat_extract(
    subset_object,
    assay = assay,
    meta1 = labels$cell_type_column,
    value_meta1 = labels$resting_label,
    meta2 = labels$hour_column,
    value_meta2 = spec$population_1$hour,
    pseudocount = pseudocount
  )
  population_2 <- SCPA::seurat_extract(
    subset_object,
    assay = assay,
    meta1 = labels$cell_type_column,
    value_meta1 = labels$activated_label,
    meta2 = labels$hour_column,
    value_meta2 = spec$population_2$hour,
    pseudocount = pseudocount
  )
  list(resting_0h = population_1, activated_24h = population_2)
}

load_phase1b_hour_results <- function(processed_dir) {
  results <- list()
  for (spec in phase1b_analysis_specs()) {
    input_path <- file.path(processed_dir, spec$output)
    if (!file.exists(input_path)) stop("Missing retained Hour-only result CSV: ", input_path)
    result <- utils::read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
    failures <- validate_scpa_result(result, pairwise = spec$pairwise)
    if (length(failures) > 0L) {
      stop("Invalid retained Hour-only result CSV ", spec$output, ": ", paste(failures, collapse = ", "))
    }
    results[[spec$id]] <- result
  }
  results
}

summarize_reference_targets <- function(result) {
  targets <- c(
    "REACTOME_ARACHIDONIC_ACID_METABOLISM",
    "KEGG_ARACHIDONIC_ACID_METABOLISM"
  )
  ranked <- rank_scpa_result(result)
  summaries <- setNames(lapply(targets, function(target) {
    row <- ranked[ranked$Pathway == target, , drop = FALSE]
    if (nrow(row) == 0L) {
      return(list(present = FALSE, rank = NULL, qval = NULL, FC = NULL))
    }
    list(
      present = TRUE,
      rank = as.integer(row$Rank[[1]]),
      qval = unname(row$qval[[1]]),
      FC = unname(row$FC[[1]])
    )
  }), targets)
  high_qval_cutoff <- unname(stats::quantile(result$qval, 0.75, na.rm = TRUE))
  high_qval_modest_fc <- sum(result$qval >= high_qval_cutoff & abs(result$FC) <= 5, na.rm = TRUE)
  reactome <- summaries[["REACTOME_ARACHIDONIC_ACID_METABOLISM"]]
  comparison <- if (!all(vapply(summaries, function(x) isTRUE(x$present), logical(1)))) {
    "TARGET_PATHWAY_MISSING"
  } else if (isTRUE(reactome$qval > 0) && isTRUE(abs(reactome$FC) <= 5) && high_qval_modest_fc > 0L) {
    "QUALITATIVELY_CONSISTENT"
  } else {
    "QUALITATIVE_DIFFERENCE_RECORDED"
  }
  list(
    targets = summaries,
    high_qval_definition = "qval at or above the result-set 75th percentile",
    high_qval_cutoff = high_qval_cutoff,
    high_qval_modest_fc_count = high_qval_modest_fc,
    high_qval_does_not_require_large_fc = high_qval_modest_fc > 0L,
    qualitative_agreement = comparison,
    comparison_recorded = TRUE,
    numerical_identity_required = FALSE
  )
}
