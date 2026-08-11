`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

phase1b_analysis_specs <- function() {
  list(
    global = list(
      id = "global_0_12_24",
      hours = c(0L, 12L, 24L),
      output = "01_global_0_12_24.csv",
      pairwise = FALSE
    ),
    `0_vs_12` = list(
      id = "pairwise_0_vs_12",
      hours = c(0L, 12L),
      output = "02_pairwise_0_vs_12.csv",
      pairwise = TRUE
    ),
    `12_vs_24` = list(
      id = "pairwise_12_vs_24",
      hours = c(12L, 24L),
      output = "03_pairwise_12_vs_24.csv",
      pairwise = TRUE
    ),
    `0_vs_24` = list(
      id = "pairwise_0_vs_24",
      hours = c(0L, 24L),
      output = "04_pairwise_0_vs_24.csv",
      pairwise = TRUE
    )
  )
}

parse_phase1b_args <- function(args) {
  analysis <- "all"
  if (length(args) > 0L) {
    equals_arg <- grep("^--analysis=", args, value = TRUE)
    plain_index <- which(args == "--analysis")
    if (length(equals_arg) == 1L) {
      analysis <- sub("^--analysis=", "", equals_arg)
    } else if (length(plain_index) == 1L && plain_index < length(args)) {
      analysis <- args[[plain_index + 1L]]
    } else {
      stop("Usage: --analysis all|global|0_vs_12|12_vs_24|0_vs_24")
    }
  }
  valid <- c("all", names(phase1b_analysis_specs()))
  if (!analysis %in% valid) {
    stop("Invalid --analysis value: ", analysis, ". Valid values: ", paste(valid, collapse = ", "))
  }
  analysis
}

selected_analysis_specs <- function(selection) {
  specs <- phase1b_analysis_specs()
  if (selection == "all") specs else specs[selection]
}

require_phase1b_packages <- function() {
  packages <- c("SCPA", "Seurat", "SeuratObject", "yaml", "jsonlite", "ggplot2", "patchwork")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing Phase 1B R packages: ", paste(missing, collapse = ", "))
  }
}

load_phase1b_object <- function(path) {
  object <- readRDS(path)
  if (!inherits(object, "Seurat")) {
    stop("Phase 1B input is not a Seurat object: ", paste(class(object), collapse = ", "))
  }
  serialized_version <- tryCatch(as.character(attributes(object)$version), error = function(e) NA_character_)
  update_required <- !"images" %in% names(attributes(object)) ||
    (!is.na(serialized_version) && package_version(serialized_version) < package_version("5.0.0"))
  messages <- character()
  warnings <- character()
  if (update_required) {
    object <- withCallingHandlers(
      SeuratObject::UpdateSeuratObject(object),
      message = function(condition) {
        messages <<- c(messages, conditionMessage(condition))
        invokeRestart("muffleMessage")
      },
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    )
  }
  list(
    object = object,
    compatibility = list(
      serialized_seurat_version = serialized_version,
      update_required = update_required,
      update_applied_in_memory = update_required,
      current_object_version = tryCatch(as.character(object@version), error = function(e) NA_character_),
      source_rds_modified = FALSE,
      messages = unname(unique(messages)),
      warnings = unname(unique(warnings))
    )
  )
}

select_timepoint_cells <- function(metadata, cell_ids, hours, downsample, seed) {
  if (!"Hour" %in% colnames(metadata)) stop("Missing required metadata column: Hour")
  if (nrow(metadata) != length(cell_ids)) stop("Metadata and cell IDs are not aligned")
  set.seed(seed)
  selected <- setNames(vector("list", length(hours)), paste0(hours, "h"))
  input_counts <- setNames(integer(length(hours)), paste0(hours, "h"))
  for (index in seq_along(hours)) {
    hour <- hours[[index]]
    candidates <- cell_ids[as.character(metadata$Hour) == as.character(hour)]
    input_counts[[index]] <- length(candidates)
    if (length(candidates) == 0L) stop("No cells found for Hour == ", hour)
    selected[[index]] <- sample(candidates, min(length(candidates), downsample), replace = FALSE)
  }
  list(selected = selected, input_counts = as.list(input_counts))
}

write_sampling_files <- function(selected_cells, specs, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  files <- list()
  for (hour_name in names(selected_cells)) {
    canonical <- file.path(directory, paste0("canonical_", hour_name, "_cells.txt"))
    writeLines(selected_cells[[hour_name]], canonical)
    files[[paste0("canonical_", hour_name)]] <- canonical
  }
  for (spec in specs) {
    prefix <- spec$id
    for (hour in spec$hours) {
      hour_name <- paste0(hour, "h")
      destination <- file.path(directory, paste0(prefix, "_", hour_name, "_cells.txt"))
      writeLines(selected_cells[[hour_name]], destination)
      files[[paste0(prefix, "_", hour_name)]] <- destination
    }
  }
  files
}

extract_phase1b_matrices <- function(object, selected_cells, hours, assay, pseudocount) {
  needed_names <- paste0(hours, "h")
  union_cells <- unique(unlist(selected_cells[needed_names], use.names = FALSE))
  subset_object <- object[, union_cells]
  matrices <- setNames(vector("list", length(hours)), needed_names)
  for (index in seq_along(hours)) {
    matrices[[index]] <- SCPA::seurat_extract(
      subset_object,
      assay = assay,
      meta1 = "Hour",
      value_meta1 = hours[[index]],
      pseudocount = pseudocount
    )
  }
  matrices
}

read_pathways_official <- function(path) {
  getFromNamespace("get_paths", "SCPA")(path)
}

summarize_numeric <- function(values) {
  stats <- stats::quantile(values, probs = c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
  list(
    min = unname(stats[[1]]),
    q1 = unname(stats[[2]]),
    median = unname(stats[[3]]),
    mean = unname(mean(values)),
    q3 = unname(stats[[4]]),
    max = unname(stats[[5]])
  )
}

inspect_pathways <- function(pathways, dataset_genes, min_genes, max_genes) {
  pathway_names <- vapply(pathways, function(x) as.character(unique(x$Pathway)[[1]]), character(1))
  pathway_genes <- lapply(pathways, function(x) unique(as.character(x$Genes[!is.na(x$Genes) & x$Genes != ""])))
  input_sizes <- vapply(pathway_genes, length, integer(1))
  matched_sizes <- vapply(pathway_genes, function(genes) sum(genes %in% dataset_genes), integer(1))
  analyzed <- matched_sizes >= min_genes & matched_sizes <= max_genes
  all_pathway_genes <- unique(unlist(pathway_genes, use.names = FALSE))
  matched_genes <- intersect(all_pathway_genes, dataset_genes)
  list(
    collection = "SCPA combined metabolic pathways (Hallmark, KEGG, Reactome)",
    input_pathway_count = length(pathways),
    analyzed_pathway_count = sum(analyzed),
    excluded_pathway_count = sum(!analyzed),
    input_size_distribution = summarize_numeric(input_sizes),
    matched_size_distribution = summarize_numeric(matched_sizes),
    unique_gene_set_genes = length(all_pathway_genes),
    matched_unique_genes = length(matched_genes),
    unmatched_unique_genes = length(setdiff(all_pathway_genes, dataset_genes)),
    excluded_pathways = unname(pathway_names[!analyzed]),
    analyzed_pathways = unname(pathway_names[analyzed])
  )
}

validate_scpa_result <- function(result, pairwise) {
  expected <- if (pairwise) {
    c("Pathway", "Pval", "adjPval", "qval", "FC")
  } else {
    c("Pathway", "Pval", "adjPval", "qval")
  }
  failures <- character()
  if (!identical(colnames(result), expected)) failures <- c(failures, "official_output_schema")
  if (nrow(result) == 0L) failures <- c(failures, "nonempty_result")
  for (column in intersect(c("Pval", "adjPval", "qval"), colnames(result))) {
    if (any(!is.finite(result[[column]]))) failures <- c(failures, paste0("finite_", column))
  }
  if (pairwise && "FC" %in% colnames(result) && any(!is.finite(result$FC))) {
    failures <- c(failures, "finite_FC")
  }
  unique(failures)
}

run_one_scpa_analysis <- function(matrices, pathway_file, scpa_config, seed, pairwise) {
  set.seed(seed)
  started <- Sys.time()
  runtime_warnings <- character()
  result <- withCallingHandlers(
    SCPA::compare_pathways(
      samples = matrices,
      pathways = pathway_file,
      downsample = scpa_config$downsample,
      min_genes = scpa_config$min_genes,
      max_genes = scpa_config$max_genes,
      parallel = isTRUE(scpa_config$parallel),
      cores = scpa_config$cores %||% NULL
    ),
    warning = function(condition) {
      runtime_warnings <<- c(runtime_warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(
    result = result,
    failures = validate_scpa_result(result, pairwise),
    warnings = unname(unique(runtime_warnings)),
    elapsed_seconds = unname(as.numeric(difftime(Sys.time(), started, units = "secs")))
  )
}

write_csv_atomic <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(data, temporary, row.names = FALSE, quote = TRUE)
  if (!file.rename(temporary, path)) stop("Could not atomically write CSV: ", path)
}

rank_scpa_result <- function(result) {
  if (!"qval" %in% colnames(result)) stop("SCPA result is missing qval")
  ordered <- result[order(result$qval, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  ordered$Rank <- seq_len(nrow(ordered))
  ordered[, c("Rank", setdiff(colnames(ordered), "Rank")), drop = FALSE]
}

markdown_result_table <- function(result, n = 10L) {
  if (is.null(result) || nrow(result) == 0L) return("No result available.")
  ordered <- if ("Rank" %in% colnames(result)) {
    result[order(result$qval, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  } else {
    rank_scpa_result(result)
  }
  shown <- utils::head(ordered, n)
  columns <- intersect(c("Rank", "Pathway", "qval", "FC"), colnames(shown))
  header <- paste0("| ", paste(columns, collapse = " | "), " |")
  divider <- paste0("| ", paste(rep("---", length(columns)), collapse = " | "), " |")
  rows <- apply(shown[, columns, drop = FALSE], 1, function(row) {
    paste0("| ", paste(row, collapse = " | "), " |")
  })
  paste(c(header, divider, rows), collapse = "\n")
}

write_phase1b_summary <- function(results, qc, path) {
  hour_names <- c("global_0_12_24", "pairwise_0_vs_12", "pairwise_12_vs_24", "pairwise_0_vs_24")
  reference_name <- "reference_resting0_vs_activated24"
  lines <- c(
    "# Phase 1B Vanilla SCPA reproduction summary",
    "",
    paste0("Gate status: `", qc$gate$status, "`"),
    "",
    "## SCPA qval convention",
    "",
    "SCPA qval is interpreted using the package convention: larger qval = stronger multivariate pathway difference. All rankings below are descending; rank 1 has the largest qval, and qval=0 belongs at the weakest end.",
    "",
    "## Official workflow followed",
    "",
    "The run uses SCPA::seurat_extract on the RNA assay's existing log1p-normalized data layer (pseudocount 0.001), then SCPA::compare_pathways with the official combined metabolic pathway collection.",
    "",
    "## Compatibility adaptation",
    "",
    paste0("Serialized Seurat ", qc$compatibility$serialized_seurat_version,
           " was updated in memory to ", qc$compatibility$current_object_version,
           "; the source RDS was not modified."),
    "",
    "## Expression and pathways",
    "",
    paste0("- Assay/layer: `", qc$expression$assay, "/", qc$expression$layer_or_slot, "`"),
    paste0("- Matrix dimensions: ", paste(names(qc$expression$matrix_dimensions),
      vapply(qc$expression$matrix_dimensions, function(x) paste0(x$genes, "x", x$cells), character(1)),
      sep = "=", collapse = ", ")),
    paste0("- Collection: ", qc$pathways$collection),
    paste0("- Input/analyzed/excluded pathways: ", qc$pathways$input_pathway_count, "/",
      qc$pathways$analyzed_pathway_count, "/", qc$pathways$excluded_pathway_count),
    "",
    "## A. Hour-only analyses",
    "",
    "These four retained analyses compare all naïve CD4 cells grouped only by observed Hour. They are valid time-point grouped analyses, not an exact paper Figure 4 or two-population tutorial reproduction."
  )
  for (analysis_name in hour_names[hour_names %in% names(results)]) {
    lines <- c(lines, "", paste0("### ", analysis_name, " — top pathways by descending qval"), "", markdown_result_table(results[[analysis_name]]))
  }

  reference <- results[[reference_name]]
  reference_qc <- qc$analyses[[reference_name]]
  lines <- c(lines, "", "## B. Official two-population reference reproduction", "")
  if (is.null(reference)) {
    lines <- c(
      lines,
      "Reference result not yet available. The required comparison is `Cell_Type=Resting AND Hour=0` versus `Cell_Type=Activated AND Hour=24`.",
      if (!is.null(reference_qc$error)) paste0("- Recorded error: ", reference_qc$error) else character()
    )
  } else {
    ranked_reference <- rank_scpa_result(reference)
    arachidonic_reference <- ranked_reference[grepl("ARACHIDONIC_ACID_METABOLISM", ranked_reference$Pathway), , drop = FALSE]
    target_summary <- reference_qc$target_summary
    lines <- c(
      lines,
      paste0("- Population 1: `", reference_qc$population_1_definition %||% "Cell_Type=Resting AND Hour=0", "`"),
      paste0("- Population 2: `", reference_qc$population_2_definition %||% "Cell_Type=Activated AND Hour=24", "`"),
      paste0("- Full/sampled cells: ", reference_qc$full_cell_count_population_1, "/", reference_qc$actual_sampled_cells_population_1,
             " vs ", reference_qc$full_cell_count_population_2, "/", reference_qc$actual_sampled_cells_population_2),
      paste0("- Seed/assay/layer: `", reference_qc$seed, "`, `", reference_qc$assay, "/", reference_qc$layer, "`"),
      paste0("- Analyzed pathways and finite qval/FC: ", reference_qc$pathway_count, "/",
             reference_qc$finite_qval_count, "/", reference_qc$finite_fc_count),
      "",
      "### Reference top pathways by descending qval",
      "",
      markdown_result_table(reference),
      "",
      "### Arachidonic-acid reference rows",
      "",
      markdown_result_table(arachidonic_reference, n = 20L),
      "",
      paste0("- Qualitative agreement: `", reference_qc$qualitative_agreement %||% "REQUIRES_REVIEW", "`"),
      if (!is.null(target_summary)) {
        paste0(
          "- High-qval pathways with modest |FC| (<=5): ",
          target_summary$high_qval_modest_fc_count,
          "; high qval does not require large FC: `",
          target_summary$high_qval_does_not_require_large_fc, "`."
        )
      } else character(),
      "- The tutorial comparison is qualitative; numerical identity is not required. High qval with modest FC supports the SCPA multivariate interpretation.",
      "- Any numerical difference may reflect the frozen package version, seed/downsampling, and the in-memory Seurat compatibility update; parameters were not retuned."
    )
  }

  if (!is.null(qc$visualizations)) {
    lines <- c(
      lines,
      "",
      "## Paper/tutorial comparison figures",
      "",
      paste0("- Visualization status: `", qc$visualizations$status, "`"),
      if (length(qc$visualizations$files) > 0L) {
        paste0("- Files: ", paste(unlist(qc$visualizations$files), collapse = ", "))
      } else {
        "- Files: none"
      },
      "- The 01–03/composite figures describe Hour-only analyses; the separate 05 figure describes the official Resting 0 h versus Activated 24 h reference. All comparisons are qualitative only."
    )
  }

  pair_024 <- results$pairwise_0_vs_24
  arachidonic <- if (!is.null(pair_024)) {
    ranked_pair_024 <- rank_scpa_result(pair_024)
    ranked_pair_024[grepl("ARACHI", ranked_pair_024$Pathway, ignore.case = TRUE), , drop = FALSE]
  } else NULL
  global <- results$global_0_12_24
  glycolysis <- if (!is.null(global)) {
    ranked_global <- rank_scpa_result(global)
    ranked_global[grepl("GLYCOLYSIS", ranked_global$Pathway, ignore.case = TRUE), , drop = FALSE]
  } else NULL
  lines <- c(
    lines,
    "",
    "## Hour-only qualitative comparison targets",
    "",
    "The official tutorial highlights arachidonic-acid metabolism as a large multivariate change that need not have large mean enrichment. The multisample tutorial highlights glycolysis-related pathways. These are review targets, not tuning targets.",
    "",
    "### Arachidonic-related rows in 0 vs 24",
    "",
    if (is.null(arachidonic) || nrow(arachidonic) == 0L) {
      "No arachidonic-related result row was found; this is an explicit qualitative mismatch/review item."
    } else {
      paste0(nrow(arachidonic), " arachidonic-related row(s) were found; significance and FC direction require GPT review.")
    },
    "",
    markdown_result_table(arachidonic, n = 20L),
    "",
    "### Glycolysis-related rows in global 0/12/24",
    "",
    if (is.null(glycolysis) || nrow(glycolysis) == 0L) {
      "No glycolysis-related result row was found; this is an explicit qualitative mismatch/review item."
    } else {
      paste0(nrow(glycolysis), " glycolysis-related row(s) were found; qval strength requires GPT review.")
    },
    "",
    markdown_result_table(glycolysis, n = 20L),
    "",
    "## Agreement, uncertainty, and review",
    "",
    "- This run preserves official SCPA expression extraction, pathway collection, statistics, and thresholds.",
    "- It is not a numerical replication of paper Figure 4: this protocol groups all cells by real Hour, whereas the paper/tutorial also use Cell_Type-specific or pseudotime-milestone populations.",
    "- Global qval alone does not identify timing or direction; review it with all three pairwise outputs.",
    "- Phase 1B PASS additionally requires the separate Resting 0 h versus Activated 24 h reference run and its qualitative review.",
    "- No parameter was tuned after viewing results.",
    paste0("- Final gate basis: `", qc$gate$pass_basis %||% "reference run pending; numerical identity is not a PASS requirement", "`."),
    paste0("- Runtime versions: SCPA ", qc$reproducibility$scpa_version,
           "; Seurat ", qc$reproducibility$Seurat_version,
           "; SeuratObject ", qc$reproducibility$SeuratObject_version, ".")
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".phase1b_summary_", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  writeLines(lines, temporary)
  if (!file.rename(temporary, path)) stop("Could not atomically write summary: ", path)
}
