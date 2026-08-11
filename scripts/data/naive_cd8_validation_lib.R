cell_type_summary <- function(metadata) {
  matches <- which(tolower(colnames(metadata)) == "cell_type")
  if (length(matches) == 0L) {
    return(list(column_present = FALSE, column = NULL, unique_values = character(), counts = list()))
  }
  column <- colnames(metadata)[matches[[1]]]
  values <- as.character(metadata[[column]])
  values[is.na(values)] <- "<NA>"
  counts <- table(values, useNA = "no")
  list(
    column_present = TRUE,
    column = column,
    unique_values = unname(names(counts)),
    counts = as.list(stats::setNames(as.integer(counts), names(counts)))
  )
}

internal_cd8_identity_evidence <- function(object, metadata) {
  evidence <- character()
  project_name <- tryCatch(as.character(object@project.name), error = function(e) "")
  if (length(project_name) > 0L && grepl("na.?ve.*cd8|cd8.*na.?ve", project_name, ignore.case = TRUE)) {
    evidence <- c(evidence, paste0("Seurat project.name: ", project_name[[1]]))
  }
  inspect_columns <- names(metadata)[grepl("cell|type|population|sample|orig|ident", names(metadata), ignore.case = TRUE)]
  for (column in inspect_columns) {
    values <- unique(as.character(metadata[[column]]))
    values <- values[!is.na(values)]
    matches <- values[grepl("na.?ve.*cd8|cd8.*na.?ve", values, ignore.case = TRUE)]
    if (length(matches) > 0L) {
      evidence <- c(evidence, paste0("metadata ", column, ": ", paste(utils::head(matches, 5L), collapse = ", ")))
    }
  }
  unname(unique(evidence))
}

load_seurat_read_only <- function(path) {
  read_error <- NULL
  object <- tryCatch(readRDS(path), error = function(e) { read_error <<- conditionMessage(e); NULL })
  if (is.null(object) || !inherits(object, "Seurat")) {
    return(list(object = NULL, read_error = read_error, original_class = if (is.null(object)) character() else class(object)))
  }
  original_class <- class(object)
  serialized_version <- tryCatch(as.character(attributes(object)$version), error = function(e) NA_character_)
  update_required <- !"images" %in% names(attributes(object)) ||
    (!is.na(serialized_version) && package_version(serialized_version) < package_version("5.0.0"))
  messages <- character()
  warnings <- character()
  update_error <- NULL
  if (update_required) {
    object <- tryCatch(
      withCallingHandlers(
        SeuratObject::UpdateSeuratObject(object),
        message = function(c) { messages <<- c(messages, conditionMessage(c)); invokeRestart("muffleMessage") },
        warning = function(c) { warnings <<- c(warnings, conditionMessage(c)); invokeRestart("muffleWarning") }
      ),
      error = function(e) { update_error <<- conditionMessage(e); NULL }
    )
  }
  list(
    object = object,
    read_error = read_error,
    original_class = original_class,
    serialized_version = serialized_version,
    update_required = update_required,
    update_applied_in_memory = update_required && !is.null(object),
    update_error = update_error,
    messages = unname(unique(messages)),
    warnings = unname(unique(warnings))
  )
}

inspect_reference_cd4 <- function(path, qc_path) {
  result <- list(status = "NOT_AVAILABLE", rds_path = path, qc_path = qc_path)
  if (!file.exists(path)) return(result)
  loaded <- load_seurat_read_only(path)
  if (is.null(loaded$object)) {
    result$status <- "ERROR"
    result$error <- loaded$read_error %||% loaded$update_error %||% "CD4 object is not a usable Seurat object"
    return(result)
  }
  object <- loaded$object
  metadata <- tryCatch(object[[]], error = function(e) data.frame())
  assays <- tryCatch(as.character(SeuratObject::Assays(object)), error = function(e) character())
  layers <- assay_layers(object, assays)
  rna_layers <- layers$RNA %||% character()
  time <- detect_time_column(metadata, c("Hour", "hour", "time", "Time", "timepoint", "time_point"))
  time_presence <- as.list(stats::setNames(
    vapply(c("0h", "12h", "24h"), function(label) any(time$canonical == label, na.rm = TRUE), logical(1)),
    c("0h", "12h", "24h")
  ))
  qc_gate <- NULL
  if (file.exists(qc_path)) {
    qc_gate <- tryCatch(jsonlite::fromJSON(qc_path, simplifyVector = FALSE)$gate$status, error = function(e) NULL)
  }
  list(
    status = "AVAILABLE",
    rds_path = path,
    qc_path = qc_path,
    qc_gate_status = qc_gate,
    features = if ("RNA" %in% assays) rownames(object[["RNA"]]) %||% character() else character(),
    feature_count = if ("RNA" %in% assays) nrow(object[["RNA"]]) else 0L,
    feature_identifier_type = infer_feature_identifier_type(if ("RNA" %in% assays) rownames(object[["RNA"]]) %||% character() else character()),
    assays = unname(assays),
    rna_assay_present = "RNA" %in% assays,
    normalized_expression_available = any(grepl("(^|\\.)data($|\\.)|norm", rna_layers, ignore.case = TRUE)),
    time_points_present = time_presence,
    cell_type = cell_type_summary(metadata)
  )
}

validate_naive_cd8_dataset <- function(
  archive_path,
  rds_path,
  download_metadata_path,
  output_path,
  cd4_rds_path,
  cd4_qc_path,
  preferred_time_columns = c("Hour", "hour", "time", "Time", "timepoint", "time_point")
) {
  if (!requireNamespace("jsonlite", quietly = TRUE) || !requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("jsonlite and SeuratObject are required for naïve CD8 validation")
  }
  expected_filename <- "GSE212270_integrated_naive_cd8.rds.gz"
  failures <- character()
  warnings <- character()
  add_failure <- function(x) failures <<- unique(c(failures, x))
  add_warning <- function(x) warnings <<- unique(c(warnings, x))

  archive_exists <- file.exists(archive_path)
  archive_size <- if (archive_exists) unname(file.info(archive_path)$size) else NA_real_
  archive_sha <- if (archive_exists) sha256_file(archive_path) else NA_character_
  archive_gzip <- archive_exists && gzip_integrity(archive_path)
  if (!archive_exists) add_failure("archive_file_exists")
  if (archive_exists && !archive_gzip) add_failure("gzip_integrity")
  if (archive_exists && is.na(archive_sha)) add_failure("sha256_available")
  filename_evidence <- identical(basename(archive_path), expected_filename)
  if (!filename_evidence) add_failure("official_geo_filename")

  download_metadata <- NULL
  if (!file.exists(download_metadata_path)) {
    add_failure("download_metadata_exists")
  } else {
    download_metadata <- tryCatch(
      jsonlite::fromJSON(download_metadata_path, simplifyVector = FALSE),
      error = function(e) { add_failure("download_metadata_valid_json"); NULL }
    )
  }
  official_source_evidence <- FALSE
  if (!is.null(download_metadata)) {
    required <- c("geo_accession", "filename", "download_source", "file_size_bytes", "sha256", "gzip_integrity", "recorded_at")
    if (length(setdiff(required, names(download_metadata))) > 0L) add_failure("download_metadata_schema")
    official_source_evidence <- identical(download_metadata$geo_accession, "GSE212270") &&
      identical(download_metadata$filename, expected_filename) &&
      grepl("^https://ftp\\.ncbi\\.nlm\\.nih\\.gov/geo/", download_metadata$download_source %||% "")
    if (!official_source_evidence) add_failure("official_geo_source_identity_evidence")
    if (archive_exists && !identical(as.numeric(download_metadata$file_size_bytes), as.numeric(archive_size))) add_failure("download_metadata_file_size_match")
    if (archive_exists && !identical(download_metadata$sha256, archive_sha)) add_failure("download_metadata_sha256_match")
    if (!isTRUE(download_metadata$gzip_integrity)) add_failure("download_metadata_gzip_integrity")
  }

  loaded <- if (file.exists(rds_path)) load_seurat_read_only(rds_path) else list(object = NULL, read_error = "RDS file does not exist")
  if (!file.exists(rds_path)) add_failure("extracted_rds_exists")
  if (is.null(loaded$object)) {
    add_failure(if (length(loaded$original_class %||% character()) > 0L) "seurat_object_class" else "rds_read_success")
    if (!is.null(loaded$update_error)) add_failure("legacy_seurat_compatibility_update")
  }

  cells <- features_count <- NA_integer_
  assays <- character(); active_assay <- NA_character_; layers <- list()
  metadata <- data.frame(); time_column <- NULL; time_candidates <- character()
  detected_time_labels <- list(); cells_per_timepoint <- list(`0h` = 0L, `12h` = 0L, `24h` = 0L)
  features <- character(); feature_type <- "unavailable"; missing_features <- duplicate_features <- NA_integer_
  counts_available <- normalized_available <- FALSE
  internal_identity <- character(); cell_types <- cell_type_summary(metadata)
  rds_read_success <- !is.null(loaded$object)
  current_object_version <- NA_character_

  if (!is.null(loaded$object)) {
    object <- loaded$object
    current_object_version <- tryCatch(as.character(object@version), error = function(e) NA_character_)
    if (isTRUE(loaded$update_applied_in_memory)) add_warning(paste0("legacy_seurat_object_updated_in_memory:", loaded$serialized_version))
    cells <- ncol(object); features_count <- nrow(object)
    if (cells <= 0L) add_failure("nonzero_cell_count")
    if (features_count <= 0L) add_failure("nonzero_feature_count")
    assays <- tryCatch(as.character(SeuratObject::Assays(object)), error = function(e) character())
    active_assay <- tryCatch(as.character(SeuratObject::DefaultAssay(object)), error = function(e) NA_character_)
    if (!"RNA" %in% assays) add_failure("rna_assay_available")
    if (is.na(active_assay) || !active_assay %in% assays) add_failure("active_assay_available")
    layers <- assay_layers(object, assays)
    rna_layers <- layers$RNA %||% character()
    counts_available <- any(grepl("^counts($|\\.)", rna_layers, ignore.case = TRUE))
    normalized_available <- any(grepl("(^|\\.)data($|\\.)|norm", rna_layers, ignore.case = TRUE))
    if (!counts_available) add_warning("rna_counts_layer_not_detected")
    if (!normalized_available) add_failure("rna_normalized_expression_available")

    metadata <- tryCatch(object[[]], error = function(e) data.frame())
    if (nrow(metadata) != cells) add_failure("metadata_cell_alignment")
    time <- detect_time_column(metadata, preferred_time_columns)
    time_column <- time$column; time_candidates <- time$candidates
    if (is.null(time_column)) {
      add_failure("time_column_detected")
    } else {
      raw_time <- as.character(metadata[[time_column]])
      for (label in names(cells_per_timepoint)) {
        cells_per_timepoint[[label]] <- as.integer(sum(time$canonical == label, na.rm = TRUE))
        detected_time_labels[[label]] <- unname(unique(raw_time[time$canonical == label & !is.na(time$canonical)]))
        if (cells_per_timepoint[[label]] <= 0L) add_failure(paste0("nonzero_cells_", label))
      }
    }
    cell_types <- cell_type_summary(metadata)
    if (!cell_types$column_present) add_warning("cell_type_column_not_detected")
    internal_identity <- internal_cd8_identity_evidence(object, metadata)
    if (length(internal_identity) == 0L) add_warning("naive_cd8_internal_metadata_or_project_evidence_not_detected")

    features <- if ("RNA" %in% assays) rownames(object[["RNA"]]) %||% character() else character()
    missing_features <- as.integer(sum(is.na(features) | !nzchar(features)))
    duplicate_features <- as.integer(sum(duplicated(features)))
    feature_type <- infer_feature_identifier_type(features)
    if (length(features) == 0L || missing_features > 0L) add_failure("usable_feature_names")
    if (duplicate_features > 0L) add_failure("unique_feature_names")
    rm(object)
    loaded$object <- NULL
    invisible(gc())
  }

  cd4 <- inspect_reference_cd4(cd4_rds_path, cd4_qc_path)
  mismatch_reasons <- character()
  comparison <- list(status = cd4$status)
  if (identical(cd4$status, "AVAILABLE")) {
    shared <- intersect(cd4$features, features)
    if (!cd4$rna_assay_present || !"RNA" %in% assays) mismatch_reasons <- c(mismatch_reasons, "rna_assay_missing")
    if (!cd4$normalized_expression_available || !normalized_available) mismatch_reasons <- c(mismatch_reasons, "normalized_expression_missing")
    if (!identical(cd4$feature_identifier_type, feature_type)) mismatch_reasons <- c(mismatch_reasons, "feature_identifier_type_mismatch")
    if (length(shared) == 0L) mismatch_reasons <- c(mismatch_reasons, "no_exact_shared_genes")
    comparison <- list(
      status = "COMPARED",
      cd4_qc_gate_status = cd4$qc_gate_status,
      cd4_feature_count = cd4$feature_count,
      cd8_feature_count = length(features),
      exact_shared_gene_count = length(shared),
      cd4_only_gene_count = length(setdiff(cd4$features, features)),
      cd8_only_gene_count = length(setdiff(features, cd4$features)),
      gene_order_identical = identical(cd4$features, features),
      rna_assay_present = list(cd4 = cd4$rna_assay_present, cd8 = "RNA" %in% assays),
      normalized_expression_available = list(cd4 = cd4$normalized_expression_available, cd8 = normalized_available),
      time_points_present = list(cd4 = cd4$time_points_present, cd8 = as.list(stats::setNames(unlist(cells_per_timepoint) > 0L, names(cells_per_timepoint)))),
      feature_naming = list(cd4 = cd4$feature_identifier_type, cd8 = feature_type, compatible = identical(cd4$feature_identifier_type, feature_type)),
      cell_type = list(cd4 = cd4$cell_type, cd8 = cell_types),
      obvious_preprocessing_or_representation_mismatch = length(mismatch_reasons) > 0L,
      mismatch_reasons = unname(mismatch_reasons),
      intersection_matrix_created = FALSE,
      source_objects_modified = FALSE
    )
  } else {
    add_warning(paste0("naive_cd4_compatibility_comparison_", tolower(cd4$status)))
  }

  comparison_unavailable <- !identical(cd4$status, "AVAILABLE")
  gate_status <- if (length(failures) > 0L) {
    "FAIL"
  } else if (length(mismatch_reasons) > 0L || comparison_unavailable) {
    "NEEDS_REVIEW"
  } else {
    "PASS"
  }
  evidence <- c(
    if (official_source_evidence) "official GEO GSE212270 download metadata" else character(),
    if (filename_evidence) paste0("official filename: ", expected_filename) else character(),
    internal_identity
  )
  report <- list(
    phase = "Future Phase 3 preparation - naïve CD8 acquisition / validation",
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    dataset = list(
      geo_accession = "GSE212270", filename = basename(archive_path), official_geo_record = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212270",
      official_download_source = if (is.null(download_metadata)) NULL else download_metadata$download_source,
      file_size_bytes = archive_size, sha256 = archive_sha, gzip_integrity = archive_gzip, extracted_rds = basename(rds_path)
    ),
    object = list(
      rds_read_success = rds_read_success, read_error = loaded$read_error %||% loaded$update_error,
      class = unname(loaded$original_class %||% character()), serialized_seurat_version = loaded$serialized_version %||% NA_character_,
      seurat_object_version = current_object_version,
      installed_seurat_version = if (requireNamespace("Seurat", quietly = TRUE)) as.character(packageVersion("Seurat")) else NA_character_,
      installed_seuratobject_version = as.character(packageVersion("SeuratObject")),
      compatibility = list(update_required = loaded$update_required %||% FALSE, update_applied_in_memory = loaded$update_applied_in_memory %||% FALSE,
        update_success = is.null(loaded$update_error), source_rds_modified = FALSE, messages = loaded$messages %||% character(), warnings = loaded$warnings %||% character(), error = loaded$update_error %||% NULL),
      cells = cells, features = features_count, assays = unname(assays), active_assay = active_assay
    ),
    metadata = list(columns = unname(colnames(metadata)), time_column = time_column, candidate_time_columns = unname(time_candidates),
      detected_time_labels = detected_time_labels, cells_per_timepoint = cells_per_timepoint, cell_type = cell_types),
    identity = list(expected = "integrated naïve CD8 T cells", official_geo_source_evidence = official_source_evidence,
      filename_evidence = filename_evidence, internal_metadata_or_project_evidence = internal_identity,
      internal_evidence_detected = length(internal_identity) > 0L, evidence = unname(evidence), matched = official_source_evidence && filename_evidence),
    expression = list(available_layers = layers, counts_available = counts_available, normalized_expression_available = normalized_available,
      rna_feature_count = length(features), feature_identifier_type = feature_type,
      duplicate_features = duplicate_features, missing_feature_names = missing_features),
    comparison_with_naive_cd4 = comparison,
    gate = list(status = gate_status, failed_checks = unname(failures), warnings = unname(warnings))
  )
  write_qc_json(report, output_path)
  list(report = report, output_path = output_path, exit_status = if (gate_status == "PASS") 0L else if (gate_status == "NEEDS_REVIEW") 2L else 1L)
}
