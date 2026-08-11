`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

sha256_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  output <- tryCatch(
    system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE),
    error = function(e) structure(conditionMessage(e), status = 127L)
  )
  status <- attr(output, "status") %||% 0L
  if (status != 0L || length(output) == 0L) {
    return(NA_character_)
  }
  strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
}

gzip_integrity <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  output <- tryCatch(
    system2(
      "gzip",
      c("-t", "--", shQuote(path)),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) structure(conditionMessage(e), status = 127L)
  )
  identical(as.integer(attr(output, "status") %||% 0L), 0L)
}

canonicalize_time <- function(values) {
  text <- tolower(trimws(as.character(values)))
  text <- gsub("[[:space:]_-]+", " ", text)
  canonical <- rep(NA_character_, length(text))
  for (hour in c(0L, 12L, 24L)) {
    pattern <- paste0("^", hour, " ?(h|hr|hrs|hour|hours)?$")
    canonical[grepl(pattern, text)] <- paste0(hour, "h")
  }
  canonical
}

detect_time_column <- function(metadata, preferred_columns = character()) {
  required <- c("0h", "12h", "24h")
  if (ncol(metadata) == 0L) {
    return(list(column = NULL, canonical = character(), candidates = character()))
  }

  canonical_values <- lapply(metadata, canonicalize_time)
  coverage <- vapply(
    canonical_values,
    function(x) sum(required %in% unique(x[!is.na(x)])),
    integer(1)
  )
  candidates <- names(coverage)[coverage == length(required)]
  if (length(candidates) == 0L) {
    best <- names(coverage)[coverage == max(coverage)]
    selected <- if (max(coverage) > 0L) best[[1]] else NULL
  } else {
    preference <- match(tolower(candidates), tolower(preferred_columns))
    preference[is.na(preference)] <- length(preferred_columns) + seq_len(sum(is.na(preference)))
    selected <- candidates[order(preference)][[1]]
  }

  list(
    column = selected,
    canonical = if (is.null(selected)) character() else canonical_values[[selected]],
    candidates = candidates,
    coverage = as.list(coverage)
  )
}

infer_feature_identifier_type <- function(features) {
  usable <- features[!is.na(features) & nzchar(features)]
  if (length(usable) == 0L) {
    return("unavailable")
  }
  ensembl_fraction <- mean(grepl("^ENSG[0-9]+(\\.[0-9]+)?$", usable))
  symbol_fraction <- mean(grepl("^[A-Za-z][A-Za-z0-9._-]*$", usable))
  if (ensembl_fraction >= 0.9) {
    "ensembl_gene_id_like"
  } else if (symbol_fraction >= 0.9) {
    "gene_symbol_like"
  } else {
    "mixed_or_unknown"
  }
}

assay_layers <- function(object, assays) {
  result <- setNames(vector("list", length(assays)), assays)
  for (assay_name in assays) {
    assay_object <- object[[assay_name]]
    layers <- tryCatch(
      SeuratObject::Layers(assay_object),
      error = function(e) character()
    )
    if (length(layers) == 0L) {
      slots <- intersect(c("counts", "data", "scale.data"), methods::slotNames(assay_object))
      layers <- slots[vapply(slots, function(slot_name) {
        value <- tryCatch(methods::slot(assay_object, slot_name), error = function(e) NULL)
        !is.null(value) && length(value) > 0L
      }, logical(1))]
    }
    result[[assay_name]] <- unname(as.character(layers))
  }
  result
}

identity_evidence <- function(object, metadata, archive_path) {
  evidence <- character()
  if (identical(basename(archive_path), "GSE212270_integrated_naive_cd4.rds.gz")) {
    evidence <- c(evidence, "official GEO filename contains integrated_naive_cd4")
  }

  project_name <- tryCatch(as.character(object@project.name), error = function(e) "")
  if (length(project_name) > 0L && grepl("na.?ve.*cd4|cd4.*na.?ve", project_name, ignore.case = TRUE)) {
    evidence <- c(evidence, paste0("Seurat project.name: ", project_name[[1]]))
  }

  inspect_columns <- names(metadata)[grepl(
    "cell|type|population|sample|orig|ident",
    names(metadata),
    ignore.case = TRUE
  )]
  for (column in inspect_columns) {
    values <- unique(as.character(metadata[[column]]))
    values <- values[!is.na(values)]
    matches <- values[grepl("na.?ve.*cd4|cd4.*na.?ve", values, ignore.case = TRUE)]
    if (length(matches) > 0L) {
      evidence <- c(
        evidence,
        paste0("metadata ", column, ": ", paste(head(matches, 5L), collapse = ", "))
      )
    }
  }
  unique(evidence)
}

write_qc_json <- function(report, output_path) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".phase1_dataset_qc_", tmpdir = dirname(output_path))
  on.exit(unlink(temporary), add = TRUE)
  jsonlite::write_json(
    report,
    temporary,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    na = "null"
  )
  if (!file.rename(temporary, output_path)) {
    stop("Could not atomically write QC JSON: ", output_path)
  }
}

validate_phase1_dataset <- function(
  archive_path,
  rds_path,
  download_metadata_path,
  output_path,
  preferred_time_columns = c("Hour", "hour", "time", "Time", "timepoint", "time_point")
) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to write the Phase 1A QC report")
  }

  failures <- character()
  warnings <- character()
  add_failure <- function(check) failures <<- unique(c(failures, check))
  add_warning <- function(check) warnings <<- unique(c(warnings, check))

  archive_exists <- file.exists(archive_path)
  archive_size <- if (archive_exists) unname(file.info(archive_path)$size) else NA_real_
  archive_sha256 <- if (archive_exists) sha256_file(archive_path) else NA_character_
  archive_gzip_ok <- archive_exists && gzip_integrity(archive_path)

  if (!archive_exists) add_failure("archive_file_exists")
  if (archive_exists && !archive_gzip_ok) add_failure("gzip_integrity")
  if (archive_exists && is.na(archive_sha256)) add_failure("sha256_available")
  if (!identical(basename(archive_path), "GSE212270_integrated_naive_cd4.rds.gz")) {
    add_failure("official_geo_filename")
  }

  download_metadata <- NULL
  if (!file.exists(download_metadata_path)) {
    add_failure("download_metadata_exists")
  } else {
    download_metadata <- tryCatch(
      jsonlite::fromJSON(download_metadata_path, simplifyVector = FALSE),
      error = function(e) {
        add_failure("download_metadata_valid_json")
        NULL
      }
    )
  }
  if (!is.null(download_metadata)) {
    required_download_keys <- c(
      "geo_accession", "filename", "download_source", "file_size_bytes",
      "sha256", "gzip_integrity"
    )
    missing_download_keys <- setdiff(required_download_keys, names(download_metadata))
    if (length(missing_download_keys) > 0L) {
      add_failure("download_metadata_schema")
    } else {
      if (!identical(download_metadata$geo_accession, "GSE212270")) {
        add_failure("download_metadata_geo_accession")
      }
      if (!identical(download_metadata$filename, basename(archive_path))) {
        add_failure("download_metadata_filename")
      }
      if (!grepl("^https://ftp\\.ncbi\\.nlm\\.nih\\.gov/geo/", download_metadata$download_source)) {
        add_failure("download_metadata_official_geo_source")
      }
      if (archive_exists && !identical(as.numeric(download_metadata$file_size_bytes), as.numeric(archive_size))) {
        add_failure("download_metadata_file_size_match")
      }
      if (archive_exists && !identical(download_metadata$sha256, archive_sha256)) {
        add_failure("download_metadata_sha256_match")
      }
      if (!isTRUE(download_metadata$gzip_integrity)) {
        add_failure("download_metadata_gzip_integrity")
      }
    }
  }

  object <- NULL
  read_error <- NULL
  if (!file.exists(rds_path)) {
    add_failure("extracted_rds_exists")
  } else {
    object <- tryCatch(
      readRDS(rds_path),
      error = function(e) {
        read_error <<- conditionMessage(e)
        NULL
      }
    )
    if (is.null(object)) add_failure("rds_read_success")
  }

  object_class <- if (is.null(object)) character() else class(object)
  is_seurat <- !is.null(object) && inherits(object, "Seurat")
  if (!is.null(object) && !is_seurat) add_failure("seurat_object_class")

  serialized_object_version <- if (is_seurat) {
    tryCatch(as.character(attributes(object)$version), error = function(e) NA_character_)
  } else {
    NA_character_
  }
  compatibility_update_required <- is_seurat && (
    !"images" %in% names(attributes(object)) ||
      (!is.na(serialized_object_version) &&
        package_version(serialized_object_version) < package_version("5.0.0"))
  )
  compatibility_update_applied <- FALSE
  compatibility_update_success <- !compatibility_update_required
  compatibility_messages <- character()
  compatibility_warnings <- character()
  compatibility_error <- NULL

  if (compatibility_update_required) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      compatibility_update_success <- FALSE
      compatibility_error <- "SeuratObject is unavailable"
      add_failure("legacy_seurat_compatibility_update")
      is_seurat <- FALSE
    } else {
      updated_object <- tryCatch(
        withCallingHandlers(
          SeuratObject::UpdateSeuratObject(object),
          message = function(message_condition) {
            compatibility_messages <<- c(
              compatibility_messages,
              conditionMessage(message_condition)
            )
            invokeRestart("muffleMessage")
          },
          warning = function(warning_condition) {
            compatibility_warnings <<- c(
              compatibility_warnings,
              conditionMessage(warning_condition)
            )
            invokeRestart("muffleWarning")
          }
        ),
        error = function(error_condition) {
          compatibility_error <<- conditionMessage(error_condition)
          NULL
        }
      )
      if (is.null(updated_object)) {
        compatibility_update_success <- FALSE
        add_failure("legacy_seurat_compatibility_update")
        is_seurat <- FALSE
      } else {
        object <- updated_object
        compatibility_update_applied <- TRUE
        compatibility_update_success <- TRUE
        add_warning(paste0(
          "legacy_seurat_object_updated_in_memory:",
          serialized_object_version,
          "->",
          as.character(object@version)
        ))
      }
    }
  }

  cells <- NA_integer_
  features_count <- NA_integer_
  assays <- character()
  active_assay <- NA_character_
  metadata <- data.frame()
  metadata_columns <- character()
  object_version <- NA_character_
  installed_seurat <- if (requireNamespace("Seurat", quietly = TRUE)) {
    as.character(utils::packageVersion("Seurat"))
  } else {
    NA_character_
  }
  installed_seurat_object <- if (requireNamespace("SeuratObject", quietly = TRUE)) {
    as.character(utils::packageVersion("SeuratObject"))
  } else {
    NA_character_
  }

  time_column <- NULL
  detected_time_labels <- list()
  cells_per_timepoint <- list(`0h` = 0L, `12h` = 0L, `24h` = 0L)
  time_candidates <- character()
  layers <- list()
  counts_available <- FALSE
  normalized_available <- FALSE
  features <- character()
  duplicate_features <- NA_integer_
  missing_feature_names <- NA_integer_
  feature_identifier_type <- "unavailable"
  identity <- character()
  internal_identity <- character()

  if (is_seurat) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      add_failure("seuratobject_package_available")
    } else {
      object_version <- tryCatch(as.character(object@version), error = function(e) NA_character_)
      cells <- ncol(object)
      features_count <- nrow(object)
      if (cells <= 0L) add_failure("nonzero_cell_count")
      if (features_count <= 0L) add_failure("nonzero_feature_count")

      assays <- tryCatch(as.character(SeuratObject::Assays(object)), error = function(e) character())
      active_assay <- tryCatch(as.character(SeuratObject::DefaultAssay(object)), error = function(e) NA_character_)
      if (length(assays) == 0L) add_failure("assays_available")
      if (is.na(active_assay) || !active_assay %in% assays) add_failure("active_assay_available")

      metadata <- tryCatch(object[[]], error = function(e) data.frame())
      metadata_columns <- colnames(metadata) %||% character()
      if (nrow(metadata) != cells) add_failure("metadata_cell_alignment")

      time_detection <- detect_time_column(metadata, preferred_time_columns)
      time_column <- time_detection$column
      time_candidates <- time_detection$candidates
      if (is.null(time_column)) {
        add_failure("time_column_detected")
      } else {
        canonical <- time_detection$canonical
        raw_values <- as.character(metadata[[time_column]])
        for (label in names(cells_per_timepoint)) {
          cells_per_timepoint[[label]] <- as.integer(sum(canonical == label, na.rm = TRUE))
          detected_time_labels[[label]] <- unname(unique(raw_values[canonical == label & !is.na(canonical)]))
          if (cells_per_timepoint[[label]] <= 0L) {
            add_failure(paste0("nonzero_cells_", label))
          }
        }
        if (length(time_candidates) > 1L) {
          add_warning(paste0("multiple_time_columns_detected:", paste(time_candidates, collapse = ",")))
        }
      }

      identity <- identity_evidence(object, metadata, archive_path)
      if (length(identity) == 0L) add_failure("naive_cd4_identity_evidence")
      internal_identity <- identity[!grepl("^official GEO filename", identity)]
      if (length(internal_identity) == 0L) {
        add_warning("naive_cd4_internal_metadata_or_project_evidence_not_detected")
      }

      layers <- assay_layers(object, assays)
      all_layers <- unname(unlist(layers, use.names = FALSE))
      counts_available <- any(grepl("^counts($|\\.)", all_layers, ignore.case = TRUE))
      normalized_available <- any(grepl("(^|\\.)data($|\\.)|norm", all_layers, ignore.case = TRUE))
      if (!counts_available) add_warning("counts_layer_not_detected")
      if (!normalized_available) add_failure("normalized_expression_available")

      features <- rownames(object) %||% character()
      missing_feature_names <- as.integer(sum(is.na(features) | !nzchar(features)))
      duplicate_features <- as.integer(sum(duplicated(features)))
      feature_identifier_type <- infer_feature_identifier_type(features)
      if (length(features) == 0L || missing_feature_names > 0L) {
        add_failure("usable_feature_names")
      }
      if (duplicate_features > 0L) add_failure("unique_feature_names")
    }
  }

  status <- if (length(failures) == 0L) "PASS" else "FAIL"
  report <- list(
    phase = "Phase 1A - Dataset acquisition / validation gate",
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    dataset = list(
      geo_accession = "GSE212270",
      filename = basename(archive_path),
      file_size_bytes = archive_size,
      sha256 = archive_sha256,
      gzip_integrity = archive_gzip_ok,
      extracted_rds = basename(rds_path)
    ),
    object = list(
      rds_read_success = !is.null(object),
      read_error = read_error,
      class = unname(object_class),
      serialized_seurat_version = serialized_object_version,
      seurat_object_version = object_version,
      installed_seurat_version = installed_seurat,
      installed_seuratobject_version = installed_seurat_object,
      compatibility = list(
        update_required = compatibility_update_required,
        update_applied_in_memory = compatibility_update_applied,
        update_success = compatibility_update_success,
        source_rds_modified = FALSE,
        messages = unname(unique(compatibility_messages)),
        warnings = unname(unique(compatibility_warnings)),
        error = compatibility_error
      ),
      cells = cells,
      features = features_count,
      assays = unname(assays),
      active_assay = active_assay
    ),
    metadata = list(
      columns = unname(metadata_columns),
      time_column = time_column,
      candidate_time_columns = unname(time_candidates),
      detected_time_labels = detected_time_labels,
      cells_per_timepoint = cells_per_timepoint
    ),
    identity = list(
      expected = "integrated naive CD4 T cells",
      evidence = unname(identity),
      matched = length(identity) > 0L,
      internal_metadata_or_project_evidence_detected = length(internal_identity) > 0L
    ),
    expression = list(
      available_layers = layers,
      counts_available = counts_available,
      normalized_expression_available = normalized_available,
      feature_identifier_type = feature_identifier_type,
      duplicate_features = duplicate_features,
      missing_feature_names = missing_feature_names
    ),
    gate = list(
      status = status,
      failed_checks = unname(failures),
      warnings = unname(warnings)
    )
  )
  write_qc_json(report, output_path)
  list(report = report, output_path = output_path, exit_status = if (status == "PASS") 0L else 1L)
}
