load_seurat_for_genept <- function(path) {
  object <- readRDS(path)
  if (!inherits(object, "Seurat")) {
    stop("GenePT input is not a Seurat object: ", paste(class(object), collapse = ", "))
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

write_lines_atomic_r <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  writeLines(lines, temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) stop("Could not atomically write: ", path)
}

write_csv_atomic_r <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(data, temporary, row.names = FALSE, quote = TRUE)
  if (!file.rename(temporary, path)) stop("Could not atomically write: ", path)
}

write_matrix_market_atomic <- function(matrix, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  Matrix::writeMM(matrix, temporary)
  if (!file.rename(temporary, path)) stop("Could not atomically write matrix: ", path)
}

export_seurat_counts_for_genept <- function(object, output_dir, dataset, source_file, compatibility = list()) {
  if (!inherits(object, "Seurat")) stop("object must inherit from Seurat")
  if (!"RNA" %in% names(object@assays)) stop("RNA assay is missing")
  counts <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  if (!inherits(counts, "sparseMatrix")) stop("RNA counts must remain a sparse Matrix")
  if (any(!is.finite(counts@x)) || any(counts@x < 0)) stop("RNA counts contain invalid values")
  genes <- rownames(counts)
  cell_ids <- colnames(counts)
  if (anyNA(genes) || any(genes == "") || anyDuplicated(genes)) stop("Gene IDs must be unique and non-empty")
  if (anyNA(cell_ids) || any(cell_ids == "") || anyDuplicated(cell_ids)) stop("Cell IDs must be unique and non-empty")

  metadata <- object[[]]
  if (!identical(rownames(metadata), cell_ids)) stop("Metadata and counts cell order differ")
  required_metadata <- c("Hour", "Cell_Type")
  missing_metadata <- setdiff(required_metadata, colnames(metadata))
  if (length(missing_metadata) > 0L) {
    stop("Missing required metadata column(s): ", paste(missing_metadata, collapse = ", "))
  }
  metadata_output <- data.frame(
    cell_id = cell_ids,
    Hour = metadata$Hour,
    Cell_Type = metadata$Cell_Type,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- list(
    counts = file.path(output_dir, paste0(dataset, "_rna_counts_genes_by_cells.mtx")),
    genes = file.path(output_dir, paste0(dataset, "_gene_ids.txt")),
    cell_ids = file.path(output_dir, paste0(dataset, "_cell_ids.txt")),
    metadata = file.path(output_dir, paste0(dataset, "_metadata.csv")),
    manifest = file.path(output_dir, paste0(dataset, "_export_manifest.json"))
  )
  write_matrix_market_atomic(counts, paths$counts)
  write_lines_atomic_r(genes, paths$genes)
  write_lines_atomic_r(cell_ids, paths$cell_ids)
  write_csv_atomic_r(metadata_output, paths$metadata)

  manifest <- list(
    dataset = dataset,
    source_file = normalizePath(source_file, mustWork = FALSE),
    source_object_modified = FALSE,
    assay = "RNA",
    layer = "counts",
    matrix_orientation = "genes_by_cells",
    genes = nrow(counts),
    cells = ncol(counts),
    nonzero_values = length(counts@x),
    sparse_class = class(counts)[[1]],
    files = lapply(paths[names(paths) != "manifest"], normalizePath),
    sha256 = lapply(paths[names(paths) != "manifest"], sha256_file),
    compatibility = compatibility
  )
  write_qc_json(manifest, paths$manifest)
  list(manifest = manifest, paths = paths)
}
