args <- commandArgs(trailingOnly = TRUE)
root_arg <- if (length(args) >= 1L) args[[1]] else "."
project_root <- normalizePath(root_arg)
source(file.path(project_root, "scripts", "data", "phase1_validation_lib.R"))

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

if (!requireNamespace("SeuratObject", quietly = TRUE) ||
    !requireNamespace("jsonlite", quietly = TRUE)) {
  stop("SeuratObject and jsonlite are required for Phase 1A validation tests")
}

make_mock_object <- function(hours) {
  counts <- matrix(
    c(
      1, 0, 2, 1, 3, 0,
      0, 1, 1, 2, 0, 3,
      2, 2, 1, 0, 1, 1,
      0, 1, 0, 2, 2, 1,
      3, 0, 1, 1, 0, 2,
      1, 1, 1, 1, 1, 1
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(
      c("IL7R", "KLF2", "PIK3IP1", "CD4", "MAL", "CCR7"),
      paste0("cell", seq_len(6))
    )
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts, project = "naive_cd4_mock")
  object <- SeuratObject::SetAssayData(
    object,
    assay = "RNA",
    layer = "data",
    new.data = log1p(counts)
  )
  object$sample_time <- hours
  object$population <- rep("naive CD4 T cells", ncol(object))
  object
}

write_fixture <- function(directory, object = NULL, invalid_rds = FALSE) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  archive <- file.path(directory, "GSE212270_integrated_naive_cd4.rds.gz")
  rds <- file.path(directory, "GSE212270_integrated_naive_cd4.rds")
  metadata <- file.path(directory, "phase1_download_metadata.json")
  qc <- file.path(directory, "phase1_dataset_qc.json")

  if (invalid_rds) {
    writeLines("not an RDS", rds)
  } else {
    saveRDS(object, rds)
  }
  input_connection <- file(rds, open = "rb")
  on.exit(try(close(input_connection), silent = TRUE), add = TRUE)
  rds_bytes <- readBin(input_connection, what = "raw", n = file.info(rds)$size)
  close(input_connection)
  output_connection <- gzfile(archive, open = "wb")
  on.exit(try(close(output_connection), silent = TRUE), add = TRUE)
  writeBin(rds_bytes, output_connection)
  close(output_connection)

  download_record <- list(
    geo_accession = "GSE212270",
    filename = basename(archive),
    download_source = paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/",
      "GSE212270/suppl/GSE212270_integrated_naive_cd4.rds.gz"
    ),
    file_size_bytes = unname(file.info(archive)$size),
    sha256 = sha256_file(archive),
    gzip_integrity = TRUE
  )
  jsonlite::write_json(download_record, metadata, auto_unbox = TRUE, pretty = TRUE)
  list(archive = archive, rds = rds, metadata = metadata, qc = qc)
}

run_fixture <- function(paths) {
  validate_phase1_dataset(
    archive_path = paths$archive,
    rds_path = paths$rds,
    download_metadata_path = paths$metadata,
    output_path = paths$qc
  )
}

test_root <- tempfile("phase1_validation_tests_")
dir.create(test_root)
on.exit(unlink(test_root, recursive = TRUE), add = TRUE)

missing_paths <- list(
  archive = file.path(test_root, "missing", "GSE212270_integrated_naive_cd4.rds.gz"),
  rds = file.path(test_root, "missing", "GSE212270_integrated_naive_cd4.rds"),
  metadata = file.path(test_root, "missing", "phase1_download_metadata.json"),
  qc = file.path(test_root, "missing", "phase1_dataset_qc.json")
)
missing_result <- run_fixture(missing_paths)
assert_true(missing_result$exit_status != 0L, "Missing files must fail")
assert_true("archive_file_exists" %in% missing_result$report$gate$failed_checks, "Missing archive failure not recorded")
assert_true(file.exists(missing_paths$qc), "Missing-file failure must still write QC JSON")

invalid_paths <- write_fixture(file.path(test_root, "invalid"), invalid_rds = TRUE)
invalid_result <- run_fixture(invalid_paths)
assert_true(invalid_result$exit_status != 0L, "Invalid RDS must fail")
assert_true("rds_read_success" %in% invalid_result$report$gate$failed_checks, "Invalid RDS failure not recorded")

missing_time_object <- make_mock_object(c("0 h", "12 h", "0 h", "12 h", "0 h", "12 h"))
missing_time_paths <- write_fixture(file.path(test_root, "missing_time"), missing_time_object)
missing_time_result <- run_fixture(missing_time_paths)
assert_true(missing_time_result$exit_status != 0L, "Missing 24 h cells must fail")
assert_true("nonzero_cells_24h" %in% missing_time_result$report$gate$failed_checks, "Missing time point failure not recorded")

success_object <- make_mock_object(c("0 h", "12h", "24 hours", "0 h", "12h", "24 hours"))
success_paths <- write_fixture(file.path(test_root, "success"), success_object)
success_result <- run_fixture(success_paths)
assert_true(success_result$exit_status == 0L, paste(
  "Three-time-point fixture should pass:",
  paste(success_result$report$gate$failed_checks, collapse = ", ")
))
assert_true(identical(success_result$report$metadata$time_column, "sample_time"), "Time column must be detected from content")
assert_true(all(unlist(success_result$report$metadata$cells_per_timepoint) > 0L), "Every time point needs cells")

qc <- jsonlite::fromJSON(success_paths$qc, simplifyVector = FALSE)
required_sections <- c("dataset", "object", "metadata", "identity", "expression", "gate")
assert_true(all(required_sections %in% names(qc)), "QC JSON schema sections are incomplete")
assert_true(identical(qc$gate$status, "PASS"), "QC JSON must record PASS for valid fixture")

cli_qc <- file.path(test_root, "cli", "phase1_dataset_qc.json")
dir.create(dirname(cli_qc), recursive = TRUE)
cli_output <- system2(
  "Rscript",
  c(
    shQuote(file.path(project_root, "scripts", "data", "validate_phase1_data.R")),
    shQuote(success_paths$archive),
    shQuote(success_paths$rds),
    shQuote(success_paths$metadata),
    shQuote(cli_qc)
  ),
  stdout = TRUE,
  stderr = TRUE
)
cli_status <- attr(cli_output, "status") %||% 0L
assert_true(cli_status == 0L, "Validation CLI should return zero for a passing fixture")
assert_true(
  any(grepl("^PHASE1A_VALIDATION_SUMMARY status=PASS", cli_output)),
  "Validation CLI summary is missing"
)

cat("Phase 1A R validation tests: PASS\n")
