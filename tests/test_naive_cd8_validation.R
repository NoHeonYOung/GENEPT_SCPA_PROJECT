args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1L) args[[1]] else ".")
source(file.path(project_root, "scripts", "data", "phase1_validation_lib.R"))
source(file.path(project_root, "scripts", "data", "naive_cd8_validation_lib.R"))

assert_true <- function(value, message) if (!isTRUE(value)) stop(message, call. = FALSE)
if (!requireNamespace("SeuratObject", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
  stop("SeuratObject and jsonlite are required for naïve CD8 validation tests")
}

make_mock <- function(identity, hours, features, cell_types = NULL) {
  counts <- matrix(seq_len(length(features) * length(hours)) %% 4L, nrow = length(features),
    dimnames = list(features, paste0(identity, "_cell", seq_along(hours))))
  object <- SeuratObject::CreateSeuratObject(counts, project = paste0("naive_", identity, "_mock"))
  object <- SeuratObject::SetAssayData(object, assay = "RNA", layer = "data", new.data = log1p(counts))
  object$sample_time <- hours
  object$population <- paste("naive", toupper(identity), "T cells")
  if (!is.null(cell_types)) object$Cell_Type <- cell_types
  object
}

write_fixture <- function(directory, cd8_object = NULL, cd4_object = NULL, invalid_cd8 = FALSE) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  archive <- file.path(directory, "GSE212270_integrated_naive_cd8.rds.gz")
  rds <- file.path(directory, "GSE212270_integrated_naive_cd8.rds")
  metadata <- file.path(directory, "naive_cd8_download_metadata.json")
  qc <- file.path(directory, "naive_cd8_dataset_qc.json")
  cd4_rds <- file.path(directory, "GSE212270_integrated_naive_cd4.rds")
  cd4_qc <- file.path(directory, "phase1_dataset_qc.json")
  if (invalid_cd8) saveRDS(list(not = "Seurat"), rds) else saveRDS(cd8_object, rds)
  saveRDS(cd4_object, cd4_rds)
  input <- file(rds, "rb"); bytes <- readBin(input, "raw", n = file.info(rds)$size); close(input)
  output <- gzfile(archive, "wb"); writeBin(bytes, output); close(output)
  record <- list(
    geo_accession = "GSE212270", filename = basename(archive),
    download_source = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212270/suppl/GSE212270_integrated_naive_cd8.rds.gz",
    file_size_bytes = unname(file.info(archive)$size), sha256 = sha256_file(archive), gzip_integrity = TRUE,
    recorded_at = "2026-08-10T00:00:00+00:00", recorded_at_utc = "2026-08-10T00:00:00+00:00"
  )
  jsonlite::write_json(record, metadata, auto_unbox = TRUE, pretty = TRUE)
  jsonlite::write_json(list(gate = list(status = "PASS")), cd4_qc, auto_unbox = TRUE)
  list(archive = archive, rds = rds, metadata = metadata, qc = qc, cd4_rds = cd4_rds, cd4_qc = cd4_qc)
}

run_fixture <- function(x) validate_naive_cd8_dataset(x$archive, x$rds, x$metadata, x$qc, x$cd4_rds, x$cd4_qc)
test_root <- tempfile("naive_cd8_validation_"); dir.create(test_root)
on.exit(unlink(test_root, recursive = TRUE), add = TRUE)

missing <- list(
  archive = file.path(test_root, "missing", "GSE212270_integrated_naive_cd8.rds.gz"),
  rds = file.path(test_root, "missing", "GSE212270_integrated_naive_cd8.rds"),
  metadata = file.path(test_root, "missing", "naive_cd8_download_metadata.json"),
  qc = file.path(test_root, "missing", "naive_cd8_dataset_qc.json"),
  cd4_rds = file.path(test_root, "missing", "GSE212270_integrated_naive_cd4.rds"),
  cd4_qc = file.path(test_root, "missing", "phase1_dataset_qc.json")
)
missing_result <- run_fixture(missing)
assert_true(missing_result$exit_status != 0L, "Missing files must fail")
assert_true("archive_file_exists" %in% missing_result$report$gate$failed_checks, "Missing archive was not recorded")
assert_true(file.exists(missing$qc), "Missing-file QC JSON was not written")

hours <- c("0 h", "12h", "24 hours", "0 h", "12h", "24 hours")
cd4 <- make_mock("cd4", hours, c("IL7R", "CCR7", "KLF2", "CD4", "MAL", "LTB"), rep(c("Resting", "Activated"), 3L))
cd8 <- make_mock("cd8", hours, c("IL7R", "CCR7", "KLF2", "CD8A", "MAL", "LTB"), rep(c("Resting", "Activated"), 3L))

invalid <- write_fixture(file.path(test_root, "invalid"), cd8, cd4, invalid_cd8 = TRUE)
invalid_result <- run_fixture(invalid)
assert_true("seurat_object_class" %in% invalid_result$report$gate$failed_checks, "Non-Seurat input was not rejected")

missing_time_cd8 <- make_mock("cd8", rep(c("0h", "12h"), 3L), rownames(cd8), rep("Resting", 6L))
missing_time <- write_fixture(file.path(test_root, "missing_time"), missing_time_cd8, cd4)
missing_time_result <- run_fixture(missing_time)
assert_true("nonzero_cells_24h" %in% missing_time_result$report$gate$failed_checks, "Missing 24 h was not detected")

malformed_time_cd8 <- make_mock("cd8", rep("not-a-timepoint", 6L), rownames(cd8), rep("Resting", 6L))
malformed_time <- write_fixture(file.path(test_root, "malformed_time"), malformed_time_cd8, cd4)
malformed_time_result <- run_fixture(malformed_time)
assert_true("time_column_detected" %in% malformed_time_result$report$gate$failed_checks, "Malformed time metadata was not rejected")

success <- write_fixture(file.path(test_root, "success"), cd8, cd4)
success_result <- run_fixture(success)
assert_true(success_result$exit_status == 0L, paste("Valid fixture failed:", paste(success_result$report$gate$failed_checks, collapse = ",")))
assert_true(success_result$report$metadata$cell_type$column_present, "Cell_Type presence was not recorded")
assert_true(success_result$report$expression$rna_feature_count == 6L, "RNA feature count is wrong")
assert_true(sum(unlist(success_result$report$metadata$cell_type$counts)) == 6L, "Cell_Type counts are wrong")
comparison <- success_result$report$comparison_with_naive_cd4
assert_true(comparison$exact_shared_gene_count == 5L, "Shared gene count is wrong")
assert_true(comparison$cd4_only_gene_count == 1L && comparison$cd8_only_gene_count == 1L, "Group-only gene counts are wrong")
assert_true(!comparison$gene_order_identical, "Different feature lists must not be order-identical")
assert_true(!comparison$intersection_matrix_created, "Validator must not create an intersection matrix")

qc <- jsonlite::fromJSON(success$qc, simplifyVector = FALSE)
required <- c("dataset", "object", "metadata", "identity", "expression", "comparison_with_naive_cd4", "gate")
assert_true(all(required %in% names(qc)), "Naïve CD8 QC JSON schema is incomplete")
assert_true(identical(qc$gate$status, "PASS"), "Valid fixture QC must report PASS")

cli_qc <- file.path(test_root, "cli", "naive_cd8_dataset_qc.json"); dir.create(dirname(cli_qc), recursive = TRUE)
cli_output <- system2("Rscript", c(
  shQuote(file.path(project_root, "scripts", "data", "validate_naive_cd8.R")),
  shQuote(success$archive), shQuote(success$rds), shQuote(success$metadata), shQuote(cli_qc),
  shQuote(success$cd4_rds), shQuote(success$cd4_qc)
), stdout = TRUE, stderr = TRUE)
assert_true((attr(cli_output, "status") %||% 0L) == 0L, "Naïve CD8 validation CLI failed")
assert_true(any(grepl("^NAIVE_CD8_VALIDATION_SUMMARY status=PASS", cli_output)), "CLI summary is missing")

cat("Naïve CD8 R validation tests: PASS\n")
