args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1L) args[[1]] else ".")
source(file.path(project_root, "scripts", "data", "phase1_validation_lib.R"))
source(file.path(project_root, "scripts", "data", "genept_export_lib.R"))

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

counts <- Matrix::Matrix(
  matrix(c(1, 0, 2, 3, 0, 1), nrow = 3L),
  sparse = TRUE,
  dimnames = list(c("G1", "G2", "G3"), c("cell1", "cell2"))
)
object <- SeuratObject::CreateSeuratObject(counts, project = "phase2_export_mock")
object$Hour <- c(0L, 24L)
object$Cell_Type <- c("Resting", "Activated")

test_root <- tempfile("phase2_export_")
dir.create(test_root)
on.exit(unlink(test_root, recursive = TRUE), add = TRUE)
result <- export_seurat_counts_for_genept(
  object,
  output_dir = test_root,
  dataset = "naive_cd4",
  source_file = file.path(test_root, "mock.rds"),
  compatibility = list(source_rds_modified = FALSE)
)
assert_true(all(file.exists(unlist(result$paths))), "One or more export files are missing")
roundtrip <- Matrix::readMM(result$paths$counts)
assert_true(identical(dim(roundtrip), c(3L, 2L)), "Sparse Matrix Market dimensions changed")
assert_true(identical(readLines(result$paths$genes), c("G1", "G2", "G3")), "Gene order changed")
assert_true(identical(readLines(result$paths$cell_ids), c("cell1", "cell2")), "Cell order changed")
metadata <- utils::read.csv(result$paths$metadata, stringsAsFactors = FALSE)
assert_true(identical(metadata$cell_id, c("cell1", "cell2")), "Metadata cell IDs changed")
manifest <- jsonlite::fromJSON(result$paths$manifest, simplifyVector = FALSE)
assert_true(identical(manifest$assay, "RNA") && identical(manifest$layer, "counts"), "Wrong assay/layer")
assert_true(isTRUE(manifest$source_object_modified == FALSE), "Source modification flag is wrong")
cat("Phase 2 Seurat sparse export tests: PASS\n")
