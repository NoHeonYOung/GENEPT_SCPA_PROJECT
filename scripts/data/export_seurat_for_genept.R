args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve export_seurat_for_genept.R location")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
source(file.path(dirname(script_path), "phase1_validation_lib.R"))
source(file.path(dirname(script_path), "genept_export_lib.R"))

dataset <- "naive_cd4"
if (length(args) >= 2L && args[[1]] == "--dataset") dataset <- args[[2]]
equals_arg <- grep("^--dataset=", args, value = TRUE)
if (length(equals_arg) == 1L) dataset <- sub("^--dataset=", "", equals_arg[[1]])
valid_datasets <- c("naive_cd4", "naive_cd8")
if (!dataset %in% valid_datasets) stop("--dataset must be naive_cd4 or naive_cd8")

filename <- if (dataset == "naive_cd4") {
  "GSE212270_integrated_naive_cd4.rds"
} else {
  "GSE212270_integrated_naive_cd8.rds"
}
source_file <- file.path(project_root, "data", "raw", "genept_scpa", filename)
output_dir <- file.path(project_root, "data", "interim", "genept_scpa", "phase2_export", dataset)
loaded <- load_seurat_for_genept(source_file)
result <- export_seurat_counts_for_genept(
  loaded$object,
  output_dir = output_dir,
  dataset = dataset,
  source_file = source_file,
  compatibility = loaded$compatibility
)
rm(loaded)
invisible(gc())
cat(
  "GENEPT_SEURAT_EXPORT",
  "status=PASS",
  paste0("dataset=", dataset),
  paste0("genes=", result$manifest$genes),
  paste0("cells=", result$manifest$cells),
  paste0("manifest=", normalizePath(result$paths$manifest)),
  sep = " "
)
cat("\n")
