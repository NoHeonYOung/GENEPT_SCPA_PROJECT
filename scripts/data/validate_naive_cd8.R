args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve validate_naive_cd8.R location")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
source(file.path(dirname(script_path), "phase1_validation_lib.R"))
source(file.path(dirname(script_path), "naive_cd8_validation_lib.R"))

archive_path <- if (length(args) >= 1L) args[[1]] else file.path(project_root, "data", "raw", "genept_scpa", "GSE212270_integrated_naive_cd8.rds.gz")
rds_path <- if (length(args) >= 2L) args[[2]] else file.path(project_root, "data", "raw", "genept_scpa", "GSE212270_integrated_naive_cd8.rds")
metadata_path <- if (length(args) >= 3L) args[[3]] else file.path(project_root, "data", "interim", "genept_scpa", "naive_cd8_download_metadata.json")
output_path <- if (length(args) >= 4L) args[[4]] else file.path(project_root, "data", "interim", "genept_scpa", "naive_cd8_dataset_qc.json")
cd4_rds_path <- if (length(args) >= 5L) args[[5]] else file.path(project_root, "data", "raw", "genept_scpa", "GSE212270_integrated_naive_cd4.rds")
cd4_qc_path <- if (length(args) >= 6L) args[[6]] else file.path(project_root, "data", "interim", "genept_scpa", "phase1_dataset_qc.json")

result <- validate_naive_cd8_dataset(archive_path, rds_path, metadata_path, output_path, cd4_rds_path, cd4_qc_path)
report <- result$report
cat("NAIVE_CD8_VALIDATION_SUMMARY", paste0("status=", report$gate$status), paste0("cells=", report$object$cells),
  paste0("features=", report$object$features), paste0("time_column=", report$metadata$time_column %||% "NA"),
  paste0("shared_genes=", report$comparison_with_naive_cd4$exact_shared_gene_count %||% "NA"),
  paste0("failed_checks=", length(report$gate$failed_checks)), paste0("warnings=", length(report$gate$warnings)),
  paste0("qc_json=", normalizePath(result$output_path, mustWork = FALSE)), sep = " ")
cat("\n")
quit(save = "no", status = result$exit_status)
