args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve run_phase3_scpa_core.R location")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
source(file.path(dirname(script_path), "scpa_core_adapter.R"))

argument_value <- function(name) {
  equals <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(equals) == 1L) return(sub(paste0("^", name, "="), "", equals[[1]]))
  index <- which(args == name)
  if (length(index) == 1L && index < length(args)) return(args[[index + 1L]])
  stop("Missing required argument: ", name)
}

genept_cd4_path <- argument_value("--genept-cd4")
genept_cd8_path <- argument_value("--genept-cd8")
expression_cd4_path <- argument_value("--expression-cd4")
expression_cd8_path <- argument_value("--expression-cd8")
output_path <- argument_value("--output")

require_scpa_core_packages()
toy <- run_scpa_core_toy_test()
if (!isTRUE(toy$passed)) stop("SCPA-core toy wiring test failed")

read_input <- function(path) {
  if (!file.exists(path)) stop("Missing Phase 3 adapter input: ", path)
  Matrix::readMM(path)
}

genept_result <- run_scpa_core_adaptation(
  read_input(genept_cd4_path),
  read_input(genept_cd8_path)
)
expression_result <- run_scpa_core_adaptation(
  read_input(expression_cd4_path),
  read_input(expression_cd8_path)
)

payload <- list(
  toy_test = toy,
  analyses = list(
    genept_w_cd4_vs_cd8 = genept_result,
    original_expression_cd4_vs_cd8 = expression_result
  ),
  interpretation_limit = paste(
    "The two raw p/q values are not representation-quality scores and must not",
    "be compared as evidence that one representation is better."
  )
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
temporary <- tempfile(pattern = ".phase3_scpa_core_", tmpdir = dirname(output_path))
on.exit(unlink(temporary), add = TRUE)
jsonlite::write_json(payload, temporary, auto_unbox = TRUE, pretty = TRUE, null = "null")
if (!file.rename(temporary, output_path)) stop("Could not atomically write: ", output_path)
cat("PHASE3_SCPA_CORE status=PASS output=", normalizePath(output_path), "\n", sep = "")
