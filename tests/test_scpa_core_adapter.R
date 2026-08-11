args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1L) args[[1]] else ".")
source(file.path(project_root, "scripts", "scpa", "scpa_core_adapter.R"))

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

toy <- run_scpa_core_toy_test()
assert_true(toy$passed, "Shifted toy population did not have a larger SCPA-style qval")
assert_true(
  toy$shifted$finite_reporting_qval > toy$nearly_identical$finite_reporting_qval,
  "Toy ordering is incorrect"
)
assert_true(toy$shifted$standard_pathway_analysis == FALSE, "Adapter was mislabeled")
assert_true(toy$shifted$input_orientation == "cells_by_features", "Orientation is incorrect")

official <- scpa_pathway_qvalues(c(0.001, 0.2))
assert_true(isTRUE(all.equal(official$adjusted_p, c(0.002, 0.4))), "Bonferroni correction differs")
assert_true(
  isTRUE(all.equal(official$qval, sqrt(-log10(c(0.002, 0.4))))),
  "Official base-10 qval transformation differs"
)

bad <- tryCatch(
  {
    run_scpa_core_adaptation(matrix(1, 4, 2), matrix(1, 4, 3))
    NULL
  },
  error = function(e) conditionMessage(e)
)
assert_true(!is.null(bad) && grepl("aligned feature count", bad), "Misalignment was not rejected")
cat("Phase 3 SCPA-core adapter tests: PASS\n")
