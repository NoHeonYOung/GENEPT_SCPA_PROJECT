require_scpa_core_packages <- function() {
  packages <- c("SCPA", "multicross", "Matrix", "jsonlite")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing SCPA-core package(s): ", paste(missing, collapse = ", "))
  }
}

validate_scpa_core_populations <- function(populations) {
  if (!is.list(populations) || length(populations) != 2L) {
    stop("SCPA-core adapter requires exactly two population matrices")
  }
  dimensions <- lapply(populations, dim)
  if (any(vapply(dimensions, is.null, logical(1)))) stop("Every population must be a matrix")
  if (any(vapply(dimensions, function(x) x[[1]] < 2L || x[[2]] < 1L, logical(1)))) {
    stop("Every population needs at least two cells and one feature")
  }
  if (dimensions[[1]][[2]] != dimensions[[2]][[2]]) {
    stop("SCPA-core populations must have the same aligned feature count")
  }
  if ((dimensions[[1]][[1]] + dimensions[[2]][[1]]) %% 2L != 0L) {
    stop("multicross::mcm matching requires an even combined cell count")
  }
  for (population in populations) {
    values <- if (inherits(population, "sparseMatrix")) population@x else as.vector(population)
    if (any(!is.finite(values))) stop("SCPA-core input contains NaN or Inf")
  }
  invisible(dimensions)
}

scpa_style_qval <- function(p_value) {
  if (!is.finite(p_value) || p_value < 0 || p_value > 1) stop("Invalid MCM p-value")
  if (p_value == 0) return(Inf)
  sqrt(-log10(p_value))
}

scpa_adjust_pathway_pvalues <- function(p_values, pathway_count = length(p_values)) {
  if (length(p_values) == 0L || any(!is.finite(p_values)) ||
      any(p_values < 0 | p_values > 1)) {
    stop("Invalid pathway MCM p-values")
  }
  if (pathway_count != length(p_values) || pathway_count < 1L) {
    stop("pathway_count must equal the eligible pathway universe")
  }
  stats::p.adjust(p_values, method = "bonferroni", n = pathway_count)
}

scpa_pathway_qvalues <- function(p_values) {
  adjusted <- scpa_adjust_pathway_pvalues(p_values)
  qval <- vapply(adjusted, scpa_style_qval, numeric(1))
  list(adjusted_p = adjusted, qval = qval)
}

run_mcm_raw <- function(population_a, population_b, level = 0.05) {
  require_scpa_core_packages()
  populations <- list(population_a, population_b)
  validate_scpa_core_populations(populations)
  runtime_warnings <- character()
  result <- withCallingHandlers(
    multicross::mcm(lapply(populations, as.matrix), level = level),
    warning = function(condition) {
      runtime_warnings <<- c(runtime_warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(
    raw_p = as.numeric(result[[1]]),
    decision = as.character(result[[2]]),
    warnings = unname(unique(runtime_warnings))
  )
}

run_scpa_core_adaptation <- function(population_a, population_b, level = 0.05) {
  require_scpa_core_packages()
  populations <- list(population_a, population_b)
  dimensions <- validate_scpa_core_populations(populations)
  dense_populations <- lapply(populations, function(x) as.matrix(x))
  result <- run_mcm_raw(dense_populations[[1]], dense_populations[[2]], level = level)
  p_value <- result$raw_p
  decision <- result$decision
  qval_raw <- scpa_style_qval(p_value)
  qval_floor <- scpa_style_qval(max(p_value, .Machine$double.xmin))
  list(
    implementation_source = list(
      scpa_version = as.character(utils::packageVersion("SCPA")),
      scpa_function = "SCPA::single_comparison -> multicross::mcm",
      multicross_version = as.character(utils::packageVersion("multicross")),
      called_function = "multicross::mcm",
      level = level
    ),
    adapter_name = "SCPA-core multivariate framework adaptation",
    standard_pathway_analysis = FALSE,
    input_orientation = "cells_by_features",
    population_cells = c(dimensions[[1]][[1]], dimensions[[2]][[1]]),
    aligned_features = dimensions[[1]][[2]],
    distance = "Euclidean (inside multicross::mcm)",
    p_value = p_value,
    decision = decision,
    scpa_style_qval = if (is.finite(qval_raw)) qval_raw else NULL,
    p_value_underflow = p_value == 0,
    finite_reporting_qval = qval_floor,
    multiple_testing = "one global hypothesis; Bonferroni factor = 1",
    warnings = result$warnings
  )
}

run_scpa_core_toy_test <- function(seed = 20260810L) {
  set.seed(seed)
  base <- matrix(stats::rnorm(40L * 4L), nrow = 40L, ncol = 4L)
  population_a <- base + matrix(stats::rnorm(40L * 4L, sd = 0.03), nrow = 40L)
  near_identical_b <- base + matrix(stats::rnorm(40L * 4L, sd = 0.03), nrow = 40L)
  shifted_b <- near_identical_b + 2
  near_result <- run_scpa_core_adaptation(population_a, near_identical_b)
  shifted_result <- run_scpa_core_adaptation(population_a, shifted_b)
  passed <- is.finite(near_result$finite_reporting_qval) &&
    is.finite(shifted_result$finite_reporting_qval) &&
    shifted_result$finite_reporting_qval > near_result$finite_reporting_qval
  list(
    seed = seed,
    nearly_identical = near_result,
    shifted = shifted_result,
    criterion = "shifted finite_reporting_qval > nearly_identical finite_reporting_qval",
    passed = passed
  )
}
