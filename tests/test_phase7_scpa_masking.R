#!/usr/bin/env Rscript

# Cheap algebra-only regression test. No SCPA/MCM call and no production data.
expression <- matrix(
  c(1.0, 2.0, 3.0, 0.5, 0.0, 1.5),
  nrow = 2,
  byrow = TRUE
)
embedding <- matrix(
  c(1.0, 0.0, 0.0, 2.0, 1.0, 1.0),
  nrow = 3,
  byrow = TRUE
)

full_projection <- expression %*% embedding
for (gene_index in seq_len(ncol(expression))) {
  zero_masked <- expression
  zero_masked[, gene_index] <- 0
  direct_projection <- zero_masked %*% embedding
  subtraction_projection <- full_projection - tcrossprod(
    expression[, gene_index], embedding[gene_index, ]
  )
  stopifnot(max(abs(direct_projection - subtraction_projection)) < 1e-12)
  stopifnot(all(zero_masked[, gene_index] == 0))
}

message("PHASE7_MASKING_ALGEBRA status=PASS scpa_calls=0 production_data=FALSE")
