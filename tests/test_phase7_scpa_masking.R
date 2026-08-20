#!/usr/bin/env Rscript

# Cheap shared-library regression test. No SCPA/MCM call and no production data.
source("scripts/scpa/gene_masking_lib.R")
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

pair <- genept_non_l2_project_pair(expression, expression, embedding)
full_projection <- pair$a
for (gene_index in seq_len(ncol(expression))) {
  vanilla <- vanilla_zero_mask_pair(expression, expression, gene_index)
  zero_masked <- vanilla$a
  direct_projection <- zero_masked %*% embedding
  genept <- genept_non_l2_subtraction_mask_pair(
    full_projection, full_projection, expression, expression, embedding, gene_index
  )
  subtraction_projection <- genept$a
  stopifnot(max(abs(direct_projection - subtraction_projection)) < 1e-12)
  stopifnot(all(zero_masked[, gene_index] == 0))
}

stopifnot(abs(masking_score(0.01) - 2) < 1e-12)
stopifnot(abs(masking_q_score(0.01) - sqrt(2)) < 1e-12)

message("PHASE7_MASKING_ALGEBRA status=PASS scpa_calls=0 production_data=FALSE")
