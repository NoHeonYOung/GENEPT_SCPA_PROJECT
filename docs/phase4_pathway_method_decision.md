# Phase 4 pathway comparison method decision

Status: frozen before the full-data result; ready for user run.

## Research question

For the same biological pathways, canonical cells and paired genes, how does a
GenePT semantic projection change pathway-level SCPA rankings relative to Vanilla
expression-space SCPA? This phase measures agreement, disagreement and relative
reordering. It does not test whether GenePT is more accurate.

## Whole-cell GenePT-w versus pathway-specific projection

Published whole-cell GenePT-w aggregates the entire transcriptome into one
1,536-dimensional vector per cell and rowwise L2-normalizes it. Its 1,536
coordinates are semantic dimensions, not genes or pathway labels. Consequently,
whole-cell GenePT-w cannot be sliced by a conventional pathway gene list, and it
must not be passed to `compare_pathways()` as 1,536 fake genes.

Phase 4 returns to full-transcriptome-normalized original expression, selects the
genes for each pathway, aligns their official GenePT embeddings and recomputes a
pathway-specific representation:

```text
X_P: cells x p paired pathway genes
E_P: p paired pathway genes x 1536 embedding dimensions
Z_P = X_P x E_P: cells x 1536 semantic dimensions
```

The GenePT branch is named **GenePT-informed pathway projection followed by an
SCPA-core multivariate comparison**. It is not published whole-cell GenePT-w and
is not standard `SCPA::compare_pathways()` pathway analysis.

## Cohort and preprocessing policy

- Comparison: naïve CD4 0 h versus naïve CD8 0 h.
- Canonical Phase 3 samples are reused without resampling: 500 cells/group,
  seed 20260810.
- Source: sparse `RNA/counts`.
- Order: normalize each cell over its complete transcriptome to total 10,000,
  apply log1p, and only then subset pathway genes.
- Pathway-local total-count normalization is prohibited because it would erase
  the full-transcriptome expression scale.

## Primary paired-gene policy

For pathway P:

```text
paired_gene_set(P) = pathway genes
                     intersection CD4 features
                     intersection CD8 features
                     intersection official GenePT keys
```

Both branches use this exact symbol set and lexicographically sorted gene order.
This design isolates representation effects from gene-loss confounding in the
GenePT branch. It intentionally differs from a secondary standard Vanilla analysis
that might use every shared CD4/CD8 gene; no such secondary analysis is primary.

`expression_coverage_cd4/cd8` is defined per group as the total log1p-normalized
expression mass in the retained paired genes divided by the mass in all pathway
genes shared by CD4 and CD8, using the canonical cells. It measures expression
mass retained after the GenePT-mappability intersection; it is not a cell-exclusion
threshold.

The official SCPA v1.6.2 combined metabolic collection is retained with Hallmark,
KEGG and Reactome pathways and the frozen 15–500 gene thresholds. Eligibility is
recomputed from the primary paired gene count, and both branches use the same
eligible pathway universe.

The preflight count is 243 input pathways and 123 eligible primary paired pathways
(120 excluded), compared with 124 analyzed in Phase 1. No threshold was changed in
response; the one-pathway difference follows the stricter GenePT-mappable paired
gene intersection and will be recorded again by the full run.

## Vanilla and GenePT-informed branches

Vanilla compares `X_CD4,P` and `X_CD8,P` directly with the MCM core used by SCPA.
The GenePT-informed branch compares `Z_CD4,P` and `Z_CD8,P` using the same
`multicross::mcm()` core. No hidden sampling occurs inside the adapter.

For each branch, raw p-values across all eligible pathways are Bonferroni-adjusted
with `n = eligible_pathway_count`; SCPA qval is then
`sqrt(-log10(adjusted_p))` using base-10 logarithms.

## L2 versus non-L2 decision

Primary projection is non-L2:

```text
Z_P = X_P x E_P
```

This preserves pathway expression magnitude together with semantic geometry,
which is appropriate for a paired comparison against Vanilla SCPA that also uses
expression magnitude. Rowwise L2-normalized Z is implemented only as a
predeclared sensitivity option. Full-data results cannot be used to change the
primary choice.

## Effective-rank caveat

Although Z has 1,536 columns, `rank(Z) <= rank(E_P) <= p`. The output therefore
does not contain 1,536 independent biological variables. Phase 4 records pathway
gene count, embedding rank, projected effective rank, rank deficiency and
singular-value summaries. It asserts `projected_rank <= paired_gene_count`.

## Comparison and interpretation

`rank_delta = genept_rank - vanilla_rank`. A negative value means that a pathway
moves upward in the GenePT-informed ranking; a positive value means that it moves
downward. Spearman/Kendall correlation, Top-10/Top-20 overlap and Jaccard values
measure ranking agreement only.

The geometries and feature spaces differ, so smaller raw p, larger qval or a
higher rank cannot establish superiority or accuracy. Phase 4 can report which
pathways agree, disagree or are relatively reordered after projection.

## Phase 5 regeneration hook

The manifest stores original/shared/paired genes, expression feature order,
embedding keys, match types, canonical cell files and checksums, preprocessing,
projection configuration and source artifact provenance. It is sufficient to
regenerate Vanilla inputs and GenePT projections without storing every large
pathway matrix.

Phase 5 may apply matched leave-one-gene-out or gene masking to both branches, but
Phase 4 does not execute contribution analysis and will not call a qval difference
a gene-importance percentage.

## Planned Phase 6 controls

Phase 6 will freeze and compare True GenePT, gene-to-embedding Permuted GenePT and
dimension-matched Random embeddings, together with repeated sampling and
robustness. Phase 7 will use biological or synthetic ground truth for any claim
about accuracy. Neither control family is run in Phase 4.

## Post-run reporting validation extension

The initial CD4 0h-vs-CD8 0h run is retained as Phase 4A. Because many pathways
had Bonferroni adjusted p=1 and qval=0, its earlier deterministic unique ordering
must not be interpreted as evidence among tied pathways. Phase 4B uses tie-aware
average ranks, reports adjusted-p detection categories and validates CD4
0h-vs-24h against official SCPA raw p-values. Its primary production scope is
exactly three within-CD4 activation comparisons (0h-vs-12h, 12h-vs-24h and
0h-vs-24h) over the same frozen 123 paired pathways. The former nine-comparison
capability is retained but is not the current production plan. This does not alter
the primary non-L2 decision, establish generalization or begin gene contribution.
