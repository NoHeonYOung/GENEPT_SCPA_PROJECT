# Phase 6 semantic-specific controls

Gate status: `READY_FOR_GPT_REVIEW`

## Frozen method

True uses the correct official GenePT gene-row correspondence. Permuted preserves the exact within-pathway vector multiset and row norms while changing more than 90% of gene assignments. Random uses independent 1536-dimensional Gaussian directions scaled to each corresponding True-row norm. The same control realization is used for both groups and for full/masked calculations; no L2 normalization is used.

## Completed scope

- Phase 6A pathway-control MCM: 6000
- Phase 6B gene-mask control MCM: 11960
- Phase 6B control-baseline MCM: 240
- Paired resampling MCM: 900
- Pathway target/control summaries: 60
- Representative target/control summaries: 12

## Interpretation limit

These controls test whether results depend on the correct GenePT gene-to-vector correspondence relative to specified null constructions. They do not establish causal genes, biological correctness, predictive superiority, or generalization to CD8 or other datasets. Phase 7, CD8 generalization, classifiers, GenePT L2 gene masking, external networks and new pathway databases were not run.
