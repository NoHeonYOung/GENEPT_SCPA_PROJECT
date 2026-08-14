# Phase 5 pathway-internal gene masking sensitivity

Gate status: `READY_FOR_GPT_REVIEW`

## Frozen scope

Targets are all 30 Phase 4B pathway-comparison instances with discordant Vanilla/GenePT non-L2 Bonferroni detection states. The same Phase 4B cells, preprocessing, pathways, paired genes and gene order were reused. GenePT L2, Phase 6 controls, CD8 and classifiers were not run.

## Method

Vanilla masks gene g by setting `X_P[:,g]=0`. GenePT non-L2 uses `Z_P - outer(X_P[:,g], E_P[g,:])`. The primary sensitivity is `delta[-log10(raw p)] = score_full - score_masked`, with raw p clipped only for scoring at 1e-300. Positive values indicate that masking weakened the observed pathway signal; negative values indicate that masking strengthened it.

## Technical results

- Targets completed: 30
- Genes evaluated per branch: 1135
- Total masking MCM evaluations: 2270
- Failed MCM calls: 0
- Runtime warnings: 0
- Baseline maximum raw-p absolute difference: 0
- Median signed-delta Spearman across defined targets: 0.155
- Median absolute-influence Spearman across defined targets: 0.334
- Vanilla threshold-flip gene instances: 194
- GenePT threshold-flip gene instances: 399

## Interpretation limit

These results describe representation-dependent gene-masking sensitivity. They do not identify causal genes, establish biological correctness, or show that either representation is superior. Semantic specificity and accuracy require the deferred Phase 6/7 controls.
