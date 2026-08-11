# Phase 4B Naive CD4 activation validation protocol

Status: `READY_FOR_USER_RUN`

## Purpose and scope

Phase 4A (Naive CD4 0h vs Naive CD8 0h) is preserved as a completed exploratory
cross-lineage feasibility analysis. Its high Bonferroni/qval floor and the former
unique ordering of tied qval=0 pathways make rank-only interpretation unsuitable,
but they do not demonstrate an SCPA implementation failure: SCPA 1.6.2, the raw-p
adapter, the log10 qval formula and average tie ranking have all been audited.

Phase 4B is the primary controlled benchmark. It fixes lineage to Naive CD4 and
varies activation time:

1. CD4 0h vs 12h (early activation)
2. CD4 12h vs 24h (later activation)
3. CD4 0h vs 24h (full activation; positive-control-like benchmark)

The current production run does not include CD8 activation or same-time lineage
comparisons. The runner's former nine-comparison capability is retained behind
the explicit `--comparison-set all_9` option, but it is not scheduled. Potential
CD8 0h-vs-24h generalization is Phase 4C and remains unrun. No generalization or
representation-superiority claim follows from Phase 4B.

## Source and canonical cells

The official GSE212270 integrated CD4/CD8 RDS objects and reusable Phase 2 sparse
exports already contain 0h, 12h and 24h cells, RNA/counts and RNA/data. No download
is required. The source audit is
`data/interim/genept_scpa/phase4_timecourse_source_audit.json`.

Seed 20260810 and 500 cells/group are frozen. CD4 0h exactly reuses Phase 3 IDs;
CD4 12h and 24h reuse the already frozen files in
`data/interim/genept_scpa/phase4_timecourse_sampling/`. The same timepoint IDs are
used by Vanilla, GenePT non-L2 and GenePT row-L2; they are never resampled between
branches or comparisons.

## Expression, pathways and projections

Preprocessing order is immutable:

`whole-transcriptome RNA/counts -> cell total 10,000 -> log1p -> pathway subset`

Pathway-local renormalization is prohibited. All three comparisons reuse the
Phase 4A frozen 123-pathway universe, exact paired genes and gene order. For every
pathway, expression column order must equal embedding row order. Any missing
pathway or changed paired list stops the run before production results are written.

Each comparison/pathway runs:

- Vanilla: `X_P` through `multicross::mcm()`
- GenePT primary: non-L2 `X_P @ E_P` through `multicross::mcm()`
- GenePT sensitivity: rowwise-L2(`X_P @ E_P`) through `multicross::mcm()`

Non-L2 remains primary regardless of the observed sensitivity result.

## SCPA, qval, ties and reporting

For each comparison and branch, raw p is `multicross::mcm()` output, adjusted p is
Bonferroni over all 123 eligible pathways, and qval is
`sqrt(-log10(adjusted_p))`. The log base remains 10. Descending ranks use average
ties, so equal qvals receive equal ranks.

Primary reporting precedes rank reporting and contains N, raw-p<0.05,
adjusted-p<0.05, qval-positive/zero counts, medians and qval-floor fractions.
Vanilla-vs-GenePT detection states use adjusted p<0.05. Non-L2-vs-L2 reporting
includes Spearman correlation and significant/qval-positive overlaps and
method-only detections. Rank correlations are secondary and descriptive.

Representative CD4 0h-vs-24h Vanilla raw p-values are cross-checked against
official `SCPA::compare_pathways()` with tolerance 1e-12. Historical Phase 4A is
read without overwriting it and is compared by qval-zero and significant counts.

## Pre-run gate and command

Python tests, the R adapter test, official raw-p cross-check, a one-pathway by
three-comparison smoke test, same-cell/same-gene assertions, qval/tie tests and
deterministic checks must pass before the user full run. The R core writes a
checkpoint after each completed comparison and reuses validated completed
comparisons on restart. Progress includes comparison, pathway, branch, elapsed
time and ETA.

```bash
PYTHONPATH=src python scripts/phase4/run_timecourse_validation.py \
  --comparison-set cd4_activation
```

This runs exactly 3 comparisons x 123 pathways x 3 branches = 1,107 MCM calls.
Full production is user-run only. Afterward, GPT must review qval floors,
significant counts, detection states, L2 sensitivity, positive-control-like
behavior, tie ranks, warnings and all cell/gene/pathway invariants before Phase 5.

## Outputs

All new production artifacts are under
`data/processed/genept_scpa/phase4_cd4_activation/`:

- `phase4_cd4_activation_all_results.csv`
- `phase4_cd4_activation_overview.csv`
- `phase4_cd4_activation_detection_states.csv`
- `phase4_cd4_activation_qc.json`
- `phase4_cd4_activation_summary.md`
- `phase4_cd4_activation_manifest.json`
- three comparison CSVs under `comparisons/`
- four primary figures plus an optional tie-aware rank scatter under `figures/`

No Phase 5 contribution analysis, CD8 generalization, semantic controls,
classifier, or unrelated analysis is run by this command.
