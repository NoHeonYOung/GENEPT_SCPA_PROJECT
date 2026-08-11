# Phase 3 SCPA qval implementation audit

Audit date: 2026-08-11

Result: `PASS`

## Compared implementations

- Installed `SCPA` version: 1.6.2
- Installed `multicross` version: 2.1.0
- Official local function inspected: `SCPA:::single_comparison`
- Project adapter inspected: `scripts/scpa/scpa_core_adapter.R`

## Official SCPA 1.6.2 convention

For every eligible pathway, SCPA obtains `multicross::mcm(path_subset,
level=0.05)[[1]]` as the raw p-value. It then computes:

```r
adjPval <- stats::p.adjust(
  Pval,
  method = "bonferroni",
  n = length(pathways_filtered)
)
qval <- sqrt(-log10(adjPval))
```

Thus the correction universe is the number of pathways retained after the
`min_genes`/`max_genes` filter, and the logarithm is base 10. Larger qval means a
stronger multivariate difference under the SCPA display convention.

## Phase 3 conclusion

Phase 3 tested one global hypothesis for each representation. Its Bonferroni
factor is therefore one, so the adapter's `sqrt(-log10(raw_p))` is exactly the
official SCPA result for that one-hypothesis universe. There is no log-base or
multiple-testing mismatch. Historical raw p-values and derived qval fields remain
unchanged.

## Phase 4 frozen implementation

Phase 4 first obtains every raw pathway p-value for the identical eligible
pathway universe in each branch. It then applies Bonferroni correction with
`n = eligible_pathway_count` separately to the Vanilla and GenePT-informed raw
p-value vectors, followed by `sqrt(-log10(adjusted_p))`.

## Phase 4B re-audit and tie handling

The installed source was inspected again before the activation validation extension.
It still confirms base-10 `log10`, not natural-log `log`. Therefore no historical
raw p-value or derived qval is changed. The separate Phase 4B CD4 activation reports
use average descending ranks so pathways with equal qval, including qval=0, receive
equal ranks. Detection categories use Bonferroni adjusted p < 0.05 and do not depend
on arbitrary ordering within ties.
