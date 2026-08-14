# Phase 3 CD4 0h vs CD8 0h summary

Gate status: `READY_FOR_GPT_REVIEW`

## Research question

Does GenePT-w preserve a detectable multivariate difference between naïve CD4 and naïve CD8 cells, and can the SCPA core framework detect it?

## Why 0h

The primary comparison uses 0h first to reduce activation/time confounding and isolate the cell-type distinction.

## GenePT-w and coverage

- CD4: 14894 cells; 14409 matched genes; median coverage 0.931967.
- CD8: 7041 cells; 14547 matched genes; median coverage 0.924674.
- CD8-minus-CD4 median coverage: -0.00729303; recorded as a potential confounder without an exclusion threshold.

## Canonical cohort

Actual 0h counts were CD4=4428 and CD8=1048. A fixed seed (20260810) selected 500/group, and the same IDs were used for both representations.

## SCPA-core adaptation

SCPA 1.6.2 calls `multicross::mcm()` on cells-by-genes pathway matrices. This adapter calls that same core function on aligned cells-by-features matrices without pretending the 1,536 GenePT dimensions are genes or a pathway.

## Original-expression reference

RNA/counts were normalized to 10,000 over each dataset's full gene set, log1p transformed, then aligned by 17085 exact shared symbols. This reference only checks that a population difference exists in original space.

## Interpretation limits

The two representations have different dimensions and geometry. Their raw MCM p/q values are not comparable representation-quality scores. This phase does not establish that GenePT is better, measure classifier accuracy, or identify important pathways/genes.

## Next planned work

Phase 4 will separately freeze separability metrics and True/Permuted/Random controls. The core pathway/gene-level interpretation question is reserved for its own later methodological decision gate.
