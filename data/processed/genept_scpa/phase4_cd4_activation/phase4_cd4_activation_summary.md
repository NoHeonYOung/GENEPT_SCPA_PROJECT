# Phase 4B Naive CD4 activation validation summary

Gate status: `READY_FOR_GPT_REVIEW`

## Source audit

12h source status: `ALREADY_PRESENT_IN_RDS`. No download was performed. The primary run uses only the frozen CD4 0h/12h/24h groups (500 cells each).

## Qval and tie audit

Installed SCPA 1.6.2 uses `sqrt(-log10(Bonferroni-adjusted p))`; no formula correction was required. Historical Phase 4A output is read-only and preserved. Current reports use average ranks, so equal qval values receive equal ranks.

## Official SCPA cross-check

Representative CD4 0h-vs-24h Vanilla pathways passed raw-p equality against `SCPA::compare_pathways()` (3 pathways; tolerance 1e-12).

## Branch-level primary metrics

| Comparison | Branch | N | raw p<.05 | adj p<.05 | q>0 | q=0 | median raw p | median adj p | median q | q-floor fraction |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| cd4_0h_vs_12h | vanilla | 123 | 96 | 84 | 93 | 30 | 4.62766e-09 | 5.69202e-07 | 2.49895 | 0.244 |
| cd4_0h_vs_12h | genept | 123 | 90 | 75 | 85 | 38 | 9.57254e-08 | 1.17742e-05 | 2.22015 | 0.309 |
| cd4_0h_vs_12h | l2 | 123 | 69 | 60 | 64 | 59 | 0.00196887 | 0.242172 | 0.784778 | 0.480 |
| cd4_12h_vs_24h | vanilla | 123 | 76 | 50 | 65 | 58 | 0.00342232 | 0.420945 | 0.613004 | 0.472 |
| cd4_12h_vs_24h | genept | 123 | 65 | 51 | 60 | 63 | 0.00946978 | 1 | 0 | 0.512 |
| cd4_12h_vs_24h | l2 | 123 | 53 | 37 | 45 | 78 | 0.179741 | 1 | 0 | 0.634 |
| cd4_0h_vs_24h | vanilla | 123 | 119 | 116 | 119 | 4 | 3.09496e-20 | 3.80681e-18 | 4.17366 | 0.033 |
| cd4_0h_vs_24h | genept | 123 | 119 | 108 | 116 | 7 | 1.7918e-17 | 2.20391e-15 | 3.82842 | 0.057 |
| cd4_0h_vs_24h | l2 | 123 | 88 | 82 | 86 | 37 | 3.63637e-13 | 4.47273e-11 | 3.21705 | 0.301 |

## Vanilla vs GenePT detection states

| Comparison | Vanilla q>0 | GenePT q>0 | L2 q>0 | Vanilla adj<.05 | GenePT adj<.05 | L2 adj<.05 | Both | V-only | G-only | Neither |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| cd4_0h_vs_12h | 93 | 85 | 64 | 84 | 75 | 60 | 74 | 10 | 1 | 38 |
| cd4_12h_vs_24h | 65 | 60 | 45 | 50 | 51 | 37 | 46 | 4 | 5 | 68 |
| cd4_0h_vs_24h | 119 | 116 | 86 | 116 | 108 | 82 | 107 | 9 | 1 | 6 |

## GenePT non-L2 vs L2 sensitivity

| Comparison | non-L2 vs L2 Spearman | significant overlap | q-positive overlap | non-L2-only sig | L2-only sig |
| --- | ---: | ---: | ---: | ---: | ---: |
| cd4_0h_vs_12h | 0.8917 | 57 | 63 | 18 | 3 |
| cd4_12h_vs_24h | 0.7753 | 34 | 42 | 17 | 3 |
| cd4_0h_vs_24h | 0.9298 | 82 | 86 | 26 | 0 |

## Positive-control interpretation

Scenario A-compatible: the CD4 activation positive-control-like contrast has more qval-positive pathways than historical resting CD4-vs-CD8 in both primary branches. This descriptively supports a stronger activation contrast, but does not establish cause or accuracy.

## Detection-state definition

At Bonferroni adjusted p < 0.05, each pathway is classified as Both significant, Vanilla-only significant, GenePT-only significant or Neither significant. These categories measure detection agreement, not accuracy.

## Timing patterns

Early, late, persistent, endpoint-only, mixed and none labels are descriptive combinations of adjusted-p significance across 0→12, 12→24 and 0→24. Rank agreement is secondary because qval floors create large tied blocks. These are not causal biological classifications.

## Scope limits

CD8 activation/generalization, gene-level contribution, True/Permuted/Random control, classifier and Phase 5 were not run. This benchmark supports representation comparison only; it does not establish GenePT superiority or generalization. Phase 5 remains blocked until GPT review.
