# Phase 7B null-calibration sensitivity 결과

상태: **COMPLETED_WITH_WARNING**
Technical QC: **PASS** — 새 MCM 0회, 기존 Phase 7 입력은 read-only, null ranking AP 재계산 최대 오차 0

## Matched null calibration

| Scope | n | Raw GenePT−Vanilla | Matched null gap | Calibrated DiD | Bonferroni p | Rank-biserial |
|---|---:|---:|---:|---:|---:|---:|
| overall_non_null | 660 | 0.072088 | 0.034609 | 0.037479 | 0.002549 | 0.154 |
| mean_shift | 220 | 0.177564 | 0.034609 | 0.142955 | 2.661e-11 | 0.535 |
| cell_subset | 220 | 0.030353 | 0.034609 | -0.004256 | 1 | -0.014 |
| mixed_direction | 220 | 0.008346 | 0.034609 | -0.026263 | 0.5755 | -0.114 |

Overall raw method gap 중 matched-null subtraction으로 설명되는 기술적 비율은 48.0%다.
이 값은 causal decomposition이 아니라 같은 pathway/draw의 arithmetic sensitivity다.

## Fallback sensitivity

- NO_FALLBACK n=400: raw gap 0.071532, null gap 0.038694, DiD 0.032838
- FALLBACK n=260: raw gap 0.072942, null gap 0.028324, DiD 0.044618

## Null truth-selection diagnostics

- A eligible-pool, feasible 9 pathways/180 units: Vanilla AP 0.334468, GenePT AP 0.369925, method gap 0.035456
- B uniform, same 9 pathways: Vanilla AP 0.253156, GenePT AP 0.253308, method gap 0.000152
- B uniform, all 11 pathways: Vanilla AP 0.255229, GenePT AP 0.255367, method gap 0.000139
- C Phase-7 fallback-matched, all 11 pathways: Vanilla AP 0.321178, GenePT AP 0.350344, method gap 0.029167
- Primary A−B on the same 9 pathways: Vanilla 0.081312, GenePT 0.116616; method-gap A−B 0.035305
- A infeasible 2 pathways는 replacement, truth-count 축소, threshold 완화 없이 `NOT_ESTIMABLE`로 기록했다.

## Baseline-expression confounding diagnostic

- 기존 null target의 within-pathway feature-rank percentile은 non-target보다 mean expression +0.319, detection +0.321, SD +0.319 높았다.
- Null masking delta와 feature의 mean Spearman은 Vanilla에서 mean/detection/SD 각각 0.140/0.141/0.140, GenePT에서 0.173/0.180/0.164였다.
- Constant masking-delta unit은 Spearman을 억지로 정의하지 않고 `CONSTANT_INPUT_NOT_ESTIMABLE`로 기록했다.

## 해석

Frozen synthetic benchmark에서 GenePT의 raw AP gap은 양수였지만 두 방법 모두 expression-conditioned
null recovery를 보였다. Pathway/draw-matched method-specific null calibration 뒤 residual difference는
overall과 mean shift에서 남았고, cell subset과 mixed direction에서는 지지되지 않았다. A와 B의 차이 및
baseline diagnostic은 expression-conditioned truth selection이 null inflation의 주요 원인과 일치함을
보이지만 인과적 분해는 아니다. 상세치는 `phase7b_expression_confounding.csv`에서 확인한다.

이 결과는 biological superiority, causal-gene identification, general validity 또는 semantic mechanism을
입증하지 않는다. Phase 6 correspondence-specificity 결론과 Phase 7/7B synthetic recovery 결론은 분리한다.
