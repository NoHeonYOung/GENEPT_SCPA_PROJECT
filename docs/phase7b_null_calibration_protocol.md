# Phase 7B — Null-calibration sensitivity analysis protocol

Status: **FROZEN BEFORE EVALUATION**

## 목적과 범위

Phase 7의 고정 synthetic benchmark에서 두 방법 모두 theoretical chance보다 높은 null AP를
보인 원인을 점검하고, 기존 non-null AP를 동일 pathway·draw의 method-specific null AP로
보정한다. Phase 7B는 새 Phase나 새 perturbation 실험이 아니라 evaluator-level sensitivity
analysis다.

- Phase 0–6과 Phase 7 perturbation/MCM은 재실행하지 않는다.
- 기존 101,920 MCM에서 저장한 full gene ranking만 재사용한다.
- 기존 `data/processed/genept_scpa/phase7_llmfree_synthetic/`는 읽기 전용이다.
- figure와 새로운 biological claim은 만들지 않는다.

## Primary endpoint와 추론 단위

- Endpoint: Average Precision(AP)
- Unit: `pathway × scenario × draw`
- 같은 unit의 strength 0.5/1.0 AP를 먼저 평균한다.
- Null matching key: `pathway × draw × method`
- Non-null scenario: mean shift, cell subset, mixed direction

각 unit에서 다음을 계산한다.

```text
Vanilla_cal = Vanilla_nonnull_AP - Vanilla_matched_null_AP
GenePT_cal  = GenePT_nonnull_AP  - GenePT_matched_null_AP
DiD         = GenePT_cal - Vanilla_cal
```

전체 평균끼리 빼지 않고 matched unit에서 먼저 계산한다.

## Primary/sensitivity test family

다음 4개 DiD를 사전 정의한다.

1. Overall non-null
2. Mean shift
3. Cell subset
4. Mixed direction

각 scope에서 n, raw/calibrated AP의 method별 mean/median, raw method gap, null method gap,
DiD mean/median, descriptive normal-approximation 95% CI, two-sided paired Wilcoxon against 0,
matched-pairs rank-biserial을 보고한다. 4개 p-value에 Bonferroni를 적용한다.

## Truth-selection fallback sensitivity

Non-null `truth_fallback_used`로 ALL, NO_FALLBACK, FALLBACK을 나누고 overall과 세 scenario를
모두 보고한다. 최소 contrast는 raw `GenePT−Vanilla`, matched null `GenePT−Vanilla`, calibrated
DiD다. NO_FALLBACK을 해석상 중요하게 보되 다른 strata를 제거하거나 primary 결과를
재정의하지 않는다.

## Diagnostic null 정의

모든 diagnostic은 null experiment의 고정 Vanilla/GenePT ranking을 그대로 사용하고 truth
label만 다시 표본화한다. seed는 20261001, resample은 pathway×draw당 1,000회이며 replacement,
truth-count 축소, threshold 완화는 금지한다. Truth 수는 Phase 7 공식을 유지한다.

### A. ELIGIBLE_POOL_MATCHED

목적은 expression-conditioned eligibility 자체가 null AP inflation을 만드는지 진단하는 것이다.
Pooled Condition A/B baseline에서 detection fraction≥0.10이고 SD>0인 gene만 pool로 쓴다.

- eligible 수≥truth_count인 9 pathway만 계산한다.
- Primary A result: 9 pathways × 20 draws = 180 units
- A와 B의 직접 비교도 반드시 같은 9 pathways/180 units에서 한다.
- 다음 두 pathway는 A를 정의하지 않고 `NOT_ESTIMABLE`로 기록한다.
  - `KEGG_GLYCOSPHINGOLIPID_BIOSYNTHESIS_LACTO_AND_NEOLACTO_SERIES`
  - `REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES`
- 사유: `eligible pool size < truth_count, therefore eligible-pool null is not identifiable`

### B. PATHWAY_WIDE_UNIFORM

해당 pathway의 전체 analysis gene에서 truth_count만큼 균등 표본화한다.

- 전체 11 pathways 결과를 보고한다.
- A와 직접 비교할 동일한 feasible 9 pathways 결과도 별도로 보고한다.
- 이는 recovery benchmark를 대체하지 않는 uniform-truth diagnostic control이다.

### C. PHASE7_FALLBACK_MATCHED

실제 Phase 7 truth-selection mechanism을 재현하는 secondary diagnostic이다. Eligible pool이
충분하면 eligible genes를, 부족하면 전체 pathway genes를 사용한다. 전체 11 pathways에서
계산하며 A와 혼합하거나 A라고 부르지 않는다.

Primary diagnostic contrast는 같은 9-pathway subset의 `A−B`다. B all-11, C all-11 및
FALLBACK/NO_FALLBACK 결과는 secondary로 투명하게 보고한다.

## Baseline-expression confounding diagnostic

Perturbation 전 pooled Condition A/B baseline에서 gene별 mean expression, detection fraction,
SD를 계산한다. Pathway size 영향을 줄이기 위해 feature를 pathway 내부 average-rank
percentile로 변환한다.

- 기존 null pseudo-truth target과 non-target의 feature-rank percentile 차이
- Vanilla/GenePT null masking delta와 각 feature의 pathway 내부 Spearman
- Vanilla/GenePT의 `-rank`와 각 feature의 pathway 내부 Spearman

이는 pseudo-truth selection과 masking sensitivity가 같은 baseline feature를 선호하는지 보는
descriptive diagnostic이며 인과 분석이 아니다.

## 입력 무결성과 출력 격리

실행 시 config/protocol 및 Phase 7 manifest, QC, rankings, metrics, aggregate, statistics, HDF5의
SHA-256를 snapshot에 기록하고 종료 전 다시 확인한다. 입력 hash가 변하면 실행을 실패시킨다.
모든 결과는 다음 별도 경로에만 저장한다.

`data/processed/genept_scpa/phase7b_null_calibration/`

필수 파일:

- `phase7b_protocol_snapshot.json`
- `phase7b_null_calibrated_metrics.csv`
- `phase7b_did_statistics.csv`
- `phase7b_fallback_sensitivity.csv`
- `phase7b_null_truth_diagnostics.csv`
- `phase7b_expression_confounding.csv`
- `phase7b_summary.json`
- `phase7b_results.md`

## 판정과 해석 제한

- `COMPLETED_WITH_WARNING`: technical validity를 유지하며 보정 뒤 일부/전체 advantage가 남지만
  scenario 및 diagnostic 제한을 함께 명시한다.
- `COMPLETED_WITH_NULL_EXPLAINED`: raw advantage의 상당 부분이 null bias로 설명되고 보정 효과가
  제한적이거나 scenario-specific이다.
- `INVALID`: leakage, 복구 불가능한 evaluator/matching 오류, 심각한 protocol 위반에만 사용한다.

보정 효과가 사라지는 것은 valid negative result이며 `INVALID` 사유가 아니다. GenePT의 biological
superiority, causal-gene identification, general accuracy, semantic mechanism을 주장하지 않는다.
Phase 6의 gene↔embedding correspondence specificity와 Phase 7/7B의 synthetic recovery/null
calibration은 서로 다른 연구 질문으로 유지한다.
