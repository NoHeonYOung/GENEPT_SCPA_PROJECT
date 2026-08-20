# Phase 8 — Mean-shift mechanism decomposition

Status: **FROZEN BEFORE CONTROL GENERATION OR NEW MCM**

## 연구 질문

Phase 7B에서 mean-shift에 남은 null-calibrated GenePT–Vanilla AP 차이를 다음 세 요소와
일치하는 패턴으로 분해한다.

1. Correct gene↔embedding correspondence
2. GenePT vector-set geometry
3. Generic norm-matched 1,536D projection

이는 mechanism을 입증하는 인과 실험이 아니다. 동결 synthetic benchmark에서 네
representation의 recovery ordering을 비교하는 최소 decomposition이다.

## 입력과 범위

- 기존 Phase 7/7B 산출물은 read-only이며 SHA-256를 실행 전후 검증한다.
- 11 pathways × 20 draws를 유지한다.
- State는 null, mean shift 0.5, mean shift 1.0만 사용한다.
- 기존 cells, pathway genes/order, truth genes, perturbation seeds/strength, fallback state를
  그대로 사용한다. 새 truth selection은 하지 않는다.
- Cell subset과 mixed direction, 새 pathway/dataset은 제외한다.
- 기존 Vanilla와 TRUE GenePT full ranking을 재사용하며 이 두 branch의 새 MCM은 0이다.

## 네 representation

### VANILLA

Phase 7 zero-mask ranking을 그대로 재사용한다.

### TRUE GenePT

Phase 7 non-L2 `X@E_TRUE`와 exact subtraction-mask ranking을 그대로 재사용한다.

### PERMUTED GenePT

Pathway 내부에서 `E_TRUE`의 1,536D row 전체를 fixed-point 없는 derangement로 이동한다.
Column이나 vector 내부 값을 섞지 않는다. Vector multiset/value와 sorted row-norm multiset은
정확히 보존한다. 이는 Phase 6의 90% 초과 mapping-change 조건보다 엄격한 correspondence
destruction이다.

### RANDOM projection

각 expression gene에 독립 standard-Gaussian 1,536D direction을 생성하고 row-wise unit
normalization 후 같은 위치의 TRUE row norm으로 rescale한다. Phase 6 RANDOM의 분포,
dimension과 corresponding-row norm matching 정의를 유지한다. 실제 GenePT vector-set
geometry와 correct correspondence는 사용하지 않는다.

Control seed는 control type별 base에 `100000×pathway_index + draw_id`를 더한다. 각
pathway×draw에는 control당 realization 하나만 만들고 그 정확한 matrix를 null/0.5/1.0
state가 공유한다. Control HDF5를 먼저 생성해 state별 재생성 차이를 방지한다. Control
생성과 masking runner는 truth label 파일을 읽지 않는다.

## Workload

Manifest의 11 pathways는 `ΣK=353`, `Σ(K+1)=364`다.

```text
control당 20 draws × 3 states × 364 = 21,840 MCM
PERMUTED + RANDOM = 43,680 new MCM
```

예상값이 다르면 full run을 시작하지 않는다. Experiment checkpoint는 두 control의 full+K
mask 결과를 함께 저장하며 smoke checkpoint는 production에서 검증 후 재사용한다.

## Anti-leakage

- Control preparation input: Phase 7 manifest의 experiment metadata와 expression/embedding HDF5
- Masking input: Phase 8 execution manifest, expression HDF5, control HDF5
- Ground-truth/evaluation labels: evaluator에서만 Phase 7 ranking을 통해 읽음
- QC fields:
  `ground_truth_read_during_control_generation=false`,
  `ground_truth_read_during_masking=false`

## Primary endpoint와 unit

Primary endpoint는 raw AP가 아닌 representation별 matched-null calibrated AP다.

```text
NonNull_AP_M = mean(AP_mean_shift_0.5, AP_mean_shift_1.0)
Calibrated_AP_M = NonNull_AP_M - AP_matched_null_M
```

Primary paired unit은 pathway×draw, n=220이다. Phase 7B unit을 변경하지 않는다.

사전 정의 contrast:

1. `C_mapping = TRUE − PERMUTED`
2. `C_geometry = PERMUTED − RANDOM`
3. `C_projection = RANDOM − VANILLA`
4. `C_total = TRUE − VANILLA`

각 contrast에서 mean/median, descriptive normal 95% CI, two-sided Wilcoxon, matched-pairs
rank-biserial을 보고하고 4개 p-value를 Bonferroni 보정한다. `C_total`은 Phase 7B mean-shift
DiD `0.1429550762075937`을 1e-12 이내에서 재현해야 한다.
Ordering 해석에서 `>`는 contrast mean이 양수이고 Bonferroni p<0.05인 경우만
통계적으로 지지됐다고 표현한다. 그 외에는 차이가 입증되지 않았다고 쓰며 수치적 동일성을
뜻하는 등호로 해석하지 않는다.

Pathway별 calibrated AP와 세 decomposition contrast를 descriptive하게만 보고하며 pathway별
significance hunting은 금지한다. FALLBACK은 non-null Phase 7 truth selection 상태로 ALL,
NO_FALLBACK, FALLBACK을 모두 보고하고 primary를 재정의하지 않는다.

## Technical QC

- 기존 Phase 7/7B input hash 불변
- 기존 Vanilla/TRUE ranking 재사용 및 새 MCM 0
- PERMUTED/RANDOM 각각 21,840 completed MCM, failed 0, partial false
- 같은 pathway×draw의 세 state가 동일 control HDF5 path/hash 사용
- PERMUTED fixed point 0, vector multiset과 sorted row-norm multiset exact 보존
- RANDOM 1,536D, finite, corresponding TRUE row norm 오차≤1e-10
- Control 생성/masking truth-label read false
- `C_total` Phase 7B 재현

Figure는 final scientific QC 이후에도 별도 요청 전까지 생성하지 않는다.

## 해석

사전 정의 CASE A–E는 ordering을 설명하는 언어로만 사용한다. 예상 밖 ordering은 억지로
분류하지 않는다. `TRUE>PERMUTED`여도 semantics proven이라고 하지 않고 correct
correspondence contribution과 일치한다고만 표현한다. Geometry/projection contrast도 causal
또는 biological superiority 근거가 아니다. Control realization은 pathway×draw당 하나이므로
결과가 애매하거나 realization sensitivity가 우려될 때만 별도 Phase 8B를 사전 동결한다.
