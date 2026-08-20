# Phase 7 LLM-free synthetic ground-truth benchmark

Status: **PREPARATION PASS / REAL-DATA SMOKE PASS / FULL RUN NOT YET COMPLETED**

## 연구 질문

Naïve CD4 0 h의 homogeneous baseline에 알려진 gene perturbation을 주입했을 때,
Vanilla SCPA zero-mask와 GenePT-informed non-L2 subtraction-mask 중 어느 쪽이 truth gene을
더 잘 회복하는지 비교한다. 이 질문은 correct gene↔embedding correspondence가 특별한지를
검정한 Phase 6과 다르며 두 결론을 합치지 않는다.

## 고정 설계

- Source: GSE212270 Naïve CD4 0 h
- Cohort: seed 20260810, pseudo A/B 각 500개, disjoint
- 전처리: raw counts → full-transcriptome total 10,000 → log1p → pathway subset
- Perturbation: B의 normalized log1p 값에만 주입, 이후 재정규화 금지
- RNG: cohort와 독립된 base seed 20260901; 모든 experiment seed 기록
- Pathways: 정확히 11개(KEGG 6, REACTOME 5, HALLMARK 0)
- Draws: pathway × scenario마다 20회
- Scenarios: null, mean shift, 30% cell subset, mixed direction
- Strength: 비-null 0.5×/1.0× pooled SD, `scale=max(sd, 0.1)`
- Truth 수: `max(3, min(10, ceil(0.15 × pathway genes)))`
- 같은 pathway/scenario/draw의 0.5×/1.0×는 truth, subset, direction seed를 공유
- Detection≥0.1/SD>0 truth pool이 고정 truth 수보다 작으면 전체 pathway gene pool로
  fallback하고 scale floor 0.1을 적용한다. 이는 frozen 11개 pathway와 truth-count 공식을
  유지하기 위한 규칙이며 저발현/zero-baseline gene의 회복이 쉬워질 수 있는 confound다.
- Mixed negative는 detection≥0.5/median>0 strict pool을 우선 사용한다. 부족하면 실제
  Condition B에서 검출된 pool로 fallback하고, 가능한 negative 수가 절반보다 작아도 최소
  1개 negative와 나머지 positive로 구성한다. Gene별 selection rule과 0-clipped cell 수를
  기록한다. 두 fallback은 production MCM 결과를 보기 전 preparation feasibility에서 동결했다.

Null은 실제 expression을 바꾸지 않는다. 다만 Recall/AP/NDCG chance sanity check를 정의하기
위해 동일 규칙으로 uninjected evaluation-target set을 뽑는다. 산출물에서
`is_evaluation_target=true`, `is_ground_truth_perturbed=false`로 구분하며 이를 perturbed
gene이라고 부르면 안 된다.

## 방법

두 방법 모두 Phase 5와 같은 shared implementation
`scripts/scpa/gene_masking_lib.R`을 사용한다.

- Vanilla: 양 condition에서 gene column을 0으로 설정
- GenePT: `Z=X@E`, `Z_minus=Z-outer(X_g,E_g)`, non-L2
- Score: `S=-log10(max(raw_p, 1e-300))`
- Ranking signal: `delta_g=S_full-S_minus_g`, 내림차순, 동률은 gene symbol 사전순
- 모든 gene의 full ranking과 delta를 저장

## 평가와 통계

- Recall@3/5/10
- Average Precision (primary)
- NDCG@3/5/10, binary relevance
- 평균±표준편차: overall non-null, scenario, scenario×strength
- Paired Wilcoxon: GenePT AP vs Vanilla AP
- 추론 unit: pathway×scenario×draw 안에서 strength AP를 먼저 평균하여 0.5×/1.0×를
  독립 반복으로 세지 않는다.
- Scenario별 세 검정은 Bonferroni 보정하고 raw p와 adjusted p를 함께 기록한다.
- Effect size: matched-pairs rank-biserial, 양수면 GenePT 방향
- Null: exact random-ranking expectation과 AP 차이의 평균/SE를 기록하고
  `|mean difference| > 2 SE`면 방법론 경고를 낸다.

## 규모와 실행

- Experiments: `11 pathways × 20 draws × (1 null + 3 scenarios × 2 strengths) = 1,540`
- MCM: `101,920`
- GPU/LLM: 사용하지 않음
- Checkpoint: experiment 단위 atomic CSV; 재실행 시 검증 후 재사용
- Checkpoint compatibility: experiment input SHA-256와 masking protocol version이 모두
  일치할 때만 재사용

Actual-data preparation은 1,540 experiments로 PASS했다. Fallback audit은 560 experiments,
1,800 target-gene rows이며 결과에서 FALLBACK/NO_FALLBACK을 분리한다. 첫 null pathway의
42 genes, 86 MCM real-data smoke는 warning 없이 3.7분에 `SMOKE_PASS`했고 즉시 재실행에서
hash/protocol-compatible checkpoint 1개가 재사용됐다. 이 속도의 단순 환산은 single-core
약 74시간이며 12-core full run은 overhead/tail을 포함해 약 7–10시간 예상이다.

```bash
PYTHONPATH=src python scripts/phase7/prepare_llmfree_benchmark.py
PYTHONPATH=src python scripts/phase7/run_llmfree_masking.py --max-experiments 1 --cores 1
PYTHONPATH=src python scripts/phase7/run_llmfree_masking.py --cores 12
PYTHONPATH=src python scripts/phase7/evaluate_llmfree_benchmark.py
```

`--max-experiments 1`은 runtime smoke일 뿐 scientific result가 아니다. Full 결과가 모두
있을 때만 기본 evaluator가 실행된다. Partial debug는 명시적 `--allow-partial`로만 가능하다.

## 산출물

`data/processed/genept_scpa/phase7_llmfree_synthetic/`:

- `phase7_llmfree_manifest.json`
- `phase7_llmfree_cell_assignments.csv`
- `phase7_llmfree_ground_truth.csv`
- `phase7_llmfree_scpa_qc.json`
- `phase7_llmfree_scpa_smoke_qc.json` (`SMOKE_PASS`, full PASS와 구분)
- `phase7_llmfree_rankings.csv`
- `phase7_llmfree_metrics.csv` / `.json`
- `phase7_llmfree_aggregate.csv`
- `phase7_llmfree_statistics.csv`
- `phase7_llmfree_results.md`

Large expression HDF5와 checkpoints는
`data/interim/genept_scpa/phase7_llmfree_synthetic/`에 둔다.

### 고정 schema

- Cell assignment CSV: `source_cell_id`, `pseudo_condition`,
  `within_condition_index`, `cohort_seed`
- Ground-truth CSV: `experiment_id`, `draw_id`, `pathway`, full `gene` universe,
  `perturbation_type`, `perturbation_strength`, `perturbation_seed`,
  `is_evaluation_target`, `is_ground_truth_perturbed`, direction, target-cell count/fraction,
  direction-selection rule, pooled SD, applied log delta, clipped-cell count
- Manifest experiment: pathway/name/source, gene count, draw/scenario/strength/seed,
  A/B/gene/embedding HDF5 paths, post-injection renormalization flag
- Ranking CSV: experiment metadata, method, every gene, strict rank, signed delta,
  evaluation-target/injected flags
- Metric CSV: experiment metadata, method, truth/gene counts, Recall@3/5/10, AP,
  NDCG@3/5/10, truth-fallback flag 및 각 exact random-chance expectation
- Statistics CSV: pairing unit, pair count, GenePT−Vanilla mean AP difference,
  Wilcoxon statistic/raw p, scenario Bonferroni p, rank-biserial effect

## 해석 제한

- 높은 metric은 이 synthetic setup에서 injected signal을 더 잘 회복했다는 뜻뿐이다.
- Biological superiority, general validity, causality를 주장하지 않는다.
- Null warning이 있으면 비-null method comparison보다 방법론 점검이 먼저다.
- GenePT 열세는 K-gene Vanilla 공간과 1,536D projection geometry 차이의 confounding 가능성을
  함께 보고한다.
- Phase 6의 “TRUE가 PERMUTED와 구분되지 않음”과 Phase 7 recovery 성능을 하나의 주장으로
  합치지 않는다.
