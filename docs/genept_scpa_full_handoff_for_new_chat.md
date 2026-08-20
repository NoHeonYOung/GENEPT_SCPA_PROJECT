# GenePT × SCPA 새 GPT 채팅용 전체 인수인계

이 문서는 기존 GPT 채팅 기록 없이도 프로젝트를 이어갈 수 있도록 만든 단일 인수인계
문서다. 새 채팅에서는 이 파일을 먼저 첨부하거나 전체 내용을 읽게 한 뒤 작업을 요청한다.

작성 기준일: **2026-08-19 (Asia/Seoul)**
프로젝트 루트: `/home/node00/nhy_python/GenePT_SCPA`
현재 Git 기준: `main`, HEAD `273491e` (`phase7-handoff-20260814` tag)
현재 active phase: **Phase 7**
가장 중요한 상태: **Phase 6은 최종 PASS/COMPLETED이며 Phase 7은 LLM-free preparation/smoke PASS, full run pending이다.**

---

## 0. 새 GPT에게 가장 먼저 전달할 짧은 요약

이 프로젝트는 GenePT의 gene-semantic embedding을 SCPA의 multivariate pathway
comparison framework와 결합할 수 있는지 단계적으로 검증한다.

- Phase 0–6은 완료됐다.
- Phase 1은 full GSE212270 naïve CD4 데이터로 Vanilla SCPA를 재현했다.
- Phase 2는 저자 공개 ada-002 1,536D embedding으로 published GenePT-w를 재현했다.
- Phase 3은 whole-cell GenePT-w 공간에서도 CD4 0 h와 CD8 0 h 차이가 검출되는지
  feasibility를 확인했다.
- Phase 4는 동일 cells/pathways/genes에서 Vanilla와 pathway-specific
  `GenePT-informed Z_P = X_P @ E_P`를 비교했다.
- Phase 5는 동일 gene을 하나씩 masking했을 때 Vanilla와 GenePT의 sensitivity ordering이
  달라지는 것을 확인했다.
- Phase 6은 True/Permuted/Random semantic controls를 수행했다. Technical gate는 PASS지만
  True와 within-pathway Permuted가 명확히 구분되지 않았다. 따라서 correct
  gene-to-embedding correspondence에 특이적인 semantic superiority는 지지되지 않는다.
- Phase 7의 이전 LLM 설계는 폐기됐다. 현재는 synthetic ground truth에서 Vanilla와
  GenePT non-L2 masking만 비교하는 CPU-only benchmark다.
- Protocol/실행 코드와 1,540-experiment preparation은 완료됐고 real-data 86-MCM smoke가
  PASS했다. Full 101,920-MCM run은 아직 하지 않았다.
- 기존 결과만으로 GenePT superiority, biological correctness, causality를 주장하면 안 된다.

---

## 1. 연구 질문과 용어

### 전체 연구 질문

1. GenePT-w 변환 뒤에도 biologically distinct cell populations의 차이가 보존되는가?
2. SCPA가 사용하는 `multicross::mcm()` multivariate comparison이 GenePT 표현에서도
   그 차이를 검출할 수 있는가?
3. pathway-specific GenePT projection이 Vanilla SCPA와 어떤 agreement/disagreement를
   보이는가?
4. pathway 내부 gene masking sensitivity가 Vanilla와 GenePT에서 어떻게 달라지는가?
5. 그 차이가 올바른 gene ↔ embedding correspondence에 특이적인 semantic effect인가?
6. synthetic ground truth가 있을 때 Vanilla와 GenePT 중 누가 perturbed genes를 더 잘
   회복하는가?

### 반드시 구분할 표현

```text
Published whole-cell GenePT-w
  Cell × Gene expression @ Gene × 1536 embedding
  → Cell × 1536
  → final row-wise unit L2 normalization

Project-specific pathway GenePT representation
  pathway expression X_P: cells × pathway genes
  pathway embeddings E_P: pathway genes × 1536
  Z_P = X_P @ E_P
  Phase 4 이후 primary는 non-L2
```

GenePT의 1,536 dimensions는 genes나 curated pathways가 아니다. 이를 SCPA의 pathway
genes처럼 부르면 안 된다. Whole-cell 1,536D 분석은 `SCPA-core multivariate framework
adaptation`, pathway 분석은 `GenePT-informed pathway projection`이라고 부른다.

### SCPA statistic

- 설치 버전: SCPA 1.6.2, multicross 2.1.0
- pathway raw p-value: `multicross::mcm()`
- multiple testing: Bonferroni
- SCPA qval convention: `sqrt(-log10(adjusted_p))`
- qval은 클수록 stronger multivariate difference다.
- qval=0은 weakest/floor이며, 동률은 average rank로 처리한다.
- qval/p-value 크기만으로 서로 다른 dimension/geometry representation의 우열을
  주장하지 않는다.

---

## 2. 데이터와 고정 입력

### Dataset A: GSE212270 integrated naïve CD4

- 원본: `GSE212270_integrated_naive_cd4.rds.gz` / extracted RDS
- full object이며 축소 예제 데이터가 아니다.
- 14,894 cells, 17,856 features
- Hour별: 0 h 4,428 / 12 h 4,547 / 24 h 5,919
- assays: RNA, integrated; active assay RNA
- RNA layers: counts, data, scale.data
- gene symbol missing/duplicate: 0/0
- serialized Seurat 3.1.5 object를 메모리에서만 현재 SeuratObject와 호환되게 갱신한다.
- source RDS는 수정하지 않았다.
- archive SHA-256:
  `f9ad7cfbe8bee87a28cec76dd66d442ad2c5ddb243942c51d905713b9b2b7842`

### Dataset B: GSE212270 integrated naïve CD8

- acquisition/QC PASS
- 7,041 cells, 17,942 features
- Hour별: 0 h 1,048 / 12 h 2,066 / 24 h 3,927
- CD4와 같은 accession family를 사용해 cross-study/platform confounding을 줄였다.
- 동일 family가 모든 batch effect를 제거한다고 가정하지 않는다.

### Pathways

- SCPA v1.6.2 `combined_metabolic_pathways.csv`
- Hallmark/KEGG/Reactome metabolic pathways 243개
- SHA-256:
  `6bc5977da3fa60f86d5ffb59fc938740bf418fa4d976182a314d65479eb8b744`
- Phase 1에서는 124개가 분석 가능했다.
- Phase 4 paired-gene policy 이후 frozen eligible universe는 123개다.

### GenePT embedding

- 논문 저자 공개 Zenodo DOI `10.5281/zenodo.10833191`
- 파일: `GenePT_gene_embedding_ada_text.pickle`
- model: `text-embedding-ada-002`
- dimension: 1,536
- 새 OpenAI API 호출로 재생성하지 않는다.
- embedding SHA-256:
  `fd297510ddd3040744033fde0b0f2cf15a40ac8b2fd2fb02f10667295e55c862`

---

## 3. Phase별 진행 내용과 결론

## Phase 0 — 프로젝트 분리와 protocol freeze

상태: **PASS**

- 기존 LOO/frequency/hybrid 프로젝트와 논리적으로 분리했다.
- package는 `src/gene_embedding_project/genept_scpa/` 아래에 둔다.
- phase gate와 상태는 `config/genept_scpa.yaml`로 관리한다.
- 연구 결정은 `docs/genept_scpa_decision_log.md`에 기록한다.
- 결과를 보고 threshold/metric/control을 바꾸지 않는 원칙을 세웠다.

## Phase 1A — full naïve CD4 acquisition/QC

상태: **PASS**

- 2.4 GB processed RDS를 다운로드하고 gzip, size, SHA-256, Seurat 구조를 검증했다.
- 다운로드/검증 스크립트가 지정한 위치가 최종 위치이므로 파일을 따로 옮길 필요가 없다.
- legacy Seurat object에서 `images` slot 오류가 났던 문제는 source object를 바꾸지 않고
  메모리상 compatibility update를 적용해 해결했다.

주요 파일:

- `data/interim/genept_scpa/phase1_download_metadata.json`
- `data/interim/genept_scpa/phase1_dataset_qc.json`
- `docs/phase1_data_download.md`

## Phase 1B — Vanilla SCPA reproduction

상태: **PASS**

고정 조건:

- `RNA/data`의 기존 log1p-normalized layer
- `SCPA::seurat_extract()`, pseudocount 0.001
- seed 20260810, 최대 500 cells/group
- matched genes 15–500
- official/default single-core CPU 실행

분석:

- Global 0/12/24
- Pairwise 0 vs 12
- Pairwise 12 vs 24
- Pairwise 0 vs 24
- 별도 official-like reference: Resting 0 h vs Activated 24 h

핵심 결과:

- 124 pathways의 qval/FC가 finite, warning/error 없음
- global 상위에 amino-acid metabolism, polyamines, respiratory electron transport,
  heme metabolism, oxidative phosphorylation, glycolysis가 검출됐다.
- official reference에서 HALLMARK glycolysis rank 5, Reactome arachidonic-acid
  metabolism rank 56, KEGG arachidonic-acid metabolism rank 69였다.
- arachidonic pathways는 작은 mean FC에서도 명확한 qval을 보여 SCPA의 multivariate
  해석과 qualitative agreement를 보였다.
- numerical identity나 exact paper Figure 4 reproduction을 주장하지 않는다.

주요 파일:

- `data/interim/genept_scpa/phase1b_reproduction_summary.md`
- `data/interim/genept_scpa/phase1b_scpa_qc.json`
- `data/processed/genept_scpa/phase1/*.csv`
- `data/processed/genept_scpa/phase1/figures/`

## Phase 2 — published GenePT-w reproduction

상태: **PASS**

전처리:

```text
RNA/counts
→ 전체 transcriptome 기준 cell-wise total 10,000 normalization
→ log1p
→ official artifact exact key / official HGNC alias key alignment
→ expression-weighted aggregation
→ official notebook denominator 사용
→ row-wise unit L2 normalization
```

핵심 결과:

- output: 14,894 cells × 1,536, float32
- exact matches 14,087
- official alias matches 322
- total matched 14,409
- unmatched 3,447
- median raw-count mass coverage 0.931967
- finite values 22,877,184 / 22,877,184
- zero vectors 0
- final row norms approximately 1
- optimized vs direct max absolute error `2.71e-08`
- deterministic repeat max error 0

중요: Phase 2의 final L2는 published whole-cell GenePT-w 재현에 필수다. 이 L2와
Phase 4 pathway-level L2 sensitivity는 서로 다른 맥락이다.

주요 파일:

- `data/interim/genept_scpa/phase2_genept_w_qc.json`
- `data/interim/genept_scpa/phase2_genept_w_summary.md`
- `data/processed/genept_scpa/phase2/`
- `docs/phase2_genept_w_protocol.md`

## Phase 3 — whole-cell GenePT-w CD4 0 h vs CD8 0 h feasibility

상태: **PASS**

- CD4/CD8 actual 0 h cells에서 seed 20260810으로 500/group을 한 번 고정했다.
- 같은 cell IDs를 original-expression과 GenePT-w 양쪽에 재사용했다.
- original-expression reference는 각 dataset full transcriptome을 total 10,000 +
  log1p한 뒤 exact shared symbols 17,085개로 정렬했다.
- GenePT-w는 1,536D, original reference는 17,085 genes다.

결과:

- GenePT-w CD4 vs CD8: p=`4.6605e-12`, SCPA-style qval=`3.3662`
- Original expression CD4 vs CD8: p=`5.5516e-60`, qval=`7.6978`
- 두 representation 모두 population difference를 검출했다.
- dimension/geometry가 다르므로 두 p/q의 크기를 representation quality로 비교하지 않는다.
- 이 Phase는 classifier accuracy나 GenePT superiority를 검증하지 않았다.

주요 파일:

- `data/interim/genept_scpa/phase3_cd4_cd8_qc.json`
- `data/interim/genept_scpa/phase3_cd4_cd8_summary.md`
- `data/interim/genept_scpa/phase3_sampling/`
- `docs/phase3_cd4_cd8_protocol.md`
- `docs/phase3_qval_implementation_audit.md`

## Phase 4A — exploratory CD4 0 h vs CD8 0 h pathway comparison

상태: **COMPLETED / HISTORICAL PRESERVED**

방법:

- 동일 500 CD4 / 500 CD8 cells
- 동일 123 pathways와 exact paired genes/order
- full-transcriptome normalization 후 pathway subset
- Vanilla: `X_P`
- GenePT primary: non-L2 `Z_P = X_P @ E_P`
- row-L2 projection은 sensitivity branch로 별도 실행

결과:

- Vanilla vs GenePT pathway rank Spearman `0.682166`
- Kendall `0.672931`
- Top-10 overlap 6, Top-20 overlap 12
- 대부분 agreement가 있지만 일부 큰 rank shift가 있었다.
- qval=0 tied pathways에 과거 unique name-order rank를 주던 artifact를 발견했다.
- raw p/qval은 보존하고 이후 report는 average tied rank를 사용했다.
- Phase 4A는 잘못된 분석이 아니지만 cross-lineage/qval-floor 문제 때문에 primary
  validation 근거로 사용하지 않는다.

주요 파일:

- `data/processed/genept_scpa/phase4/vanilla_vs_genept_pathway_comparison.csv`
- `data/processed/genept_scpa/phase4/pathway_projection_manifest.json`
- `data/interim/genept_scpa/phase4_pathway_qc.json`
- `data/interim/genept_scpa/phase4_pathway_summary.md`
- `data/processed/genept_scpa/phase4/figures/`

## Phase 4B — primary Naïve CD4 activation validation

상태: **PASS / COMPLETED**

Primary comparisons:

- CD4 0 h vs 12 h
- CD4 12 h vs 24 h
- CD4 0 h vs 24 h

각 Hour에서 frozen 500 cells, 동일 123 pathways/paired genes를 사용했다.

Primary Vanilla vs GenePT-informed 결과:

| Comparison | Vanilla adj p<.05 | GenePT adj p<.05 | Both | Vanilla-only | GenePT-only | Neither |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 h vs 12 h | 84 | 75 | 74 | 10 | 1 | 38 |
| 12 h vs 24 h | 50 | 51 | 46 | 4 | 5 | 68 |
| 0 h vs 24 h | 116 | 108 | 107 | 9 | 1 | 6 |

- total discordant pathway-comparison targets: 30
- Vanilla-only 23, GenePT-only 7
- 0 h vs 24 h가 strongest positive-control-like activation contrast였다.
- 3 representative Vanilla pathways는 official `SCPA::compare_pathways()` raw p와
  tolerance `1e-12`에서 일치했다.
- detection categories는 agreement/disagreement이지 accuracy가 아니다.

### Phase 4 L2 처리

- primary science branch는 non-L2 `X_P @ E_P`다.
- row-L2 sensitivity도 계산되어 원 CSV/summary/QC에 보존돼 있다.
- 발표에서는 L2 결과를 사용하지 않기로 결정했다.
- 2026-08-18에 Phase 4 figure의 L2 막대, 범례, `non-L2` 표기를 제거하고
  `GenePT-informed`로 통일했다.
- 원 결과 테이블의 L2 columns는 재현성을 위해 삭제하지 않았다.

주요 파일:

- `data/processed/genept_scpa/phase4_cd4_activation/phase4_cd4_activation_all_results.csv`
- `data/processed/genept_scpa/phase4_cd4_activation/phase4_cd4_activation_detection_states.csv`
- `data/processed/genept_scpa/phase4_cd4_activation/phase4_cd4_activation_overview.csv`
- `data/processed/genept_scpa/phase4_cd4_activation/phase4_cd4_activation_qc.json`
- `data/processed/genept_scpa/phase4_cd4_activation/phase4_cd4_activation_summary.md`
- `data/processed/genept_scpa/phase4_cd4_activation/figures/`

Phase 4C CD8 generalization은 **NOT SCHEDULED / NOT RUN**이다.

## Phase 5 — pathway-internal leave-one-gene-out masking

상태: **PASS / COMPLETED**

대상:

- Phase 4B discordant 30 pathway-comparison instances
- breakdown: 0–12 h 11 / 12–24 h 9 / 0–24 h 10
- branch당 총 1,135 gene instances

masking:

```text
Vanilla:
  X_minus_g[:, g] = 0

GenePT-informed:
  Z_full = X @ E
  Z_minus_g = Z_full - outer(X[:, g], E[g, :])

score = -log10(raw_p), raw_p clip only for scoring at 1e-300
delta_g = score_full - score_minus_g
```

- positive delta: gene masking이 observed pathway signal을 약화
- negative delta: gene masking이 signal을 강화
- primary: signed delta ordering
- secondary: absolute influence ordering

결과:

- masking MCM evaluations 2,270
- failed MCM 0, warnings 0
- Phase 4B baseline raw-p max difference 0
- Vanilla zero-mask vs physical-removal difference 0
- GenePT subtraction vs direct recomputation max difference `8.8818e-16`
- median signed-delta Spearman `0.155`
- median absolute-influence Spearman `0.334`
- Vanilla threshold-flip gene instances 194
- GenePT threshold-flip gene instances 399

해석:

- Vanilla와 GenePT가 gene-level sensitivity ordering을 다르게 만든다는 결과다.
- causal gene, biological importance, correctness, superiority 증거가 아니다.

주요 파일:

- `data/processed/genept_scpa/phase5_gene_contribution/phase5_target_pathways.csv`
- `data/processed/genept_scpa/phase5_gene_contribution/phase5_gene_masking_all_results.csv`
- `data/processed/genept_scpa/phase5_gene_contribution/phase5_gene_rank_comparison.csv`
- `data/processed/genept_scpa/phase5_gene_contribution/phase5_pathway_summary.csv`
- `data/processed/genept_scpa/phase5_gene_contribution/phase5_gene_contribution_qc.json`
- `data/processed/genept_scpa/phase5_gene_contribution/phase5_gene_contribution_summary.md`
- `data/processed/genept_scpa/phase5_gene_contribution/figures/`

## Phase 6 — True/Permuted/Random semantic controls

상태: **최종 PASS / COMPLETED**
주의: QC 파일의 historical runtime gate 문자열은 `READY_FOR_GPT_REVIEW`지만,
decision log D-0030에서 scientific review 후 최종 PASS로 닫았다.

Representations:

- True: 정확한 official GenePT gene ↔ embedding row correspondence
- Permuted: pathway 내부에서 exact vector multiset/value/dimension/norm을 유지하면서
  gene assignment를 90% 초과 변경
- Random: 독립 Gaussian 1,536D direction을 대응 True row norm으로 scaling
- 같은 control realization을 두 group과 full/masked 계산에 사용
- Phase 6에는 L2 normalization을 사용하지 않았다.

실행 규모:

- Phase 6A pathway controls: 30 targets × 2 controls × 100 = 6,000 MCM
- Phase 6B gene masks: 6 representative targets × 2 controls × 20 reps = 11,960 MCM
- Phase 6B control baselines: 240 MCM
- paired resampling: 10 × 30 × 3 = 900 MCM
- total: 19,100 MCM
- failed MCM 0, warnings 0, technical criteria 19/19 true

robustness:

- median True rank stability: `-0.179933`
- True-minus-Permuted sign consistency: `0.466667`
- True-minus-Random sign consistency: `0.566667`

### Phase 6 최종 과학적 결론 — 문구를 약화/강화하지 말 것

> GenePT-informed representation은 pathway/gene sensitivity를 바꾸지만, 정확한
> gene-to-embedding correspondence는 pathway 내부 permutation control과 구분되지
> 않았다. Random과의 차이는 structured embedding-geometry effect 가능성과 양립하지만
> 이를 입증하지 않으며, resampling robustness는 약했다. 따라서 semantic-specific
> superiority는 현재 결과에서 지지되지 않는다.

다음처럼 말하면 안 된다.

- “GenePT가 무의미하다”
- “embedding geometry effect가 입증됐다”
- “GenePT가 더 생물학적으로 정확하다”
- “gene importance/causality를 찾았다”

주요 파일:

- `data/processed/genept_scpa/phase6_semantic_controls/phase6_control_targets.csv`
- `data/processed/genept_scpa/phase6_semantic_controls/phase6_pathway_control_all_results.csv`
- `data/processed/genept_scpa/phase6_semantic_controls/phase6_pathway_control_summary.csv`
- `data/processed/genept_scpa/phase6_semantic_controls/phase6_gene_control_all_results.csv`
- `data/processed/genept_scpa/phase6_semantic_controls/phase6_gene_control_summary.csv`
- `data/processed/genept_scpa/phase6_semantic_controls/phase6_resampling_results.csv`
- `data/processed/genept_scpa/phase6_semantic_controls/phase6_semantic_control_qc.json`
- `data/processed/genept_scpa/phase6_semantic_controls/phase6_semantic_control_summary.md`
- `data/processed/genept_scpa/phase6_semantic_controls/figures/`

## Phase 7 — LLM-free synthetic ground-truth recovery

상태: **PREPARATION PASS / REAL-DATA SMOKE PASS / FULL RUN PENDING**

Phase 6 negative result 때문에 연구가 실패한 것이 아니다. Phase 6은 “correct semantic
correspondence가 control보다 특이적인가?”라는 가설 검정에서 명확한 지지를 얻지 못한
결과다. Phase 7은 실제 known perturbation ground truth를 만들어 accuracy를 직접 비교하기
위해 설계됐다.

### Phase 7 research question

Vanilla SCPA zero-mask와 GenePT-informed non-L2 subtraction-mask 중 어느 방법이 known
synthetic perturbed genes의 full ranking을 더 잘 회복하는가?

### Frozen cohort

- source: GSE212270 naïve CD4 0 h
- sorted source cell IDs
- NumPy PCG64 seed 20260810
- 1,000 unique cells without replacement
- seeded permutation의 첫 500 = pseudo Condition A
- 나머지 500 = pseudo Condition B
- A/B는 disjoint

### Frozen pathways

- Phase 4 paired universe만 사용
- 기존 동결 selection 결과의 정확한 pathway 이름을 config에 고정
- KEGG 6개, REACTOME 5개
- eligible HALLMARK 0개이며 기준을 완화하지 않음
- total expected 11 pathways
- preparation은 PASS했지만 full MCM이 아직이라 scientific result table은 생성되지 않음

### Perturbations

- null
- mean shift
- B cells 중 seeded 30% cell subset shift
- mixed direction
- effect strengths: pooled baseline SD의 0.5, 1.0
- perturbation은 normalized log1p space의 B에만 주입
- 주입 후 재정규화하지 않음
- truth gene count: `max(3, min(10, ceil(0.15 * pathway_gene_count)))`
- pathway×scenario마다 20 independent draws
- cohort와 분리된 perturbation base seed 20260901
- null은 uninjected evaluation-target을 뽑아 chance metric을 정의하되 perturbed gene으로
  부르지 않음

### 비교 methods

- Vanilla SCPA: pathway expression gene column zero-mask
- GenePT non-L2: `Z_minus=Z_full-outer(X_g,E_g)`
- 두 방법 모두 `scripts/scpa/gene_masking_lib.R`의 Phase 5 shared implementation 사용
- Primary ranking: `delta_g=S_full-S_minus_g`, `S=-log10(raw_p)`, signed descending

평가:

- Recall@3/5/10, Average Precision(primary), NDCG@3/5/10
- mean±SD overall/scenario/scenario×strength
- paired Wilcoxon AP와 matched-pairs rank-biserial effect
- strength pair는 같은 truth를 공유하므로 pathway×scenario×draw에서 평균 후 추론
- scenario별 세 p-value만 Bonferroni 보정
- null은 exact random-ranking chance와 비교하고 warning을 먼저 판정

### Runtime

- GPU와 LLM을 사용하지 않음
- 11 pathways × 20 draws × 7 scenario-strength 조건 = 1,540 experiments
- 예상 101,920 MCM calls
- CPU multi-core와 experiment atomic checkpoint/resume 지원
- 1-experiment smoke 뒤 full run을 사용자가 명시적으로 실행

### Phase 7 현재 상태

- Synthetic preparation PASS: 1,540 experiments / expected 101,920 MCM
- Real-data smoke PASS: 첫 42-gene null experiment, 86 MCM, 3.7분, warning 0
- Checkpoint input hash/protocol version 검증과 재사용 PASS
- Fallback audit: 560 experiments / 1,800 target-gene rows; 별도 집계 필수
- Full SCPA/evaluation은 아직 실행하지 않음
- LLM backend를 다시 추가하지 않음
- Phase 0–6 산출물을 수정하거나 재실행하지 않음

주요 문서/코드:

- `config/phase7_llmfree_synthetic.yaml`
- `docs/phase7_llmfree_protocol.md`
- `scripts/phase7/`
- `src/gene_embedding_project/genept_scpa/phase7/synthetic_benchmark_llmfree/`
- `tests/test_phase7_llmfree.py`
- `tests/test_phase7_scpa_masking.R`

---

## 4. L2를 어디에 사용했는지 정확한 정리

| 위치 | L2 사용 여부 | 의미/현재 처리 |
| --- | --- | --- |
| Phase 2 whole-cell GenePT-w | 사용 | Published method의 final row-wise unit L2. 유지 필수 |
| Phase 3 GenePT-w | 사용된 Phase 2 representation 재사용 | 추가 L2를 하지 않음 |
| Phase 4A/4B pathway projection primary | 미사용 | Primary는 non-L2 `X_P @ E_P` |
| Phase 4A/4B sensitivity | 별도 사용 | 계산 결과/CSV는 보존, 발표 figure에서는 제외 |
| Phase 5 gene masking | 미사용 | GenePT non-L2 subtraction masking |
| Phase 6 controls | 미사용 | QC의 `no_l2=true` |
| Phase 7 planned GenePT masking | 미사용 | non-L2 primary로 frozen |

---

## 5. 실행 중 발견하고 수정한 주요 문제

1. **Legacy Seurat `images` slot error**
   - serialized Seurat 3.1.5 object와 현재 SeuratObject compatibility 문제였다.
   - source RDS를 수정하지 않고 메모리에서만 update하여 해결했다.

2. **Phase 4 `numpy.int64 is not JSON serializable`**
   - QC JSON 작성 시 NumPy scalar가 들어간 문제였다.
   - native Python scalar로 변환하도록 고쳤다.

3. **장시간 작업에 progress가 보이지 않음**
   - Phase 4/6에 진행률, task/MCM count, elapsed, ETA, active workers, 현재 gene/task를
     출력하도록 추가했다.

4. **Phase 6 후반 active workers가 12가 아니라 2로 감소**
   - worker 설정 오류가 아니라 240 tasks 중 대부분이 끝나고 긴 마지막 2 tasks만 남은
     tail imbalance였다.
   - 중단 후 checkpoint/resume가 가능하며 완료 checkpoint를 다시 계산하지 않는다.

5. **Phase 6 robustness 후처리 crash**
   - R에서 나온 lowercase representation label `true`를 Python parser가 boolean True로
     오인해 true rows가 0개가 된 문제였다.
   - parser를 수정하고 이미 완료된 900 robustness checkpoints를 재사용해 aggregation했다.
   - 900 MCM을 다시 계산하지 않았다.

6. **과거 Phase runner 실행 시 `Phase 4 must remain active`**
   - config active phase가 다음 단계로 넘어간 뒤 옛 runner를 다시 실행했기 때문이다.
   - 현재 active phase는 7이므로 Phase 1–6 production runner를 단순 재실행하지 않는다.
   - figure만 다시 만들 때는 아래 plotting-only 명령을 사용한다.

7. **GPU 관련 판단**
   - SCPA/multicross MCM과 현재 R masking loops는 CPU 작업이다.
   - GPU를 붙인다고 실질적으로 빨라지지 않는다.
   - Phase 6의 `--cores 12`는 대략 12 CPU workers를 의미한다.
   - 현재 LLM-free Phase 7도 CPU-only이며 GPU를 사용하지 않는다.

---

## 6. Figure와 발표 상태

2026-08-18에 Phase 1–6 기존 scientific figures를 시각적으로 재검수했다.

- Phase 1 composite의 과도한 빈 공간, subtitle 잘림, heatmap label 크기를 수정했다.
- Phase 4의 긴 pathway 제목/축 label과 figure ratio를 수정했다.
- Phase 4B 발표 figure에서 L2 bar/legend/text를 제거했다.
- Phase 4B detection heatmap은 123개 전체를 작게 표시하지 않고, 적어도 한 comparison에서
  discordant한 30 pathways만 표시하며 title에 selection을 명시했다.
- Phase 5 subplot 제목과 heatmap label을 wrap했다.
- Phase 5 detection-flip figure를 horizontal grouped bar로 바꿨다.
- Phase 6 control heatmap x축을 Permuted/Random으로 명시했다.
- Phase 6 gene overlap figure를 읽기 쉬운 horizontal layout으로 바꿨다.
- Phase 2/3에는 별도 dedicated scientific-result figure artifact가 없다.

Figure-only 재생성:

```bash
Rscript scripts/scpa/replot_phase1b_figures.R
PYTHONPATH=src:. python scripts/replot_phase1_to_phase6_figures.py
```

이 명령은 기존 CSV/JSON만 읽으며 SCPA/MCM을 다시 실행하지 않는다.

중간보고에서 우선 고려할 figure:

1. Phase 1 official reference/Phase 1B comparison
2. Phase 4 Vanilla vs GenePT agreement/discordance
3. Phase 5 gene-rank or detection-flip comparison
4. Phase 6 semantic control figure 1–2개

Phase 6 결론을 가장 직접 보여주는 후보:

- `01_true_vs_control_pathway_scores.png`
- `06_resampling_robustness.png`
- 필요하면 `03_control_difference_heatmap.png`

negative Phase 6 결과를 숨기지 말고 “실패”가 아닌 frozen hypothesis test 결과로 말한다.

---

## 7. 재현성, Git, 대용량 파일 경계

현재 portable handoff commit/tag:

```text
273491e Make repository portable and preserve project state through Phase 7
tag: phase7-handoff-20260814
```

새 시스템에 별도 전달해야 하는 필수 대용량 파일은 두 개다.

1. `data/interim/genept_scpa/phase2_export/naive_cd4/naive_cd4_rna_counts_genes_by_cells.mtx`
   - 371,775,504 bytes
   - SHA-256 `6ea1626a0610d701fd23ae330ae384cfcbc90d013a458568f24d686d50ad9e88`
2. `data/reference/genept_scpa/genept_ada002/GenePT_gene_embedding_ada_text.pickle`
   - 460,797,248 bytes
   - SHA-256 `fd297510ddd3040744033fde0b0f2cf15a40ac8b2fd2fb02f10667295e55c862`

두 파일의 원본 SHA-256은 Phase 2 export manifest와 Phase 2 embedding provenance/QC에
기록돼 있으며 Phase 7 preparation이 실행 전에 다시 검증한다.

대용량 raw RDS, checkpoints, prepared HDF5, generated NumPy matrices, caches와 virtual
environment는 Git에 넣지 않는다. `.env`, token, API key도 절대 commit하지 않는다.
`git add .`, `git add -A`, broad `git commit -a`를 사용하지 않는다.

현재 작업 트리는 2026-08-18 figure 개선 때문에 dirty 상태다. 주요 변경은 plotting code,
Phase 1/4/5/6 rendered figures와 다음 plotting-only scripts다.

- `scripts/scpa/replot_phase1b_figures.R`
- `scripts/replot_phase1_to_phase6_figures.py`

이 변경은 아직 이 문서 작성 시점 HEAD `273491e` 이후 별도 commit되지 않았다.

---

## 8. 현재 검증 상태와 안전한 명령

최근 figure 수정 후 실행한 전체 Python test:

```bash
PYTHONPATH=src:. MPLCONFIGDIR=/tmp/genept_scpa_plot_cache/matplotlib pytest -q
```

기록 당시 결과: **99 passed**. Phase 7 교체 뒤 현재 test 결과는 새 검증 기록을 따른다.

기본 repository test:

```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test_*.py'
Rscript tests/test_phase7_scpa_masking.R
git diff --check
```

Phase 7 cheap tests:

```bash
PYTHONPATH=src pytest -q tests/test_phase7_llmfree.py
Rscript tests/test_phase7_scpa_masking.R
```

Phase 7 실행 순서:

```bash
PYTHONPATH=src python scripts/phase7/prepare_llmfree_benchmark.py
PYTHONPATH=src python scripts/phase7/run_llmfree_masking.py --max-experiments 1 --cores 1
PYTHONPATH=src python scripts/phase7/run_llmfree_masking.py --cores 12
PYTHONPATH=src python scripts/phase7/evaluate_llmfree_benchmark.py
```

### 현재 주의할 명령/행동

- 101,920 MCM full run을 smoke 없이 시작
- `--allow-partial` 결과를 production scientific result로 보고
- LLM backend를 Phase 7에 다시 추가
- Phase 1–6 production runner 재실행
- Phase 4 L2를 다시 primary result처럼 발표

---

## 9. 새 GPT가 반드시 지켜야 할 해석 원칙

1. Phase 6은 이미 끝났다. 다시 시작하거나 Phase 6A/6B를 미완료라고 말하지 않는다.
2. Technical PASS와 scientific support는 다르다.
3. Phase 6 technical gate는 PASS지만 semantic-specific superiority는 지지되지 않았다.
4. Random과 다르다는 사실만으로 semantic correctness를 주장하지 않는다.
5. GenePT가 sensitivity ordering을 바꾼다는 사실만으로 biological accuracy를 주장하지 않는다.
6. CD4 activation은 positive-control-like benchmark이지 완전한 pathway ground truth가 아니다.
7. Real-data reference에 없는 pathway를 즉시 false positive로 부르지 않는다.
8. SCPA qval이 클수록 stronger difference다.
9. exact same cells/pathways/paired genes/order를 유지한다.
10. source RDS와 official embedding original을 수정하지 않는다.
11. 새 분석 전 metric, threshold, seed, control, stopping rule을 먼저 고정한다.
12. Phase 7 null evaluation-target은 실제 perturbed gene이라고 부르지 않는다.
13. Phase 7 partial smoke 결과를 scientific comparison으로 해석하지 않는다.
14. 기존 Phase 0–6 결과를 덮어쓰지 않는다.

---

## 10. 프로젝트의 현재 정직한 결론

- Vanilla SCPA reproduction은 성공했다.
- Published GenePT-w 구현은 기술적으로 정확하게 재현됐다.
- Whole-cell GenePT-w는 CD4/CD8 multivariate difference를 보존했다.
- Pathway-specific GenePT projection은 Vanilla와 전반적 agreement를 보이면서 일부
  pathway ranking/detection을 바꿨다.
- Gene masking sensitivity ordering도 Vanilla와 GenePT에서 달라졌다.
- 그러나 correct GenePT correspondence는 within-pathway permutation control과 명확히
  구분되지 않았다.
- 따라서 현재 자료는 GenePT의 semantic-specific superiority를 지지하지 않는다.
- 연구가 망한 것이 아니라, “표현이 결과를 바꾸는 것”과 “올바른 semantics 때문에 더
  정확한 것”을 분리해 낸 결과다.
- 다음 타당한 질문은 synthetic ground truth에서 실제 perturbed-gene recovery accuracy를
  비교하는 Phase 7이다.

---

## 11. 새 GPT 채팅에 보낼 추천 시작 프롬프트

아래 문장을 이 파일과 함께 보낸다.

```text
프로젝트 루트는 /home/node00/nhy_python/GenePT_SCPA 이다.

먼저 docs/genept_scpa_full_handoff_for_new_chat.md를 처음부터 끝까지 읽고,
genept_scpa_experiment_plan.md, docs/genept_scpa_decision_log.md,
config/genept_scpa.yaml, config/phase7_llmfree_synthetic.yaml과 대조해 현재 상태를
파악해라.

중요:
- Phase 0–6은 완료됐고 Phase 6은 최종 PASS/COMPLETED다.
- Phase 6의 frozen scientific conclusion을 바꾸지 마라.
- 현재 active phase는 Phase 7이다.
- Phase 7의 기존 LLM 설계는 폐기됐고, 현재는 CPU-only Vanilla vs GenePT synthetic
  recovery protocol/code가 구현됐다. Full run은 아직 하지 않았다.
- 기존 Phase 0–6 결과를 덮어쓰지 마라.
- 기존 자료만으로 GenePT superiority, biological correctness 또는 causality를
  주장하지 마라.
- 작업 전에 git status와 관련 문서/코드를 확인하고, 현재 사용자 변경을 보존해라.

문서를 읽은 뒤 먼저 현재 상태, 완료된 것, 잠긴 것, 다음 decision gate를 짧게
요약하고 내가 요청한 범위만 진행해라.
```

---

## 12. 원본 authoritative 문서 우선순위

이 인수인계와 원본이 충돌하면 다음 순서로 확인한다.

1. `docs/genept_scpa_decision_log.md`의 가장 최신 Accepted decision
2. `config/genept_scpa.yaml`과 `config/phase7_llmfree_synthetic.yaml`의 frozen design
3. `genept_scpa_experiment_plan.md`
4. Phase별 QC JSON과 summary MD
5. 이 인수인계 문서

이 문서는 새 채팅 context 복원을 위한 요약이며 historical raw artifacts 자체를
대체하지 않는다.
