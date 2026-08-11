# GenePT × SCPA Project Experimental Plan

목적:

- 연구 질문, 단계별 목적, 구현 범위, 평가 지표와 중단 조건을 일관되게 관리한다.
- 각 Phase는 이전 Phase의 PASS와 필요한 methodological decision 이후에만 진행한다.
- 결과를 보기 전에 주요 분석 설계를 고정한다.

프로젝트 루트: `/home/node00/nhy_python/GenePT_SCPA`

## Current project hypothesis

Primary hypothesis to be tested:

> GenePT-w의 gene-semantic embedding을 이용한 cell representation이 original
> gene-expression space의 biologically meaningful cell-population differences를
> 보존하거나 더 명확하게 표현할 수 있는지 검증한다.

이는 검증할 가설이며 GenePT-w의 우수성을 전제하지 않는다.

- Initial feasibility benchmark: GSE212270 naïve CD4 0h vs naïve CD8 0h
- Current primary validation benchmark: naïve CD4 activation 0h-vs-12h,
  12h-vs-24h, 0h-vs-24h
- Current primary evaluation: lineage를 CD4로 고정하고 동일 cells/pathways/paired
  genes에서 Vanilla pathway SCPA와 GenePT-informed pathway projection + SCPA-core의
  detection state, qval floor와 adjusted-p signal을 비교
- Downstream question: GenePT semantic information이 pathway 내부 gene의
  differential contribution을 더 잘 설명하는지는 별도 ground truth/control 아래 검증

## Current project status

| Item | Status |
| --- | --- |
| Phase 0 — project reset / protocol freeze | **PASS** |
| Phase 1A — CD4 dataset acquisition / validation | **PASS** |
| Phase 1B — Vanilla SCPA reproduction | **PASS** |
| Phase 2 — Published GenePT-w reproduction | **PASS** |
| Phase 3 — Whole-cell GenePT-w CD4/CD8 SCPA-core feasibility | **PASS** |
| Phase 4A — Exploratory CD4 0h vs CD8 0h lineage comparison | **COMPLETED / HISTORICAL PRESERVED** |
| Phase 4B — Naive CD4 activation primary validation (3 comparisons) | **PASS / COMPLETED** |
| Phase 4C — Optional CD8 generalization | **NOT SCHEDULED / NOT RUN** |
| Phase 5 — Pathway-internal gene masking sensitivity | **IN PROGRESS / METHOD FROZEN** |
| Phase 6 이후 | **NOT STARTED** |
| Naïve CD8 acquisition / QC | **PASS** |

문서 변경만으로 어떤 Phase도 새로 PASS 처리하지 않는다.

# 1. Primary research question

> Single-cell gene-expression profiles를 GenePT-w의 1536-dimensional cell
> representation으로 변환한 뒤에도 biologically distinct cell populations을
> 구별하는 정보가 보존되는가? 그리고 SCPA가 사용하는 multivariate comparison
> framework가 GenePT-w 공간에서 이러한 population difference를 검출할 수 있는가?

GenePT-w는 population당 하나의 평균 vector를 만드는 방법이 아니다.

```text
Cell × Gene expression + Gene × 1536 GenePT embedding
→ Cell × 1536 GenePT-w representation
```

CD4와 CD8의 각 individual cell이 하나의 1536D vector를 가진다.

구분할 질문:

```text
A. GenePT-w 이후에도 CD4/CD8 known distinction이 보존되는가?
B. SCPA-core multivariate framework가 그 population difference를 검출하는가?
C. True GenePT가 Permuted/Random controls보다 나아 semantic mapping의
   추가 효과를 지지하는가?
```

A나 B만으로 semantic improvement를 주장하지 않는다.

Supervisor clarification에 따라 프로젝트는 cell-type separability에 한정되지 않는다.
또 하나의 core downstream question은 GenePT semantic information을 이용했을 때
pathway 내부 gene들의 differential contribution을 기존 방식보다 더 정확하게
표현·해석할 수 있는가이다. 여기서 "더 정확하게"의 metric/algorithm은 아직
확정하지 않으며 Phase 6 methodological decision gate에서 ground truth와 평가 기준을
먼저 고정한다.

# 2. Established methods and project-specific adaptation

## 2.1 SCPA

- pathway gene들의 single-cell multivariate expression distribution 변화를 평가
- 공식 `compare_pathways()`는 genes × cells expression matrices와 pathway gene
  sets를 입력으로 받음
- `qval`은 primary statistic, FC는 secondary statistic
- 공식/default reproduction: 최대 500 cells/group, matched genes 15–500

Reference: Bibby JA et al. *Systematic single-cell pathway analysis to
characterize early T cell activation.* Cell Reports. 2022;41(8):111697.
DOI: 10.1016/j.celrep.2022.111697

## 2.2 GenePT and GenePT-w

- NCBI Gene summary를 text-embedding-ada-002로 embedding
- published gene representation dimension = 1,536
- normalized expression을 weight로 gene embeddings를 cell별 aggregate하고 final
  vector를 L2-normalize
- GenePT 논문은 embedding cosine-similarity graph 기반 gene programs도 분석

Reference: Chen Y, Zou J. *Simple and effective embedding model for single-cell
biology built from ChatGPT.* Nature Biomedical Engineering. 2025;9:483–493.
DOI: 10.1038/s41551-024-01284-6

## 2.3 Project-specific adaptation

GenePT-w의 1536 dimensions는 genes나 curated pathways가 아니다. 이를 일반
`compare_pathways()`에 pathway genes처럼 직접 넣는 것은 original SCPA pathway
analysis가 아니다. Phase 3에서 검토할 방법은 다음처럼 구분한다.

```text
adaptation of the SCPA multivariate comparison framework
to GenePT-w cell representations
```

입력 타당성, dimension effect, null behavior와 해석 범위는 구현 전 별도 gate에서
고정한다.

# 3. Dataset roles

## Dataset A — GSE212270 integrated naïve CD4

역할:

1. Phase 1 original SCPA reproduction
2. Phase 2 GenePT-w reproduction/QC
3. Phase 3 main experiment의 CD4 side

확보 파일: `GSE212270_integrated_naive_cd4.rds.gz`

Phase 1A 검증 결과:

- 14,894 cells, 17,856 features
- time metadata: `Hour`
- 0 h: 4,428; 12 h: 4,547; 24 h: 5,919 cells
- assays: RNA, integrated; active assay: RNA
- RNA layers: counts, data, scale.data
- gene identifiers: gene-symbol-like; duplicate/missing: 0/0
- legacy object는 메모리에서만 호환 갱신, source RDS 수정 없음
- Phase 1A: PASS

## Dataset B — GSE212270 integrated naïve CD8

역할:

- main CD4-vs-CD8 experiment를 위해 이후 acquisition/validation
- Dataset A와 동일한 QC protocol 적용
- CD8 0 h cells가 primary comparison의 CD8 side

현재 상태: **ACQUIRED / QC PASS**. 7,041 cells, 17,942 features이며 RNA counts가
존재한다. 이 사실은 향후 Phase 3 준비 상태일 뿐, Phase 2에서는 CD8 projection이나
CD4-vs-CD8 비교를 실행하지 않는다.

## Same-family rationale

CD4와 CD8 모두 동일 GSE212270 family를 사용해 서로 다른 연구·플랫폼·처리
pipeline에서 오는 confounding을 줄인다. 동일 accession family가 모든 batch effect를
제거한다고 가정하지는 않는다.

# 4. Reference labels and biological ground truth

```text
CD4/CD8 annotation
→ known/reference class label
→ classification/separability benchmark에 사용 가능

Biological pathway relevance
→ reference/expected pathway
→ 완전한 pathway ground truth로 가정하지 않음
```

Real data에서 reference에 없는 pathway가 상위에 왔다고 즉시 false positive라고
부르지 않는다.

# 5. Phase overview

```text
Phase 0  Project reset / protocol freeze                         PASS
Phase 1  Original SCPA reproduction
  1A     Dataset acquisition / validation                       PASS
  1B     Vanilla SCPA global + pairwise + reference             PASS
Phase 2  Published GenePT-w reproduction and QC                  PASS
Phase 3  Whole-cell GenePT-w CD4/CD8 SCPA-core feasibility       PASS
Phase 4A Exploratory lineage pathway comparison                  COMPLETED / HISTORICAL PRESERVED
Phase 4B Naive CD4 activation primary validation                 PASS / COMPLETED
Phase 4C Optional CD8 0h-vs-24h generalization                   NOT SCHEDULED / NOT RUN
Phase 5  Pathway-internal gene masking sensitivity               IN PROGRESS / METHOD FROZEN
Phase 6  Semantic-specific controls and robustness               NOT STARTED
Phase 7  External/biological/synthetic reference validation      NOT STARTED
Phase 8  Optional classifier/separability analysis               NOT STARTED
Phase 9  Additional timepoint analyses                            NOT STARTED
Phase 10 Final interpretation / presentation / report            NOT STARTED
```

# 6. Common principles and methodological cautions

1. 결과를 본 뒤 metric, classifier, threshold 또는 control을 바꾸지 않는다.
2. 비교 시 same cells, labels, preprocessing/splits와 paired seeds를 사용한다.
3. dataset/version/seed/package/gene coverage와 sampling IDs를 기록한다.
4. GenePT-w 1536 dimensions는 genes가 아니다.
5. GenePT-w에 일반 `compare_pathways()`를 직접 적용하는 것은 original SCPA
   pathway analysis와 동일하지 않다.
6. SCPA-core statistic을 GenePT-w에 쓰면 adapted method이며 별도 validation이
   필요하다.
7. 낮은 p-value는 높은 classification accuracy와 동일하지 않다.
8. Original expression과 GenePT-w의 dimension이 다르므로 p-value/qval 크기만으로
   어느 representation이 더 좋다고 주장하지 않는다.
9. Multivariate statistic, separability metrics, negative controls와 robustness를
   함께 해석한다.
10. CD4/CD8 comparison은 time point를 먼저 고정한다. Primary는 CD4 0 h vs
    CD8 0 h다.
11. Gene mapping coverage가 group과 연관되면 artificial separation을 만들 수 있어
    group별 coverage를 확인한다.
12. Source RDS와 official embedding 원본은 수정하지 않는다.

---

# Phase 0. Project reset and protocol freeze

Status: **PASS**

기존 pathway-recovery 프로젝트와 논리적으로 분리하고 package/script/config/docs,
phase gate와 decision log를 만들었다. 기존 LOO/frequency/hybrid 코드는 삭제하지
않았다. 초기 tests/static checks가 통과했으며 Phase 0 중 full analysis는 실행하지
않았다. 이 완료 상태를 되돌리지 않는다.

---

# Phase 1. Original SCPA reproduction

## 목적

GenePT 없이 공식 SCPA workflow가 현재 환경과 full naïve CD4 object에서 정상
재현되는지 확인한다.

## Phase 1A. Dataset acquisition / validation

Status: **PASS**

GEO GSE212270 full integrated naïve CD4 object를 확보했고 Dataset A에 기록한 QC를
통과했다.

## Phase 1B. Vanilla SCPA reproduction

Status: **PASS**

`Hour`별 전체 cells로 네 분석을 수행한다.

```text
Global:    0 h + 12 h + 24 h
Pairwise:  0 h vs 12 h
Pairwise: 12 h vs 24 h
Pairwise:  0 h vs 24 h
```

Workflow와 고정 조건:

```text
existing Seurat RNA/data log1p-normalized expression
→ SCPA::seurat_extract()
→ SCPA::compare_pathways()
```

- SCPA v1.6.2 combined metabolic pathways 243개
- seed 20260810; downsample 500/group; matched genes 15–500
- 동일 Hour별 sampled cells를 네 분석에 재사용
- new normalization/ScaleData/DEG filtering 없음
- official default single-core, GPU 없음

Hour-only 설계는 tutorial의 Resting@0 vs Activated@24 및 paper pseudotime
milestones와 다르므로 qualitative comparison만 한다.

예상 결과:

```text
01_global_0_12_24.csv
02_pairwise_0_vs_12.csv
03_pairwise_12_vs_24.csv
04_pairwise_0_vs_24.csv
phase1b_scpa_qc.json
phase1b_reproduction_summary.md
sampling IDs and paper/tutorial comparison figures
```

Hour-only 네 분석과 별도 official reference(`Resting 0 h` vs `Activated 24 h`)가
모두 완료됐다. Reference는 4,341/1,697 full cells에서 각각 500 cells를 사용했고,
124 pathways의 qval/FC가 모두 finite였다. Reactome/KEGG arachidonic-acid pathways가
각각 rank 56/69, qval 6.6355/5.8762, FC -0.8939/-0.2206으로 확인되어 작은 mean
shift와 multivariate difference가 공존하는 qualitative example을 재현했다.
고정된 9개 gate criterion이 모두 true이고 warning/error 및 parameter tuning이
없으므로 Phase 1B는 PASS다. Phase 1에서 GenePT는 사용하지 않았다.

---

# Phase 2. Published GenePT-w reproduction and QC

Status: **PASS**

SCPA와 결합하기 전에 published GenePT embedding과 cell-wise GenePT-w pipeline을
독립적으로 재현한다. 먼저 validated naïve CD4 object에서 검증한다.

## Embedding and mapping QC

- 저자 공개 embedding 우선; 원 논문 비교 시 새 model로 재생성하지 않음
- source/version/checksum, gene count, dimension=1536, dtype, NaN/Inf
- gene identifiers/aliases/duplicates
- exact/alias matches, unmatched genes와 group별 coverage
- mapping rule과 gene order를 결과 전에 고정

## Published preprocessing gate

최종 peer-reviewed paper, 저자 methods와 official aorta notebook을 확인해 다음을
고정했다.

- source expression: sparse `RNA/counts`
- all dataset genes를 포함한 cell-wise 10,000 transcript normalization
- 이어서 element-wise `log1p`
- official artifact key에 exact match; artifact가 포함하는 HGNC alias key만 허용
- unmatched gene은 library-size normalization에는 포함하고 projection 기여는 zero
- official notebook과 같이 weighted sum을 전체 dataset gene 수로 나눔
- final row-wise unit L2 normalization

Official aorta notebook은 이미 전처리된 `.X`를 읽으므로 upstream normalization을
코드에서 반복하지 않지만, 저자 methods는 counts → 10,000 normalization → log1p와
final L2를 명시한다. 이번 pipeline은 이를 raw Seurat counts에서 명시적으로 재현한다.

```text
X: cells × matched genes
E: matched genes × 1536
→ gene-aligned expression-weighted aggregation
→ cell × 1536 vectors
→ row-wise L2 normalization
```

필수 QC/tests:

```text
matched/unmatched counts; input/output dimensions; finite values
weighted aggregation correctness; vector norm distribution
zero-expression handling; deterministic output; gene-order invariance
reproducibility metadata
```

Phase 2에서는 SCPA와 결합하지 않는다.

상세 protocol과 실행 절차는 `docs/phase2_genept_w_protocol.md`에 고정한다.

실제 full naïve CD4 결과는 14,894 × 1,536이며 17,856 dataset genes 중 14,409가
official lookup과 매칭됐다(exact 14,087; official alias 322). Median raw-count mass
coverage는 0.931967, zero vector는 0이고 모든 값이 finite이며 post-L2 norm은 1에
수치적으로 일치했다. Synthetic, optimized/direct, determinism과 source-object
read-only gate가 모두 통과했고 warning/failed check가 없으므로 Phase 2를 PASS로
확정한다. 기존 QC의 historical gate value `READY_FOR_GPT_REVIEW`는 변경하지 않는다.

---

# Phase 3. Primary CD4 vs CD8 GenePT-w benchmark

Status: **PASS**

Primary comparison: **CD4 0 h vs CD8 0 h**.

목적은 GenePT-w representation에 CD4/CD8 biological distinction이 보존되는지,
그리고 SCPA의 multivariate comparison core가 두 population 차이를 검출하는지
평가하는 것이다. 0 h를 먼저 사용해 activation/time effect를 최소화한다. 이
Phase는 어느 representation이 더 우수한지 평가하지 않는다.

## 3.1 Naïve CD8 GenePT-w

CD8 acquisition/QC는 이미 PASS다(7,041 cells; 17,942 RNA features; 0/12/24 h =
1,048/2,066/3,927). Phase 2에서 검증한 동일 ada-002 artifact와 동일한
RNA/counts → total 10,000 → log1p → exact official mapping → weighted aggregation →
row-wise L2 pipeline을 재사용한다. 새 artifact를 다운로드하거나 CD4와 다른
preprocessing을 쓰지 않는다. CD4/CD8 median mapping coverage 차이를 potential
confounder로 기록하되 사후 exclusion threshold를 만들지 않는다.

## 3.2 Representations

동일한 CD4/CD8 0 h cells에서 비교한다.

```text
A. Original expression representation
B. Published-workflow GenePT-w representation (cells × 1536)
```

같은 cell IDs, labels와 paired sampling을 사용한다. 실제 metadata의 0 h count를
재확인한 뒤 SCPA default와 Phase 1 convention에 근거해 seed 20260810으로
500 cells/group을 고정한다. 모든 0 h IDs와 canonical sampled IDs를 모두 저장한다.

Original-expression reference는 각 dataset의 모든 RNA/counts genes로 cell total을
10,000에 맞추고 log1p한 다음 exact shared gene symbols로 정렬한다. 현재 QC의 shared
count 17,085를 gate로 확인한다. 이 policy는 gene position concatenate를 금지하며,
GenePT projection이 반드시 같은 17,085 genes만 사용해야 한다는 뜻은 아니다.

## 3.3 SCPA-core methodological decision gate

다음을 명시적으로 검토한다.

```text
CD4 GenePT-w cells × 1536 vs CD8 GenePT-w cells × 1536
```

- SCPA 1.6.2 `single_comparison()`의 core는 `multicross::mcm()`이다.
- Official pathway workflow는 pathway subset을 transpose해 cells × genes로 전달한다.
- Adapter는 cells × features의 두 population을 같은 `multicross::mcm()`에 직접
  전달하며 GenePT dimensions를 genes/pathway로 위장하지 않는다.
- MCM 내부 distance는 Euclidean이고 adapter 내부 hidden downsampling은 없다.
- 한 global hypothesis의 p-value와 SCPA-style `sqrt(-log10(p))`를 기록한다.
- GenePT-w에는 추가 normalization을 적용하지 않는다. 이미 row-wise unit L2다.

이를 `SCPA-core multivariate framework adaptation`이라고 부르며 standard pathway
analysis라고 부르지 않는다. Nearly-identical vs shifted toy population에서 shifted
case의 statistic이 더 큰지를 full run 전에 검사한다. Original/GenePT의 차원과
geometry가 다르므로 raw MCM p/q 크기로 어느 representation이 더 좋다고 결론 내리지
않는다.

실제 결과는 canonical 500 CD4/500 CD8 cells에서 GenePT-w p=4.6605e-12,
original expression p=5.5516e-60이었고 두 공간 모두 multivariate difference를
검출했다. QC 14개 criterion이 모두 true이고 warning/failed check가 없어 Phase 3를
PASS로 닫았다. 이는 classifier accuracy나 GenePT 우수성을 의미하지 않는다.

---

# Current Phase 4–10 roadmap (frozen 2026-08-11)

## Phase 4. Pathway-specific Vanilla vs GenePT-informed SCPA comparison

Status: **IN PROGRESS — CD4 ACTIVATION VALIDATION READY FOR USER RUN**

### Phase 4A. Initial lineage comparison

Status: **COMPLETED, HISTORICAL OUTPUT PRESERVED**

Primary cohort는 Phase 3 canonical naïve CD4 0 h 500 cells와 naïve CD8 0 h 500
cells를 그대로 재사용한다. RNA/counts를 전체 transcriptome 기준 total 10,000으로
정규화하고 log1p한 뒤 pathway gene을 subset한다. Pathway 내부 재정규화는 하지 않는다.

Primary paired gene set은 다음과 같다.

```text
pathway genes ∩ CD4 genes ∩ CD8 genes ∩ official GenePT keys
```

Vanilla와 GenePT branch가 동일한 cells, pathway universe, paired genes와 gene order를
사용한다. Vanilla는 cells × paired genes를, GenePT-informed branch는
`Z_P = X_P × E_P`로 만든 cells × 1,536 representation을 `multicross::mcm()`으로
비교한다. Primary Z는 non-L2이며 rowwise L2는 사전 선언된 sensitivity option이다.
각 pathway의 embedding/projected effective rank를 기록한다.

Official SCPA 1.6.2 convention에 따라 동일 eligible pathway universe에서 raw p를
Bonferroni 보정하고 `qval = sqrt(-log10(adjusted_p))`를 계산한다. Spearman, Kendall,
Top-10/20 overlap/Jaccard와 `rank_delta = genept_rank - vanilla_rank`를 agreement 및
reordering metric으로 사용한다. 이를 accuracy 또는 superiority metric으로 해석하지
않는다. Full production run은 사용자가 수행하고 gate는 `READY_FOR_GPT_REVIEW`에서
멈춘다.

실제 Phase 4A full run은 123 pathways에서 완료됐으나 Bonferroni adjusted p=1,
qval=0 tie가 많았다. Historical raw p와 기존 산출물은 보존한다. 기존 unique rank는
qval=0 ties에 임의 순서를 부여하므로 validation extension의 primary reporting에는
사용하지 않는다.

### Phase 4B. Primary Naive CD4 activation validation

Status: **PASS / COMPLETED**

Lineage를 Naïve CD4로 고정하고 activation time만 바꾼 다음 세 comparison을 현재
primary production scope로 사용한다.

```text
CD4 0h vs 12h   early activation
CD4 12h vs 24h  later activation
CD4 0h vs 24h   full activation / positive-control-like benchmark
```

이미 frozen된 각 500 cell IDs, 동일 123 paired pathways와 Vanilla/non-L2/L2 branch를
사용한다. 0h-vs-24h 대표 pathway의 Vanilla raw p를 동일 입력의 official
`SCPA::compare_pathways()`와 직접 비교하고 tolerance 1e-12를 요구한다.

설치된 SCPA 1.6.2의 실제 공식 qval은
`sqrt(-log10(Bonferroni-adjusted p))`이므로 formula mismatch는 없다. Regenerated
reporting은 `ties.method="average"`를 사전 고정하고 adjusted p<0.05 기준으로 Both,
Vanilla-only, GenePT-only, Neither detection states를 기록한다.

Primary reporting은 rank보다 N, raw/adjusted-p significant counts, qval positive/zero,
medians와 qval floor fraction을 먼저 사용한다. Historical Phase 4A의 qval-zero와
significant counts를 비교하되 기존 파일을 덮어쓰지 않는다. 이 비교는 representation
차이를 기술하지만 superiority, accuracy 또는 generalization을 입증하지 않는다.

Production은 세 comparison의 123 pathways와 Vanilla/non-L2/L2에서 모두 완료됐다.
QC failed check와 runtime warning은 0이고 official raw-p cross-check, average-tie rank,
0h-vs-24h positive-control-like signal을 확인했으므로 Phase 4B를 PASS로 닫는다.

### Phase 4C. Optional CD8 generalization

Status: **NOT SCHEDULED / NOT RUN**

Phase 4B review 후 필요할 때 CD8 0h-vs-24h를 generalization 후보로 별도 사전 고정할
수 있다. 현재 CD8 activation과 12h/24h lineage comparison은 실행하지 않는다. 이전
9-comparison runner 기능은 `--comparison-set all_9` explicit option으로만 보존하며
default/primary production에 포함하지 않는다.

## Phase 5. Pathway-internal gene contribution analysis

Status: **IN PROGRESS — METHOD FROZEN**

Phase 4B에서 adjusted p<0.05 기준 Vanilla-only 또는 GenePT-only였던 모든
pathway-comparison pair를 target으로 고정한다. Expected target은 0h-vs-12h 11개,
12h-vs-24h 9개, 0h-vs-24h 10개의 총 30개다. 결과를 보고 target을 바꾸지 않는다.

같은 frozen cells, preprocessing, 123-pathway paired genes/order를 재사용하고 두
branch에서 같은 gene을 masking한다. Vanilla는 `X_P[:,g]=0`, GenePT non-L2는
`Z_P - outer(X_P[:,g], E_P[g,:])`를 사용한다. Primary sensitivity는
`delta[-log10(raw p)] = score_full - score_masked`이며 raw p의 score 변환에만
1e-300 clipping을 사용한다. Signed/absolute delta를 average-tie rank로 보고한다.

이는 gene importance 또는 biological correctness가 아니라 gene masking sensitivity다.
GenePT L2 gene-level LOO, Phase 6 True/Permuted/Random, CD8 generalization과 classifier는
실행하지 않는다. Pathway-comparison 단위 atomic checkpoint/resume와 baseline Phase 4B
raw-p reproduction, masking equivalence, deterministic ranking을 production gate로 둔다.

## Phase 6. Semantic-specific controls

Status: **NOT STARTED**

True GenePT, gene↔embedding assignment만 바꾼 Permuted GenePT, dimension-matched
Random embedding을 비교한다. Repeated sampling, seed/sample-size robustness와
coverage/dimension controls를 함께 사전 고정한다.

## Phase 7. External, biological and synthetic reference validation

Status: **NOT STARTED**

Known biological reference 또는 injected synthetic perturbation ground truth와 평가
기준을 결과 전에 고정한다. `more accurate`라는 결론은 이 Phase의 external criterion을
통과한 경우에만 검토한다.

## Phase 8. Optional classifier/separability analysis

Status: **NOT STARTED / OPTIONAL**

Known CD4/CD8 labels를 이용한 classifier/CV는 pathway comparison과 분리한다.
Classifier, folds, hyperparameters와 primary metric을 구현 전에 고정한다.

## Phase 9. Additional time-point analyses beyond Phase 4C

Status: **NOT STARTED**

Phase 4B/4C 이외에 추가 timepoint task가 필요한 경우에만 연다.
Cell-type comparison과 within-cell-type activation comparison을 하나의 task로 섞지
않는다.

## Phase 10. Final interpretation, presentation and report

Status: **NOT STARTED**

Phase 1–9의 reproduction, feasibility, pathway reordering, contribution, semantic
controls, reference validation과 robustness를 각 주장 수준에 맞춰 통합한다.

---

# Historical pre-Phase-4 roadmap (retained, superseded 2026-08-11)

아래 Phase 4–9 내용은 2026-08-11 이전 계획 history 보존용이며 현재 실행 순서에는
사용하지 않는다.

## Former Phase 4. Separability metrics and negative controls

Status: **NOT STARTED**

Known CD4/CD8 class labels로 representation의 cell-type information을 평가한다.

Representations:

1. Original expression
2. True GenePT-w
3. Permuted GenePT-w
4. Random dimension-matched embedding projection

Permuted GenePT는 embedding vector 집합을 유지하고 gene↔embedding mapping만
shuffle한다. Random baseline은 GenePT와 dimension을 맞추고 distribution/scale/norm
규칙을 미리 고정한다.

Classifier 후보는 logistic regression 또는 k-nearest-neighbour다. 최종 classifier,
hyperparameters, preprocessing과 split/CV는 구현 전에 고정하며 모든 representation에
동일 folds/cells를 사용한다. 후보 metrics는 AUROC, Accuracy, F1, balanced accuracy다.
Primary metric 역시 결과 전에 고정한다.

SCPA-core statistic은 classification 성능이 아니며 supervised score는 known-label
separability를 측정할 뿐 biological mechanism을 직접 증명하지 않는다.

---

## Former Phase 5. Time-course extension

Status: **NOT STARTED**

Primary 0 h benchmark 검토 후에만 확장한다.

Question A — cell-type separability:

```text
0 h: CD4 vs CD8
12 h: CD4 vs CD8
24 h: CD4 vs CD8
```

Question B — activation-state separability:

```text
within CD4: 0 / 12 / 24 h
within CD8: 0 / 12 / 24 h
```

두 질문을 하나의 task로 섞지 않는다.

---

## Former Phase 6. Core pathway/gene-level GenePT-SCPA interpretation

Status: **NOT STARTED**

Supervisor clarification에 따라 이 Phase는 optional extension이 아니라 핵심 연구
질문이다.

> GenePT semantic information을 이용했을 때 pathway 내부 gene의 differential
> contribution을 기존 방식보다 더 잘 설명할 수 있는가?

비교 축은 Vanilla SCPA pathway analysis, GenePT-informed pathway representation,
gene-level contribution/sensitivity, 그리고 SCPA 사용/비사용 시 pathway 내부 gene
ranking 차이다. 단 "더 정확하게"는 다른 ranking이 나왔다는 사실로 정의하지 않는다.
별도 methodological decision gate에서 contribution 정의, ground truth/reference와
평가 규칙을 먼저 고정한다.

표준으로 확정하지 않은 후보:

- leave-one-gene-out sensitivity
- gene perturbation 또는 pathway-score sensitivity
- GenePT semantic similarity-based grouping
- known markers/reference genes와 curated reference pathways
- synthetic ground-truth pathway benchmark

기존 GenePT-derived gene programs, pathway-specific GenePT-w와 synthetic pathway
benchmark 설계는 이 Phase의 후보로 보존한다.

## 6A. GenePT-derived gene programs → original SCPA

```text
GenePT embeddings → cosine similarity → gene graph
→ Leiden/Louvain programs → gene sets → original expression-based SCPA
```

## 6B. Pathway-specific GenePT-w projection

```text
X_P: cells × pathway genes
E_P: pathway genes × 1536
Z_P = X_P × E_P followed by frozen GenePT-w normalization
```

## 6C. Synthetic and real-reference evaluation candidates

기존 null/mean-shift/partial-gene/cell-subset/mixed-direction synthetic scenarios와
real pathway reference 설계를 삭제하지 않는다. Metric, replicate 수, reference
grading은 결과 전에 별도로 고정한다.

---

## Former Phase 7. Robustness and sensitivity

Status: **NOT STARTED**

- repeated matched cell subsampling; 모든 method에서 동일 sampled IDs
- sample-size/seed sensitivity
- label permutation과 rank stability
- group별 GenePT gene-coverage confounding
- feature/dimension and random-projection sensitivity

실제 sensitivity sample sizes와 반복 수, CI 방법은 결과 전에 별도로 고정한다. Phase 1의
official 500-cell reproduction과 Phase 7 sensitivity를 구분한다. Gene mapping
permutation과 label permutation도 서로 다른 null이다.

---

## Former Phase 8. Optional external embedding baseline

Status: **NOT STARTED / SECONDARY**

Primary experiment 후 첫 후보로 Gene2vec를 고려한다.

```text
literature-derived embedding vs expression/co-expression-derived embedding
```

처음부터 여러 external model을 모두 구현하지 않는다. Coverage/dimension/
preprocessing을 별도 QC하고 동일 cells/splits를 사용한다.

---

## Former Phase 9. Final interpretation, figures, and report

Status: **NOT STARTED**

최소 질문:

1. Original SCPA reproduction이 성공했는가?
2. Published GenePT-w를 재현했는가?
3. GenePT-w가 naïve CD4/CD8 distinction을 보존하는가?
4. Original expression 대비 separability는 어떠한가?
5. True GenePT가 Permuted/Random보다 나은가?
6. 결과가 0/12/24 h에서 안정적인가?
7. sampling/coverage/dimension artifact 가능성은 없는가?
8. Optional pathway extension이 추가 insight를 제공하는가?

해석 예:

```text
GenePT-w >= Original and True > controls
→ semantic representation이 distinction을 보존/강화할 가능성

True ~= Permuted/Random
→ semantic correspondence보다 projection effect 가능성

GenePT-w < Original
→ transformation이 cell-type signal을 약화할 가능성
```

모두 정당한 결과다. GenePT의 보편적 우수성, p-value만으로 classification 성능,
완전한 pathway truth 또는 Hour-only Figure 4 수치 복제를 주장하지 않는다.

# 8. Expected outputs and testing policy

Phase별 interim QC와 processed output을 `data/interim/genept_scpa/` 및
`data/processed/genept_scpa/phaseN/`에 저장한다. Large raw data와 generated
embeddings는 Git에 commit하지 않는다. 이후 Phase의 exact schema/filename은 해당
Phase protocol을 고정할 때 결정한다.

Tests:

- Phase 1: three/two-sample SCPA, Hour filtering, schema/finite values, seed/QC/figures
- Phase 2: embedding provenance/dimension, mapping, weighted aggregation, L2 norm,
  gene-order consistency
- Phase 3: same cells, matrix orientation, SCPA-core null/calibration, gate enforcement
- Phase 4: exact canonical cells, paired gene identity/order, hand-computed projection,
  no pathway renormalization, official qval correction, effective rank, rank metrics
- Phase 5: matched gene removal/contribution-method wiring and regeneration manifest
- Phase 6: True/Permuted/Random mapping, paired sampling, seed determinism, rank stability
- Phase 7: injected pathway/null records and locked biological/synthetic references
- Phase 8: identical classifier folds, leakage prevention and frozen supervised metrics
- Phase 9: time point 고정 및 cell-type/activation task 분리

Full dataset을 unit-test fixture로 사용하지 않는다.

# 9. Codex execution policy

Codex는 사용자가 요청한 Phase만 수행하고 자동 진행하지 않는다. 별도 승인 없이
dataset/embedding download, full analysis, gate PASS 변경, classifier/metric 확정 또는
optional extension을 시작하지 않는다.

Phase 보고에는 status, files, implemented/not run, QC, tests, assumptions, pending
decisions, exact user command와 검토할 산출물을 포함한다.

# 10. Revised execution order

```text
0. Project/protocol setup                                      PASS
1A. Full naïve CD4 acquisition/QC                             PASS
1B. Vanilla SCPA global + pairwise execution                  PASS
2. Official GenePT embedding + cell-wise GenePT-w QC          PASS
3. Whole-cell GenePT-w CD4/CD8 SCPA-core feasibility          PASS
4. Pathway-specific Vanilla vs GenePT-informed comparison    NEXT
5. Pathway-internal gene contribution
6. True/Permuted/Random semantic controls + robustness
7. External/biological/synthetic reference validation
8. Optional classifier/separability analysis
9. 12 h / 24 h extension
10. Final interpretation/presentation/report
```

# 11. Most important principle

```text
먼저 Vanilla SCPA를 재현한다.
그 다음 published GenePT-w를 각 cell에 대해 재현한다.
그 다음 동일 GSE212270 family의 CD4 0 h와 CD8 0 h에서
known cell-type distinction 보존 여부를 검증한다.
그 다음 동일 pathways/cells/paired genes에서 pathway-specific representation 차이를
평가한다. 그 후에만 contribution, semantic controls, reference validation,
classifier와 time-course extension을 수행한다.
```

최종 목표는 “GenePT를 넣었더니 결과가 달라졌다”가 아니라 다음을 검증하는 것이다.

> GenePT-w가 known CD4/CD8 population distinction을 보존하는가, 그 차이가
> semantic mapping에 의한 것인가, 그리고 SCPA multivariate framework의 adapted
> use가 타당한가? 또한 별도 ground truth와 contribution 정의 아래 GenePT semantic
> information이 pathway 내부 gene 차이의 해석을 개선하는가?
