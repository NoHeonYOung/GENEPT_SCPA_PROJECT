# Phase 3 — CD4 0 h vs CD8 0 h GenePT-w benchmark protocol

Status: `READY_FOR_USER_RUN`.

## Research question and scope

Primary question은 GenePT-w에 naïve CD4/CD8 biological distinction이 보존되고,
SCPA가 사용하는 multivariate core가 그 distributional difference를 검출하는가이다.
0 h를 먼저 비교해 activation/time effect를 섞지 않는다. 이 Phase에서는 classifier,
representation 우열, pathway 또는 gene importance를 평가하지 않는다.

## Naïve CD8 GenePT-w

Phase 2에서 checksum 검증한 `GenePT_gene_embedding_ada_text.pickle`
(`text-embedding-ada-002`, 1,536D)을 그대로 사용한다. 새 artifact를 다운로드하지
않는다.

```text
sparse RNA/counts
→ per-cell total 10,000 over all dataset genes
→ log1p
→ identical official exact/alias lookup policy
→ expression-weighted GenePT aggregation
→ row-wise unit L2 normalization
```

CD8 output은 `data/processed/genept_scpa/phase3/`에 별도로 기록하며 기존 CD4 Phase 2
산출물은 이동하거나 수정하지 않는다. CD8 QC는 exact/alias/unmatched/duplicates,
expression coverage, finite/zero-vector/L2, synthetic/direct/determinism 및 CD4 median
coverage와의 차이를 기록한다. Coverage threshold로 dataset을 제외하지 않는다.

## Canonical cohort

- Primary populations: naïve CD4 `Hour=0` vs naïve CD8 `Hour=0`
- Expected full metadata counts: 4,428 vs 1,048; 실행 시 다시 검증
- Sampling: seed 20260810, without replacement, 500 cells/group
- Basis: SCPA 1.6.2 default downsample과 Phase 1 frozen convention
- Same canonical IDs: GenePT-w와 original expression에서 동일하게 재사용

모든 0 h ID와 sampled canonical ID를 `phase3_sampling/`에 저장한다.

## SCPA-core methodology

Source: SCPA v1.6.2 official repository
<https://github.com/jackbibby1/SCPA/tree/v1.6.2> and the locally installed SCPA 1.6.2
namespace used by this project.

로컬 설치 source를 확인한 결과 SCPA 1.6.2의 `compare_pathways()`는
`single_comparison()`으로 연결되고, 각 pathway matrix를 genes × cells에서
cells × genes로 transpose/sort한 뒤 `multicross::mcm(..., level=0.05)`을 호출한다.
SCPA는 pathway별 p-value에 Bonferroni correction을 적용하고
`qval = sqrt(-log10(adjPval))`을 계산한다.

Phase 3 adapter는 1,536 dimensions를 genes/pathway로 포장하지 않는다. Explicitly
sampled cells × features matrices 두 개를 동일한 `multicross::mcm()`에 직접 전달한다.
한 global hypothesis이므로 Bonferroni factor는 1이고 결과를 SCPA-style qval과 함께
기록한다. `multicross::mcm()`은 combined cells의 Euclidean distance를 사용한다.
Adapter 내부 hidden downsampling이나 추가 GenePT normalization은 없다.
Non-bipartite matching을 위해 combined cell count가 even인지도 입력 gate에서 확인한다.
MCM p-value가 machine precision 아래로 underflow하면 raw p=0 사실을 보존하고,
finite display용 qval은 `.Machine$double.xmin`에서 cap한 별도 필드로 기록한다.

SCPA wrapper의 `min_genes`/`max_genes`는 pathway filtering parameter이지 MCM 자체의
feature limit가 아니다. 이번 global adapter에는 pathway filter를 적용하지 않고 각
test에서 GenePT 1,536 dimensions 또는 aligned original 17,085 genes를 모두 사용한다.
이 큰 dimension 차이는 distance geometry와 p-value에 영향을 줄 수 있으므로 raw
statistics를 서로 representation-quality comparison에 사용하지 않는다.

따라서 명칭은 `SCPA-core multivariate framework adaptation`이며 standard SCPA
pathway analysis가 아니다. Nearly-identical/shifted toy populations의 qval ordering을
full run 전에 검사한다.

## Original-expression reference

CD4/CD8 feature position은 다르므로 concatenate하지 않는다. Selected cells의 각
RNA/counts matrix를 해당 dataset의 모든 genes 기준 total 10,000으로 normalize하고
log1p한 후 exact shared gene symbols를 lexicographic order로 정렬한다. Expected shared
count는 17,085다. GenePT-w 자체를 이 intersection으로 다시 만들지는 않는다.

Original reference는 original space에서도 population difference가 존재하는지를
확인할 뿐이다. 17,085D original과 1,536D unit-L2 GenePT의 raw p/q 값은 representation
quality score로 직접 비교하지 않는다.

## Commands

```bash
PYTHONPATH=src python scripts/genept/build_genept_w.py --dataset naive_cd8
PYTHONPATH=src python scripts/phase3/run_cd4_cd8_benchmark.py
```

첫 명령은 full CD8 GenePT-w와 독립 QC를 생성한다. 두 번째 명령은 canonical 0 h
cohort와 두 representation을 준비하고 R SCPA-core adapter를 호출해 Phase 3 QC와
summary를 만든다. 성공 gate는 `READY_FOR_GPT_REVIEW`이며 Phase 4로 자동 진행하지
않는다.
