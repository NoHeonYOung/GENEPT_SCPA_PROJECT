# Phase 1B — Vanilla SCPA reproduction protocol

Status: `PASS`. Hour-only 네 분석과 official two-population reference가 완료되었고,
고정 gate criterion 9개를 모두 통과했다. 이 문서는 완료된 Phase 1B protocol을
규정한다.
GenePT는 사용하지 않는다.

## Official sources reviewed

- SCPA paper: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10704209/>
- Two-population naïve CD4 tutorial:
  <https://jackbibby1.github.io/SCPA/articles/comparing_two_populations.html>
- Multisample/pseudotime tutorial:
  <https://jackbibby1.github.io/SCPA/articles/pseudotime.html>
- `compare_pathways()` reference:
  <https://jackbibby1.github.io/SCPA/reference/compare_pathways.html>
- `seurat_extract()` reference:
  <https://jackbibby1.github.io/SCPA/reference/seurat_extract.html>
- Paper analysis code:
  <https://github.com/jackbibby1/scpa_paper/blob/main/figure_4.R>
- SCPA source/tag: <https://github.com/jackbibby1/SCPA/tree/v1.6.2>

## Frozen decisions

1. **Assay:** `RNA`. 공식 `seurat_extract()` 기본값이며 paper code도
   `RNA@data`를 사용한다. `integrated` assay는 사용하지 않는다.
2. **Expression:** 이미 객체에 존재하는 log1p-normalized `data` layer/slot.
   counts, scale.data 또는 새 normalization을 사용하지 않는다.
3. **Extraction:** `SCPA::seurat_extract(..., assay="RNA", meta1="Hour")`.
   함수의 공식 기본 pseudocount `0.001`을 유지한다.
4. **Pathways:** SCPA v1.6.2의
   `combined_metabolic_pathways.csv`. Hallmark 5, KEGG 54, Reactome 184의
   수동 curated metabolic collection, 총 243 pathways다.
5. **Pathway provenance:** Git tag `v1.6.2`, Git blob
   `4b794712ae7e411753c5199afb72322ae58d52b5`, downloaded file SHA-256
   `6bc5977da3fa60f86d5ffb59fc938740bf418fa4d976182a314d65479eb8b744`.
6. **Input pathway sizes:** min 5, median 22, mean 44.11, max 739 genes.
7. **SCPA filtering:** `min_genes=15`, `max_genes=500`, dataset에 실제 매칭된
   genes 기준으로 inclusive filtering한다. 이는 현재 SCPA 공식 defaults다.
8. **Downsampling:** condition당 500 cells. 500보다 적으면 전부 사용한다.
9. **Seed:** `20260810`. SCPA의 `random_cells()`가 base R `sample()`을 쓰므로
   같은 방식으로 timepoint별 500 cell IDs를 명시적으로 먼저 선택하고 저장한다.
10. **Sampling reuse:** 선택된 각 timepoint의 동일한 500 cells를 global과 모든
    pairwise 분석에서 재사용한다. 이후 fair comparison에서 같은 cells를 재사용할
    수 있도록 ID를 저장한다.
11. **Global API:** `compare_pathways(samples=list(h0,h12,h24), ...)`.
12. **Pairwise API:** 같은 함수에 두 matrices를 순서대로 전달한다.
13. **Statistics:** 모든 결과는 `Pathway`, `Pval`, `adjPval`, `qval`을 가진다.
    두 sample일 때만 `FC`가 추가된다. FC는 population1 minus population2이며
    음수이면 population2에서 더 높다는 뜻이다. SCPA qval은 package convention상
    값이 클수록 multivariate pathway difference가 강하다. 모든 표와 rank는 qval
    내림차순이며 rank 1은 최대 qval, qval=0은 weakest end다. qval을 primary
    statistic으로 사용하고 FC는 secondary로만 해석한다.
14. **Parallelism:** 최초 reproduction은 official default인 single-core
    (`parallel=false`)로 고정한다. GPU는 사용하지 않는다.
15. **Figures:** 실제 결과 생성 후 paper/tutorial 비교용 global qval rank,
    0-vs-24 `-FC`/qval scatter, 네 분석 qval heatmap을 자동 생성한다. heatmap은
    분석 간 timing pattern 검토를 위해 추가한 진단 그림이며 paper 원본 panel은
    아니다.

## Population definition and deviation from the tutorial

완료된 네 secondary 분석은 metadata `Hour`의 전체 cell counts(0 h 4,428;
12 h 4,547; 24 h 5,919)를 사용한다. 따라서 `Hour`만으로 subset하며
`Cell_Type`을 추가 filtering하지 않는다.

이는 공식 two-population tutorial의 `Resting@0` 대 `Activated@24` 비교와
완전히 동일하지 않다. 실제 객체에는 각 시간대에 Resting, Intermediate,
Activated, Treg labels가 섞여 있다. 또한 global 0/12/24 real-time comparison은
paper의 pseudotime milestone comparison과 다른, 요청된 multisample extension이다.
이 차이는 결과 해석에서 반드시 유지하며 결과를 paper Figure 4의 수치적 재현이라
부르지 않는다.

Phase 1B final reference는 이와 별개로 실제 metadata spelling/case를 확인한 뒤
다음 두 population을 사용한다.

```text
population 1: Cell_Type=Resting   AND Hour=0
population 2: Cell_Type=Activated AND Hour=24
```

각 population은 기존과 동일한 seed와 base R `sample()` 방식으로 최대 500 cells를
선택한다. `RNA/data`, `SCPA::seurat_extract()`, pseudocount 0.001, official pathway
file, min/max genes 15/500 및 single-core 설정도 동일하게 유지하며 결과를 맞추기
위한 tuning은 하지 않는다.

## Four analyses

```text
global_0_12_24       list(h0, h12, h24)   qval, no FC
pairwise_0_vs_12     list(h0, h12)        qval + FC (0 h - 12 h)
pairwise_12_vs_24    list(h12, h24)       qval + FC (12 h - 24 h)
pairwise_0_vs_24     list(h0, h24)        qval + FC (0 h - 24 h)
```

Global qval은 세 distribution 중 동일하지 않은 것이 있는지를 나타낼 뿐 변화
시점이나 방향을 나타내지 않는다. pairwise qval과 FC를 함께 확인해야 한다.

## Paper/tutorial comparison figures

전체 실행은 `data/processed/genept_scpa/phase1/figures/` 아래에 다음 파일을
추가로 생성한다.

```text
01_global_qval_rank.png
02_0_vs_24_enrichment_qval.png
03_four_analysis_qval_heatmap.png
phase1b_paper_comparison.png
phase1b_paper_comparison.pdf
phase1b_figure_notes.md
```

첫 패널은 공식 global qval/rank 표현을 따르고 glycolysis-related pathway를
강조한다. 두 번째 패널은 공식 two-population tutorial처럼 `-FC` 대 `qval`을
그리고 arachidonic-related pathway를 강조한다. 세 번째 heatmap은 네 분석을
같이 비교하기 위해 추가한 것으로, 논문 그림을 그대로 복제한 panel이 아니다.

이 그림들은 Hour-only 결과와 논문의 Cell_Type-specific/pseudotime 결과 사이의
qualitative comparison만 지원한다. 정확한 좌표나 순위가 일치해야 한다는 뜻은
아니다. CSV가 이미 있다면 SCPA 재계산 없이 아래 명령으로 그림만 다시 만든다.

```bash
Rscript scripts/scpa/plot_phase1b_results.R
```

## Official two-population reference 실행

기존 01–04 CSV를 보존한 상태에서 아래 명령 하나만 실행한다.

```bash
Rscript scripts/scpa/reproduce_scpa_phase1b_reference.R
```

예상 산출물:

```text
data/processed/genept_scpa/phase1/05_reference_resting0_vs_activated24.csv
data/processed/genept_scpa/phase1/figures/05_reference_resting0_vs_activated24_qval_fc.png
data/interim/genept_scpa/phase1b_sampling/reference_resting0_cells.txt
data/interim/genept_scpa/phase1b_sampling/reference_activated24_cells.txt
data/interim/genept_scpa/phase1b_scpa_qc.json
data/interim/genept_scpa/phase1b_reproduction_summary.md
```

05 figure는 x축에 `-FC`를 사용하므로 양수는 Activated 24 h enrichment를 뜻한다.
Reactome/KEGG arachidonic-acid pathways를 강조하며 numerical identity가 아닌
qualitative reproduction만 평가한다.

## Biological reproduction targets

- 공식 tutorial/paper는 arachidonic-acid metabolism이 큰 multivariate change를
  보이지만 반드시 큰 mean enrichment를 보이는 것은 아니라고 강조한다.
- multisample pseudotime tutorial은 glycolysis-related pathways가 상위에 있음을
  예시로 든다.
- 이번 real-time/Hour-only 설계는 이 신호를 qualitative review 대상으로만
  사용한다. 결과를 맞추기 위한 parameter 변경은 금지한다.

## PASS candidate requirements

기존 Hour-only 네 분석 PASS, qval 내림차순 해석 수정, reference 실행 성공,
finite qval/FC, official workflow 보존, 두 arachidonic target 확인, qualitative
comparison 기록, critical runtime warning/error 없음, parameter tuning 없음이 모두
충족되면 `PASS`다. Reference 실행 또는 population 재현이 실패하면
`NEEDS_REVIEW`다. Numerical identity는 PASS 조건이 아니다.
