# GenePT × SCPA Decision Log

결과를 확인한 뒤 분석 조건을 바꾸지 않도록, Phase를 여는 결정과 프로토콜
변경을 여기에 기록한다.

## D-0001 — Phase 0 프로젝트 분리

- 날짜: 2026-08-10
- 상태: Accepted
- 결정: GenePT × SCPA 코드를 `src/gene_embedding_project/genept_scpa/` 아래에
  분리하고 모든 실행 조건은 `config/genept_scpa.yaml`에서 관리한다.
- 근거: 기존 LOO/frequency/hybrid 분석과 논리적으로 분리하고 Phase gate를
  코드에서도 확인하기 위함이다.

## D-0002 — Phase 1 전체 데이터 원본

- 날짜: 2026-08-10
- 상태: Accepted
- 결정: Phase 1에는 GEO GSE212270의
  `GSE212270_integrated_naive_cd4.rds.gz`를 사용한다.
- 범위: naïve CD4 T cell의 0 h, 12 h, 24 h가 통합된 저자 공개 Seurat
  객체 전체이며, 축소 예제 객체를 사용하지 않는다.
- 근거: SCPA 논문 코드 저장소가 GSE212270을 raw/processed sequencing data의
  공식 위치로 지정하며, GEO가 해당 파일을 2.4 GB processed RDS로 제공한다.
- 제외: Phase 1은 Cell Ranger부터의 원시 재처리 재현이 아니므로
  `GSE212270_RAW.tar` 및 SRA FASTQ는 primary input으로 사용하지 않는다.

## D-0003 — Phase 1은 아직 잠금 상태

- 날짜: 2026-08-10
- 상태: Accepted
- 결정: Phase 0 중에는 데이터 다운로드 및 SCPA 실행을 하지 않는다.
  다운로드 스크립트와 검증 절차만 준비한다.
- 다음 결정: Phase 0 PASS 확인 후 `active_phase`와 `max_phase_allowed`를 1로
  바꾸는 별도 기록을 남겨 Phase 1을 연다.

## D-0004 — Phase 0 PASS

- 날짜: 2026-08-10
- 상태: Accepted
- PASS 근거:
  - 새 프로젝트 package/script/config/docs 구조가 존재한다.
  - YAML config loader와 phase gate unit test 3개가 통과했다.
  - decision log가 존재하고 dataset 및 phase-lock 결정을 기록한다.
  - R 스크립트 parse check와 Bash syntax check가 통과했다.
  - Phase 1 다운로드 스크립트의 lock 동작을 확인했으며 데이터 파일은 생성되지
    않았다.
- 결정: Phase 0은 PASS로 종료한다. 단, Phase 1은 별도 승인 기록 전까지
  잠금 상태를 유지한다.

## D-0005 — Phase 1A dataset acquisition / validation gate 시작

- 날짜: 2026-08-10
- 상태: Accepted
- 결정: Phase 1을 활성화하고 상태를 `in_progress`로 둔다. 이번 범위는
  Phase 1A acquisition/inspection gate까지이며 Phase 1 전체 PASS가 아니다.
- Dataset: GEO GSE212270의 full processed
  `GSE212270_integrated_naive_cd4.rds.gz`.
- 준비된 산출물:
  - 다운로드 후 실제 byte size, SHA-256, gzip integrity를 기록하는
    `phase1_download_metadata.json` 생성 절차
  - Seurat object, time points, identity evidence, assay/layer, feature identifiers를
    검사하는 `phase1_dataset_qc.json` 생성 절차
- 현재 사실: 실제 2.4 GB 파일은 아직 다운로드하지 않았으며 실제 Seurat
  validation 및 SCPA 분석도 실행하지 않았다.
- Gate: 실제 QC의 `gate.status`가 `PASS`가 되기 전에는 Phase 1B SCPA
  reproduction으로 진행하지 않는다.

## D-0006 — Phase 1A dataset gate PASS

- 날짜: 2026-08-10
- 상태: Accepted
- 실제 입력: GEO GSE212270
  `GSE212270_integrated_naive_cd4.rds.gz`.
- Acquisition QC:
  - file size: 2,530,955,413 bytes
  - SHA-256: `f9ad7cfbe8bee87a28cec76dd66d442ad2c5ddb243942c51d905713b9b2b7842`
  - gzip integrity: PASS
- Seurat QC:
  - serialized version 3.1.5; SeuratObject 5.4.0으로 메모리상 호환 갱신 성공
  - source RDS modification: false
  - cells: 14,894; features: 17,856
  - assays: RNA, integrated; active assay: RNA
  - RNA layers: counts, data, scale.data
  - 0 h: 4,428; 12 h: 4,547; 24 h: 5,919 cells
  - gene identifiers: gene-symbol-like, missing 0, duplicates 0
- Gate result: failed checks 0, warnings 2, `gate.status = PASS`.
- Warning 해석:
  - legacy Seurat compatibility update는 메모리상 수행됐고 원본은 변경하지 않았다.
  - object 내부에 naïve CD4라는 명시적 metadata/project label은 없어서 공식 GEO
    filename/accession을 identity 근거로 사용한다.
- 결정: Phase 1A는 PASS. Phase 1 전체는 계속 `in_progress`이며 Phase 1B SCPA
  reproduction은 별도 요청 전까지 진행하지 않는다.

## D-0007 — Phase 1B Vanilla SCPA protocol freeze

- 날짜: 2026-08-10
- 상태: Accepted
- 목적: GenePT 없이 global 0/12/24 및 세 pairwise Vanilla SCPA analyses를
  재현 가능하게 실행한다.
- 공식 representation: `RNA` assay의 기존 log1p-normalized `data`,
  `SCPA::seurat_extract()` 사용, pseudocount 0.001.
- Pathways: SCPA v1.6.2 `combined_metabolic_pathways.csv`, 243 curated metabolic
  pathways (Hallmark/KEGG/Reactome), SHA-256 고정.
- Parameters: seed 20260810, downsample 500, min/max matched genes 15/500,
  single-core official default.
- Sampling: SCPA `random_cells()`와 같은 base R `sample()`로 각 timepoint의
  500 IDs를 먼저 고정하고 모든 분석에 재사용한다.
- Population decision: 이번 요청의 전체 Hour별 cell counts를 따르므로
  `Hour`만 subset하고 `Cell_Type`은 filtering하지 않는다.
- 공식 tutorial과의 차이: tutorial의 직접 0/24 example은 Resting@0 대
  Activated@24이며, paper multisample example은 real-time 세 group이 아니라
  pseudotime milestones다. 이번 네 분석은 공식 API/representation/pathways를
  유지한 real-time Hour-only extension으로 해석한다.
- 상태: 실행 코드 준비 후 `READY_FOR_USER_RUN`; 실제 full SCPA 결과 없이는
  Phase 1B PASS로 바꾸지 않는다.

## D-0008 — Phase 1B paper/tutorial comparison figures

- 날짜: 2026-08-10
- 상태: Accepted
- 결정: 결과 CSV 생성 후 global qval rank, 0-vs-24 enrichment/qval scatter,
  네 분석 qval heatmap 및 세 panel 조합 PNG/PDF를 자동 생성한다.
- 공식 대응: global rank/heat 색상과 two-population `-FC`/qval 표현은 SCPA 공식
  Figure 4 reproduction tutorial을 따른다. glycolysis와 arachidonic-related
  pathways는 사전에 정한 qualitative target으로만 강조한다.
- 추가 진단: 네 분석 heatmap은 논문 원본 panel이 아니라 real Hour 분석의 변화
  시점을 비교하기 위한 project-specific visualization이다.
- 해석 제한: 현재 분석은 전체 Hour population을 사용하므로 논문의
  Cell_Type-specific 또는 pseudotime milestone 결과와 수치적으로 동일하다고
  주장하지 않는다.

## D-0009 — Primary biological benchmark refinement

- 날짜: 2026-08-10
- 상태: Accepted
- 배경: 의도한 biological benchmark를 명확히 한 뒤 primary GenePT-SCPA
  feasibility question을 재정의했다.
- 기존 강조점: pathway-specific GenePT/SCPA integration과 synthetic pathway
  recovery를 main route로 우선 검토.
- 새 primary emphasis: published GenePT-w가 GSE212270 naïve CD4-vs-CD8
  cell-population separability를 보존하는지, 그리고 SCPA multivariate comparison
  framework가 이 차이를 검출하도록 adaptation될 수 있는지 검토한다.
- 첫 benchmark: activation confounding을 줄이기 위해 CD4 0 h vs CD8 0 h를
  primary comparison으로 둔다. CD4/CD8 annotation은 known/reference class label로
  supervised separability 평가에 사용한다.
- Dataset 역할: validated integrated naïve CD4 object는 Phase 1 reproduction,
  Phase 2 GenePT-w validation 및 main experiment의 CD4 side에 사용한다. 동일
  GSE212270 family의 naïve CD8 object는 이후 별도 acquisition/QC하며 아직 수치나
  object 구조를 가정하지 않는다.
- Method boundary: GenePT-w 1536 dimensions는 genes가 아니므로 일반
  `compare_pathways()` 입력으로 취급하지 않는다. SCPA core statistic을 사용한다면
  `adaptation of the SCPA multivariate comparison framework`로 구분하고 별도
  methodological gate를 거친다.
- Secondary scope: GenePT-derived gene programs, pathway-specific GenePT-w,
  synthetic pathway benchmark와 external embedding은 primary CD4/CD8 experiment
  이후의 exploratory/optional extension으로 유지한다.
- 근거: known cell-type labels를 이용해 초기 biological question을 직접 검증하고,
  semantic mapping 효과를 Permuted/Random controls와 구분할 수 있다.
- 상태 영향: Phase 0과 Phase 1A PASS 유지. Phase 1B는 NOT RUN/NEXT, Phase 2
  이후와 naïve CD8 acquisition은 NOT STARTED이며 자동 Phase 진행은 없다.

## D-0010 — Naïve CD8 acquisition pipeline prepared

- 날짜: 2026-08-10
- 상태: Accepted
- 결정: 향후 Phase 3 CD4-vs-CD8 benchmark를 위해 official GSE212270
  `GSE212270_integrated_naive_cd8.rds.gz`의 resumable download, provenance metadata,
  dynamic Seurat validation 및 read-only CD4 compatibility pipeline을 준비한다.
- 공식 근거: GSE212270 Series는 naïve CD8 0/12/24 h samples와 약 1.1 GB의
  integrated naïve CD8 processed RDS를 제공한다.
- 제한: 실제 CD8 파일은 아직 다운로드/검증하지 않았고 cell/feature/assay 구조를
  가정하지 않는다. CD4/CD8 merge, intersection matrix, GenePT-w, SCPA 또는
  classifier를 실행하지 않는다.
- 상태 영향: Phase 0/1A/1B 및 Phase 2 이후의 gate/status를 변경하지 않는다.

## D-0011 — Phase 1B qval correction and official reference finalization

- 날짜: 2026-08-10
- 상태: Accepted
- 해석 수정: SCPA qval은 클수록 multivariate pathway difference가 강하므로 모든
  pathway ranking을 내림차순으로 고쳤다. rank 1은 최대 qval이며 qval=0은
  weakest end다. 기존 CSV statistic은 수정하지 않는다.
- 분석 구분: 완료된 네 Hour-only 분석은 유효한 secondary time-point grouped
  analysis로 보존하지만 exact paper Figure 4 또는 tutorial population reproduction으로
  부르지 않는다.
- reference 결정: official two-population workflow에 맞춰 `Resting 0 h`와
  `Activated 24 h`를 비교하는 별도 reference pipeline과 05 output을 추가했다.
- 고정 설정: RNA/data, `SCPA::seurat_extract()`, pseudocount 0.001, SCPA v1.6.2
  combined metabolic pathways, seed 20260810, downsample 500, min/max 15/500,
  single-core를 유지하고 result-driven tuning을 금지한다.
- Gate: 실제 full reference 실행 전 상태는 `REFERENCE_READY_FOR_USER_RUN`이다.
  Phase 1B PASS는 성공적인 pipeline/QC와 qualitative reference reproduction에
  근거하며 tutorial과 numerical identity를 요구하지 않는다. 실패하거나 population
  definition을 재현하지 못하면 `NEEDS_REVIEW`다.
- 범위: GenePT 및 Phase 2/3는 시작하지 않는다.

## D-0012 — Phase 1B final PASS

- 날짜: 2026-08-10
- 상태: Accepted
- 실제 reference: `Cell_Type=Resting AND Hour=0` 4,341 cells와
  `Cell_Type=Activated AND Hour=24` 1,697 cells에서 고정 seed로 각 500 cells를
  사용했다.
- 결과 QC: 124 pathways의 qval/FC가 모두 finite이고 runtime warning/error가 없다.
  Reactome/KEGG arachidonic-acid targets가 존재하며 작은 FC 대비 명확한 qval을
  보여 qualitative agreement가 기록됐다.
- Gate: Hour-only 네 분석, qval 해석, official workflow, reference, qualitative
  comparison 및 no-tuning criterion이 모두 true이므로 Phase 1B를 PASS로 닫는다.
- 제한: Numerical identity를 요구하지 않으며 Hour-only 분석은 exact Figure 4로
  부르지 않는다.

## D-0013 — Phase 2 Published GenePT-w protocol freeze

- 날짜: 2026-08-10
- 상태: Accepted
- 목적: GSE212270 naïve CD4 전체에서 published GenePT-w representation 생성의
  correctness만 검증한다. SCPA와 downstream separability 분석은 제외한다.
- Primary embedding: author-provided Zenodo DOI `10.5281/zenodo.10833191`의
  `GenePT_gene_embedding_ada_text.pickle`, model `text-embedding-ada-002`, 1,536D.
  OpenAI API로 재생성하지 않는다.
- Expression/preprocessing: sparse Seurat `RNA/counts` → 모든 dataset genes를 포함한
  cell-wise 10,000 normalization → log1p → official artifact exact-key alignment →
  expression-weighted aggregation → row-wise unit L2 normalization.
- Mapping: fuzzy/case conversion은 금지한다. Primary NCBI key가 아닌 artifact key는
  저자 methods의 official HGNC alias lookup으로 구분하고 coverage를 기록한다.
- Correctness: hand-calculable synthetic test, sparse-vs-direct real-cell check,
  determinism, finite/zero-vector/L2/cell-ID gate를 통과해야 한다.
- 실행 상태: acquisition과 full CD4 calculation은 사용자 실행 대상으로 두고 Phase 2는
  `READY_FOR_USER_RUN`이다. Full 결과는 `READY_FOR_GPT_REVIEW`에서 멈춘다.
- 범위 제한: CD8 projection, CD4-vs-CD8, classifier, SCPA, Phase 3는 시작하지 않는다.

## D-0014 — Phase 2 Published GenePT-w PASS

- 날짜: 2026-08-11
- 상태: Accepted
- 실제 결과: GSE212270 naïve CD4 14,894 cells × 1,536 dimensions를 생성했다.
  17,856 RNA features 중 official GenePT lookup과 14,409가 매칭됐고(exact 14,087;
  official alias 322), median raw-count mass coverage는 0.931967이었다.
- QC: all finite, zero vectors 0, post-L2 norms approximately 1, synthetic/direct/
  determinism tests PASS, source Seurat object unmodified, warning 및 failed check 0.
- 결정: 기존 `phase2_genept_w_qc.json`의 `READY_FOR_GPT_REVIEW` 결과를 검토해
  Phase 2를 PASS로 닫는다. Historical QC gate 문자열은 실행 provenance로 보존한다.
- 범위: 이 PASS는 published GenePT-w correctness에 대한 것이며 CD4/CD8 difference,
  classifier 성능 또는 pathway/gene-level 우수성을 의미하지 않는다.

## D-0015 — Supervisor clarification on primary and pathway-level research questions

- 날짜: 2026-08-11
- 상태: Accepted
- 결정: CD4/CD8 expression을 GenePT-w로 변환하고 SCPA multivariate framework로
  population difference를 검출하는 해석 방향을 확인했다. CD4 0 h vs CD8 0 h를 첫
  direct feasibility benchmark로 유지한다.
- 추가 핵심 질문: 프로젝트는 cell-type separability에 한정되지 않는다. GenePT의
  semantic gene information이 pathway 내부 gene 차이/contribution의 표현과 해석을
  개선할 수 있는지도 core downstream question으로 평가한다.
- 해석 제한: "더 정확한" gene-level interpretation은 다른 ranking이 나왔다는
  사실만으로 성립하지 않는다. Contribution 정의, reference/ground truth와 평가
  criterion을 별도 methodological gate에서 먼저 고정한다.
- 보존 후보: leave-one-gene-out, perturbation, pathway-score sensitivity, GenePT
  semantic grouping, known markers/reference pathways, GenePT-derived programs,
  pathway-specific GenePT-w와 synthetic ground-truth benchmark를 후보로 유지하되
  어떤 것도 현재 standard method로 확정하지 않는다.
- Phase 영향: Phase 6을 optional extension에서 core pathway/gene-level interpretation
  phase로 승격하고 robustness를 Phase 7로 이동한다.

## D-0016 — Phase 3 CD4 0 h vs CD8 0 h method freeze

- 날짜: 2026-08-11
- 상태: Accepted
- CD8 GenePT-w: Phase 2의 동일 ada-002 1,536D artifact와 동일 RNA/counts → total
  10,000 → log1p → official exact/alias lookup → weighted aggregation → row L2
  pipeline을 재사용한다. 새 embedding download는 금지한다.
- Cohort: 실제 metadata의 0 h를 재확인하고 SCPA default 및 Phase 1 convention에
  근거해 seed 20260810, 500 cells/group을 한 번 명시적으로 sampling한다. 같은 IDs를
  GenePT-w와 original-expression reference에 재사용한다.
- SCPA-core: SCPA 1.6.2 `single_comparison()`이 호출하는 `multicross::mcm()`을
  cells × features 두 population에 직접 호출한다. 1,536 dimensions를 genes 또는
  pathway로 취급하지 않으며 결과를 `SCPA-core multivariate framework adaptation`으로
  명명한다.
- Original reference: 각 dataset의 모든 RNA/counts genes로 total 10,000 normalization
  후 log1p하고 exact shared symbols(현재 QC 기대 17,085)로 정렬한다.
- 해석 제한: Original과 GenePT-w의 dimension/geometry가 다르므로 raw MCM p/q 값의
  크기로 representation 우열을 주장하지 않는다. Classifier와 controls는 Phase 4다.
- 실행 상태: code/small tests 이후 `READY_FOR_USER_RUN`; full 결과 후
  `READY_FOR_GPT_REVIEW`에서 멈춘다.

## D-0017 — Phase 3 final PASS

- 날짜: 2026-08-11
- 상태: Accepted
- 근거: `phase3_cd4_cd8_qc.json`의 gate가 `READY_FOR_GPT_REVIEW`이고 14개
  criterion이 모두 true이며 failed check와 warning이 없다. 동일한 canonical 500
  CD4/500 CD8 세포에서 GenePT-w와 original-expression 양쪽 모두 CD4/CD8
  multivariate difference를 검출했다.
- 결정: Phase 3를 PASS로 닫고 Phase 4를 연다. 이 PASS는 whole-cell GenePT-w
  feasibility에 대한 것이며 classifier accuracy나 GenePT 우수성을 뜻하지 않는다.

## D-0018 — Phase 3 qval implementation audit PASS

- 날짜: 2026-08-11
- 상태: Accepted
- 직접 확인: 로컬 SCPA 1.6.2 `single_comparison()`과 Phase 3
  `scpa_core_adapter.R`를 비교했다.
- 공식 계산: pathway별 `multicross::mcm()` raw p-value를 얻고
  `stats::p.adjust(Pval, method="bonferroni", n=eligible_pathway_count)` 후
  `qval = sqrt(-log10(adjPval))`을 사용한다. 로그 밑은 10이다.
- Phase 3 판단: 단일 global hypothesis이므로 Bonferroni factor가 1이며 기존
  `sqrt(-log10(raw_p))`와 공식 convention이 정확히 같다. mismatch가 없으므로 raw
  p-value와 historical output을 수정하지 않는다.
- Phase 4 결정: 두 branch 모두 동일한 eligible pathway universe를 사용하고 각
  branch의 전체 raw p-value vector에 위 공식 Bonferroni/qval convention을 적용한다.

## D-0019 — Phase 4 pathway-specific paired comparison method freeze

- 날짜: 2026-08-11
- 상태: Accepted
- Primary analysis: Vanilla pathway expression SCPA와 pathway-specific GenePT
  semantic projection 뒤의 SCPA-core multivariate comparison을 비교한다.
- Whole-cell 1,536D GenePT-w를 pathway별로 slice하지 않는다. 1,536 dimensions는
  pathway label이나 gene이 아니므로 각 pathway의 original expression으로 돌아가
  `Z_P = X_P E_P`를 새로 계산한다.
- Primary paired gene policy: `pathway genes ∩ CD4 genes ∩ CD8 genes ∩ official
  GenePT keys`. 두 branch는 동일한 canonical cells, eligible pathways, gene symbols,
  gene order를 사용한다. 이는 GenePT branch의 mapping loss가 representation 차이로
  혼입되는 것을 막기 위한 것이다.
- Preprocessing: RNA/counts를 cell별 전체 transcriptome total 10,000으로 정규화하고
  log1p한 뒤 pathway를 subset한다. pathway 내부 재정규화는 금지한다.
- Projection: primary는 non-L2 `X_P E_P`다. Vanilla가 가진 pathway expression
  magnitude를 primary paired comparison에서 제거하지 않기 위함이다. Rowwise L2는
  사전 선언된 sensitivity option이며 결과를 보고 primary를 바꾸지 않는다. 이 방법은
  published whole-cell GenePT-w가 아니라 `GenePT-informed pathway projection`이다.
- Effective rank: 1,536 output columns가 독립 변수 1,536개를 뜻하지 않는다.
  embedding/projected rank와 singular-value summary를 pathway별로 기록하고
  `projected_rank <= paired_gene_count`를 확인한다.
- Interpretation: rank correlation, overlap과 rank shift는 agreement/disagreement
  metric이다. qval/p-value/rank 크기로 GenePT의 우수성이나 정확도를 주장하지 않는다.
- 후속 범위: Phase 5는 동일 gene-removal policy의 pathway 내부 contribution,
  Phase 6은 True/Permuted/Random embedding과 repeated-sampling controls, Phase 7은
  biological/synthetic ground truth 기반 accuracy validation이다. 이번 Phase에서는
  어느 것도 실행하지 않는다.

## D-0020 — Phase 4A completed with reporting validation required

- 날짜: 2026-08-11
- 상태: Accepted
- 결과: 123 paired pathways의 naïve CD4 0h-vs-CD8 0h Vanilla/non-L2/L2 실행은
  완료됐고 raw p와 QC는 보존한다. 다수 adjusted p=1/qval=0 pathway가 존재한다.
- 문제: 기존 구현은 동일 qval=0에도 pathway 이름 순으로 서로 다른 unique rank를
  부여했다. 이 순서는 통계적 차이를 나타내지 않으므로 후속 primary reporting에서
  사용하지 않는다.
- 결정: Phase 4A는 계산 완료 상태로 보존하되 Phase 4B/C validation 검토 전 Phase 5
  PASS로 연결하지 않는다. Regenerated output은 별도 directory에 저장한다.

## D-0021 — SCPA qval re-audit and tie-aware reporting

- 날짜: 2026-08-11
- 상태: Accepted
- 재감사: 설치된 SCPA 1.6.2 `single_comparison()`은 Bonferroni adjusted p에
  `sqrt(-log10(adjPval))`을 적용한다. R natural `log()`이 아니므로 Phase 3/4 raw p와
  qval formula mismatch는 없다.
- 결정: historical raw p/qval을 덮어쓰지 않는다. Phase 4A/4B/4C regenerated
  reporting은 descending qval의 average tied rank를 사용한다. Equal qval은 equal
  rank를 받는다.
- Primary categories: adjusted p<0.05를 기준으로 Both significant, Vanilla-only,
  GenePT-only, Neither significant를 기록한다. qval positive/zero와 raw/adjusted
  p<0.05 counts도 함께 기록한다.

## D-0022 — Phase 4B/C activation validation freeze

- 날짜: 2026-08-11
- 상태: Superseded by D-0023 (9-comparison capability는 보존)
- Source audit: CD4/CD8 integrated RDS와 기존 Phase 2 export 모두 0h/12h/24h를
  포함한다. CD4 counts는 4,428/4,547/5,919, CD8 counts는
  1,048/2,066/3,927이다. RNA/counts와 RNA/data가 존재하므로 download하지 않는다.
- Canonical groups: seed 20260810, 각 500 cells. CD4/CD8 0h는 Phase 3 IDs를 그대로
  재사용하고 12h/24h는 한 번 생성한 ID와 SHA256을 이후 모든 comparison에서 재사용한다.
- Pathways: Phase 4A의 동일 123 pathways와 exact paired gene lists를 여섯 group
  전체에서 고정한다. Timepoint별 universe 변경이나 결과 기반 threshold 변경은 금지한다.
- Phase 4B: CD4 0h-vs-24h를 positive control로 두고 representative pathways의
  adapter Vanilla raw p를 official `SCPA::compare_pathways()`와 tolerance 1e-12로
  cross-check한다.
- Phase 4C: CD4 activation 3, CD8 activation 3, same-time lineage 3의 9 comparisons를
  Vanilla, GenePT non-L2 primary, GenePT rowwise-L2 sensitivity로 실행한다.
- Gate: code/toy/smoke 뒤 `READY_FOR_USER_RUN`, full output 뒤
  `READY_FOR_GPT_REVIEW`. Phase 5는 별도 검토 전 시작하지 않는다.

## D-0023 — Phase 4B primary scope를 Naive CD4 activation으로 제한

- 날짜: 2026-08-11
- 상태: Accepted
- 연구 판단: Phase 4A CD4 0h-vs-CD8 0h는 잘못된 분석이 아니라 exploratory
  cross-lineage feasibility 분석으로 보존한다. 다만 큰 qval floor와 과거 qval=0
  unique-rank artifact 때문에 통제된 representation benchmark의 primary 근거로는
  사용하지 않는다. SCPA 1.6.2 raw-p adapter, Bonferroni와
  `sqrt(-log10(adjusted_p))`, average ties audit는 PASS 상태를 유지한다.
- Primary Phase 4B: lineage를 Naive CD4로 고정하고 activation만 바꾼
  `0h-vs-12h`, `12h-vs-24h`, `0h-vs-24h` 세 comparison만 production으로
  실행한다. 마지막 comparison은 original SCPA activation study와 방향이 맞는
  positive-control-like benchmark이지 ground-truth accuracy proof가 아니다.
- Cells/pathways: 이미 frozen된 CD4 0h/12h/24h 각 500 IDs와 Phase 4A의 동일
  123 pathways/exact paired genes/order를 재사용한다. Whole-transcriptome
  RNA/counts total=10,000 후 log1p하고 pathway를 subset하며 pathway-local
  normalization은 금지한다.
- Branch/reporting: Vanilla, GenePT non-L2 primary, row-L2 sensitivity를 모두
  실행한다. Primary report는 qval floor, raw/adjusted-p counts와 medians,
  detection states이며 average-tie rank는 secondary다. Phase 4A historical
  qval-zero/significant counts를 함께 비교하되 기존 파일을 덮어쓰지 않는다.
- Deferred: CD8 activation 및 12h/24h lineage comparisons는 current production에
  포함하지 않는다. 이전 `all_9` runner 기능만 explicit option으로 보존한다.
  가능한 CD8 0h-vs-24h generalization은 Phase 4C 후보이며 현재 실행하지 않는다.
- Gate: 정확히 3 comparisons x 123 pathways x 3 branches = 1,107 MCM calls.
  Codex는 production을 실행하지 않고 `READY_FOR_USER_RUN`에서 멈춘다. Full run 후
  `READY_FOR_GPT_REVIEW` 검토 전 Phase 5로 진행하지 않는다.

## D-0024 — Phase 4B final PASS

- 날짜: 2026-08-11
- 상태: Accepted
- 결과: CD4 0h-vs-12h, 12h-vs-24h, 0h-vs-24h의 각 123 pathways에서
  Vanilla, GenePT non-L2, row-L2가 모두 완료됐다. QC failed checks=0,
  runtime warnings=0이며 official SCPA raw-p cross-check와 average-tie rank가 PASS했다.
- 결정: 0h-vs-24h positive-control-like activation signal을 포함한 technical gate를
  충족했으므로 Phase 4B를 PASS/COMPLETED로 닫는다. 이는 GenePT superiority 또는
  generalization 근거가 아니다. Phase 4A historical output은 그대로 보존한다.

## D-0025 — Phase 5 paired gene-masking sensitivity method freeze

- 날짜: 2026-08-11
- 상태: Accepted
- Target: Phase 4B adjusted p<0.05 기준 Vanilla-only 및 GenePT-only인 모든
  pathway-comparison pair. Frozen expected N=30(0h-vs-12h 11,
  12h-vs-24h 9, 0h-vs-24h 10)이며 결과 기반 추가/삭제는 금지한다.
- Paired masking: 동일 frozen cells, preprocessing, 123-pathway exact paired genes와
  order를 유지한다. Vanilla는 `X_P[:,g]=0`, GenePT non-L2는
  `Z_P - outer(X_P[:,g], E_P[g,:])`로 같은 gene을 masking한다. Zero-mask/physical
  removal과 subtraction/direct recomputation equivalence를 production 전에 검증한다.
- Metric: primary는 `delta[-log10(raw p)] = score_full - score_masked`; score 계산의
  raw p clipping은 1e-300으로 고정한다. Signed supporting rank와 absolute influence
  rank는 average ties를 사용하며 qval/rank shift를 primary contribution metric으로
  사용하지 않는다.
- Interpretation: gene importance/causality/biological correctness가 아니라
  representation-dependent gene masking sensitivity로만 부른다.
- Deferred: GenePT L2 gene-level analysis, Phase 6 True/Permuted/Random, CD8
  generalization, classifier를 실행하지 않는다. Full run은 baseline Phase 4B raw-p,
  same cells/genes/order, equivalence, deterministic 및 checkpoint/resume gate PASS 후에만
  허용하고 `READY_FOR_GPT_REVIEW`에서 멈춘다.

## Decision template

```text
## D-XXXX — 제목
- 날짜:
- 상태: Proposed | Accepted | Rejected | Superseded
- 결정:
- 근거:
- 고정되는 설정:
- 영향을 받는 Phase:
```
