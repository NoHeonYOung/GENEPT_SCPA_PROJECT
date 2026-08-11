# Phase 2 — Published GenePT-w reproduction protocol

Status: `PASS`. 이 단계는 naïve CD4에서 published GenePT-w matrix를
정확히 생성하는지만 검증한다. SCPA, CD4-vs-CD8 비교와 classifier는 실행하지 않는다.

## Official sources

- Final peer-reviewed paper: Chen Y, Zou J. *Simple and effective embedding model
  for single-cell biology built from ChatGPT.* Nature Biomedical Engineering 9,
  483–493 (2025). DOI: <https://doi.org/10.1038/s41551-024-01284-6>
- Official repository: <https://github.com/yiqunchen/GenePT>
- Official cell-level example:
  <https://github.com/yiqunchen/GenePT/blob/main/aorta_data_analysis.ipynb>
- Author-provided embeddings: <https://zenodo.org/records/10833191>

Final paper를 primary citation으로 사용했다. Subscription page의 공개 abstract와
supplement, author methods, repository code를 함께 확인해 구현 세부사항을 고정했다.

## Primary embedding artifact

```text
Zenodo DOI:       10.5281/zenodo.10833191
Archive:          GenePT_emebdding_v2.zip
Archive size:     574,395,233 bytes
Zenodo MD5:       3f6ce4317e3a0091978ae5cb8fbf05a3
Primary file:     GenePT_gene_embedding_ada_text.pickle
Model:            text-embedding-ada-002
Dimension:        1,536
Schema:           Python dictionary; uppercase gene/alias key -> numpy vector
Summary keys:     NCBI_summary_of_genes.json
```

새로운 `text-embedding-3-large` artifact는 primary reproduction에 사용하지 않는다.
OpenAI API를 호출하지 않고 저자 공개 precomputed artifact를 사용한다. Archive와
추출 파일의 SHA-256은 다운로드 후 provenance JSON에 기록한다.

## Published GenePT-w preprocessing

저자 methods는 다음 순서를 명시한다.

```text
sparse raw RNA counts
→ each cell total normalized to 10,000 transcripts using all dataset genes
→ element-wise log(1 + x)
→ align genes to the official GenePT lookup
→ expression-weighted gene-embedding aggregation
→ unit row-wise L2 normalization
```

Official aorta notebook은 이미 전처리된 AnnData `.X`를 읽고, dataset gene 순서의
lookup matrix를 만든 뒤 unmatched gene을 zero vector로 두고
`X @ lookup / number_of_dataset_genes`를 계산한다. Notebook 자체에는 upstream
normalize/log1p와 final L2 코드가 보이지 않지만 paper methods와 Figure 1은 이를
명시한다. 이번 pipeline은 raw counts에서 이 두 단계를 명시적으로 수행하면서
notebook의 lookup 및 전체 dataset gene 수 denominator를 유지한다. 마지막 L2 이후
공통 scalar denominator는 방향에 영향을 주지 않는다.

## Gene matching policy

- Dataset feature string과 artifact key의 exact match만 허용한다.
- 대소문자 변경, fuzzy matching 또는 임의 alias substitution은 하지 않는다.
- 저자 methods가 HGNC aliases를 lookup에 추가했으므로, artifact exact key 중 primary
  `NCBI_summary_of_genes.json` key가 아닌 항목을 `official_alias`로 기록한다.
- Unmatched genes는 cell library-size normalization mass에는 포함하고 projection에는
  zero contribution을 갖는다.
- Dataset/embedding duplicate와 exact/alias/unmatched count, raw-count mass coverage
  six-number summary를 QC에 기록한다.

## Sparse interoperability and outputs

R exporter는 Seurat object를 메모리에서만 compatibility-update하고 `RNA/counts`를
genes × cells Matrix Market으로 sparse export한다. Python은 이를 cells × genes CSR로
읽고 batch-wise projection한다. 전체 expression matrix를 dense로 만들지 않는다.

```text
data/processed/genept_scpa/phase2/naive_cd4_genept_w.npy
data/processed/genept_scpa/phase2/naive_cd4_genept_w_cell_ids.txt
data/processed/genept_scpa/phase2/naive_cd4_genept_w_metadata.csv
data/processed/genept_scpa/phase2/genept_gene_mapping.csv
data/interim/genept_scpa/phase2_genept_w_qc.json
data/interim/genept_scpa/phase2_genept_w_summary.md
```

Matrix는 14,894 × 1,536 float32 `.npy`로 저장한다. Metadata에는 `cell_id`, `Hour`,
`Cell_Type`을 같은 순서로 보존한다. Source RDS는 수정하지 않는다.

## Correctness and gate

- 3 genes × 2 dimensions × 2 cells toy example을 hand calculation과 비교한다.
- Sparse optimized calculation을 explicit Python loop와 비교한다.
- Gene-order invariance와 batch-size determinism을 검사한다.
- Full run에서는 실제 3 cells의 optimized/direct 결과도 비교한다.
- All cells/IDs, 1,536 dimensions, finite values, no unexpected zero vector와
  `abs(L2 norm - 1) < 1e-5`를 확인한다.

모든 조건이 맞으면 `READY_FOR_GPT_REVIEW`, 아니면 `NEEDS_REVIEW`다. Phase 2는
결과를 검토한 뒤에도 자동으로 Phase 3로 넘어가지 않는다.

## User commands

```bash
PYTHONPATH=src python scripts/genept/prepare_genept_embeddings.py --download
PYTHONPATH=src python scripts/genept/build_genept_w.py --dataset naive_cd4
```

첫 명령은 약 574 MB official archive를 resumable download하고 checksum/schema를
검사한다. 두 번째 명령이 sparse Seurat export와 full GenePT-w/QC를 연속 수행한다.

실제 full run은 14,894 × 1,536 output, 14,409 matched genes, median expression
coverage 0.931967, zero vectors 0, all-finite output과 unit L2 norms를 확인했다. 모든
gate criterion이 true이고 warning/failed check가 없어 2026-08-11 review에서 PASS로
확정했다. Historical QC JSON의 `READY_FOR_GPT_REVIEW` 값은 실행 당시 상태로 보존한다.
