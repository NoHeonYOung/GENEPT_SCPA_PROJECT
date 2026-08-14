# Phase 3 preparation - naïve CD8 GenePT-w reproduction summary

Gate status: `READY_FOR_GPT_REVIEW`

## Official source and primary embedding

- Final paper: Chen & Zou, Nature Biomedical Engineering 9, 483–493 (2025).
- Official code: `yiqunchen/GenePT`, including `aorta_data_analysis.ipynb`.
- Artifact: Zenodo DOI `10.5281/zenodo.10833191`.
- Primary model: `text-embedding-ada-002`, 1,536 dimensions.
- OpenAI API calls: none; the author-provided precomputed artifact was used.

## Input and published preprocessing

- Source: `/home/node00/nhy_python/GenePT_SCPA/data/raw/genept_scpa/GSE212270_integrated_naive_cd8.rds` (`RNA/counts`, 7041 cells x 17942 genes).
- Raw counts -> cell-wise total-count normalization to 10,000 -> log1p.
- Dataset genes are then aligned by exact official artifact key; no fuzzy or case conversion.
- Expression-weighted GenePT aggregation uses the official notebook denominator (all dataset genes), followed by row-wise unit L2 normalization.
- Unmatched genes participate in library-size normalization and then contribute a zero embedding vector.

## Gene matching and expression coverage

- Embedding keys: 93800
- Exact primary-symbol matches: 14226
- Official HGNC-alias-key matches: 321
- Unmatched dataset genes: 3395
- Duplicate mappings: 0
- Raw-count mass coverage summary: `{'min': 0.897535139539621, 'q1': 0.9160095989029825, 'median': 0.924673784104389, 'mean': 0.9251688761737864, 'q3': 0.9336362161285089, 'max': 0.9951311874492832}`
- CD8-minus-CD4 median coverage: -0.00729303 (recorded as a potential confounder; no exclusion threshold applied).

## Output and correctness

- Matrix: `/home/node00/nhy_python/GenePT_SCPA/data/processed/genept_scpa/phase3/naive_cd8_genept_w.npy` (7041 x 1536, float32).
- Finite values: 10814976/10814976
- Zero vectors: 0
- Post-L2 norm summary: `{'min': 0.9999999999999997, 'q1': 0.9999999999999999, 'median': 1.0, 'mean': 1.0, 'q3': 1.0, 'max': 1.0000000000000002}`
- Optimized-vs-direct maximum absolute error: 1.67e-08
- Deterministic repeat maximum absolute error: 0

## Compatibility and scope

The official notebook loads an already-preprocessed AnnData `.X`; this pipeline makes the paper's explicit normalization/log1p and final L2 steps reproducible from Seurat raw counts. The notebook's lookup/aggregation rule is retained.

No SCPA, CD4-vs-CD8 comparison, classifier, separability metric, or pathway/gene-level interpretation was run by this builder.
