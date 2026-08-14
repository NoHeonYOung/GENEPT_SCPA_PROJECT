# Phase 7 environment snapshot

Generated: `2026-08-14T03:14:10.053130+00:00`

No model, scientific result, production synthetic data, or production SCPA was accessed.

## System

- OS: `Linux-5.4.0-216-generic-x86_64-with-glibc2.31`
- Distribution: `Ubuntu 20.04.6 LTS`
- Python: `3.9.18` (`/home/node00/anaconda3/bin/python`)
- R: `R version 4.5.2 (2025-10-31) `
- GPU: `NVIDIA GeForce RTX 2080 SUPER`
- Compute capability: `[7, 5]`
- NVIDIA driver: `535.230.02`
- CUDA runtime: `12.8`
- System RAM bytes: `135056318464`
- Disk free bytes: `4236722176`

## Python packages

- `PyYAML`: `6.0.1`
- `numpy`: `1.26.4`
- `scipy`: `1.12.0`
- `h5py`: `3.1.0`
- `matplotlib`: `3.7.5`
- `torch`: `2.8.0`
- `transformers`: `4.57.6`
- `triton`: `3.4.0`
- `accelerate`: `NOT_INSTALLED`
- `kernels`: `NOT_INSTALLED`
- `openai-harmony`: `NOT_INSTALLED`

## R packages

- `SCPA`: `1.6.2`
- `multicross`: `2.1.0`
- `Matrix`: `1.7.6`
- `jsonlite`: `2.0.0`
- `rhdf5`: `2.54.1`
- `Seurat`: `5.5.1`
- `SeuratObject`: `5.4.0`

## Runtime gate

- Status: **UNSUPPORTED_PRIMARY**
- Primary supported: `False`
- Reasons: `['VRAM 8.00GiB is below primary minimum 16GiB', 'free disk 3.95GiB is below 30GiB', 'required package is missing: accelerate', 'required package is missing: kernels']`
- Expected runtime: `transformers_mxfp4_single_cuda_no_cpu_offload`
- Model/network access performed: `False`

## Git

- Commit at snapshot: `5f771c359e8ac47e688fd14ee01effffd4e6065c`
- Branch: `main`
- Dirty: `True`
- Porcelain status:

```text
M README.md
 M config/genept_scpa.yaml
 M docs/genept_scpa_decision_log.md
 M genept_scpa_experiment_plan.md
 M src/gene_embedding_project/genept_scpa/config.py
 M tests/test_config.py
?? artifacts/
?? config/phase7_gpt_oss_synthetic.yaml
?? data/processed/genept_scpa/phase5_gene_contribution/comparisons/
?? data/processed/genept_scpa/phase5_gene_contribution/figures/
?? data/processed/genept_scpa/phase5_gene_contribution/phase5_gene_contribution_manifest.json
?? data/processed/genept_scpa/phase5_gene_contribution/phase5_gene_contribution_qc.json
?? data/processed/genept_scpa/phase5_gene_contribution/phase5_gene_contribution_summary.md
?? data/processed/genept_scpa/phase5_gene_contribution/phase5_gene_masking_all_results.csv
?? data/processed/genept_scpa/phase5_gene_contribution/phase5_gene_rank_comparison.csv
?? data/processed/genept_scpa/phase5_gene_contribution/phase5_pathway_summary.csv
?? data/processed/genept_scpa/phase6_semantic_controls/
?? docs/phase7_environment_snapshot.md
?? docs/phase7_git_boundaries.md
?? docs/phase7_gpt_oss_synthetic_protocol.md
?? docs/phase7_home_server_handoff.md
?? docs/phase7_transfer_manifest.md
?? scripts/phase6/
?? scripts/phase7/
?? scripts/scpa/run_phase6_semantic_controls_core.R
?? scripts/scpa/run_phase7_synthetic_masking_core.R
?? src/gene_embedding_project/genept_scpa/phase7/
?? tests/test_phase6_semantic_controls.py
?? tests/test_phase7_cohort.py
?? tests/test_phase7_execution_gates.py
?? tests/test_phase7_gpt_oss_backend.py
?? tests/test_phase7_gpt_oss_runtime.py
?? tests/test_phase7_llm_pipeline.py
?? tests/test_phase7_ranking_metrics.py
?? tests/test_phase7_scpa_masking.R
?? tests/test_phase7_synthetic.py
?? tests/test_phase7_toy_smoke.py
?? tests/test_phase7_transfer_integrity.py
```
