# Phase 7 Git boundaries

This file records the original selective Phase 7 handoff boundary. The later
repository-portability audit intentionally extends the Git boundary to include
small Phase 1--6 results, provenance, frozen IDs and export sidecars. Do not use
`git add .`, `git add -A`, or a broad `git commit -a`.

## Phase 7 files safe to commit

- `artifacts/phase7_environment_snapshot.json`
- `artifacts/phase7_transfer_manifest.json`
- `config/phase7_gpt_oss_synthetic.yaml`
- `docs/phase7_environment_snapshot.md`
- `docs/phase7_git_boundaries.md`
- `docs/phase7_gpt_oss_synthetic_protocol.md`
- `docs/phase7_home_server_handoff.md`
- `docs/phase7_transfer_manifest.md`
- `scripts/phase7/*.py`
- `scripts/scpa/run_phase7_synthetic_masking_core.R`
- `src/gene_embedding_project/genept_scpa/phase7/*.py`
- `tests/test_phase7*.py`
- `tests/test_phase7_scpa_masking.R`

## Shared files intentionally modified for Phase 7

- `config/genept_scpa.yaml`
- `docs/genept_scpa_decision_log.md`
- `genept_scpa_experiment_plan.md`
- `src/gene_embedding_project/genept_scpa/config.py`
- `tests/test_config.py`

These shared changes record Phase 6 closure, unlock/freeze Phase 7 and validate the
new active phase. They are part of the Phase 7 transition and handoff.

## Files reviewed by the repository-portability audit

- `README.md`
- `data/processed/genept_scpa/phase5_gene_contribution/**`
- `data/processed/genept_scpa/phase6_semantic_controls/**`
- `scripts/phase6/**`
- `scripts/scpa/run_phase6_semantic_controls_core.R`
- `tests/test_phase6_semantic_controls.py`

These were previously left out only because they pre-dated the Phase 7 selective
commit. They are small, non-secret project state and are included by the later
portability commit after fresh-clone validation.

The CD4 counts matrix and GenePT embedding remain separate. The small CD4 export
sidecars are Git-tracked. Exact hashes and destination paths are recorded in
`artifacts/phase7_transfer_manifest.json`.
