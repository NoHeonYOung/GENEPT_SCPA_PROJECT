# Phase 7 home/server handoff

This checklist prepares a destination machine and **stops before model download,
model loading, inference, production synthetic generation, or production SCPA**.

The frozen primary scientific runtime remains:

- `openai/gpt-oss-20b`
- Hugging Face Transformers with pretrained MXFP4
- one CUDA GPU with at least 16GiB VRAM
- no CPU offload
- no automatic backend or precision fallback

## A. Research-lab machine — completed

- Phase 6 is PASS/COMPLETED with its reviewed negative semantic-specific result.
- Phase 7 biological protocol, 11-pathway universe, synthetic scenarios, metrics,
  masking definitions and LLM controls are frozen.
- Transformers backend is lazy, local-files-only and blocked by execution gates.
- Transfer manifest/checksums and environment snapshot are recorded.
- Production split count is intentionally not selected.
- No gpt-oss weights or scientific Phase 7 results were generated.

Stop work on the research-lab machine after the final test, transfer-integrity,
runtime-gate and Git-status checks. Do not make a data copy as part of the commit.

## B. Home machine setup — allowed next steps only

### 1. Restore the repository

Clone or pull the Phase 7 handoff commit and enter the repository root.

```bash
git clone <repository-url> GenePT_SCPA
cd GenePT_SCPA
```

If the repository already exists, use the normal reviewed pull workflow. Do not
merge local result directories into Git.

### 2. Restore separately transferred data

Place every separately transferred file at the exact repository-relative path
listed in `artifacts/phase7_transfer_manifest.json`. The required non-Git inputs
are the CD4 sparse export files and GenePT embedding artifact. Do not rename files,
edit the export manifest, or guess alternate locations.

The optional source RDS is provenance only. Phase 7 consumes the validated sparse
RNA/counts export.

### 3. Verify checksums before any analysis

```bash
PYTHONPATH=src python scripts/phase7/check_transfer_integrity.py \
  --manifest artifacts/phase7_transfer_manifest.json \
  --root "$PWD"
```

Required result: `PHASE7_TRANSFER_INTEGRITY status=PASS`. If it fails, stop and
restore the exact missing/corrupt file. The checker never downloads, writes, or
searches alternate paths.

### 4. Create and activate the Python environment

Use the exact versions recorded in `artifacts/phase7_environment_snapshot.json`
as the starting reproducibility target. A minimal isolated environment may be
created as follows:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

Install the CUDA-compatible PyTorch/Transformers/MXFP4 runtime dependencies only
as a reviewed environment setup step. Do not install or download model weights.
Do not enable CPU offload or another quantization/backend to make the gate pass.

### 5. Verify R dependencies

```bash
Rscript -e 'pkgs <- c("SCPA","multicross","Matrix","jsonlite","rhdf5"); for (p in pkgs) cat(p, if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else "NOT_INSTALLED", "\n")'
```

Compare versions with `artifacts/phase7_environment_snapshot.json`. Do not run
SCPA at this stage.

### 6. Run the full unit tests

```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test_*.py'
Rscript tests/test_phase7_scpa_masking.R
git diff --check
```

The R test is algebra-only and makes no SCPA/MCM call.

### 7. Run the model-free runtime gate

```bash
PYTHONPATH=src python scripts/phase7/check_gpt_oss_runtime.py
```

This command inspects hardware/software only. It must report
`network_or_model_access_performed=false`.

### 8. STOP and review the gate

Do not download gpt-oss and do not run `run_gpt_oss_inference.py`.

The planned home machine has RTX 5070 12GB, 32GB RAM and about 300GB free SSD.
It is expected to fail the frozen primary **>=16GiB VRAM** requirement. This is
not a protocol error and the 16GiB requirement must not be relaxed.

After reviewing the actual home-machine gate, make one new explicit decision:

- A. Preserve the primary protocol and use a >=16–24GB GPU for production; or
- B. Define a separate secondary runtime-only pilot backend for RTX 5070 12GB,
  without changing or replacing the frozen primary scientific backend.

## C. Future runtime-only pilot — not authorized yet

Only after a separate decision may a pilot test one split, two frozen pathways,
stats-only and one candidate order. It may expose load success, VRAM/RAM, latency,
tokens/sec, JSON validity and retry/parser behavior only.

It must not expose Recall, AP, NDCG, truth-gene positions, method comparisons or
scientific ranking interpretation. Production split seeds/count must be frozen
from runtime measurements before scientific metrics are inspected.

## D. Future scientific production — not authorized yet

Scientific production requires all of the following:

- an explicitly approved runtime decision;
- a `SUPPORTED_PRIMARY` gate for the frozen primary backend;
- a pinned local model and tokenizer revision;
- production pseudo-split seeds/count frozen before metrics;
- separately unlocked production synthetic, SCPA and real-inference gates.

Until then, do not change any execution-gate boolean and do not invoke model or
production runners.
