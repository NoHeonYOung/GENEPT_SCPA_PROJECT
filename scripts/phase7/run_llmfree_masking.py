#!/usr/bin/env python3
"""Launch resumable Phase 7 Vanilla/GenePT masking; no LLM or GPU is used."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path,
                        default=PROJECT_ROOT / "config/phase7_llmfree_synthetic.yaml")
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--max-experiments", type=int, default=0,
                        help="Use 1 for a pre-production smoke; 0 runs all experiments")
    parser.add_argument("--progress-every-genes", type=int, default=10)
    args = parser.parse_args()
    if args.cores < 1 or args.max_experiments < 0:
        parser.error("--cores must be positive and --max-experiments non-negative")
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if not config["execution"]["llm_backend_forbidden"]:
        raise RuntimeError("Refusing a Phase 7 configuration that permits an LLM backend")
    artifacts = config["artifacts"]
    interim = PROJECT_ROOT / artifacts["interim_directory"]
    processed = PROJECT_ROOT / artifacts["processed_directory"]
    qc_name = artifacts["smoke_qc"] if args.max_experiments > 0 else artifacts["core_qc"]
    command = [
        "Rscript", str(PROJECT_ROOT / "scripts/scpa/run_phase7_llmfree_masking_core.R"),
        "--input-h5", str(interim / artifacts["expression_h5"]),
        "--manifest", str(processed / artifacts["manifest"]),
        "--checkpoint-dir", str(interim / artifacts["checkpoints"]),
        "--output-json", str(processed / qc_name),
        "--cores", str(args.cores),
        "--max-experiments", str(args.max_experiments),
        "--progress-every-genes", str(args.progress_every_genes),
    ]
    print("[Phase 7] Starting CPU-only LLM-free SCPA masking", flush=True)
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
