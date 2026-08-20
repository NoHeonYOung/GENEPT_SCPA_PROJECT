#!/usr/bin/env python3
"""Run resumable Phase 8 PERMUTED/RANDOM SCPA masking controls."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path,
                        default=PROJECT_ROOT / "config/phase8_mean_shift_mechanism.yaml")
    parser.add_argument("--cores", type=int, default=12)
    parser.add_argument("--max-experiments", type=int, default=0,
                        help="Use 1 for smoke; 0 runs all 660 experiments")
    parser.add_argument("--progress-every-genes", type=int, default=10)
    args = parser.parse_args()
    if args.cores < 1 or args.max_experiments < 0:
        parser.error("--cores must be positive and --max-experiments non-negative")
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    artifacts = config["artifacts"]
    interim = PROJECT_ROOT / artifacts["interim_directory"]
    processed = PROJECT_ROOT / artifacts["processed_directory"]
    qc_name = artifacts["smoke_qc"] if args.max_experiments else artifacts["control_qc"]
    command = [
        "Rscript", str(PROJECT_ROOT / "scripts/scpa/run_phase8_mean_shift_controls_core.R"),
        "--expression-h5", str(PROJECT_ROOT / config["source"]["phase7_expression_h5"]),
        "--controls-h5", str(interim / artifacts["controls_h5"]),
        "--manifest", str(processed / artifacts["manifest"]),
        "--checkpoint-dir", str(interim / artifacts["checkpoints"]),
        "--output-json", str(processed / qc_name),
        "--cores", str(args.cores),
        "--max-experiments", str(args.max_experiments),
        "--progress-every-genes", str(args.progress_every_genes),
    ]
    print(
        f"[Phase 8] Starting controls: expected PERMUTED={config['workload']['expected_permuted_mcm']} "
        f"RANDOM={config['workload']['expected_random_mcm']} cores={args.cores}",
        flush=True,
    )
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
