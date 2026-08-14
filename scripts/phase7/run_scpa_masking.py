#!/usr/bin/env python3
"""Run Phase 7 SCPA masking only after the production gate is explicitly unlocked."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def run(
    config_path: Path, input_h5: Path, manifest: Path,
    checkpoint_dir: Path, output_json: Path,
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not config["execution_gate"]["production_scpa_allowed"]:
        raise RuntimeError("Phase 7 production SCPA remains locked")
    command = [
        "Rscript", str(PROJECT_ROOT / "scripts/scpa/run_phase7_synthetic_masking_core.R"),
        "--input-h5", str(input_h5), "--manifest", str(manifest),
        "--checkpoint-dir", str(checkpoint_dir), "--output-json", str(output_json),
    ]
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=PROJECT_ROOT / "config/phase7_gpt_oss_synthetic.yaml")
    parser.add_argument("--input-h5", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--checkpoint-dir", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()
    run(args.config, args.input_h5, args.manifest, args.checkpoint_dir, args.output_json)
