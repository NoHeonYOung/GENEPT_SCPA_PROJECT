#!/usr/bin/env python3
"""Report the Phase 7 gpt-oss runtime gate without model/network access."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from gene_embedding_project.genept_scpa.io import write_json_atomic  # noqa: E402
from gene_embedding_project.genept_scpa.phase7.runtime import (  # noqa: E402
    build_runtime_report, report_json,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path,
        default=PROJECT_ROOT / "config/phase7_gpt_oss_synthetic.yaml",
    )
    parser.add_argument("--cache-path", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    report = build_runtime_report(config, args.cache_path)
    if args.output:
        write_json_atomic(report, args.output)
    print(report_json(report))
    print(f"PHASE7_GPT_OSS_RUNTIME status={report['classification']['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
