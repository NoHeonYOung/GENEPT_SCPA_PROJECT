#!/usr/bin/env python3
"""Verify the exact Phase 7 transfer manifest without writes or downloads."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from gene_embedding_project.genept_scpa.phase7.transfer import verify_transfer_manifest  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path,
        default=PROJECT_ROOT / "artifacts/phase7_transfer_manifest.json",
    )
    parser.add_argument("--root", type=Path, default=PROJECT_ROOT)
    args = parser.parse_args()
    report = verify_transfer_manifest(args.manifest, args.root)
    for row in report["results"]:
        print(
            f"[{row['status']}] required={str(row['required']).lower()} "
            f"path={row['repository_relative_path']}"
        )
    print(
        f"PHASE7_TRANSFER_INTEGRITY status={report['status']} "
        f"required={report['required_passed']}/{report['required_count']}"
    )
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
