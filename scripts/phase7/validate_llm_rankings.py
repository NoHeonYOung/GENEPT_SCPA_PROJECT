#!/usr/bin/env python3
"""Validate Phase 7 mock/future LLM ranking artifacts without ground truth."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from gene_embedding_project.genept_scpa.phase7.schemas import validate_llm_response  # noqa: E402


def validate_directory(directory: Path) -> int:
    count = 0
    for path in sorted(directory.glob("*.json")):
        artifact = json.loads(path.read_text(encoding="utf-8"))
        response = artifact.get("response", artifact)
        validate_llm_response(response)
        if artifact.get("backend") and artifact["backend"] != response["backend"]:
            raise ValueError(f"Backend mismatch in {path}")
        count += 1
    if count == 0:
        raise ValueError("No ranking artifacts found")
    return count


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--response-dir", type=Path, required=True)
    args = parser.parse_args()
    total = validate_directory(args.response_dir)
    print(f"PHASE7_LLM_VALIDATION status=PASS rankings={total}")
