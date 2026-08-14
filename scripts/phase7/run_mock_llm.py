#!/usr/bin/env python3
"""Run only the deterministic Phase 7 mock backend over prepared requests."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from gene_embedding_project.genept_scpa.io import write_json_atomic  # noqa: E402
from gene_embedding_project.genept_scpa.phase7.llm_backend import MockLLMBackend  # noqa: E402


def run_mock(request_dir: Path, output_dir: Path) -> int:
    backend = MockLLMBackend()
    output_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    for request_path in sorted(request_dir.glob("*.json")):
        request = json.loads(request_path.read_text(encoding="utf-8"))
        response = backend.rank(request)
        artifact = {
            "backend": backend.name,
            "scientific_evaluation_allowed": backend.scientific_evaluation_allowed,
            "request_file": str(request_path),
            "response": response,
        }
        write_json_atomic(artifact, output_dir / request_path.name)
        count += 1
    return count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    total = run_mock(args.request_dir, args.output_dir)
    print(f"PHASE7_MOCK_LLM status=PASS backend={MockLLMBackend.name} responses={total} scientific=false")
