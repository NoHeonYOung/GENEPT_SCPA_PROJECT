#!/usr/bin/env python3
"""Run the frozen Transformers backend only after all explicit gates are unlocked."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from gene_embedding_project.genept_scpa.io import write_json_atomic  # noqa: E402
from gene_embedding_project.genept_scpa.phase7.gpt_oss_backend import (  # noqa: E402
    TransformersGPTOSSBackend,
)
from gene_embedding_project.genept_scpa.phase7.runtime import (  # noqa: E402
    SUPPORTED_PRIMARY, build_runtime_report,
)


def run(
    config_path: Path,
    request_dir: Path,
    trace_dir: Path,
    invalid_raw_dir: Path,
    runtime_report_path: Path,
    *,
    model_path: Path | None = None,
    cache_path: Path | None = None,
) -> int:
    protocol = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not protocol["execution_gate"]["real_llm_inference_allowed"]:
        raise RuntimeError("Real gpt-oss inference remains locked")
    report = build_runtime_report(protocol, cache_path)
    write_json_atomic(report, runtime_report_path)
    if report["classification"]["status"] != SUPPORTED_PRIMARY:
        raise RuntimeError(
            f"Primary runtime gate failed: {report['classification']['primary_failure_reasons']}"
        )
    backend = TransformersGPTOSSBackend(
        protocol,
        model_path=model_path,
        cache_path=cache_path,
        execution_authorized=True,
        runtime_report=report,
    )
    trace_dir.mkdir(parents=True, exist_ok=True)
    invalid_raw_dir.mkdir(parents=True, exist_ok=True)
    completed = 0
    invalid = 0
    for request_path in sorted(request_dir.glob("*.json")):
        request = json.loads(request_path.read_text(encoding="utf-8"))
        trace = backend.rank_with_trace(request, invalid_raw_directory=invalid_raw_dir)
        write_json_atomic(trace, trace_dir / request_path.name)
        if trace["status"] == "PASS":
            completed += 1
        else:
            invalid += 1
        print(
            f"[Phase7 GPT-OSS] completed={completed} invalid={invalid} "
            f"run_id={request['run_id']} attempts={trace['attempt_count']}"
        )
    if completed + invalid == 0:
        raise ValueError("No Phase 7 LLM requests found")
    if invalid:
        raise RuntimeError(f"Strict gpt-oss output failed for {invalid} request(s)")
    return completed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path,
        default=PROJECT_ROOT / "config/phase7_gpt_oss_synthetic.yaml",
    )
    parser.add_argument("--request-dir", type=Path, required=True)
    parser.add_argument("--trace-dir", type=Path, required=True)
    parser.add_argument("--invalid-raw-dir", type=Path, required=True)
    parser.add_argument("--runtime-report", type=Path, required=True)
    parser.add_argument("--model-path", type=Path)
    parser.add_argument("--cache-path", type=Path)
    args = parser.parse_args()
    completed = run(
        args.config, args.request_dir, args.trace_dir, args.invalid_raw_dir,
        args.runtime_report, model_path=args.model_path, cache_path=args.cache_path,
    )
    print(f"PHASE7_GPT_OSS_INFERENCE status=PASS calls={completed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
