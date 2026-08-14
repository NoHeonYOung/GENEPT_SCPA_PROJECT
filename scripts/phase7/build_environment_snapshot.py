#!/usr/bin/env python3
"""Write a read-only Phase 7 software/hardware/Git environment snapshot."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from gene_embedding_project.genept_scpa.io import write_json_atomic  # noqa: E402
from gene_embedding_project.genept_scpa.phase7.runtime import build_runtime_report  # noqa: E402


PYTHON_PACKAGES = (
    "PyYAML", "numpy", "scipy", "h5py", "matplotlib", "torch",
    "transformers", "triton", "accelerate", "kernels", "openai-harmony",
)
R_PACKAGES = ("SCPA", "multicross", "Matrix", "jsonlite", "rhdf5", "Seurat", "SeuratObject")


def _version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=PROJECT_ROOT, check=True, capture_output=True, text=True
    ).stdout.strip()


def _r_snapshot() -> dict[str, Any]:
    expression = (
        "pkgs <- c(" + ",".join(f'\"{name}\"' for name in R_PACKAGES) + ");"
        "cat(R.version.string, '\\n');"
        "for (p in pkgs) cat(p, if (requireNamespace(p, quietly=TRUE)) "
        "as.character(packageVersion(p)) else 'NOT_INSTALLED', '\\n')"
    )
    completed = subprocess.run(
        ["Rscript", "-e", expression], check=True, capture_output=True, text=True
    )
    lines = completed.stdout.splitlines()
    versions: dict[str, str | None] = {}
    for line in lines[1:]:
        name, value = line.split(maxsplit=1)
        value = value.strip()
        versions[name] = None if value == "NOT_INSTALLED" else value
    return {"version": lines[0], "packages": versions}


def _os_release() -> dict[str, str]:
    path = Path("/etc/os-release")
    if not path.is_file():
        return {}
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", maxsplit=1)
            values[key] = value.strip().strip('"')
    return values


def build_snapshot(config_path: Path) -> dict[str, Any]:
    protocol = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    runtime = build_runtime_report(protocol)
    dirty_lines = _git("status", "--porcelain=v1").splitlines()
    return {
        "schema_version": "phase7.environment.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "repository_root": str(PROJECT_ROOT),
        "os": {"platform": platform.platform(), "release": _os_release()},
        "python": {
            "version": platform.python_version(),
            "executable": sys.executable,
            "packages": {name: _version(name) for name in PYTHON_PACKAGES},
        },
        "r": _r_snapshot(),
        "runtime_gate": runtime,
        "git": {
            "head_commit": _git("rev-parse", "HEAD"),
            "branch": _git("branch", "--show-current"),
            "dirty": bool(dirty_lines),
            "status_porcelain": dirty_lines,
        },
        "constraints": {
            "primary_min_vram_gib": 16,
            "cpu_offload": False,
            "automatic_backend_or_precision_fallback": False,
            "model_download_performed": False,
            "model_loaded": False,
            "real_inference_performed": False,
        },
    }


def markdown(snapshot: dict[str, Any]) -> str:
    runtime = snapshot["runtime_gate"]
    observed = runtime["snapshot"]
    classification = runtime["classification"]
    lines = [
        "# Phase 7 environment snapshot",
        "",
        f"Generated: `{snapshot['generated_at_utc']}`",
        "",
        "No model, scientific result, production synthetic data, or production SCPA was accessed.",
        "",
        "## System",
        "",
        f"- OS: `{snapshot['os']['platform']}`",
        f"- Distribution: `{snapshot['os']['release'].get('PRETTY_NAME')}`",
        f"- Python: `{snapshot['python']['version']}` (`{snapshot['python']['executable']}`)",
        f"- R: `{snapshot['r']['version']}`",
        f"- GPU: `{(observed.get('gpu') or {}).get('name')}`",
        f"- Compute capability: `{(observed.get('gpu') or {}).get('compute_capability')}`",
        f"- NVIDIA driver: `{(observed.get('gpu') or {}).get('driver_version')}`",
        f"- CUDA runtime: `{observed.get('cuda', {}).get('runtime')}`",
        f"- System RAM bytes: `{observed.get('system_ram_bytes')}`",
        f"- Disk free bytes: `{observed.get('disk', {}).get('free_bytes')}`",
        "",
        "## Python packages",
        "",
    ]
    lines.extend(
        f"- `{name}`: `{version if version is not None else 'NOT_INSTALLED'}`"
        for name, version in snapshot["python"]["packages"].items()
    )
    lines.extend(["", "## R packages", ""])
    lines.extend(
        f"- `{name}`: `{version if version is not None else 'NOT_INSTALLED'}`"
        for name, version in snapshot["r"]["packages"].items()
    )
    lines.extend([
        "",
        "## Runtime gate",
        "",
        f"- Status: **{classification['status']}**",
        f"- Primary supported: `{classification['primary_supported']}`",
        f"- Reasons: `{classification['primary_failure_reasons']}`",
        f"- Expected runtime: `{observed['expected_runtime_mode']}`",
        f"- Model/network access performed: `{observed['network_or_model_access_performed']}`",
        "",
        "## Git",
        "",
        f"- Commit at snapshot: `{snapshot['git']['head_commit']}`",
        f"- Branch: `{snapshot['git']['branch']}`",
        f"- Dirty: `{snapshot['git']['dirty']}`",
        "- Porcelain status:",
        "",
        "```text",
        *snapshot["git"]["status_porcelain"],
        "```",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path,
        default=PROJECT_ROOT / "config/phase7_gpt_oss_synthetic.yaml",
    )
    parser.add_argument(
        "--json-output", type=Path,
        default=PROJECT_ROOT / "artifacts/phase7_environment_snapshot.json",
    )
    parser.add_argument(
        "--markdown-output", type=Path,
        default=PROJECT_ROOT / "docs/phase7_environment_snapshot.md",
    )
    args = parser.parse_args()
    snapshot = build_snapshot(args.config)
    write_json_atomic(snapshot, args.json_output)
    args.markdown_output.parent.mkdir(parents=True, exist_ok=True)
    args.markdown_output.write_text(markdown(snapshot), encoding="utf-8")
    print(
        f"PHASE7_ENVIRONMENT_SNAPSHOT status=PASS "
        f"runtime={snapshot['runtime_gate']['classification']['status']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
