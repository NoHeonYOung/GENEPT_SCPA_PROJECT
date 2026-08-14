"""Read-only gpt-oss runtime capability inspection and fail-closed gating."""

from __future__ import annotations

import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
from typing import Any, Mapping


SUPPORTED_PRIMARY = "SUPPORTED_PRIMARY"
UNSUPPORTED_PRIMARY = "UNSUPPORTED_PRIMARY"
PILOT_ONLY = "PILOT_ONLY"


def default_hf_cache_path() -> Path:
    configured = os.environ.get("HF_HOME")
    return Path(configured).expanduser() if configured else Path.home() / ".cache" / "huggingface"


def _package_version(distribution: str) -> str | None:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None


def _version_tuple(value: str | None) -> tuple[int, ...]:
    if not value:
        return tuple()
    return tuple(int(part) for part in re.findall(r"\d+", value)[:3])


def _existing_parent(path: Path) -> Path:
    candidate = path.expanduser().resolve(strict=False)
    while not candidate.exists() and candidate != candidate.parent:
        candidate = candidate.parent
    return candidate


def _gpu_snapshot() -> tuple[dict[str, Any] | None, str | None]:
    command = [
        "nvidia-smi",
        "--query-gpu=name,memory.total,memory.free,driver_version,compute_cap",
        "--format=csv,noheader,nounits",
    ]
    try:
        completed = subprocess.run(command, check=True, capture_output=True, text=True, timeout=15)
    except (FileNotFoundError, subprocess.SubprocessError) as error:
        return None, f"{type(error).__name__}: {error}"
    rows = [row.strip() for row in completed.stdout.splitlines() if row.strip()]
    if len(rows) != 1:
        return None, f"Primary protocol requires exactly one visible GPU; found {len(rows)}"
    fields = [field.strip() for field in rows[0].split(",")]
    if len(fields) != 5:
        return None, f"Unexpected nvidia-smi row: {rows[0]}"
    major, minor = (int(value) for value in fields[4].split(".", maxsplit=1))
    return {
        "name": fields[0],
        "vram_total_bytes": int(fields[1]) * 1024 * 1024,
        "vram_free_bytes": int(fields[2]) * 1024 * 1024,
        "driver_version": fields[3],
        "compute_capability": [major, minor],
    }, None


def collect_runtime_snapshot(cache_path: str | Path | None = None) -> dict[str, Any]:
    """Inspect environment metadata without touching Hugging Face network APIs."""

    resolved_cache = Path(cache_path).expanduser() if cache_path else default_hf_cache_path()
    disk_anchor = _existing_parent(resolved_cache)
    disk = shutil.disk_usage(disk_anchor)
    packages = {
        name: _package_version(name)
        for name in ("torch", "transformers", "accelerate", "openai-harmony", "kernels", "triton")
    }
    cuda: dict[str, Any] = {"available": False, "runtime": None, "device_count": 0}
    torch_error: str | None = None
    if packages["torch"] is not None:
        try:
            torch = importlib.import_module("torch")
            cuda = {
                "available": bool(torch.cuda.is_available()),
                "runtime": torch.version.cuda,
                "device_count": int(torch.cuda.device_count()),
            }
        except Exception as error:  # pragma: no cover - hardware/driver dependent
            torch_error = f"{type(error).__name__}: {error}"
    gpu, gpu_error = _gpu_snapshot()
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "gpu": gpu,
        "gpu_probe_error": gpu_error,
        "system_ram_bytes": int(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")),
        "disk": {
            "probe_path": str(disk_anchor),
            "total_bytes": int(disk.total),
            "free_bytes": int(disk.free),
        },
        "cuda": cuda,
        "torch_probe_error": torch_error,
        "packages": packages,
        "model_cache_path": str(resolved_cache.resolve(strict=False)),
        "model_cache_exists": resolved_cache.exists(),
        "network_or_model_access_performed": False,
        "expected_runtime_mode": "transformers_mxfp4_single_cuda_no_cpu_offload",
    }


def classify_runtime(
    snapshot: Mapping[str, Any], requirements: Mapping[str, Any]
) -> dict[str, Any]:
    """Classify the frozen primary and pilot classes from a supplied snapshot."""

    gib = 1024**3
    gpu = snapshot.get("gpu")
    packages = snapshot.get("packages", {})
    primary_reasons: list[str] = []
    pilot_reasons: list[str] = []
    if not isinstance(gpu, Mapping):
        primary_reasons.append("exactly one NVIDIA GPU was not detected")
        pilot_reasons.append("pilot CUDA GPU was not detected")
        capability = (0, 0)
        vram_gib = 0.0
    else:
        capability = tuple(int(value) for value in gpu.get("compute_capability", (0, 0)))
        vram_gib = float(gpu.get("vram_total_bytes", 0)) / gib
    minimum_capability = tuple(int(value) for value in requirements["min_compute_capability"])
    if capability < minimum_capability:
        message = f"compute capability {capability} is below {minimum_capability}"
        primary_reasons.append(message)
        pilot_reasons.append(message)
    if vram_gib < float(requirements["min_vram_gib_primary"]):
        primary_reasons.append(
            f"VRAM {vram_gib:.2f}GiB is below primary minimum "
            f"{requirements['min_vram_gib_primary']}GiB"
        )
    if vram_gib < float(requirements["min_vram_gib_pilot_class"]):
        pilot_reasons.append(
            f"VRAM {vram_gib:.2f}GiB is below pilot-class minimum "
            f"{requirements['min_vram_gib_pilot_class']}GiB"
        )
    ram_gib = float(snapshot.get("system_ram_bytes", 0)) / gib
    if ram_gib < float(requirements["min_system_ram_gib"]):
        message = f"RAM {ram_gib:.2f}GiB is below {requirements['min_system_ram_gib']}GiB"
        primary_reasons.append(message)
        pilot_reasons.append(message)
    free_gib = float(snapshot.get("disk", {}).get("free_bytes", 0)) / gib
    if free_gib < float(requirements["min_free_disk_gib"]):
        message = f"free disk {free_gib:.2f}GiB is below {requirements['min_free_disk_gib']}GiB"
        primary_reasons.append(message)
        pilot_reasons.append(message)
    if not bool(snapshot.get("cuda", {}).get("available")):
        primary_reasons.append("torch CUDA is unavailable")
        pilot_reasons.append("torch CUDA is unavailable")
    for package in requirements["required_packages"]:
        if not packages.get(package):
            message = f"required package is missing: {package}"
            primary_reasons.append(message)
            pilot_reasons.append(message)
    if _version_tuple(packages.get("triton")) < _version_tuple(requirements["min_triton_version"]):
        message = f"triton must be >= {requirements['min_triton_version']}"
        primary_reasons.append(message)
        pilot_reasons.append(message)

    if not primary_reasons:
        status = SUPPORTED_PRIMARY
    elif not pilot_reasons:
        status = PILOT_ONLY
    else:
        status = UNSUPPORTED_PRIMARY
    return {
        "status": status,
        "primary_supported": status == SUPPORTED_PRIMARY,
        "pilot_class_supported": status in {SUPPORTED_PRIMARY, PILOT_ONLY},
        "primary_failure_reasons": primary_reasons,
        "pilot_failure_reasons": pilot_reasons,
        "observed_vram_gib": vram_gib,
        "observed_ram_gib": ram_gib,
        "observed_free_disk_gib": free_gib,
    }


def build_runtime_report(config: Mapping[str, Any], cache_path: str | Path | None = None) -> dict[str, Any]:
    snapshot = collect_runtime_snapshot(cache_path)
    classification = classify_runtime(snapshot, config["inference"]["runtime_requirements"])
    return {
        "phase": 7,
        "model": config["inference"]["model"],
        "backend": config["inference"]["backend_name"],
        "snapshot": snapshot,
        "classification": classification,
        "production_inference_authorized": bool(
            config["execution_gate"]["real_llm_inference_allowed"]
        ),
        "model_download_authorized": bool(
            config["execution_gate"]["gpt_oss_download_allowed"]
        ),
    }


def report_json(report: Mapping[str, Any]) -> str:
    return json.dumps(report, indent=2, ensure_ascii=False)
