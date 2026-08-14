"""Read-only verification of exact Phase 7 transfer-manifest paths."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from gene_embedding_project.genept_scpa.io import sha256_file


TRANSFER_SCHEMA_VERSION = "phase7.transfer.v1"


def verify_transfer_manifest(
    manifest_path: str | Path, repository_root: str | Path
) -> dict[str, Any]:
    """Verify existence, size and SHA256 without modifying or relocating anything."""

    manifest_file = Path(manifest_path)
    root = Path(repository_root).resolve()
    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != TRANSFER_SCHEMA_VERSION:
        raise ValueError("Unexpected Phase 7 transfer manifest schema")
    resources = manifest.get("resources")
    if not isinstance(resources, list) or not resources:
        raise ValueError("Transfer manifest has no resources")
    results: list[dict[str, Any]] = []
    failures: list[str] = []
    required_count = 0
    for resource in resources:
        relative = resource.get("repository_relative_path")
        if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
            raise ValueError("Every transfer resource needs one exact repository-relative path")
        required = resource.get("required_status") == "required"
        required_count += int(required)
        path = root / relative
        row: dict[str, Any] = {
            "id": resource.get("id"),
            "repository_relative_path": relative,
            "required": required,
            "exists": path.is_file(),
            "size_matches": False,
            "sha256_matches": False,
            "status": "FAIL" if required else "OPTIONAL_MISSING",
        }
        if path.is_file():
            row["observed_size_bytes"] = path.stat().st_size
            row["size_matches"] = path.stat().st_size == int(resource["size_bytes"])
            if row["size_matches"]:
                observed_hash = sha256_file(path)
                row["observed_sha256"] = observed_hash
                row["sha256_matches"] = observed_hash == resource["sha256"]
            row["status"] = (
                "PASS" if row["size_matches"] and row["sha256_matches"] else "FAIL"
            )
        if required and row["status"] != "PASS":
            failures.append(f"{relative}: {row['status']}")
        results.append(row)
    return {
        "status": "PASS" if not failures else "FAIL",
        "manifest": str(manifest_file.resolve()),
        "repository_root": str(root),
        "resource_count": len(results),
        "required_count": required_count,
        "required_passed": sum(row["required"] and row["status"] == "PASS" for row in results),
        "failures": failures,
        "results": results,
        "files_modified": False,
        "downloads_performed": False,
        "alternate_paths_inferred": False,
    }
