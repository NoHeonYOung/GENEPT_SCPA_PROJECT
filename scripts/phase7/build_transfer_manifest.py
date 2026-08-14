#!/usr/bin/env python3
"""Build the small Phase 7 transfer inventory; never copies source resources."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from gene_embedding_project.genept_scpa.io import sha256_file, write_json_atomic  # noqa: E402
from gene_embedding_project.genept_scpa.phase7.transfer import TRANSFER_SCHEMA_VERSION  # noqa: E402


DATA_RESOURCES = [
    ("cd4_export_manifest", "data/interim/genept_scpa/phase2_export/naive_cd4/naive_cd4_export_manifest.json", "Frozen RNA/counts export dimensions and component hashes", "required"),
    ("cd4_counts", "data/interim/genept_scpa/phase2_export/naive_cd4/naive_cd4_rna_counts_genes_by_cells.mtx", "GSE212270 naïve CD4 RNA/counts sparse matrix", "required"),
    ("cd4_genes", "data/interim/genept_scpa/phase2_export/naive_cd4/naive_cd4_gene_ids.txt", "RNA/counts gene axis", "required"),
    ("cd4_cells", "data/interim/genept_scpa/phase2_export/naive_cd4/naive_cd4_cell_ids.txt", "Original source cell IDs", "required"),
    ("cd4_metadata", "data/interim/genept_scpa/phase2_export/naive_cd4/naive_cd4_metadata.csv", "Cell metadata including Hour and Cell_Type", "required"),
    ("phase4_pathway_manifest", "data/processed/genept_scpa/phase4/pathway_projection_manifest.json", "Frozen Phase 4 paired pathway universe", "required"),
    ("metabolic_pathways", "data/reference/genept_scpa/combined_metabolic_pathways.csv", "Source pathway/gene collection and provenance", "required"),
    ("genept_embeddings", "data/reference/genept_scpa/genept_ada002/GenePT_gene_embedding_ada_text.pickle", "Official GenePT ada-002 gene embeddings", "required"),
    ("ncbi_descriptions", "data/reference/genept_scpa/genept_ada002/NCBI_summary_of_genes.json", "Primary gene keys and sanitized LLM descriptions", "required"),
    ("source_cd4_rds", "data/raw/genept_scpa/GSE212270_integrated_naive_cd4.rds", "Optional source-object provenance; validated export is sufficient for Phase 7", "optional"),
    ("phase1_download_metadata", "data/interim/genept_scpa/phase1_download_metadata.json", "Optional GEO acquisition provenance", "optional"),
    ("phase1_dataset_qc", "data/interim/genept_scpa/phase1_dataset_qc.json", "Optional source-dataset validation provenance", "optional"),
]


EXPLICIT_CODE = [
    "pyproject.toml",
    "config/genept_scpa.yaml",
    "config/phase7_gpt_oss_synthetic.yaml",
    "genept_scpa_experiment_plan.md",
    "docs/genept_scpa_decision_log.md",
    "docs/phase7_gpt_oss_synthetic_protocol.md",
    "scripts/genept/build_genept_w.py",
    "scripts/phase3/run_cd4_cd8_benchmark.py",
    "scripts/phase4/run_pathway_comparison.py",
    "scripts/phase4/run_timecourse_validation.py",
    "scripts/scpa/scpa_core_adapter.R",
    "scripts/scpa/run_phase7_synthetic_masking_core.R",
    "src/gene_embedding_project/genept_scpa/config.py",
    "src/gene_embedding_project/genept_scpa/gene_mapping.py",
    "src/gene_embedding_project/genept_scpa/genept_projection.py",
    "src/gene_embedding_project/genept_scpa/io.py",
    "src/gene_embedding_project/genept_scpa/pathway_projection.py",
    "tests/test_config.py",
]


def _git_tracked(relative: str) -> bool:
    result = subprocess.run(
        ["git", "ls-files", "--error-unmatch", relative],
        cwd=PROJECT_ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    return result.returncode == 0


def _resource(
    identifier: str, relative: str, purpose: str, required_status: str, category: str
) -> dict[str, Any]:
    path = PROJECT_ROOT / relative
    if not path.is_file():
        if required_status == "required":
            raise FileNotFoundError(f"Required Phase 7 transfer resource is missing: {path}")
        raise FileNotFoundError(f"Optional resource should be filtered before inventory: {path}")
    tracked = _git_tracked(relative)
    return {
        "id": identifier,
        "category": category,
        "current_absolute_path": str(path.resolve()),
        "repository_relative_path": relative,
        "purpose": purpose,
        "tracked_by_git_at_snapshot": tracked,
        "must_transfer_separately": not tracked,
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "required_status": required_status,
    }


def build_resources() -> list[dict[str, Any]]:
    resources: list[dict[str, Any]] = []
    for identifier, relative, purpose, status in DATA_RESOURCES:
        if status == "optional" and not (PROJECT_ROOT / relative).is_file():
            continue
        resources.append(_resource(identifier, relative, purpose, status, "data"))
    code_paths = set(EXPLICIT_CODE)
    code_paths.update(
        str(path.relative_to(PROJECT_ROOT))
        for base in (PROJECT_ROOT / "scripts/phase7", PROJECT_ROOT / "src/gene_embedding_project/genept_scpa/phase7")
        for path in base.glob("*.py")
    )
    code_paths.update(
        str(path.relative_to(PROJECT_ROOT))
        for path in (PROJECT_ROOT / "tests").glob("test_phase7*")
        if path.is_file()
    )
    for relative in sorted(code_paths):
        resources.append(_resource(
            "code_" + relative.replace("/", "_").replace(".", "_"),
            relative,
            "Phase 7 implementation, protocol dependency, or regression test",
            "required",
            "code_or_protocol",
        ))
    resources.sort(key=lambda row: (row["category"], row["repository_relative_path"]))
    return resources


def markdown(manifest: dict[str, Any]) -> str:
    separate = manifest["summary"]["separately_transferable_required"]
    lines = [
        "# Phase 7 transfer manifest",
        "",
        f"Generated: `{manifest['generated_at_utc']}`",
        "",
        f"Repository root at snapshot: `{manifest['repository_root_at_snapshot']}`",
        "",
        f"Schema: `{manifest['schema_version']}`",
        "",
        "This inventory does not copy or download any resource. On the destination machine,",
        "restore each file at its exact repository-relative path and run the integrity checker.",
        "",
        "## Summary",
        "",
        f"- Resources: {manifest['summary']['resource_count']}",
        f"- Required: {manifest['summary']['required_count']}",
        f"- Required separate-transfer files: {separate['file_count']}",
        f"- Required separate-transfer bytes: {separate['total_size_bytes']}",
        "- Optional source RDS is provenance only; the checksum-validated sparse export is the Phase 7 input.",
        "- The historical export manifest contains lab absolute paths. Phase 7 uses the exact",
        "  repo-relative `source.counts_export_files` map and verifies each original manifest hash.",
        "",
        "## Resources",
        "",
        "| Required | Separate | Git | Size (bytes) | Repository-relative path | SHA256 | Purpose |",
        "|---|---:|---:|---:|---|---|---|",
    ]
    for row in manifest["resources"]:
        lines.append(
            f"| {row['required_status']} | {str(row['must_transfer_separately']).lower()} | "
            f"{str(row['tracked_by_git_at_snapshot']).lower()} | {row['size_bytes']} | "
            f"`{row['repository_relative_path']}` | `{row['sha256']}` | {row['purpose']} |"
        )
    lines.extend([
        "",
        "## Verification",
        "",
        "```bash",
        "PYTHONPATH=src python scripts/phase7/check_transfer_integrity.py \\",
        "  --manifest artifacts/phase7_transfer_manifest.json \\",
        "  --root \"$PWD\"",
        "```",
        "",
        "The checker never writes, downloads, or searches alternate paths.",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-output", type=Path, default=PROJECT_ROOT / "artifacts/phase7_transfer_manifest.json")
    parser.add_argument("--markdown-output", type=Path, default=PROJECT_ROOT / "docs/phase7_transfer_manifest.md")
    args = parser.parse_args()
    resources = build_resources()
    separate_required = [
        row for row in resources
        if row["required_status"] == "required" and row["must_transfer_separately"]
    ]
    manifest = {
        "schema_version": TRANSFER_SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "repository_root_at_snapshot": str(PROJECT_ROOT),
        "resources": resources,
        "summary": {
            "resource_count": len(resources),
            "required_count": sum(row["required_status"] == "required" for row in resources),
            "optional_count": sum(row["required_status"] == "optional" for row in resources),
            "separately_transferable_required": {
                "file_count": len(separate_required),
                "total_size_bytes": sum(row["size_bytes"] for row in separate_required),
                "paths": [row["repository_relative_path"] for row in separate_required],
            },
        },
        "files_copied": False,
        "downloads_performed": False,
    }
    write_json_atomic(manifest, args.json_output)
    args.markdown_output.parent.mkdir(parents=True, exist_ok=True)
    args.markdown_output.write_text(markdown(manifest), encoding="utf-8")
    print(
        f"PHASE7_TRANSFER_MANIFEST status=PASS resources={len(resources)} "
        f"separate_required_bytes={manifest['summary']['separately_transferable_required']['total_size_bytes']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
