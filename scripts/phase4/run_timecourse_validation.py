#!/usr/bin/env python3
"""Prepare or run Phase 4B CD4 activation validation.

The frozen primary scope is the three Naive CD4 activation comparisons.  The
previous nine-comparison runner remains available only through the explicit
``--comparison-set all_9`` option.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import time
from typing import Any, Sequence

import h5py
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from scripts.genept.build_genept_w import (  # noqa: E402
    load_and_validate_export,
    require_export,
    write_lines_atomic,
)
from scripts.phase3.run_cd4_cd8_benchmark import (  # noqa: E402
    canonical_hour,
    rows_for_ids,
)
from scripts.phase4.run_pathway_comparison import (  # noqa: E402
    read_json,
    write_csv_atomic,
    write_markdown_atomic,
)
from gene_embedding_project.genept_scpa.config import load_config  # noqa: E402
from gene_embedding_project.genept_scpa.gene_mapping import (  # noqa: E402
    load_official_genept_embeddings,
)
from gene_embedding_project.genept_scpa.genept_projection import (  # noqa: E402
    normalize_log1p_sparse,
)
from gene_embedding_project.genept_scpa.io import (  # noqa: E402
    sha256_file,
    write_json_atomic,
)
from gene_embedding_project.genept_scpa.pathway_projection import (  # noqa: E402
    average_rank_descending,
    build_paired_pathways,
    filter_eligible_pathways,
    read_wide_pathway_csv,
    significance_state,
)


GROUPS = ("cd4_0h", "cd4_12h", "cd4_24h", "cd8_0h", "cd8_12h", "cd8_24h")
PRIMARY_GROUPS = ("cd4_0h", "cd4_12h", "cd4_24h")
PRIMARY_COMPARISONS = (
    {"id": "cd4_0h_vs_12h", "comparison_class": "CD4 activation", "group_a": "cd4_0h", "group_b": "cd4_12h"},
    {"id": "cd4_12h_vs_24h", "comparison_class": "CD4 activation", "group_a": "cd4_12h", "group_b": "cd4_24h"},
    {"id": "cd4_0h_vs_24h", "comparison_class": "CD4 activation positive-control-like benchmark", "group_a": "cd4_0h", "group_b": "cd4_24h"},
)
OPTIONAL_COMPARISONS = (
    {"id": "cd8_0h_vs_12h", "comparison_class": "CD8 activation", "group_a": "cd8_0h", "group_b": "cd8_12h"},
    {"id": "cd8_12h_vs_24h", "comparison_class": "CD8 activation", "group_a": "cd8_12h", "group_b": "cd8_24h"},
    {"id": "cd8_0h_vs_24h", "comparison_class": "CD8 activation", "group_a": "cd8_0h", "group_b": "cd8_24h"},
    {"id": "cd4_vs_cd8_0h", "comparison_class": "same-time lineage", "group_a": "cd4_0h", "group_b": "cd8_0h"},
    {"id": "cd4_vs_cd8_12h", "comparison_class": "same-time lineage", "group_a": "cd4_12h", "group_b": "cd8_12h"},
    {"id": "cd4_vs_cd8_24h", "comparison_class": "same-time lineage", "group_a": "cd4_24h", "group_b": "cd8_24h"},
)
ALL_COMPARISONS = PRIMARY_COMPARISONS + OPTIONAL_COMPARISONS
# Backward-compatible public constant for callers that inspect the preserved
# nine-comparison capability. Production selection must use comparison_set().
COMPARISONS = ALL_COMPARISONS


def comparison_set(name: str) -> tuple[dict[str, str], ...]:
    if name == "cd4_activation":
        return PRIMARY_COMPARISONS
    if name == "all_9":
        return ALL_COMPARISONS
    raise ValueError(f"Unknown comparison set: {name}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--comparison-set",
        choices=("cd4_activation", "all_9"),
        default="cd4_activation",
        help="Primary default is exactly the three CD4 activation contrasts.",
    )
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--smoke-test", action="store_true")
    return parser.parse_args(argv)


def read_metadata(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not {"cell_id", "Hour"}.issubset(reader.fieldnames):
            raise ValueError(f"Metadata lacks cell_id/Hour: {path}")
        rows = list(reader)
    if len({row["cell_id"] for row in rows}) != len(rows):
        raise ValueError(f"Duplicate cell IDs in {path}")
    return rows


def metadata_hour_counts(rows: Sequence[dict[str, str]]) -> dict[str, int]:
    counts = {"0h": 0, "12h": 0, "24h": 0}
    for row in rows:
        hour = canonical_hour(row["Hour"])
        if hour in counts:
            counts[hour] += 1
    return counts


def choose_cells(
    metadata: Sequence[dict[str, str]], hour: str, *, sample_size: int, seed: int
) -> list[str]:
    candidates = sorted(
        row["cell_id"] for row in metadata if canonical_hour(row["Hour"]) == hour
    )
    if len(candidates) < sample_size:
        raise ValueError(f"{hour} has only {len(candidates)} cells; need {sample_size}")
    generator = np.random.default_rng(seed)
    indices = generator.choice(len(candidates), size=sample_size, replace=False)
    return [candidates[int(index)] for index in indices]


def create_source_audit() -> dict[str, Any]:
    interim = PROJECT_ROOT / "data/interim/genept_scpa"
    raw = PROJECT_ROOT / "data/raw/genept_scpa"
    definitions = {
        "naive_cd4": ("phase1_dataset_qc.json", "phase1_download_metadata.json", "GSE212270_integrated_naive_cd4.rds"),
        "naive_cd8": ("naive_cd8_dataset_qc.json", "naive_cd8_download_metadata.json", "GSE212270_integrated_naive_cd8.rds"),
    }
    datasets: dict[str, Any] = {}
    for dataset, (qc_name, download_name, rds_name) in definitions.items():
        qc = read_json(interim / qc_name)
        download = read_json(interim / download_name)
        export = read_json(
            interim / "phase2_export" / dataset / f"{dataset}_export_manifest.json"
        )
        rds = raw / rds_name
        archive = raw / f"{rds_name}.gz"
        hour_counts = qc["metadata"]["cells_per_timepoint"]
        required_layers = qc["expression"]["available_layers"]["RNA"]
        checks = {
            "source_rds_exists": rds.is_file(),
            "source_archive_exists": archive.is_file(),
            "rds_read_success": qc["object"]["rds_read_success"] is True,
            "all_hours_present": all(int(hour_counts[key]) >= 500 for key in ("0h", "12h", "24h")),
            "rna_counts_present": "counts" in required_layers,
            "rna_data_present": "data" in required_layers,
            "export_includes_all_cells": int(export["cells"]) == int(qc["object"]["cells"]),
            "source_object_unmodified": export["source_object_modified"] is False,
        }
        datasets[dataset] = {
            "source_rds": str(rds),
            "source_rds_size_bytes": rds.stat().st_size if rds.exists() else None,
            "source_archive": str(archive),
            "source_archive_size_bytes": archive.stat().st_size if archive.exists() else None,
            "archive_sha256": download["sha256"],
            "gzip_integrity": download["gzip_integrity"],
            "rds_read_success": qc["object"]["rds_read_success"],
            "cells": qc["object"]["cells"],
            "genes": qc["object"]["features"],
            "hour_counts": hour_counts,
            "rna_layers": required_layers,
            "export_manifest": str(interim / "phase2_export" / dataset / f"{dataset}_export_manifest.json"),
            "checks": checks,
        }
        if not all(checks.values()):
            raise RuntimeError(f"12h source audit failed for {dataset}: {checks}")
    audit = {
        "phase": "Phase 4B CD4 activation source audit",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_status": "ALREADY_PRESENT_IN_RDS",
        "download_performed": False,
        "official_source": "NCBI GEO GSE212270",
        "datasets": datasets,
        "existing_source_objects_modified": False,
        "gate": {"status": "PASS", "failed_checks": [], "warnings": []},
    }
    write_json_atomic(audit, interim / "phase4_timecourse_source_audit.json")
    return audit


def prepare_canonical_sampling(
    seed: int = 20260810,
    sample_size: int = 500,
    groups: Sequence[str] = GROUPS,
) -> dict[str, Any]:
    interim = PROJECT_ROOT / "data/interim/genept_scpa"
    sampling_dir = interim / "phase4_timecourse_sampling"
    sampling_dir.mkdir(parents=True, exist_ok=True)
    metadata = {
        dataset: read_metadata(
            interim / "phase2_export" / dataset / f"{dataset}_metadata.csv"
        )
        for dataset in ("naive_cd4", "naive_cd8")
    }
    selected: dict[str, list[str]] = {}
    phase3_dir = interim / "phase3_sampling"
    for lineage in ("cd4", "cd8"):
        dataset = f"naive_{lineage}"
        for hour in ("0h", "12h", "24h"):
            group = f"{lineage}_{hour}"
            if group not in groups:
                continue
            destination = sampling_dir / f"{group}_cells.txt"
            if hour == "0h":
                cells = (phase3_dir / f"{lineage}_0h_cells.txt").read_text(encoding="utf-8").splitlines()
            elif destination.exists():
                cells = destination.read_text(encoding="utf-8").splitlines()
            else:
                cells = choose_cells(metadata[dataset], hour, sample_size=sample_size, seed=seed)
            candidates = {
                row["cell_id"] for row in metadata[dataset]
                if canonical_hour(row["Hour"]) == hour
            }
            if len(cells) != sample_size or len(set(cells)) != sample_size or not set(cells) <= candidates:
                raise RuntimeError(f"Invalid frozen canonical cells for {group}")
            write_lines_atomic(cells, destination)
            selected[group] = cells
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "seed": seed,
        "sample_size_per_group": sample_size,
        "sampling": "numpy Generator(PCG64), without replacement, from sorted IDs",
        "phase3_zero_hour_reused": True,
        "groups": {
            group: {
                "cell_count": len(cells),
                "file": str(sampling_dir / f"{group}_cells.txt"),
                "sha256": sha256_file(sampling_dir / f"{group}_cells.txt"),
            }
            for group, cells in selected.items()
        },
    }
    write_json_atomic(manifest, sampling_dir / "sampling_manifest.json")
    return manifest


def prepare_core_inputs(
    config: dict[str, Any],
    work: Path,
    sampling: dict[str, Any],
    comparisons: Sequence[dict[str, str]],
    *,
    comparison_set_name: str,
) -> tuple[dict[str, Any], Path, Path]:
    active_groups = tuple(
        group
        for group in GROUPS
        if any(group in (item["group_a"], item["group_b"]) for item in comparisons)
    )
    print("[Phase 4B] Loading and validating reusable sparse exports...", flush=True)
    exports: dict[str, tuple[dict, Any, list[str], list[str]]] = {}
    for dataset in ("naive_cd4", "naive_cd8"):
        exports[dataset] = load_and_validate_export(require_export(dataset, force=False))
    cd4_manifest, cd4_counts, cd4_genes, cd4_ids = exports["naive_cd4"]
    cd8_manifest, cd8_counts, cd8_genes, cd8_ids = exports["naive_cd8"]

    sampling_dir = PROJECT_ROOT / "data/interim/genept_scpa/phase4_timecourse_sampling"
    selected = {
        group: (sampling_dir / f"{group}_cells.txt").read_text(encoding="utf-8").splitlines()
        for group in active_groups
    }
    print(
        f"[Phase 4B] Normalizing {len(active_groups)} active groups over their full transcriptomes...",
        flush=True,
    )
    group_log: dict[str, Any] = {}
    for lineage, counts, ids in (
        ("cd4", cd4_counts, cd4_ids), ("cd8", cd8_counts, cd8_ids)
    ):
        lineage_groups = [group for group in active_groups if group.startswith(lineage)]
        if not lineage_groups:
            continue
        union_ids = [cell for group in lineage_groups for cell in selected[group]]
        union_rows = rows_for_ids(ids, union_ids)
        normalized = normalize_log1p_sparse(counts[union_rows], normalization_target=10_000.0)
        sample_size = int(sampling["sample_size_per_group"])
        for offset, group in enumerate(lineage_groups):
            start = offset * sample_size
            group_log[group] = normalized[start:start + sample_size]

    phase4a_manifest = read_json(
        PROJECT_ROOT / "data/processed/genept_scpa/phase4/pathway_projection_manifest.json"
    )
    embedding_path = PROJECT_ROOT / "data/reference/genept_scpa/genept_ada002/GenePT_gene_embedding_ada_text.pickle"
    embeddings = load_official_genept_embeddings(embedding_path, expected_dimension=1536)
    pathways = read_wide_pathway_csv(
        PROJECT_ROOT / config["phase4"]["pathways"]["file"]
    )
    paired = build_paired_pathways(pathways, cd4_genes, cd8_genes, set(embeddings))
    eligible = filter_eligible_pathways(
        paired,
        min_genes=int(config["phase4"]["pathways"]["min_genes"]),
        max_genes=int(config["phase4"]["pathways"]["max_genes"]),
    )
    frozen_by_name = {pathway["pathway"]: pathway for pathway in phase4a_manifest["pathways"]}
    current_by_name = {pathway.definition.name: pathway for pathway in eligible}
    if set(frozen_by_name) != set(current_by_name) or len(eligible) != 123:
        raise RuntimeError("Six-group pathway universe differs from frozen Phase 4A 123 pathways")
    for name, frozen in frozen_by_name.items():
        if list(current_by_name[name].paired_genes) != frozen["paired_genes"]:
            raise RuntimeError(f"Frozen paired genes changed for {name}")

    global_genes = sorted({gene for pathway in eligible for gene in pathway.paired_genes})
    global_index = {gene: index for index, gene in enumerate(global_genes)}
    cd4_index = {gene: index for index, gene in enumerate(cd4_genes)}
    cd8_index = {gene: index for index, gene in enumerate(cd8_genes)}
    group_dense: dict[str, np.ndarray] = {}
    for group in active_groups:
        lookup = cd4_index if group.startswith("cd4") else cd8_index
        indices = np.asarray([lookup[gene] for gene in global_genes], dtype=np.int64)
        group_dense[group] = group_log[group][:, indices].toarray().astype(np.float64)
    embedding_global = np.stack([embeddings[gene] for gene in global_genes]).astype(np.float64)

    pathway_records = [
        {
            "pathway": pathway.definition.name,
            "source_database": pathway.definition.source_database,
            "n_primary_paired_genes": len(pathway.paired_genes),
            "paired_genes": list(pathway.paired_genes),
            "global_gene_indices": [global_index[gene] for gene in pathway.paired_genes],
        }
        for pathway in eligible
    ]
    manifest = {
        "phase": "Phase 4B primary Naive CD4 activation validation",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "comparison_set": comparison_set_name,
        "seed": int(sampling["seed"]),
        "sample_size_per_group": 500,
        "groups": {group: sampling["groups"][group] for group in active_groups},
        "comparisons": list(comparisons),
        "preprocessing": {
            "source": "RNA/counts",
            "normalization": "full-transcriptome cell-wise total=10000 then log1p",
            "pathway_renormalization": False,
        },
        "pathway_universe": {
            "source": "Phase 4A frozen primary paired universe",
            "paired_gene_policy": "exact Phase 4A frozen pathway/CD4/CD8/GenePT intersection",
            "input_count": 243,
            "eligible_count": 123,
            "dropped_from_phase4a": [],
            "identical_paired_genes": True,
            "same_across_comparisons": True,
            "gene_order_identical_across_comparisons": True,
        },
        "projection": {
            "vanilla": "X_P",
            "genept_primary": "non-L2 X_P @ E_P",
            "genept_sensitivity": "rowwise-L2(X_P @ E_P)",
        },
        "qval": {
            "official_scpa_1_6_2": "sqrt(-log10(Bonferroni-adjusted p))",
            "log_base": 10,
            "tie_method": "average",
            "significance_threshold": "adjusted p < 0.05",
        },
        "embedding": {"model": "text-embedding-ada-002", "dimension": 1536, "sha256": sha256_file(embedding_path)},
        "source_exports": {
            "cd4_source_object_modified": cd4_manifest["source_object_modified"],
            "cd8_source_object_modified": cd8_manifest["source_object_modified"],
        },
        "same_cells_across_branches": True,
        "global_gene_order": global_genes,
        "pathways": pathway_records,
        "scope": {"gene_contribution_run": False, "semantic_controls_run": False, "phase5_run": False},
    }
    h5_path = work / "phase4_timecourse_inputs.h5"
    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(h5_path, "w") as handle:
        handle.create_dataset("gene_names", data=np.asarray(global_genes, dtype=object), dtype=string_dtype)
        handle.create_dataset("embeddings", data=embedding_global, compression="gzip", shuffle=True)
        for group, matrix in group_dense.items():
            handle.create_dataset(f"expression/{group}", data=matrix, compression="gzip", shuffle=True)
    manifest_path = work / "timecourse_validation_manifest.json"
    write_json_atomic(manifest, manifest_path)
    print(
        f"[Phase 4B] Core input ready: {len(active_groups)} groups, "
        f"{len(comparisons)} comparisons, 123 pathways, {len(global_genes)} genes.",
        flush=True,
    )
    return manifest, h5_path, manifest_path


def run_r_core(
    h5_path: Path, manifest_path: Path, output_dir: Path, *, smoke_test: bool
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    result_path = output_dir / "timecourse_core_results.csv"
    qc_path = output_dir / "timecourse_core_qc.json"
    command = [
        "Rscript", str(PROJECT_ROOT / "scripts/scpa/run_phase4_timecourse_core.R"),
        "--input-h5", str(h5_path), "--manifest", str(manifest_path),
        "--output-csv", str(result_path), "--output-json", str(qc_path),
    ]
    if smoke_test:
        command += ["--max-pathways", "1", "--max-comparisons", "3"]
    print("[Phase 4B] Starting R core; comparison/pathway/branch progress and ETA follow...", flush=True)
    started = time.monotonic()
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
    print(f"[Phase 4B] R core completed in {(time.monotonic() - started) / 60:.1f} min.", flush=True)
    with result_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return rows, read_json(qc_path)


def parse_core_rows(rows: Sequence[dict[str, str]]) -> list[dict[str, Any]]:
    integer = {"n_primary_paired_genes"}
    numeric = {
        f"{method}_{field}"
        for method in ("vanilla", "genept", "l2")
        for field in ("raw_p", "adjusted_p", "qval", "rank")
    } | {"rank_delta"}
    parsed: list[dict[str, Any]] = []
    for source in rows:
        row: dict[str, Any] = {}
        for key, value in source.items():
            if key in integer:
                row[key] = int(value)
            elif key in numeric:
                row[key] = float(value)
            else:
                row[key] = value
        parsed.append(row)
    return parsed


def method_counts(rows: Sequence[dict[str, Any]], method: str) -> dict[str, Any]:
    raw = np.asarray([row[f"{method}_raw_p"] for row in rows], dtype=float)
    adjusted = np.asarray([row[f"{method}_adjusted_p"] for row in rows], dtype=float)
    qval = np.asarray([row[f"{method}_qval"] for row in rows], dtype=float)
    return {
        f"{method}_n_pathways": int(len(rows)),
        f"{method}_qval_positive_count": sum(row[f"{method}_qval"] > 0 for row in rows),
        f"{method}_raw_p_lt_0_05": sum(row[f"{method}_raw_p"] < 0.05 for row in rows),
        f"{method}_adj_p_lt_0_05": sum(row[f"{method}_adjusted_p"] < 0.05 for row in rows),
        f"{method}_qval_zero_count": sum(row[f"{method}_qval"] == 0 for row in rows),
        f"{method}_median_raw_p": float(np.median(raw)),
        f"{method}_median_adjusted_p": float(np.median(adjusted)),
        f"{method}_median_qval": float(np.median(qval)),
        f"{method}_qval_floor_fraction": float(np.mean(qval == 0)),
    }


def rank_correlation(values_a: Sequence[float], values_b: Sequence[float]) -> float | None:
    ranks_a = average_rank_descending(values_a)
    ranks_b = average_rank_descending(values_b)
    if np.std(ranks_a) == 0 or np.std(ranks_b) == 0:
        return None
    return float(np.corrcoef(ranks_a, ranks_b)[0, 1])


def l2_sensitivity_metrics(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    non_l2_significant = {
        row["pathway"] for row in rows if row["genept_adjusted_p"] < 0.05
    }
    l2_significant = {
        row["pathway"] for row in rows if row["l2_adjusted_p"] < 0.05
    }
    non_l2_positive = {row["pathway"] for row in rows if row["genept_qval"] > 0}
    l2_positive = {row["pathway"] for row in rows if row["l2_qval"] > 0}
    return {
        "genept_non_l2_vs_l2_spearman": rank_correlation(
            [row["genept_qval"] for row in rows],
            [row["l2_qval"] for row in rows],
        ),
        "genept_non_l2_l2_significant_overlap": len(non_l2_significant & l2_significant),
        "genept_non_l2_l2_qval_positive_overlap": len(non_l2_positive & l2_positive),
        "genept_non_l2_only_significant": len(non_l2_significant - l2_significant),
        "genept_l2_only_significant": len(l2_significant - non_l2_significant),
        "genept_non_l2_only_qval_positive": len(non_l2_positive - l2_positive),
        "genept_l2_only_qval_positive": len(l2_positive - non_l2_positive),
    }


def build_reporting(
    rows: Sequence[dict[str, Any]],
    comparisons: Sequence[dict[str, str]] = PRIMARY_COMPARISONS,
    *,
    expected_pathways: int = 123,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    by_comparison = {
        comparison["id"]: [row for row in rows if row["comparison"] == comparison["id"]]
        for comparison in comparisons
    }
    overview: list[dict[str, Any]] = []
    detection: list[dict[str, Any]] = []
    for comparison in comparisons:
        comparison_rows = by_comparison[comparison["id"]]
        if len(comparison_rows) != expected_pathways:
            raise RuntimeError(f"Incomplete result for {comparison['id']}: {len(comparison_rows)}")
        state_counts = {
            "Both significant": 0,
            "Vanilla-only significant": 0,
            "GenePT-only significant": 0,
            "Neither significant": 0,
        }
        for row in comparison_rows:
            state = significance_state(row["vanilla_adjusted_p"], row["genept_adjusted_p"])
            state_counts[state] += 1
            detection.append({
                "comparison": comparison["id"],
                "comparison_class": comparison["comparison_class"],
                "pathway": row["pathway"],
                "vanilla_adjusted_p": row["vanilla_adjusted_p"],
                "genept_adjusted_p": row["genept_adjusted_p"],
                "l2_adjusted_p": row["l2_adjusted_p"],
                "detection_state": state,
                "l2_significant": row["l2_adjusted_p"] < 0.05,
            })
        summary: dict[str, Any] = {
            "comparison": comparison["id"],
            "comparison_class": comparison["comparison_class"],
        }
        for method in ("vanilla", "genept", "l2"):
            summary.update(method_counts(comparison_rows, method))
        summary.update(l2_sensitivity_metrics(comparison_rows))
        summary.update({
            "both_significant": state_counts["Both significant"],
            "vanilla_only_significant": state_counts["Vanilla-only significant"],
            "genept_only_significant": state_counts["GenePT-only significant"],
            "neither_significant": state_counts["Neither significant"],
        })
        overview.append(summary)

    timing: list[dict[str, Any]] = []
    for lineage in ("cd4", "cd8"):
        contrast_ids = (
            f"{lineage}_0h_vs_12h", f"{lineage}_12h_vs_24h", f"{lineage}_0h_vs_24h"
        )
        if not all(comparison in by_comparison for comparison in contrast_ids):
            continue
        lookup = {
            comparison: {row["pathway"]: row for row in by_comparison[comparison]}
            for comparison in contrast_ids
        }
        for pathway in sorted(lookup[contrast_ids[0]]):
            for method in ("vanilla", "genept", "l2"):
                flags = [lookup[comparison][pathway][f"{method}_adjusted_p"] < 0.05 for comparison in contrast_ids]
                if all(flags):
                    pattern = "persistent"
                elif flags[0] and not flags[1]:
                    pattern = "early"
                elif flags[1] and not flags[0]:
                    pattern = "late"
                elif flags[2] and not flags[0] and not flags[1]:
                    pattern = "endpoint_only"
                elif not any(flags):
                    pattern = "none"
                else:
                    pattern = "mixed"
                timing.append({
                    "lineage": lineage.upper(), "pathway": pathway, "method": method,
                    "significant_0h_vs_12h": flags[0],
                    "significant_12h_vs_24h": flags[1],
                    "significant_0h_vs_24h": flags[2],
                    "descriptive_timing_pattern": pattern,
                })
    return overview, detection, timing


def ties_are_valid(rows: Sequence[dict[str, Any]]) -> bool:
    for comparison in {row["comparison"] for row in rows}:
        subset = [row for row in rows if row["comparison"] == comparison]
        for method in ("vanilla", "genept", "l2"):
            observed: dict[float, set[float]] = {}
            for row in subset:
                observed.setdefault(row[f"{method}_qval"], set()).add(row[f"{method}_rank"])
            if any(len(ranks) != 1 for ranks in observed.values()):
                return False
    return True


def load_phase4a_historical_reference() -> dict[str, Any]:
    historical_path = PROJECT_ROOT / "data/processed/genept_scpa/phase4/vanilla_vs_genept_pathway_comparison.csv"
    with historical_path.open("r", encoding="utf-8", newline="") as handle:
        source_rows = list(csv.DictReader(handle))
    rows: list[dict[str, Any]] = []
    numeric_fields = {
        f"{method}_{field}"
        for method in ("vanilla", "genept", "l2")
        for field in ("raw_p", "adjusted_p", "qval", "rank")
    }
    for source in source_rows:
        row: dict[str, Any] = {"comparison": "cd4_vs_cd8_0h"}
        row.update(source)
        for field in numeric_fields:
            row[field] = float(row[field])
        rows.append(row)
    reference: dict[str, Any] = {
        "comparison": "cd4_vs_cd8_0h",
        "comparison_class": "historical Phase 4A resting lineage",
        "path": str(historical_path),
        "sha256": sha256_file(historical_path),
        "pathway_count": len(rows),
    }
    for method in ("vanilla", "genept", "l2"):
        reference.update(method_counts(rows, method))
    return reference


def create_validation_figures(
    overview: Sequence[dict[str, Any]],
    detection: Sequence[dict[str, Any]],
    rows: Sequence[dict[str, Any]],
    output_dir: Path,
    historical_reference: dict[str, Any] | None = None,
) -> list[Path]:
    cache = Path("/tmp/genept_scpa_plot_cache")
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache / "xdg"))
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415

    output_dir.mkdir(parents=True, exist_ok=True)
    count_rows = list(overview)
    if historical_reference is not None:
        count_rows = [historical_reference, *count_rows]
    def comparison_label(value: str) -> str:
        labels = {
            "cd4_vs_cd8_0h": "CD4 vs CD8\n(resting, 0 h)",
            "cd4_0h_vs_12h": "CD4 activation\n0 h vs 12 h",
            "cd4_12h_vs_24h": "CD4 activation\n12 h vs 24 h",
            "cd4_0h_vs_24h": "CD4 activation\n0 h vs 24 h",
        }
        return labels.get(value, value.replace("_", " "))

    def pathway_label(value: str, width: int = 34) -> str:
        return "\n".join(textwrap.wrap(value.replace("_", " "), width=width))

    def save(path: Path, fig: Any) -> None:
        fig.savefig(path, dpi=180, bbox_inches="tight", pad_inches=0.15)
        plt.close(fig)
        files.append(path)

    labels = [comparison_label(row["comparison"]) for row in count_rows]
    x = np.arange(len(labels))
    width = 0.36
    files: list[Path] = []

    path = output_dir / "01_qval_zero_counts.png"
    fig, ax = plt.subplots(figsize=(11, 6), layout="constrained")
    for index, (method, label) in enumerate((("vanilla", "Vanilla"), ("genept", "GenePT-informed"))):
        ax.bar(x + (index - .5) * width, [row[f"{method}_qval_zero_count"] for row in count_rows], width, label=label)
    ax.set(xticks=x, xticklabels=labels, ylabel="Pathways with qval = 0", title="Qval-floor counts")
    ax.legend(); save(path, fig)

    path = output_dir / "02_bonferroni_significant_counts.png"
    fig, ax = plt.subplots(figsize=(11, 6), layout="constrained")
    for index, (method, label) in enumerate((("vanilla", "Vanilla"), ("genept", "GenePT-informed"))):
        ax.bar(x + (index - .5) * width, [row[f"{method}_adj_p_lt_0_05"] for row in count_rows], width, label=label)
    ax.set(xticks=x, xticklabels=labels, ylabel="Pathways with adjusted p < 0.05", title="Bonferroni-significant pathway counts")
    ax.legend(); save(path, fig)

    states = ("Neither significant", "GenePT-only significant", "Vanilla-only significant", "Both significant")
    state_code = {state: index for index, state in enumerate(states)}
    primary_labels = [row["comparison"] for row in overview]
    all_pathways = sorted({row["pathway"] for row in detection})
    detection_lookup = {
        (row["pathway"], row["comparison"]): row["detection_state"] for row in detection
    }
    pathways = [
        pathway for pathway in all_pathways
        if any(
            detection_lookup[(pathway, comparison)] in {"GenePT-only significant", "Vanilla-only significant"}
            for comparison in primary_labels
        )
    ]
    heatmap_title = "Detection states for pathways discordant in at least one comparison"
    if not pathways:
        pathways = all_pathways
        heatmap_title = "Pathway detection states"
    state_matrix = np.asarray([
        [state_code[detection_lookup[(pathway, comparison)]] for comparison in primary_labels]
        for pathway in pathways
    ])
    path = output_dir / "03_detection_state_heatmap.png"
    fig, ax = plt.subplots(figsize=(10, max(9, len(pathways) * 0.32)), layout="constrained")
    image = ax.imshow(state_matrix, aspect="auto", cmap="viridis", vmin=0, vmax=3)
    ax.set(
        yticks=np.arange(len(pathways)), yticklabels=[pathway_label(item, 42) for item in pathways],
        xticks=np.arange(len(primary_labels)), xticklabels=[comparison_label(item) for item in primary_labels],
        title=heatmap_title,
    )
    ax.tick_params(axis="y", labelsize=7)
    colorbar = fig.colorbar(image, ax=ax, ticks=range(4))
    colorbar.ax.set_yticklabels(states)
    save(path, fig)

    by_key = {(row["comparison"], row["pathway"]): row for row in rows}
    contrasts = tuple(item["id"] for item in PRIMARY_COMPARISONS)
    if all(any(row["comparison"] == contrast for row in rows) for contrast in contrasts):
        representative = sorted(
            all_pathways,
            key=lambda pathway: -max(
                max(by_key[(contrast, pathway)]["vanilla_qval"], by_key[(contrast, pathway)]["genept_qval"])
                for contrast in contrasts
            ),
        )[:8]
        path = output_dir / "04_representative_pathways.png"
        fig, axes = plt.subplots(2, 4, figsize=(18, 10), sharex=True, layout="constrained")
        for ax, pathway in zip(axes.flat, representative):
            ax.plot(range(3), [by_key[(contrast, pathway)]["vanilla_qval"] for contrast in contrasts], marker="o", label="Vanilla")
            ax.plot(range(3), [by_key[(contrast, pathway)]["genept_qval"] for contrast in contrasts], marker="o", label="GenePT-informed")
            ax.set_title(pathway_label(pathway, 28), fontsize=9)
            ax.set_xticks(range(3), ("0→12", "12→24", "0→24"))
        axes.flat[0].legend(fontsize=7)
        fig.supylabel("qval")
        fig.suptitle("Representative CD4 activation pathways (descriptive)")
        save(path, fig)

        path = output_dir / "optional_05_rank_scatter.png"
        fig, axes = plt.subplots(1, 3, figsize=(15, 5), layout="constrained")
        for ax, contrast in zip(axes, contrasts):
            subset = [row for row in rows if row["comparison"] == contrast]
            ax.scatter(
                [row["vanilla_rank"] for row in subset],
                [row["genept_rank"] for row in subset],
                s=12, alpha=0.65,
            )
            ax.set(title=comparison_label(contrast), xlabel="Vanilla average rank", ylabel="GenePT average rank")
        fig.suptitle("Tie-aware ranks (secondary/descriptive)")
        save(path, fig)
    return files


def build_qc(
    audit: dict[str, Any],
    sampling: dict[str, Any],
    manifest: dict[str, Any],
    rows: Sequence[dict[str, Any]],
    overview: Sequence[dict[str, Any]],
    core_qc: dict[str, Any],
    figures: Sequence[Path],
    historical_reference: dict[str, Any],
    comparisons: Sequence[dict[str, str]],
    *,
    comparison_set_name: str,
) -> dict[str, Any]:
    comparison_ids = {item["id"] for item in comparisons}
    observed_ids = {row["comparison"] for row in rows}
    expected_rows = len(comparisons) * 123
    active_groups = manifest["groups"]
    criteria = {
        "source_audit_pass": audit["gate"]["status"] == "PASS",
        "all_active_groups_500_cells": all(group["cell_count"] == 500 for group in active_groups.values()),
        "phase3_cd4_zero_hour_reused": sampling["phase3_zero_hour_reused"] is True,
        "requested_comparisons_complete": observed_ids == comparison_ids,
        "frozen_123_pathways_each": len(rows) == expected_rows and all(
            sum(row["comparison"] == item["id"] for row in rows) == 123
            for item in comparisons
        ),
        "paired_universe_unchanged": manifest["pathway_universe"]["identical_paired_genes"] is True,
        "same_pathways_and_genes_across_comparisons": (
            manifest["pathway_universe"]["same_across_comparisons"] is True
            and manifest["pathway_universe"]["gene_order_identical_across_comparisons"] is True
        ),
        "same_cells_across_branches": manifest["same_cells_across_branches"] is True,
        "official_qval_log10_verified": core_qc["log_base"] == 10,
        "tie_aware_average_ranking": core_qc["tie_method"] == "average" and ties_are_valid(rows),
        "official_scpa_raw_p_crosscheck": all(item["passed"] for item in core_qc["official_scpa_crosscheck"]),
        "historical_phase4a_preserved": (
            historical_reference["pathway_count"] == 123
            and sha256_file(Path(historical_reference["path"])) == historical_reference["sha256"]
        ),
        "runtime_warnings_absent": not core_qc["warnings"],
        "positive_control_evaluated": any(row["comparison"] == "cd4_0h_vs_24h" for row in rows),
        "source_objects_unmodified": manifest["source_exports"]["cd4_source_object_modified"] is False and manifest["source_exports"]["cd8_source_object_modified"] is False,
        "phase5_not_run": manifest["scope"]["phase5_run"] is False,
    }
    if comparison_set_name == "cd4_activation":
        criteria["primary_scope_exactly_three_cd4_comparisons"] = (
            tuple(item["id"] for item in comparisons)
            == tuple(item["id"] for item in PRIMARY_COMPARISONS)
            and all(not group.startswith("cd8") for group in active_groups)
        )
    failed = [name for name, passed in criteria.items() if not passed]
    status = "READY_FOR_GPT_REVIEW" if not failed else "NEEDS_REVIEW"
    overview_by_id = {row["comparison"]: row for row in overview}
    return {
        "phase": "Phase 4B primary Naive CD4 activation validation",
        "comparison_set": comparison_set_name,
        "status": status,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_audit": audit,
        "canonical_sampling": sampling,
        "pathway_universe": manifest["pathway_universe"],
        "preprocessing": manifest["preprocessing"],
        "methods": manifest["projection"],
        "qval_audit": {
            "status": "PASS_NO_FORMULA_CHANGE",
            "scpa_version": core_qc["scpa_version"],
            "multicross_version": core_qc["multicross_version"],
            "raw_p_source": core_qc["raw_p_source"],
            "multiple_testing": core_qc["multiple_testing"],
            "official_formula": core_qc["qval_formula"],
            "log_base": core_qc["log_base"],
            "historical_phase4a_sha256": historical_reference["sha256"],
            "historical_output_overwritten": False,
            "regenerated_reporting": "separate Phase 4B outputs with average tie ranks",
        },
        "tie_ranking": {"method": "average", "validated": ties_are_valid(rows)},
        "official_scpa_crosscheck": core_qc["official_scpa_crosscheck"],
        "comparisons": list(overview),
        "positive_control": overview_by_id["cd4_0h_vs_24h"],
        "resting_lineage_reference": historical_reference,
        "figures": [str(path) for path in figures],
        "scope": manifest["scope"],
        "gate": {"status": status, "criteria": criteria, "failed_checks": failed, "warnings": core_qc["warnings"]},
    }


def summary_lines(qc: dict[str, Any], overview: Sequence[dict[str, Any]]) -> list[str]:
    header = (
        "| Comparison | Vanilla q>0 | GenePT q>0 | L2 q>0 | "
        "Vanilla adj<.05 | GenePT adj<.05 | L2 adj<.05 | Both | V-only | G-only | Neither |"
    )
    divider = "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
    table = [header, divider]
    for row in overview:
        table.append(
            f"| {row['comparison']} | {row['vanilla_qval_positive_count']} | "
            f"{row['genept_qval_positive_count']} | {row['l2_qval_positive_count']} | "
            f"{row['vanilla_adj_p_lt_0_05']} | {row['genept_adj_p_lt_0_05']} | "
            f"{row['l2_adj_p_lt_0_05']} | {row['both_significant']} | "
            f"{row['vanilla_only_significant']} | {row['genept_only_significant']} | "
            f"{row['neither_significant']} |"
        )
    metrics_header = (
        "| Comparison | Branch | N | raw p<.05 | adj p<.05 | q>0 | q=0 | "
        "median raw p | median adj p | median q | q-floor fraction |"
    )
    metrics = [metrics_header, "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in overview:
        for method in ("vanilla", "genept", "l2"):
            metrics.append(
                f"| {row['comparison']} | {method} | {row[f'{method}_n_pathways']} | "
                f"{row[f'{method}_raw_p_lt_0_05']} | {row[f'{method}_adj_p_lt_0_05']} | "
                f"{row[f'{method}_qval_positive_count']} | {row[f'{method}_qval_zero_count']} | "
                f"{row[f'{method}_median_raw_p']:.6g} | {row[f'{method}_median_adjusted_p']:.6g} | "
                f"{row[f'{method}_median_qval']:.6g} | {row[f'{method}_qval_floor_fraction']:.3f} |"
            )
    l2_table = [
        "| Comparison | non-L2 vs L2 Spearman | significant overlap | q-positive overlap | non-L2-only sig | L2-only sig |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in overview:
        correlation = row["genept_non_l2_vs_l2_spearman"]
        correlation_text = "undefined (constant ranks)" if correlation is None else f"{correlation:.4f}"
        l2_table.append(
            f"| {row['comparison']} | {correlation_text} | "
            f"{row['genept_non_l2_l2_significant_overlap']} | "
            f"{row['genept_non_l2_l2_qval_positive_overlap']} | "
            f"{row['genept_non_l2_only_significant']} | {row['genept_l2_only_significant']} |"
        )
    positive = qc["positive_control"]
    resting = qc["resting_lineage_reference"]
    if (
        positive["vanilla_qval_positive_count"] > resting["vanilla_qval_positive_count"]
        and positive["genept_qval_positive_count"] > resting["genept_qval_positive_count"]
    ):
        scenario = "Scenario A-compatible: the CD4 activation positive-control-like contrast has more qval-positive pathways than historical resting CD4-vs-CD8 in both primary branches. This descriptively supports a stronger activation contrast, but does not establish cause or accuracy."
    else:
        scenario = "Scenario B or mixed pattern requires review: the CD4 activation positive control is not uniformly stronger than resting CD4-vs-CD8 across both primary branches. Review the official cross-check and preprocessing before Phase 5."
    return [
        "# Phase 4B Naive CD4 activation validation summary", "",
        f"Gate status: `{qc['gate']['status']}`", "",
        "## Source audit", "",
        f"12h source status: `{qc['source_audit']['source_status']}`. No download was performed. The primary run uses only the frozen CD4 0h/12h/24h groups (500 cells each).", "",
        "## Qval and tie audit", "",
        "Installed SCPA 1.6.2 uses `sqrt(-log10(Bonferroni-adjusted p))`; no formula correction was required. Historical Phase 4A output is read-only and preserved. Current reports use average ranks, so equal qval values receive equal ranks.", "",
        "## Official SCPA cross-check", "",
        f"Representative CD4 0h-vs-24h Vanilla pathways passed raw-p equality against `SCPA::compare_pathways()` ({len(qc['official_scpa_crosscheck'])} pathways; tolerance 1e-12).", "",
        "## Branch-level primary metrics", "", *metrics, "",
        "## Vanilla vs GenePT detection states", "", *table, "",
        "## GenePT non-L2 vs L2 sensitivity", "", *l2_table, "",
        "## Positive-control interpretation", "", scenario, "",
        "## Detection-state definition", "",
        "At Bonferroni adjusted p < 0.05, each pathway is classified as Both significant, Vanilla-only significant, GenePT-only significant or Neither significant. These categories measure detection agreement, not accuracy.", "",
        "## Timing patterns", "",
        "Early, late, persistent, endpoint-only, mixed and none labels are descriptive combinations of adjusted-p significance across 0→12, 12→24 and 0→24. Rank agreement is secondary because qval floors create large tied blocks. These are not causal biological classifications.", "",
        "## Scope limits", "",
        "CD8 activation/generalization, gene-level contribution, True/Permuted/Random control, classifier and Phase 5 were not run. This benchmark supports representation comparison only; it does not establish GenePT superiority or generalization. Phase 5 remains blocked until GPT review.",
    ]


def output_profile(comparison_set_name: str) -> dict[str, Any]:
    if comparison_set_name == "cd4_activation":
        directory = PROJECT_ROOT / "data/processed/genept_scpa/phase4_cd4_activation"
        return {
            "directory": directory,
            "all_results": "phase4_cd4_activation_all_results.csv",
            "overview": "phase4_cd4_activation_overview.csv",
            "detection": "phase4_cd4_activation_detection_states.csv",
            "timing": "phase4_cd4_activation_timing_patterns.csv",
            "manifest": "phase4_cd4_activation_manifest.json",
            "qc": "phase4_cd4_activation_qc.json",
            "summary": "phase4_cd4_activation_summary.md",
            "comparison_filenames": {
                "cd4_0h_vs_12h": "cd4_0_vs_12.csv",
                "cd4_12h_vs_24h": "cd4_12_vs_24.csv",
                "cd4_0h_vs_24h": "cd4_0_vs_24.csv",
            },
        }
    directory = PROJECT_ROOT / "data/processed/genept_scpa/phase4_validation"
    return {
        "directory": directory,
        "all_results": "timecourse_all_pathway_results.csv",
        "overview": "comparison_overview.csv",
        "detection": "detection_state_by_pathway.csv",
        "timing": "activation_timing_patterns.csv",
        "manifest": "timecourse_validation_manifest.json",
        "qc": "timecourse_validation_qc.json",
        "summary": "timecourse_validation_summary.md",
        "comparison_filenames": {
            item["id"]: f"{item['id']}.csv" for item in ALL_COMPARISONS
        },
    }


def run(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    comparisons = comparison_set(args.comparison_set)
    protocol = load_config(PROJECT_ROOT / "config/genept_scpa.yaml")
    protocol.require_phase(4)
    if protocol.active_phase != 4:
        raise RuntimeError("Phase 4 must remain active")
    config = dict(protocol.values)
    audit = create_source_audit()
    requested_groups = tuple(
        group
        for group in GROUPS
        if any(group in (item["group_a"], item["group_b"]) for item in comparisons)
    )
    sampling = prepare_canonical_sampling(
        seed=20260810, sample_size=500, groups=requested_groups
    )
    print(
        "[Phase 4B] Source audit PASS: 12h already present; canonical sampling files reused.",
        flush=True,
    )
    if args.prepare_only:
        print(
            f"PHASE4_CD4_ACTIVATION_PREP status=PASS comparison_set={args.comparison_set} "
            f"comparisons={len(comparisons)} download=false"
        )
        return 0

    with tempfile.TemporaryDirectory(prefix="genept_scpa_phase4_timecourse_") as directory:
        work = Path(directory)
        manifest, h5_path, manifest_path = prepare_core_inputs(
            config,
            work,
            sampling,
            comparisons,
            comparison_set_name=args.comparison_set,
        )
        checkpoint = (
            work
            if args.smoke_test
            else PROJECT_ROOT
            / "data/interim/genept_scpa"
            / (
                "phase4_cd4_activation_core_checkpoint"
                if args.comparison_set == "cd4_activation"
                else "phase4_timecourse_core_checkpoint"
            )
        )
        core_rows, core_qc = run_r_core(h5_path, manifest_path, checkpoint, smoke_test=args.smoke_test)
        if args.smoke_test:
            if core_qc["status"] != "PASS" or core_qc["warnings"]:
                raise RuntimeError(f"Time-course smoke failed: {core_qc}")
            if core_qc["comparison_count"] != min(3, len(comparisons)) or core_qc["pathway_count"] != 1:
                raise RuntimeError("Smoke scope was not one pathway across three comparisons")
            print(
                f"PHASE4_CD4_ACTIVATION_SMOKE status=PASS comparisons={core_qc['comparison_count']} "
                f"pathways=1 branches=3 official_crosschecks={len(core_qc['official_scpa_crosscheck'])} "
                "production_outputs_written=false"
            )
            return 0
        rows = parse_core_rows(core_rows)
        overview, detection, timing = build_reporting(rows, comparisons)
        profile = output_profile(args.comparison_set)
        processed = profile["directory"]
        processed.mkdir(parents=True, exist_ok=True)
        columns = list(rows[0])
        write_csv_atomic(rows, columns, processed / profile["all_results"])
        comparison_dir = processed / "comparisons"
        for comparison in comparisons:
            subset = [row for row in rows if row["comparison"] == comparison["id"]]
            write_csv_atomic(
                subset,
                columns,
                comparison_dir / profile["comparison_filenames"][comparison["id"]],
            )
        write_csv_atomic(overview, list(overview[0]), processed / profile["overview"])
        write_csv_atomic(detection, list(detection[0]), processed / profile["detection"])
        if timing:
            write_csv_atomic(timing, list(timing[0]), processed / profile["timing"])
        write_json_atomic(manifest, processed / profile["manifest"])
        historical = load_phase4a_historical_reference()
        figures = create_validation_figures(
            overview, detection, rows, processed / "figures", historical
        )
        qc = build_qc(
            audit,
            sampling,
            manifest,
            rows,
            overview,
            core_qc,
            figures,
            historical,
            comparisons,
            comparison_set_name=args.comparison_set,
        )
        write_json_atomic(qc, processed / profile["qc"])
        write_markdown_atomic(summary_lines(qc, overview), processed / profile["summary"])
    print(
        f"PHASE4_CD4_ACTIVATION status={qc['gate']['status']} "
        f"comparison_set={args.comparison_set} comparisons={len(comparisons)} pathways=123 "
        f"qc={processed / profile['qc']}"
    )
    return 0 if qc["gate"]["status"] == "READY_FOR_GPT_REVIEW" else 1


if __name__ == "__main__":
    raise SystemExit(run())
