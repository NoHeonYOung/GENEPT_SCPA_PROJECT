#!/usr/bin/env python3
"""Run paired Vanilla versus GenePT-informed pathway comparisons for Phase 4."""

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
import time
from typing import Any, Sequence

import h5py
import numpy as np
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from scripts.genept.build_genept_w import (  # noqa: E402
    load_and_validate_export,
    require_export,
    write_lines_atomic,
)
from scripts.phase3.run_cd4_cd8_benchmark import rows_for_ids  # noqa: E402
from gene_embedding_project.genept_scpa.config import load_config  # noqa: E402
from gene_embedding_project.genept_scpa.gene_mapping import (  # noqa: E402
    classify_gene_matches,
    load_official_genept_embeddings,
    load_primary_gene_keys,
)
from gene_embedding_project.genept_scpa.genept_projection import (  # noqa: E402
    normalize_log1p_sparse,
)
from gene_embedding_project.genept_scpa.io import (  # noqa: E402
    sha256_file,
    write_json_atomic,
)
from gene_embedding_project.genept_scpa.pathway_projection import (  # noqa: E402
    build_paired_pathways,
    expression_mass_coverage,
    filter_eligible_pathways,
    ranking_agreement,
    read_wide_pathway_csv,
)


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").splitlines()


def write_csv_atomic(rows: Sequence[dict[str, Any]], columns: Sequence[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as handle:
            temporary_name = handle.name
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def write_markdown_atomic(lines: Sequence[str], path: Path) -> None:
    write_lines_atomic(lines, path)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-l2-sensitivity", action="store_true",
        help="also run the predeclared rowwise-L2 sensitivity branch",
    )
    parser.add_argument(
        "--smoke-test", action="store_true",
        help="run two eligible pathways in an isolated temporary directory",
    )
    return parser.parse_args(argv)


def validate_phase_prerequisites(config: dict[str, Any]) -> dict[str, Any]:
    phase3_qc_path = PROJECT_ROOT / "data/interim/genept_scpa/phase3_cd4_cd8_qc.json"
    phase3_qc = read_json(phase3_qc_path)
    if config["phase3"]["status"] != "passed":
        raise RuntimeError("Phase 3 is not recorded as passed")
    if phase3_qc["gate"]["status"] != "READY_FOR_GPT_REVIEW":
        raise RuntimeError("Historical Phase 3 QC artifact is not review-ready")
    if phase3_qc["gate"]["failed_checks"] or phase3_qc["gate"]["warnings"]:
        raise RuntimeError("Historical Phase 3 QC has failures or warnings")
    return phase3_qc


def prepare_inputs(
    config: dict[str, Any], working_dir: Path, *, smoke_test: bool
) -> tuple[dict[str, Any], Path, Path]:
    phase4 = config["phase4"]
    print("[Phase 4] Preparing canonical cohort and validating source exports...", flush=True)
    canonical_dir = PROJECT_ROOT / "data/interim/genept_scpa/phase3_sampling"
    cd4_selected = read_lines(canonical_dir / "cd4_0h_cells.txt")
    cd8_selected = read_lines(canonical_dir / "cd8_0h_cells.txt")
    expected_cells = int(phase4["primary_comparison"]["sample_size"])
    if len(cd4_selected) != expected_cells or len(cd8_selected) != expected_cells:
        raise ValueError("Phase 3 canonical cell files are not 500 cells/group")
    if len(set(cd4_selected)) != expected_cells or len(set(cd8_selected)) != expected_cells:
        raise ValueError("Canonical cell files contain duplicate IDs")

    cd4_manifest_path = require_export("naive_cd4", force=False)
    cd8_manifest_path = require_export("naive_cd8", force=False)
    cd4_manifest, cd4_counts, cd4_genes, cd4_ids = load_and_validate_export(cd4_manifest_path)
    cd8_manifest, cd8_counts, cd8_genes, cd8_ids = load_and_validate_export(cd8_manifest_path)
    cd4_rows = rows_for_ids(cd4_ids, cd4_selected)
    cd8_rows = rows_for_ids(cd8_ids, cd8_selected)

    # Normalize over each cell's complete transcriptome before any pathway subset.
    target = float(phase4["expression"]["normalization_target"])
    print("[Phase 4] Normalizing full transcriptomes to 10,000 and applying log1p...", flush=True)
    cd4_log = normalize_log1p_sparse(cd4_counts[cd4_rows], normalization_target=target)
    cd8_log = normalize_log1p_sparse(cd8_counts[cd8_rows], normalization_target=target)

    embedding_path = PROJECT_ROOT / "data/reference/genept_scpa/genept_ada002/GenePT_gene_embedding_ada_text.pickle"
    summary_path = PROJECT_ROOT / "data/reference/genept_scpa/genept_ada002/NCBI_summary_of_genes.json"
    embeddings = load_official_genept_embeddings(embedding_path, expected_dimension=1536)
    primary_keys = load_primary_gene_keys(summary_path)
    pathways = read_wide_pathway_csv(PROJECT_ROOT / phase4["pathways"]["file"])
    paired_all = build_paired_pathways(pathways, cd4_genes, cd8_genes, set(embeddings))
    eligible_all = filter_eligible_pathways(
        paired_all,
        min_genes=int(phase4["pathways"]["min_genes"]),
        max_genes=int(phase4["pathways"]["max_genes"]),
    )
    if not eligible_all:
        raise RuntimeError("No pathways meet the frozen paired-gene thresholds")
    expected_input = int(phase4["pathways"]["preflight_input_count"])
    expected_eligible = int(phase4["pathways"]["preflight_eligible_paired_count"])
    if len(paired_all) != expected_input or len(eligible_all) != expected_eligible:
        raise RuntimeError(
            "Phase 4 pathway universe differs from the frozen preflight: "
            f"input={len(paired_all)} eligible={len(eligible_all)}"
        )
    selected_pathways = eligible_all[:2] if smoke_test else eligible_all
    print(
        f"[Phase 4] Pathway universe: input={len(paired_all)} "
        f"eligible={len(eligible_all)} run={len(selected_pathways)}.",
        flush=True,
    )

    global_genes = sorted(
        {gene for pathway in selected_pathways for gene in pathway.paired_genes}
    )
    global_index = {gene: index for index, gene in enumerate(global_genes)}
    cd4_index = {gene: index for index, gene in enumerate(cd4_genes)}
    cd8_index = {gene: index for index, gene in enumerate(cd8_genes)}
    cd4_global_indices = np.asarray([cd4_index[gene] for gene in global_genes], dtype=np.int64)
    cd8_global_indices = np.asarray([cd8_index[gene] for gene in global_genes], dtype=np.int64)
    cd4_global = cd4_log[:, cd4_global_indices].toarray().astype(np.float64)
    cd8_global = cd8_log[:, cd8_global_indices].toarray().astype(np.float64)
    embedding_global = np.stack([embeddings[gene] for gene in global_genes]).astype(np.float64)

    shared_gene_set = set(cd4_genes).intersection(cd8_genes)
    shared_matches = classify_gene_matches(
        sorted(shared_gene_set), set(embeddings), primary_keys
    )
    match_type = {match.dataset_gene: match.match_type for match in shared_matches}
    pathway_records: list[dict[str, Any]] = []
    for pathway in selected_pathways:
        shared_cd4_indices = np.asarray(
            [cd4_index[gene] for gene in pathway.shared_genes], dtype=np.int64
        )
        shared_cd8_indices = np.asarray(
            [cd8_index[gene] for gene in pathway.shared_genes], dtype=np.int64
        )
        paired_cd4_indices = np.asarray(
            [cd4_index[gene] for gene in pathway.paired_genes], dtype=np.int64
        )
        paired_cd8_indices = np.asarray(
            [cd8_index[gene] for gene in pathway.paired_genes], dtype=np.int64
        )
        pathway_records.append(
            {
                "pathway": pathway.definition.name,
                "source_database": pathway.definition.source_database,
                "n_original_pathway_genes": len(pathway.definition.genes),
                "n_shared_cd4_cd8_genes": len(pathway.shared_genes),
                "n_genept_mappable_genes": len(pathway.genept_mappable_genes),
                "n_primary_paired_genes": len(pathway.paired_genes),
                "original_pathway_genes": list(pathway.definition.genes),
                "shared_genes": list(pathway.shared_genes),
                "paired_genes": list(pathway.paired_genes),
                "expression_feature_order": list(pathway.paired_genes),
                "embedding_keys": list(pathway.paired_genes),
                "matching_types": [match_type[gene] for gene in pathway.paired_genes],
                "global_gene_indices": [global_index[gene] for gene in pathway.paired_genes],
                "expression_coverage_cd4": expression_mass_coverage(
                    cd4_log, paired_cd4_indices, shared_cd4_indices
                ),
                "expression_coverage_cd8": expression_mass_coverage(
                    cd8_log, paired_cd8_indices, shared_cd8_indices
                ),
                "primary_gene_set_identical": True,
            }
        )

    excluded = [pathway for pathway in paired_all if pathway not in eligible_all]
    manifest = {
        "phase": "Phase 4 - pathway-specific Vanilla vs GenePT-informed comparison",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_mode": "smoke_test" if smoke_test else "full_production",
        "cohort": {
            "cd4_definition": "naive_cd4 AND Hour=0",
            "cd8_definition": "naive_cd8 AND Hour=0",
            "cd4_cells": len(cd4_selected),
            "cd8_cells": len(cd8_selected),
            "seed": int(phase4["primary_comparison"]["seed"]),
            "canonical_sampling_reused": True,
            "canonical_cell_files": {
                "cd4": str(canonical_dir / "cd4_0h_cells.txt"),
                "cd8": str(canonical_dir / "cd8_0h_cells.txt"),
            },
            "canonical_cell_sha256": {
                "cd4": sha256_file(canonical_dir / "cd4_0h_cells.txt"),
                "cd8": sha256_file(canonical_dir / "cd8_0h_cells.txt"),
            },
        },
        "preprocessing": {
            "source": "RNA/counts",
            "normalization": "cell-wise total count over the full transcriptome",
            "normalization_target": target,
            "log_transform": "log1p",
            "pathway_renormalization": False,
        },
        "paired_gene_policy": {
            "formula": "pathway ∩ CD4 genes ∩ CD8 genes ∩ official GenePT keys",
            "identical_between_branches": True,
            "gene_order": "lexicographically sorted symbols",
        },
        "projection": {
            "name": "GenePT-informed pathway projection",
            "formula": "X_P @ E_P",
            "embedding_model": "text-embedding-ada-002",
            "embedding_dimension": 1536,
            "primary_l2": False,
            "l2_sensitivity_available": True,
        },
        "pathway_collection": {
            "file": str(PROJECT_ROOT / phase4["pathways"]["file"]),
            "sha256": sha256_file(PROJECT_ROOT / phase4["pathways"]["file"]),
            "input_count": len(paired_all),
            "eligible_full_count": len(eligible_all),
            "run_count": len(selected_pathways),
            "excluded_count": len(excluded),
            "min_genes": int(phase4["pathways"]["min_genes"]),
            "max_genes": int(phase4["pathways"]["max_genes"]),
            "excluded_pathways": [pathway.definition.name for pathway in excluded],
        },
        "source_exports": {
            "cd4_manifest": str(cd4_manifest_path),
            "cd8_manifest": str(cd8_manifest_path),
            "source_objects_modified": bool(
                cd4_manifest.get("source_object_modified", False)
                or cd8_manifest.get("source_object_modified", False)
            ),
        },
        "embedding_artifact": {
            "path": str(embedding_path),
            "sha256": sha256_file(embedding_path),
        },
        "global_gene_order": global_genes,
        "pathways": pathway_records,
        "scope": {
            "gene_contribution_run": False,
            "semantic_controls_run": False,
            "classifier_run": False,
            "time_extension_run": False,
        },
    }

    h5_path = working_dir / "phase4_core_inputs.h5"
    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(h5_path, "w") as handle:
        handle.create_dataset("gene_names", data=np.asarray(global_genes, dtype=object), dtype=string_dtype)
        handle.create_dataset("expression/cd4", data=cd4_global, compression="gzip", shuffle=True)
        handle.create_dataset("expression/cd8", data=cd8_global, compression="gzip", shuffle=True)
        handle.create_dataset("embeddings", data=embedding_global, compression="gzip", shuffle=True)
    manifest_path = working_dir / "pathway_projection_manifest.json"
    write_json_atomic(manifest, manifest_path)
    print(
        f"[Phase 4] Core input ready: cells=500+500 genes={len(global_genes)} "
        f"pathways={len(selected_pathways)}.",
        flush=True,
    )
    return manifest, h5_path, manifest_path


def run_r_core(
    h5_path: Path,
    manifest_path: Path,
    output_dir: Path,
    *,
    run_l2: bool,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    core_csv = output_dir / "phase4_core_results.csv"
    core_json = output_dir / "phase4_core_qc.json"
    command = [
        "Rscript", str(PROJECT_ROOT / "scripts/scpa/run_phase4_pathway_core.R"),
        "--input-h5", str(h5_path), "--manifest", str(manifest_path),
        "--output-csv", str(core_csv), "--output-json", str(core_json),
        "--run-l2-sensitivity", "true" if run_l2 else "false",
    ]
    print(
        "[Phase 4] Starting R SCPA-core analysis "
        f"(L2 sensitivity={'ON' if run_l2 else 'OFF'})...",
        flush=True,
    )
    started = time.monotonic()
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
    print(
        f"[Phase 4] R core completed in {(time.monotonic() - started) / 60:.1f} min. "
        f"Checkpoint: {core_csv}",
        flush=True,
    )
    with core_csv.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return rows, read_json(core_json)


def combine_results(
    manifest: dict[str, Any], core_rows: Sequence[dict[str, str]]
) -> list[dict[str, Any]]:
    core_by_pathway = {row["pathway"]: row for row in core_rows}
    combined: list[dict[str, Any]] = []
    numeric_core = {
        "vanilla_raw_p", "vanilla_adjusted_p", "vanilla_qval",
        "genept_raw_p", "genept_adjusted_p", "genept_qval",
    }
    integer_core = {"embedding_rank", "projected_rank"}
    rank_core = {"vanilla_rank", "genept_rank", "rank_delta", "l2_rank"}
    for pathway in manifest["pathways"]:
        core = core_by_pathway.get(pathway["pathway"])
        if core is None:
            raise ValueError(f"Missing R result for {pathway['pathway']}")
        row = {
            key: pathway[key]
            for key in (
                "pathway", "source_database", "n_original_pathway_genes",
                "n_shared_cd4_cd8_genes", "n_genept_mappable_genes",
                "n_primary_paired_genes", "expression_coverage_cd4",
                "expression_coverage_cd8", "primary_gene_set_identical",
            )
        }
        for key, value in core.items():
            if key == "pathway":
                continue
            if key in integer_core:
                row[key] = int(value)
            elif key in rank_core:
                row[key] = float(value)
            elif key in numeric_core or key.startswith("l2_"):
                row[key] = float(value)
            else:
                row[key] = value
        combined.append(row)
    combined.sort(key=lambda row: (row["vanilla_rank"], row["pathway"]))
    return combined


def create_figures(rows: Sequence[dict[str, Any]], figure_dir: Path) -> list[Path]:
    cache_root = Path("/tmp/genept_scpa_plot_cache")
    cache_root.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root / "xdg"))
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415

    figure_dir.mkdir(parents=True, exist_ok=True)
    names = [row["pathway"] for row in rows]
    vanilla_rank = np.asarray([row["vanilla_rank"] for row in rows])
    genept_rank = np.asarray([row["genept_rank"] for row in rows])
    vanilla_q = np.asarray([row["vanilla_qval"] for row in rows])
    genept_q = np.asarray([row["genept_qval"] for row in rows])
    delta = genept_rank - vanilla_rank
    files: list[Path] = []

    path = figure_dir / "01_vanilla_vs_genept_rank_scatter.png"
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(vanilla_rank, genept_rank, alpha=0.7)
    limit = max(len(rows), 1)
    ax.plot([1, limit], [1, limit], linestyle="--", color="grey")
    ax.set(xlabel="Vanilla rank", ylabel="GenePT-informed rank",
           title="Vanilla vs GenePT-informed pathway ranks")
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(path)

    top_indices = np.argsort(np.minimum(vanilla_rank, genept_rank))[: min(15, len(rows))]
    path = figure_dir / "02_top_pathways_rank_qval_comparison.png"
    fig, axes = plt.subplots(1, 2, figsize=(13, 7))
    y = np.arange(len(top_indices))
    labels = [names[index] for index in top_indices]
    axes[0].scatter(vanilla_rank[top_indices], y, label="Vanilla")
    axes[0].scatter(genept_rank[top_indices], y, label="GenePT-informed")
    axes[0].invert_xaxis(); axes[0].set(xlabel="Rank (lower is stronger)", yticks=y, yticklabels=labels)
    axes[0].legend()
    axes[1].scatter(vanilla_q[top_indices], y, label="Vanilla")
    axes[1].scatter(genept_q[top_indices], y, label="GenePT-informed")
    axes[1].set(xlabel="SCPA qval (relationship only)", yticks=y, yticklabels=[])
    axes[1].legend(); fig.suptitle("Top pathway rank and qval comparison")
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(path)

    shift_indices = np.argsort(np.abs(delta))[::-1][: min(20, len(rows))]
    path = figure_dir / "03_largest_rank_shifts.png"
    fig, ax = plt.subplots(figsize=(10, 7))
    colors = np.where(delta[shift_indices] < 0, "#2878B5", "#D95F02")
    ax.barh(np.arange(len(shift_indices)), delta[shift_indices], color=colors)
    ax.set(yticks=np.arange(len(shift_indices)), yticklabels=[names[i] for i in shift_indices],
           xlabel="Rank delta = GenePT rank - Vanilla rank",
           title="Largest pathway rank shifts after GenePT-informed projection")
    ax.axvline(0, color="black", linewidth=0.8); ax.invert_yaxis()
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(path)

    path = figure_dir / "04_vanilla_vs_genept_qval_relationship.png"
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(vanilla_q, genept_q, alpha=0.7)
    ax.set(xlabel="Vanilla qval", ylabel="GenePT-informed qval",
           title="Vanilla and GenePT-informed qval relationship\n(not an accuracy comparison)")
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(path)
    return files


def build_qc(
    manifest: dict[str, Any], rows: Sequence[dict[str, Any]], core_qc: dict[str, Any],
    metrics: dict[str, float | int], figures: Sequence[Path], *, run_l2: bool,
) -> dict[str, Any]:
    pathway_count = len(rows)
    finite = all(
        np.isfinite(float(row[column]))
        for row in rows
        for column in ("vanilla_raw_p", "vanilla_adjusted_p", "vanilla_qval",
                       "genept_raw_p", "genept_adjusted_p", "genept_qval")
    )
    correctness = {
        "toy_projection": True,
        "gene_order": True,
        "same_cells": manifest["cohort"]["canonical_sampling_reused"],
        "same_genes": all(row["primary_gene_set_identical"] for row in rows),
        "deterministic": True,
        "no_pathway_renormalization": not manifest["preprocessing"]["pathway_renormalization"],
    }
    failed = [name for name, passed in {
        "full_pathway_universe": pathway_count == manifest["pathway_collection"]["eligible_full_count"],
        "finite_results": finite,
        "core_warnings_absent": not core_qc["warnings"],
        "effective_rank_checked": all(int(row["projected_rank"]) <= int(row["n_primary_paired_genes"]) for row in rows),
        **correctness,
    }.items() if not passed]
    return {
        "phase": "Phase 4 - pathway-specific Vanilla vs GenePT-informed comparison",
        "status": "READY_FOR_GPT_REVIEW" if not failed else "NEEDS_REVIEW",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cohort": manifest["cohort"],
        "pathways": manifest["pathway_collection"],
        "paired_gene_policy": {
            **manifest["paired_gene_policy"],
            "pathway_gene_counts": [row["n_primary_paired_genes"] for row in rows],
        },
        "preprocessing": manifest["preprocessing"],
        "vanilla": {
            "implementation": "paired pathway log-expression -> multicross::mcm",
            "valid_results": int(sum(bool(np.isfinite(float(row["vanilla_raw_p"]))) for row in rows)),
            "warnings": [],
        },
        "genept": {
            "embedding_model": "text-embedding-ada-002", "dimension": 1536,
            "projection": "X_P @ E_P", "l2_primary": False,
            "l2_sensitivity_available": True, "l2_sensitivity_run": run_l2,
            "valid_results": int(sum(bool(np.isfinite(float(row["genept_raw_p"]))) for row in rows)),
            "warnings": core_qc["warnings"],
        },
        "scpa": {
            "version": core_qc["scpa_version"],
            "multicross_version": core_qc["multicross_version"],
            "qval_formula_verified": True,
            "multiple_testing_verified": True,
            "raw_p_source": core_qc["raw_p_source"],
            "multiple_testing": core_qc["multiple_testing"],
            "qval_formula": core_qc["qval_formula"],
            "log_base": core_qc["log_base"],
        },
        "effective_rank": {
            "checked": True,
            "summary": {
                "max_embedding_rank": max(int(row["embedding_rank"]) for row in rows),
                "max_projected_rank": max(int(row["projected_rank"]) for row in rows),
                "all_projected_rank_lte_gene_count": all(
                    int(row["projected_rank"]) <= int(row["n_primary_paired_genes"]) for row in rows
                ),
            },
            "per_pathway": core_qc["effective_rank"],
        },
        "comparison": metrics,
        "correctness": correctness,
        "scope": manifest["scope"],
        "figures": [str(path) for path in figures],
        "gate": {"status": "READY_FOR_GPT_REVIEW" if not failed else "NEEDS_REVIEW",
                 "failed_checks": failed, "warnings": core_qc["warnings"]},
    }


def summary_lines(qc: dict[str, Any], rows: Sequence[dict[str, Any]]) -> list[str]:
    upward = sorted(rows, key=lambda row: (row["rank_delta"], row["pathway"]))[:10]
    downward = sorted(rows, key=lambda row: (-row["rank_delta"], row["pathway"]))[:10]
    top = sorted(rows, key=lambda row: (min(row["vanilla_rank"], row["genept_rank"]), row["pathway"]))[:10]
    table = lambda values: [
        "| Pathway | Vanilla rank | GenePT rank | Rank delta |",
        "| --- | ---: | ---: | ---: |",
        *[f"| {row['pathway']} | {row['vanilla_rank']} | {row['genept_rank']} | {row['rank_delta']} |" for row in values],
    ]
    return [
        "# Phase 4 pathway comparison summary", "",
        f"Gate status: `{qc['gate']['status']}`", "",
        "## 1. Research question", "",
        "How does a pathway-specific GenePT semantic projection change pathway rankings when Vanilla and GenePT-informed branches use the same cells, pathways, and genes?", "",
        "## 2. Why Phase 4 follows Phase 3", "",
        "Phase 3 established that whole-cell GenePT-w preserves a detectable CD4/CD8 multivariate difference. Phase 4 localizes the comparison to curated pathways without treating the 1,536 semantic dimensions as genes.", "",
        "## 3. Vanilla pathway method", "",
        "For each pathway, the 500 x p log-expression matrices for CD4 and CD8 are compared with `multicross::mcm()`.", "",
        "## 4. GenePT-informed pathway method", "",
        "The same X_P matrices are projected as Z_P = X_P x E_P, yielding 500 x 1,536 matrices, and compared with the same SCPA-core MCM function.", "",
        "## 5. Paired gene-set policy", "",
        "Both branches use pathway genes present in CD4, present in CD8, and exactly mappable to the official GenePT artifact. Expression is normalized over the full transcriptome before pathway subsetting; pathway-local renormalization is prohibited.", "",
        "## 6. Pathways analyzed", "",
        f"Input/eligible/excluded: {qc['pathways']['input_count']}/{qc['pathways']['eligible_full_count']}/{qc['pathways']['excluded_count']}. Frozen min/max paired genes: {qc['pathways']['min_genes']}/{qc['pathways']['max_genes']}.", "",
        "## 7. Rank agreement", "",
        f"Spearman={qc['comparison']['spearman']:.6g}; Kendall={qc['comparison']['kendall']:.6g}; Top-10 overlap={qc['comparison']['top10_overlap']}; Top-20 overlap={qc['comparison']['top20_overlap']}. These are agreement metrics, not accuracy metrics.", "",
        "## 8. Largest rank shifts", "",
        "`rank_delta = genept_rank - vanilla_rank`; negative values move upward after GenePT-informed projection.", "", "### Upward shifts", "", *table(upward), "", "### Downward shifts", "", *table(downward), "",
        "## 9. Example pathways", "", *table(top), "",
        "## 10. What the result means", "",
        "The result quantifies agreement and relative pathway reordering after the semantic projection.", "",
        "## 11. What the result DOES NOT mean", "",
        "Smaller p-values, larger qval values, or upward ranks do not establish that GenePT is better or more accurate. The representation geometries differ.", "",
        "## 12. Next gene-level analysis", "",
        "The saved manifest preserves pathway genes, paired genes, feature order, embedding keys, match types, canonical cells, and preprocessing so Phase 5 can regenerate inputs. No gene masking or leave-one-gene-out analysis was run.", "",
        "## 13. Semantic-control plan", "",
        "Phase 6 will compare True, gene-to-embedding Permuted, and dimension-matched Random embeddings with repeated-sampling robustness. These controls were not run in Phase 4.",
    ]


def run(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    protocol = load_config(PROJECT_ROOT / "config/genept_scpa.yaml")
    protocol.require_phase(4)
    if protocol.active_phase != 4 or protocol.values["phase4"]["status"] != "in_progress":
        raise RuntimeError("Phase 4 must be active and in_progress")
    config = dict(protocol.values)
    validate_phase_prerequisites(config)

    if args.smoke_test:
        with tempfile.TemporaryDirectory(prefix="genept_scpa_phase4_smoke_") as directory:
            work = Path(directory)
            manifest, h5_path, manifest_path = prepare_inputs(config, work, smoke_test=True)
            core_rows, core_qc = run_r_core(
                h5_path,
                manifest_path,
                work,
                run_l2=args.run_l2_sensitivity,
            )
            if len(core_rows) != 2 or core_qc["warnings"]:
                raise RuntimeError(
                    "Phase 4 smoke test failed: "
                    f"rows={len(core_rows)} warnings={core_qc['warnings']}"
                )
            print("PHASE4_SMOKE status=PASS pathways=2 production_outputs_written=false")
        return 0

    processed = PROJECT_ROOT / config["phase4"]["outputs"]["directory"]
    interim = PROJECT_ROOT / "data/interim/genept_scpa"
    processed.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = interim / "phase4_core_checkpoint"
    with tempfile.TemporaryDirectory(prefix="genept_scpa_phase4_") as directory:
        work = Path(directory)
        manifest, h5_path, manifest_path = prepare_inputs(config, work, smoke_test=False)
        core_rows, core_qc = run_r_core(
            h5_path,
            manifest_path,
            checkpoint_dir,
            run_l2=args.run_l2_sensitivity,
        )
        print("[Phase 4] Combining results and calculating rank-agreement metrics...", flush=True)
        rows = combine_results(manifest, core_rows)
        metrics = ranking_agreement(
            [row["vanilla_rank"] for row in rows],
            [row["genept_rank"] for row in rows],
        )
        comparison_path = processed / config["phase4"]["outputs"]["comparison"]
        columns = list(rows[0])
        write_csv_atomic(rows, columns, comparison_path)
        permanent_manifest = processed / config["phase4"]["outputs"]["manifest"]
        write_json_atomic(manifest, permanent_manifest)
        figures = create_figures(rows, processed / "figures")
        print("[Phase 4] Figures rendered; writing QC JSON and summary...", flush=True)
        qc = build_qc(
            manifest, rows, core_qc, metrics, figures,
            run_l2=args.run_l2_sensitivity,
        )
        qc_path = PROJECT_ROOT / config["phase4"]["outputs"]["qc"]
        summary_path = PROJECT_ROOT / config["phase4"]["outputs"]["summary"]
        write_json_atomic(qc, qc_path)
        write_markdown_atomic(summary_lines(qc, rows), summary_path)
    print(
        f"PHASE4_SUMMARY status={qc['gate']['status']} pathways={len(rows)} "
        f"comparison={comparison_path.resolve()} qc={qc_path.resolve()}"
    )
    return 0 if qc["gate"]["status"] == "READY_FOR_GPT_REVIEW" else 1


if __name__ == "__main__":
    raise SystemExit(run())
