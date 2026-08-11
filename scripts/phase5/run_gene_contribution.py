#!/usr/bin/env python3
"""Run Phase 5 paired Vanilla/GenePT pathway gene-masking sensitivity."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any, Sequence

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from scripts.phase4.run_pathway_comparison import (  # noqa: E402
    read_json,
    write_csv_atomic,
    write_markdown_atomic,
)
from scripts.phase4.run_timecourse_validation import (  # noqa: E402
    PRIMARY_COMPARISONS,
    create_source_audit,
    prepare_canonical_sampling,
    prepare_core_inputs,
)
from gene_embedding_project.genept_scpa.config import load_config  # noqa: E402
from gene_embedding_project.genept_scpa.io import sha256_file, write_json_atomic  # noqa: E402
from gene_embedding_project.genept_scpa.pathway_projection import (  # noqa: E402
    average_rank_descending,
)


PHASE4_DIR = PROJECT_ROOT / "data/processed/genept_scpa/phase4_cd4_activation"
OUTPUT_DIR = PROJECT_ROOT / "data/processed/genept_scpa/phase5_gene_contribution"
CHECKPOINT_DIR = PROJECT_ROOT / "data/interim/genept_scpa/phase5_gene_contribution_checkpoint"
TARGET_COUNT = 30
PATHWAY_COUNT = 123
RAW_P_CLIP = 1e-300
SIGNIFICANCE_THRESHOLD = 0.05


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--smoke-test", action="store_true")
    return parser.parse_args(argv)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def phase4b_gate() -> tuple[dict[str, Any], list[dict[str, str]], dict[str, str]]:
    qc_path = PHASE4_DIR / "phase4_cd4_activation_qc.json"
    results_path = PHASE4_DIR / "phase4_cd4_activation_all_results.csv"
    manifest_path = PHASE4_DIR / "phase4_cd4_activation_manifest.json"
    qc = read_json(qc_path)
    if qc.get("gate", {}).get("status") != "READY_FOR_GPT_REVIEW":
        raise RuntimeError("Phase 4B QC is not READY_FOR_GPT_REVIEW")
    if qc["gate"].get("failed_checks") or qc["gate"].get("warnings"):
        raise RuntimeError("Phase 4B has failed checks or runtime warnings")
    if not all(qc["gate"]["criteria"].values()):
        raise RuntimeError("Phase 4B QC criteria are not all true")
    rows = read_csv(results_path)
    if len(rows) != 3 * PATHWAY_COUNT:
        raise RuntimeError(f"Phase 4B expected 369 rows, found {len(rows)}")
    hashes = {
        "qc": sha256_file(qc_path),
        "results": sha256_file(results_path),
        "manifest": sha256_file(manifest_path),
        "historical_phase4a": sha256_file(
            PROJECT_ROOT / "data/processed/genept_scpa/phase4/vanilla_vs_genept_pathway_comparison.csv"
        ),
    }
    return qc, rows, hashes


def select_targets(rows: Sequence[dict[str, str]]) -> list[dict[str, Any]]:
    comparison_lookup = {item["id"]: item for item in PRIMARY_COMPARISONS}
    targets: list[dict[str, Any]] = []
    for row in rows:
        vanilla_significant = float(row["vanilla_adjusted_p"]) < SIGNIFICANCE_THRESHOLD
        genept_significant = float(row["genept_adjusted_p"]) < SIGNIFICANCE_THRESHOLD
        if vanilla_significant == genept_significant:
            continue
        state = "Vanilla-only significant" if vanilla_significant else "GenePT-only significant"
        comparison = comparison_lookup[row["comparison"]]
        targets.append({
            "comparison": row["comparison"],
            "group_a": comparison["group_a"],
            "group_b": comparison["group_b"],
            "pathway": row["pathway"],
            "detection_state": state,
            "n_paired_genes": int(row["n_primary_paired_genes"]),
            "vanilla_raw_p": float(row["vanilla_raw_p"]),
            "vanilla_adjusted_p": float(row["vanilla_adjusted_p"]),
            "vanilla_qval": float(row["vanilla_qval"]),
            "genept_raw_p": float(row["genept_raw_p"]),
            "genept_adjusted_p": float(row["genept_adjusted_p"]),
            "genept_qval": float(row["genept_qval"]),
        })
    order = {item["id"]: index for index, item in enumerate(PRIMARY_COMPARISONS)}
    targets.sort(key=lambda row: (order[row["comparison"]], row["pathway"]))
    observed = {
        comparison: sum(target["comparison"] == comparison for target in targets)
        for comparison in order
    }
    expected = {"cd4_0h_vs_12h": 11, "cd4_12h_vs_24h": 9, "cd4_0h_vs_24h": 10}
    if len(targets) != TARGET_COUNT or observed != expected:
        raise RuntimeError(f"Frozen Phase 5 target mismatch: total={len(targets)}, by comparison={observed}")
    return targets


def target_columns() -> list[str]:
    return [
        "comparison", "pathway", "detection_state", "n_paired_genes",
        "vanilla_raw_p", "vanilla_adjusted_p", "vanilla_qval",
        "genept_raw_p", "genept_adjusted_p", "genept_qval",
    ]


def representative_targets(targets: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    for comparison in (item["id"] for item in PRIMARY_COMPARISONS):
        subset = [target for target in targets if target["comparison"] == comparison]
        for state in ("GenePT-only significant", "Vanilla-only significant"):
            candidates = [target for target in subset if target["detection_state"] == state]
            if not candidates:
                continue
            if state.startswith("GenePT"):
                difference = lambda row: -math.log10(max(row["genept_raw_p"], RAW_P_CLIP)) + math.log10(max(row["vanilla_raw_p"], RAW_P_CLIP))
            else:
                difference = lambda row: -math.log10(max(row["vanilla_raw_p"], RAW_P_CLIP)) + math.log10(max(row["genept_raw_p"], RAW_P_CLIP))
            winner = sorted(candidates, key=lambda row: (-difference(row), row["pathway"]))[0]
            selected.append({
                "comparison": comparison,
                "detection_state": state,
                "pathway": winner["pathway"],
                "baseline_directional_score_difference": difference(winner),
            })
    return selected


def build_manifest(
    base: dict[str, Any], targets: Sequence[dict[str, Any]], phase4_hashes: dict[str, str]
) -> dict[str, Any]:
    base["phase"] = "Phase 5 pathway-internal gene masking sensitivity"
    base["phase5"] = {
        "target_rule": "all Phase 4B Vanilla-only and GenePT-only pathway-comparison pairs",
        "target_count": len(targets),
        "targets": list(targets),
        "eligible_pathway_count": PATHWAY_COUNT,
        "raw_p_clip": RAW_P_CLIP,
        "primary_metric": "delta[-log10(raw p)] = full score - masked score",
        "masking": {
            "vanilla": "set the same gene expression column to zero in both populations",
            "genept_non_l2": "Z - outer(X[:,g], E[g,:]) in both populations",
            "same_gene_between_branches": True,
            "genept_l2_gene_level_run": False,
        },
        "rank": {"supporting": "descending signed delta", "absolute": "descending abs(delta)", "ties": "average"},
        "representative_selection_rule": "per comparison, maximum baseline directional score difference within each discordant state",
        "representative_targets": representative_targets(targets),
        "phase4b_artifact_sha256": phase4_hashes,
        "scope": {
            "phase6_controls_run": False,
            "cd8_generalization_run": False,
            "classifier_run": False,
            "genept_l2_gene_masking_run": False,
        },
    }
    return base


def run_r(
    h5_path: Path,
    manifest_path: Path,
    checkpoint_dir: Path,
    qc_path: Path,
    *,
    preflight_only: bool = False,
    max_targets: int = 0,
    max_genes: int = 0,
) -> dict[str, Any]:
    command = [
        "Rscript", str(PROJECT_ROOT / "scripts/scpa/run_phase5_gene_masking_core.R"),
        "--input-h5", str(h5_path), "--manifest", str(manifest_path),
        "--checkpoint-dir", str(checkpoint_dir), "--output-json", str(qc_path),
    ]
    if preflight_only:
        command.append("--preflight-only")
    if max_targets:
        command += ["--max-targets", str(max_targets)]
    if max_genes:
        command += ["--max-genes", str(max_genes)]
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
    return read_json(qc_path)


def parse_masking_rows(rows: Sequence[dict[str, str]]) -> list[dict[str, Any]]:
    integer = {"n_paired_genes", "gene_index"}
    boolean = {
        "vanilla_significant_full", "vanilla_significant_masked",
        "genept_significant_full", "genept_significant_masked",
        "vanilla_detection_flip", "genept_detection_flip",
    }
    text_fields = {"comparison", "group_a", "group_b", "pathway", "detection_state", "gene"}
    parsed: list[dict[str, Any]] = []
    for source in rows:
        row: dict[str, Any] = {}
        for key, value in source.items():
            if key in integer:
                row[key] = int(value)
            elif key in boolean:
                row[key] = value.upper() == "TRUE"
            elif key in text_fields:
                row[key] = value
            else:
                row[key] = float(value)
        parsed.append(row)
    return parsed


def collect_checkpoints(targets: Sequence[dict[str, Any]], directory: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for target in targets:
        path = directory / target["comparison"] / f"{target['pathway']}_gene_masking.csv"
        if not path.is_file():
            raise RuntimeError(f"Missing Phase 5 checkpoint: {path}")
        parsed = parse_masking_rows(read_csv(path))
        if len(parsed) != target["n_paired_genes"]:
            raise RuntimeError(f"Incomplete Phase 5 checkpoint: {path}")
        rows.extend(parsed)
    return rows


def safe_spearman(values_a: Sequence[float], values_b: Sequence[float]) -> float | None:
    rank_a = average_rank_descending(values_a)
    rank_b = average_rank_descending(values_b)
    if np.std(rank_a) == 0 or np.std(rank_b) == 0:
        return None
    return float(np.corrcoef(rank_a, rank_b)[0, 1])


def add_ranks(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    output = [dict(row) for row in rows]
    keys = sorted({(row["comparison"], row["pathway"]) for row in output})
    for key in keys:
        subset = [row for row in output if (row["comparison"], row["pathway"]) == key]
        n = len(subset)
        for method in ("vanilla", "genept"):
            signed = average_rank_descending([row[f"{method}_delta_score"] for row in subset])
            absolute = average_rank_descending([abs(row[f"{method}_delta_score"]) for row in subset])
            for row, signed_rank, absolute_rank in zip(subset, signed, absolute):
                row[f"{method}_supporting_rank"] = float(signed_rank)
                row[f"{method}_absolute_rank"] = float(absolute_rank)
                row[f"{method}_absolute_percentile"] = 1.0 if n == 1 else float((n - absolute_rank) / (n - 1))
        for row in subset:
            row["absolute_percentile_difference_genept_minus_vanilla"] = (
                row["genept_absolute_percentile"] - row["vanilla_absolute_percentile"]
            )
    return output


def top_set(rows: Sequence[dict[str, Any]], method: str, k: int) -> set[str]:
    cutoff = min(k, len(rows))
    return {row["gene"] for row in rows if row[f"{method}_absolute_rank"] <= cutoff}


def build_pathway_summary(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    summary: list[dict[str, Any]] = []
    keys = sorted({(row["comparison"], row["pathway"]) for row in rows})
    for comparison, pathway in keys:
        subset = [row for row in rows if row["comparison"] == comparison and row["pathway"] == pathway]
        vanilla_support = min(subset, key=lambda row: (row["vanilla_supporting_rank"], row["gene"]))
        genept_support = min(subset, key=lambda row: (row["genept_supporting_rank"], row["gene"]))
        vanilla_absolute = min(subset, key=lambda row: (row["vanilla_absolute_rank"], row["gene"]))
        genept_absolute = min(subset, key=lambda row: (row["genept_absolute_rank"], row["gene"]))
        genept_dominant = max(subset, key=lambda row: (row["absolute_percentile_difference_genept_minus_vanilla"], row["gene"]))
        vanilla_dominant = min(subset, key=lambda row: (row["absolute_percentile_difference_genept_minus_vanilla"], row["gene"]))
        top20 = max(1, math.ceil(0.2 * len(subset)))
        summary.append({
            "comparison": comparison,
            "pathway": pathway,
            "detection_state": subset[0]["detection_state"],
            "n_genes": len(subset),
            "vanilla_top_supporting_gene": vanilla_support["gene"],
            "genept_top_supporting_gene": genept_support["gene"],
            "vanilla_top_absolute_gene": vanilla_absolute["gene"],
            "genept_top_absolute_gene": genept_absolute["gene"],
            "signed_delta_spearman": safe_spearman(
                [row["vanilla_delta_score"] for row in subset],
                [row["genept_delta_score"] for row in subset],
            ),
            "absolute_rank_spearman": safe_spearman(
                [abs(row["vanilla_delta_score"]) for row in subset],
                [abs(row["genept_delta_score"]) for row in subset],
            ),
            "top5_overlap": len(top_set(subset, "vanilla", 5) & top_set(subset, "genept", 5)),
            "top10_overlap": len(top_set(subset, "vanilla", 10) & top_set(subset, "genept", 10)),
            "top20pct_overlap": len(top_set(subset, "vanilla", top20) & top_set(subset, "genept", top20)),
            "n_vanilla_detection_flip_genes": sum(row["vanilla_detection_flip"] for row in subset),
            "n_genept_detection_flip_genes": sum(row["genept_detection_flip"] for row in subset),
            "largest_genept_vs_vanilla_rank_difference_gene": genept_dominant["gene"],
            "largest_genept_vs_vanilla_percentile_difference": genept_dominant["absolute_percentile_difference_genept_minus_vanilla"],
            "largest_vanilla_vs_genept_rank_difference_gene": vanilla_dominant["gene"],
            "largest_vanilla_vs_genept_percentile_difference": vanilla_dominant["absolute_percentile_difference_genept_minus_vanilla"],
        })
    return summary


def rank_comparison_rows(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    fields = [
        "comparison", "pathway", "detection_state", "gene",
        "vanilla_delta_score", "genept_delta_score",
        "vanilla_supporting_rank", "genept_supporting_rank",
        "vanilla_absolute_rank", "genept_absolute_rank",
        "vanilla_absolute_percentile", "genept_absolute_percentile",
        "absolute_percentile_difference_genept_minus_vanilla",
    ]
    return [{field: row[field] for field in fields} for row in rows]


def create_figures(
    targets: Sequence[dict[str, Any]],
    rows: Sequence[dict[str, Any]],
    summary: Sequence[dict[str, Any]],
    representatives: Sequence[dict[str, Any]],
    output_dir: Path,
) -> list[Path]:
    cache = Path("/tmp/genept_scpa_plot_cache")
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache / "xdg"))
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415

    output_dir.mkdir(parents=True, exist_ok=True)
    files: list[Path] = []
    comparisons = [item["id"] for item in PRIMARY_COMPARISONS]
    labels = [item.replace("cd4_", "CD4 ").replace("h_vs_", "-").replace("h", "") for item in comparisons]

    path = output_dir / "01_target_overview.png"
    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(3)
    vanilla = [sum(t["comparison"] == c and t["detection_state"].startswith("Vanilla") for t in targets) for c in comparisons]
    genept = [sum(t["comparison"] == c and t["detection_state"].startswith("GenePT") for t in targets) for c in comparisons]
    ax.bar(x - 0.2, vanilla, 0.4, label="Vanilla-only")
    ax.bar(x + 0.2, genept, 0.4, label="GenePT-only")
    ax.set(xticks=x, xticklabels=labels, ylabel="Target pathway-comparison count", title="Phase 5 frozen discordant targets")
    ax.legend(); fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(path)

    rep_keys = [(item["comparison"], item["pathway"]) for item in representatives]
    path = output_dir / "02_gene_contribution_scatter.png"
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    for ax, key in zip(axes.flat, rep_keys):
        subset = [r for r in rows if (r["comparison"], r["pathway"]) == key]
        ax.scatter([r["vanilla_delta_score"] for r in subset], [r["genept_delta_score"] for r in subset], s=18, alpha=.65)
        divergent = sorted(subset, key=lambda r: -abs(r["absolute_percentile_difference_genept_minus_vanilla"]))[:3]
        for row in divergent:
            ax.annotate(row["gene"], (row["vanilla_delta_score"], row["genept_delta_score"]), fontsize=7)
        ax.axhline(0, color="grey", lw=.5); ax.axvline(0, color="grey", lw=.5)
        ax.set(title=f"{key[0]}\n{key[1]}", xlabel="Vanilla delta score", ylabel="GenePT delta score")
    fig.suptitle("Representative within-pathway masking sensitivity")
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(path)

    path = output_dir / "03_top_gene_contributions.png"
    fig, axes = plt.subplots(2, 3, figsize=(17, 11))
    for ax, key in zip(axes.flat, rep_keys):
        subset = [r for r in rows if (r["comparison"], r["pathway"]) == key]
        chosen = sorted(subset, key=lambda r: -max(abs(r["vanilla_delta_score"]), abs(r["genept_delta_score"])))[:8]
        y = np.arange(len(chosen))
        ax.barh(y - .18, [r["vanilla_delta_score"] for r in chosen], .36, label="Vanilla")
        ax.barh(y + .18, [r["genept_delta_score"] for r in chosen], .36, label="GenePT")
        ax.set(yticks=y, yticklabels=[r["gene"] for r in chosen], title=f"{key[0]}\n{key[1]}")
        ax.invert_yaxis()
    axes.flat[0].legend(); fig.suptitle("Top absolute gene-masking effects")
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(path)

    divergent_rows: list[dict[str, Any]] = []
    for key in rep_keys:
        subset = [r for r in rows if (r["comparison"], r["pathway"]) == key]
        divergent_rows.extend(sorted(subset, key=lambda r: -abs(r["absolute_percentile_difference_genept_minus_vanilla"]))[:4])
    matrix = np.asarray([[r["vanilla_absolute_percentile"], r["genept_absolute_percentile"]] for r in divergent_rows])
    path = output_dir / "04_rank_divergence_heatmap.png"
    fig, ax = plt.subplots(figsize=(8, max(7, len(divergent_rows) * .3)))
    image = ax.imshow(matrix, aspect="auto", cmap="coolwarm", vmin=0, vmax=1)
    ax.set(xticks=(0, 1), xticklabels=("Vanilla", "GenePT"), yticks=np.arange(len(divergent_rows)),
           yticklabels=[f"{r['comparison']} | {r['pathway']} | {r['gene']}" for r in divergent_rows],
           title="Absolute-influence percentile divergence")
    ax.tick_params(axis="y", labelsize=6); fig.colorbar(image, ax=ax)
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(path)

    correlations = [row["absolute_rank_spearman"] for row in summary if row["absolute_rank_spearman"] is not None]
    path = output_dir / "05_pathway_rank_agreement.png"
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(correlations, bins=min(12, max(3, len(correlations))), edgecolor="black")
    ax.set(xlabel="Vanilla vs GenePT absolute-rank Spearman", ylabel="Target count", title="Pathway-level gene-rank agreement")
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(path)

    path = output_dir / "06_detection_flip_counts.png"
    fig, ax = plt.subplots(figsize=(15, 6))
    x = np.arange(len(summary))
    ax.bar(x - .2, [r["n_vanilla_detection_flip_genes"] for r in summary], .4, label="Vanilla")
    ax.bar(x + .2, [r["n_genept_detection_flip_genes"] for r in summary], .4, label="GenePT")
    ax.set(xticks=x, xticklabels=[f"{r['comparison']} | {r['pathway']}" for r in summary], ylabel="Genes causing threshold flip", title="Detection-state sensitivity to single-gene masking")
    ax.tick_params(axis="x", rotation=90, labelsize=5); ax.legend()
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig); files.append(path)
    return files


def build_qc(
    phase4_qc: dict[str, Any],
    phase4_hashes: dict[str, str],
    manifest: dict[str, Any],
    targets: Sequence[dict[str, Any]],
    rows: Sequence[dict[str, Any]],
    summary: Sequence[dict[str, Any]],
    preflight: dict[str, Any],
    core: dict[str, Any],
    checkpoint_smoke: bool,
    figures: Sequence[Path],
) -> dict[str, Any]:
    current_hashes = {
        "qc": sha256_file(PHASE4_DIR / "phase4_cd4_activation_qc.json"),
        "results": sha256_file(PHASE4_DIR / "phase4_cd4_activation_all_results.csv"),
        "manifest": sha256_file(PHASE4_DIR / "phase4_cd4_activation_manifest.json"),
        "historical_phase4a": sha256_file(PROJECT_ROOT / "data/processed/genept_scpa/phase4/vanilla_vs_genept_pathway_comparison.csv"),
    }
    target_counts = {
        comparison: sum(target["comparison"] == comparison for target in targets)
        for comparison in (item["id"] for item in PRIMARY_COMPARISONS)
    }
    expected_evaluations = 2 * sum(target["n_paired_genes"] for target in targets)
    finite_fields = [key for key in rows[0] if key not in {
        "comparison", "group_a", "group_b", "pathway", "detection_state", "gene"
    }]
    all_finite = all(
        isinstance(row[key], (int, float, bool)) and math.isfinite(float(row[key]))
        for row in rows for key in finite_fields
    )
    baseline_items = list(core.get("baseline_reproduction", [])) + list(preflight.get("baseline_reproduction", []))
    baseline_differences = [
        max(float(item["vanilla_absolute_difference"]), float(item["genept_absolute_difference"]))
        for item in baseline_items
    ]
    criteria = {
        "phase4b_pass": phase4_qc["gate"]["status"] == "READY_FOR_GPT_REVIEW",
        "exact_30_targets": len(targets) == TARGET_COUNT,
        "target_counts_by_comparison": target_counts == {"cd4_0h_vs_12h": 11, "cd4_12h_vs_24h": 9, "cd4_0h_vs_24h": 10},
        "same_cells_across_branches": manifest["same_cells_across_branches"] is True,
        "frozen_pathways_genes_order": manifest["pathway_universe"]["identical_paired_genes"] is True,
        "vanilla_masking_equivalence": preflight["toy"]["vanilla_zero_vs_removal_pass"] is True,
        "genept_subtraction_equivalence": preflight["toy"]["genept_subtraction_vs_direct_pass"] is True,
        "baseline_phase4b_reproduced": bool(baseline_differences) and max(baseline_differences) <= 1e-12,
        "all_gene_mask_outputs_finite": all_finite,
        "all_targets_completed": len(summary) == TARGET_COUNT and len(rows) == sum(target["n_paired_genes"] for target in targets),
        "all_gene_mask_evaluations_completed": core["gene_mask_evaluation_count"] == expected_evaluations,
        "failed_mcm_calls_zero": core["failed_mcm_calls"] == 0,
        "runtime_warnings_absent": not core["warnings"],
        "checkpoint_resume_pass": checkpoint_smoke and core["checkpoint_resume"] is True,
        "average_tie_ranking": manifest["phase5"]["rank"]["ties"] == "average",
        "deterministic_ranking": list(rows) == add_ranks(rows),
        "raw_p_clip_frozen": core["raw_p_clip"] == RAW_P_CLIP,
        "historical_outputs_unmodified": current_hashes == phase4_hashes,
        "no_phase6_cd8_classifier_l2_gene_run": not any(manifest["phase5"]["scope"].values()),
    }
    failed = [name for name, passed in criteria.items() if not passed]
    status = "READY_FOR_GPT_REVIEW" if not failed else "NEEDS_REVIEW"
    return {
        "phase": "Phase 5 pathway-internal gene contribution sensitivity",
        "status": status,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "target_rule": manifest["phase5"]["target_rule"],
        "target_count": len(targets),
        "target_counts_by_comparison": target_counts,
        "total_genes": len(rows),
        "gene_mask_branch_evaluations": expected_evaluations,
        "masking": manifest["phase5"]["masking"],
        "primary_metric": manifest["phase5"]["primary_metric"],
        "raw_p_clip": RAW_P_CLIP,
        "preflight": preflight,
        "core": core,
        "baseline_max_absolute_difference": max(baseline_differences) if baseline_differences else None,
        "phase4_artifact_sha256_before": phase4_hashes,
        "phase4_artifact_sha256_after": current_hashes,
        "figures": [str(path) for path in figures],
        "scope": manifest["phase5"]["scope"],
        "gate": {"status": status, "criteria": criteria, "failed_checks": failed, "warnings": core["warnings"]},
    }


def summary_lines(qc: dict[str, Any], summary: Sequence[dict[str, Any]]) -> list[str]:
    correlations = [row["absolute_rank_spearman"] for row in summary if row["absolute_rank_spearman"] is not None]
    signed = [row["signed_delta_spearman"] for row in summary if row["signed_delta_spearman"] is not None]
    vanilla_flips = sum(row["n_vanilla_detection_flip_genes"] for row in summary)
    genept_flips = sum(row["n_genept_detection_flip_genes"] for row in summary)
    return [
        "# Phase 5 pathway-internal gene masking sensitivity", "",
        f"Gate status: `{qc['gate']['status']}`", "",
        "## Frozen scope", "",
        "Targets are all 30 Phase 4B pathway-comparison instances with discordant Vanilla/GenePT non-L2 Bonferroni detection states. The same Phase 4B cells, preprocessing, pathways, paired genes and gene order were reused. GenePT L2, Phase 6 controls, CD8 and classifiers were not run.", "",
        "## Method", "",
        "Vanilla masks gene g by setting `X_P[:,g]=0`. GenePT non-L2 uses `Z_P - outer(X_P[:,g], E_P[g,:])`. The primary sensitivity is `delta[-log10(raw p)] = score_full - score_masked`, with raw p clipped only for scoring at 1e-300. Positive values indicate that masking weakened the observed pathway signal; negative values indicate that masking strengthened it.", "",
        "## Technical results", "",
        f"- Targets completed: {qc['target_count']}",
        f"- Genes evaluated per branch: {qc['total_genes']}",
        f"- Total masking MCM evaluations: {qc['gene_mask_branch_evaluations']}",
        f"- Failed MCM calls: {qc['core']['failed_mcm_calls']}",
        f"- Runtime warnings: {len(qc['core']['warnings'])}",
        f"- Baseline maximum raw-p absolute difference: {qc['baseline_max_absolute_difference']:.3g}",
        f"- Median signed-delta Spearman across defined targets: {float(np.median(signed)):.3f}" if signed else "- Median signed-delta Spearman: undefined",
        f"- Median absolute-influence Spearman across defined targets: {float(np.median(correlations)):.3f}" if correlations else "- Median absolute-influence Spearman: undefined",
        f"- Vanilla threshold-flip gene instances: {vanilla_flips}",
        f"- GenePT threshold-flip gene instances: {genept_flips}", "",
        "## Interpretation limit", "",
        "These results describe representation-dependent gene-masking sensitivity. They do not identify causal genes, establish biological correctness, or show that either representation is superior. Semantic specificity and accuracy require the deferred Phase 6/7 controls.",
    ]


def copy_comparison_outputs(targets: Sequence[dict[str, Any]]) -> None:
    for target in targets:
        source = CHECKPOINT_DIR / target["comparison"] / f"{target['pathway']}_gene_masking.csv"
        destination = OUTPUT_DIR / "comparisons" / target["comparison"] / source.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def run(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    protocol = load_config(PROJECT_ROOT / "config/genept_scpa.yaml")
    protocol.require_phase(5)
    if protocol.active_phase != 5:
        raise RuntimeError("Phase 5 must be the active phase")
    phase4_qc, phase4_rows, phase4_hashes = phase4b_gate()
    targets = select_targets(phase4_rows)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    columns = target_columns()
    write_csv_atomic(
        [{column: target[column] for column in columns} for target in targets],
        columns,
        OUTPUT_DIR / "phase5_target_pathways.csv",
    )
    print("[Phase 5] Phase 4B PASS; frozen targets=30 (11/9/10).", flush=True)

    audit = create_source_audit()
    if audit["gate"]["status"] != "PASS":
        raise RuntimeError("Source audit failed")
    sampling = prepare_canonical_sampling(groups=("cd4_0h", "cd4_12h", "cd4_24h"))
    if args.prepare_only:
        print("PHASE5_PREP status=PASS targets=30 production_run=false")
        return 0

    with tempfile.TemporaryDirectory(prefix="genept_scpa_phase5_") as directory:
        work = Path(directory)
        base_manifest, h5_path, _ = prepare_core_inputs(
            dict(protocol.values), work, sampling, PRIMARY_COMPARISONS,
            comparison_set_name="cd4_activation",
        )
        manifest = build_manifest(base_manifest, targets, phase4_hashes)
        manifest_path = work / "phase5_gene_contribution_manifest.json"
        write_json_atomic(manifest, manifest_path)

        preflight_path = work / "phase5_preflight_qc.json"
        print("[Phase 5] Running masking equivalence and Phase 4B baseline preflight...", flush=True)
        preflight = run_r(h5_path, manifest_path, work / "preflight_checkpoint", preflight_path, preflight_only=True)
        if preflight["status"] != "PASS":
            raise RuntimeError("Phase 5 preflight did not pass")

        smoke_checkpoint = work / "resume_smoke"
        smoke_first_path = work / "resume_smoke_first.json"
        smoke_second_path = work / "resume_smoke_second.json"
        print("[Phase 5] Running two-gene checkpoint/resume smoke...", flush=True)
        first = run_r(h5_path, manifest_path, smoke_checkpoint, smoke_first_path, max_targets=1, max_genes=2)
        second = run_r(h5_path, manifest_path, smoke_checkpoint, smoke_second_path, max_targets=1, max_genes=2)
        checkpoint_smoke = (
            first["status"] == "PASS" and second["status"] == "PASS"
            and first["reused_checkpoint_count"] == 0
            and second["reused_checkpoint_count"] == 1
        )
        if not checkpoint_smoke:
            raise RuntimeError("Phase 5 checkpoint/resume smoke failed")
        if args.smoke_test:
            print("PHASE5_SMOKE status=PASS targets=30 preflight=PASS checkpoint_resume=PASS production_run=false")
            return 0

        core_qc_path = CHECKPOINT_DIR / "phase5_core_qc.json"
        print("[Phase 5] Pre-run gates PASS. Starting full production with checkpoint/resume...", flush=True)
        started = time.monotonic()
        core = run_r(h5_path, manifest_path, CHECKPOINT_DIR, core_qc_path)
        print(f"[Phase 5] Full R core completed in {(time.monotonic() - started) / 60:.1f} min.", flush=True)
        rows = add_ranks(collect_checkpoints(targets, CHECKPOINT_DIR))
        summary = build_pathway_summary(rows)
        rank_rows = rank_comparison_rows(rows)
        write_csv_atomic(rows, list(rows[0]), OUTPUT_DIR / "phase5_gene_masking_all_results.csv")
        write_csv_atomic(summary, list(summary[0]), OUTPUT_DIR / "phase5_pathway_summary.csv")
        write_csv_atomic(rank_rows, list(rank_rows[0]), OUTPUT_DIR / "phase5_gene_rank_comparison.csv")
        write_json_atomic(manifest, OUTPUT_DIR / "phase5_gene_contribution_manifest.json")
        copy_comparison_outputs(targets)
        figures = create_figures(targets, rows, summary, manifest["phase5"]["representative_targets"], OUTPUT_DIR / "figures")
        qc = build_qc(
            phase4_qc, phase4_hashes, manifest, targets, rows, summary,
            preflight, core, checkpoint_smoke, figures,
        )
        write_json_atomic(qc, OUTPUT_DIR / "phase5_gene_contribution_qc.json")
        write_markdown_atomic(summary_lines(qc, summary), OUTPUT_DIR / "phase5_gene_contribution_summary.md")
    print(
        f"PHASE5 status={qc['gate']['status']} targets={len(targets)} genes={len(rows)} "
        f"mask_evaluations={qc['gene_mask_branch_evaluations']}"
    )
    return 0 if qc["gate"]["status"] == "READY_FOR_GPT_REVIEW" else 1


if __name__ == "__main__":
    raise SystemExit(run())
