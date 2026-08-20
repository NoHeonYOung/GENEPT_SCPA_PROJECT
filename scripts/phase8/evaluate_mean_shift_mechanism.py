#!/usr/bin/env python3
"""Evaluate the frozen Phase 8 mean-shift mechanism decomposition."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable

import numpy as np
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from gene_embedding_project.genept_scpa.phase7.synthetic_benchmark_llmfree.metrics import average_precision  # noqa: E402
from scripts.phase7b.evaluate_null_calibration import (  # noqa: E402
    as_bool,
    descriptive,
    rank_biserial,
    read_csv,
    sha256_file,
    signed_rank,
    write_csv_atomic,
    write_json_atomic,
    write_text_atomic,
)


REPRESENTATIONS = ("vanilla", "true_genept", "permuted_genept", "random_projection")
CONTRASTS = (
    ("C_mapping", "true_genept", "permuted_genept"),
    ("C_geometry", "permuted_genept", "random_projection"),
    ("C_projection", "random_projection", "vanilla"),
    ("C_total", "true_genept", "vanilla"),
)


def _load_phase7_selected_rankings(path: Path) -> dict[tuple[str, str], dict[str, Any]]:
    groups: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in read_csv(path):
        if row["perturbation_type"] in {"null", "mean_shift"}:
            groups.setdefault((row["experiment_id"], row["method"]), []).append(row)
    output: dict[tuple[str, str], dict[str, Any]] = {}
    for key, rows in groups.items():
        ordered = sorted(rows, key=lambda row: int(row["rank"]))
        genes = [row["gene"] for row in ordered]
        if len(genes) != len(set(genes)) or [int(row["rank"]) for row in ordered] != list(range(1, len(rows) + 1)):
            raise ValueError(f"Invalid Phase 7 strict ranking: {key}")
        first = ordered[0]
        output[key] = {
            "experiment_id": first["experiment_id"],
            "draw_id": int(first["draw_id"]),
            "pathway": first["pathway"],
            "scenario": first["perturbation_type"],
            "strength": float(first["perturbation_strength"]),
            "fallback": as_bool(first["truth_fallback_used"]),
            "genes": genes,
            "truth": {row["gene"] for row in ordered if as_bool(row["is_evaluation_target"])},
            "delta": {row["gene"]: float(row["delta_score"]) for row in rows},
        }
    return output


def _phase7_metric_map(path: Path) -> dict[tuple[str, str], float]:
    return {
        (row["experiment_id"], row["method"]): float(row["average_precision"])
        for row in read_csv(path)
        if row["perturbation_type"] in {"null", "mean_shift"}
    }


def _metric_row(meta: dict[str, Any], representation: str, ranking: list[str]) -> dict[str, Any]:
    return {
        "experiment_id": meta["experiment_id"],
        "pathway": meta["pathway"],
        "draw_id": meta["draw_id"],
        "perturbation_type": meta["scenario"],
        "perturbation_strength": meta["strength"],
        "truth_fallback_used": meta["fallback"],
        "representation": representation,
        "gene_count": len(ranking),
        "truth_count": len(meta["truth"]),
        "average_precision": average_precision(ranking, meta["truth"]),
    }


def build_metrics(
    phase7_rankings: dict[tuple[str, str], dict[str, Any]],
    checkpoint_paths: list[Path], phase7_metrics: dict[tuple[str, str], float],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], float]:
    metrics: list[dict[str, Any]] = []
    control_rankings: list[dict[str, Any]] = []
    max_phase7_error = 0.0
    experiment_ids = sorted({key[0] for key in phase7_rankings})
    for experiment_id in experiment_ids:
        vanilla = phase7_rankings[(experiment_id, "vanilla_scpa")]
        true = phase7_rankings[(experiment_id, "genept_scpa")]
        if vanilla["truth"] != true["truth"] or set(vanilla["genes"]) != set(true["genes"]):
            raise ValueError(f"Phase 7 method alignment failure: {experiment_id}")
        for representation, method, meta in (
            ("vanilla", "vanilla_scpa", vanilla),
            ("true_genept", "genept_scpa", true),
        ):
            row = _metric_row(meta, representation, meta["genes"])
            max_phase7_error = max(
                max_phase7_error,
                abs(float(row["average_precision"]) - phase7_metrics[(experiment_id, method)]),
            )
            metrics.append(row)

    checkpoint_by_id = {path.name.removesuffix("_controls.csv"): path for path in checkpoint_paths}
    if set(checkpoint_by_id) != set(experiment_ids):
        missing = set(experiment_ids) - set(checkpoint_by_id)
        extra = set(checkpoint_by_id) - set(experiment_ids)
        raise RuntimeError(f"Phase 8 checkpoint mismatch: missing={len(missing)} extra={len(extra)}")
    for experiment_id in experiment_ids:
        rows = read_csv(checkpoint_by_id[experiment_id])
        meta = phase7_rankings[(experiment_id, "vanilla_scpa")]
        if {row["gene"] for row in rows} != set(meta["genes"]):
            raise ValueError(f"Control/Phase 7 gene universe mismatch: {experiment_id}")
        for short, representation in (("permuted", "permuted_genept"),
                                      ("random", "random_projection")):
            ordered = sorted(rows, key=lambda row: (-float(row[f"{short}_delta_score"]), row["gene"]))
            ranking = [row["gene"] for row in ordered]
            metrics.append(_metric_row(meta, representation, ranking))
            for rank, row in enumerate(ordered, 1):
                control_rankings.append({
                    "experiment_id": experiment_id,
                    "pathway": meta["pathway"],
                    "draw_id": meta["draw_id"],
                    "perturbation_type": meta["scenario"],
                    "perturbation_strength": meta["strength"],
                    "truth_fallback_used": meta["fallback"],
                    "representation": representation,
                    "gene": row["gene"],
                    "rank": rank,
                    "delta_score": float(row[f"{short}_delta_score"]),
                    "is_evaluation_target": row["gene"] in meta["truth"],
                    "is_ground_truth_perturbed": (
                        row["gene"] in meta["truth"] and meta["scenario"] == "mean_shift"
                    ),
                })
    return metrics, control_rankings, max_phase7_error


def build_calibrated_units(metrics: list[dict[str, Any]]) -> list[dict[str, Any]]:
    index = {
        (row["pathway"], int(row["draw_id"]), row["perturbation_type"],
         float(row["perturbation_strength"]), row["representation"]): row
        for row in metrics
    }
    pathways = sorted({row["pathway"] for row in metrics})
    draws = sorted({int(row["draw_id"]) for row in metrics})
    output: list[dict[str, Any]] = []
    for pathway in pathways:
        for draw_id in draws:
            values: dict[str, dict[str, float]] = {}
            fallback_values = set()
            for representation in REPRESENTATIONS:
                null = index[(pathway, draw_id, "null", 0.0, representation)]
                weak = index[(pathway, draw_id, "mean_shift", 0.5, representation)]
                strong = index[(pathway, draw_id, "mean_shift", 1.0, representation)]
                fallback_values.update([bool(weak["truth_fallback_used"]), bool(strong["truth_fallback_used"])])
                nonnull = float(np.mean([weak["average_precision"], strong["average_precision"]]))
                null_ap = float(null["average_precision"])
                values[representation] = {
                    "nonnull": nonnull, "null": null_ap, "calibrated": nonnull - null_ap,
                }
            if len(fallback_values) != 1:
                raise ValueError(f"Mean-shift fallback mismatch: {pathway}/draw {draw_id}")
            row: dict[str, Any] = {
                "pathway": pathway, "draw_id": draw_id,
                "inference_unit": "pathway_draw_mean_across_mean_shift_strengths",
                "truth_fallback_used": fallback_values.pop(),
            }
            for representation in REPRESENTATIONS:
                row[f"{representation}_nonnull_ap"] = values[representation]["nonnull"]
                row[f"{representation}_matched_null_ap"] = values[representation]["null"]
                row[f"{representation}_calibrated_ap"] = values[representation]["calibrated"]
            for name, left, right in CONTRASTS:
                row[name] = values[left]["calibrated"] - values[right]["calibrated"]
            output.append(row)
    if len(output) != 220:
        raise ValueError(f"Expected 220 calibrated units, got {len(output)}")
    return output


def contrast_statistics(units: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output = []
    for name, left, right in CONTRASTS:
        values = [float(row[name]) for row in units]
        summary = descriptive(values)
        statistic, p_value = signed_rank(values)
        output.append({
            "contrast": name,
            "left_representation": left,
            "right_representation": right,
            "pairing_unit": "pathway_draw_mean_across_mean_shift_strengths",
            "n": len(values),
            "mean_contrast": summary["mean"],
            "median_contrast": summary["median"],
            "contrast_sd": summary["sd"],
            "descriptive_ci95_low": summary["ci95_low"],
            "descriptive_ci95_high": summary["ci95_high"],
            "wilcoxon_statistic": statistic,
            "wilcoxon_raw_p_value": p_value,
            "wilcoxon_bonferroni_4_p_value": min(1.0, 4.0 * p_value),
            "rank_biserial_positive_is_left": rank_biserial(values),
        })
    return output


def fallback_sensitivity(units: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output = []
    for stratum in ("ALL", "NO_FALLBACK", "FALLBACK"):
        selected = units if stratum == "ALL" else [
            row for row in units
            if bool(row["truth_fallback_used"]) == (stratum == "FALLBACK")
        ]
        for name, left, right in CONTRASTS:
            summary = descriptive(float(row[name]) for row in selected)
            output.append({
                "fallback_stratum": stratum,
                "contrast": name,
                "left_representation": left,
                "right_representation": right,
                "n": len(selected),
                "mean_contrast": summary["mean"],
                "median_contrast": summary["median"],
                "descriptive_ci95_low": summary["ci95_low"],
                "descriptive_ci95_high": summary["ci95_high"],
            })
    return output


def pathway_summary(units: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output = []
    for pathway in sorted({row["pathway"] for row in units}):
        selected = [row for row in units if row["pathway"] == pathway]
        row: dict[str, Any] = {"pathway": pathway, "n_draws": len(selected)}
        for representation in REPRESENTATIONS:
            summary = descriptive(float(item[f"{representation}_calibrated_ap"]) for item in selected)
            row[f"{representation}_calibrated_ap_mean"] = summary["mean"]
            row[f"{representation}_calibrated_ap_median"] = summary["median"]
        for name, _, _ in CONTRASTS[:3]:
            summary = descriptive(float(item[name]) for item in selected)
            row[f"{name}_mean"] = summary["mean"]
            row[f"{name}_median"] = summary["median"]
        output.append(row)
    return output


def _means(units: list[dict[str, Any]]) -> dict[str, float]:
    return {
        representation: float(np.mean([
            row[f"{representation}_calibrated_ap"] for row in units
        ]))
        for representation in REPRESENTATIONS
    }


def _scientific_pattern(statistics: list[dict[str, Any]], alpha: float) -> str:
    supported = {
        row["contrast"]: row["mean_contrast"] > 0
        and row["wilcoxon_bonferroni_4_p_value"] < alpha
        for row in statistics
    }
    if supported["C_mapping"] and not supported["C_geometry"]:
        return "CASE_A_COMPATIBLE_MAPPING_CONTRIBUTION"
    if not supported["C_mapping"] and supported["C_geometry"]:
        return "CASE_B_COMPATIBLE_GEOMETRY_EFFECT"
    if (not supported["C_mapping"] and not supported["C_geometry"]
            and supported["C_projection"]):
        return "CASE_C_COMPATIBLE_GENERIC_PROJECTION_EFFECT"
    if not any(supported.values()):
        return "CASE_D_NO_PRIMARY_POSITIVE_CONTRAST_SUPPORTED"
    if supported["C_mapping"] and supported["C_geometry"] and supported["C_projection"]:
        return "CASE_E_COMPATIBLE_MULTI_COMPONENT_PATTERN"
    return "MIXED_OR_UNEXPECTED_ORDERING"


def evaluate(config_path: Path) -> dict[str, Any]:
    config_path = config_path.resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    artifacts = config["artifacts"]
    interim = PROJECT_ROOT / artifacts["interim_directory"]
    processed = PROJECT_ROOT / artifacts["processed_directory"]
    source = {name: PROJECT_ROOT / path for name, path in config["source"].items()}
    manifest_path = processed / artifacts["manifest"]
    qc_path = processed / artifacts["control_qc"]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    if not (qc.get("status") == "PASS" and qc.get("partial_run") is False
            and int(qc.get("mcm_count_permuted", -1)) == int(config["workload"]["expected_permuted_mcm"])
            and int(qc.get("mcm_count_random", -1)) == int(config["workload"]["expected_random_mcm"])
            and int(qc.get("failed_mcm_calls", -1)) == 0):
        raise RuntimeError("Phase 8 full control MCM QC is not PASS")

    snapshot = json.loads((processed / artifacts["input_snapshot"]).read_text(encoding="utf-8"))
    for name, record in snapshot["input_sha256"].items():
        actual = sha256_file(PROJECT_ROOT / record["path"])
        if actual != record["sha256"]:
            raise RuntimeError(f"Frozen input hash changed: {name}")
    if sha256_file(manifest_path) != snapshot["generated"]["manifest"]["sha256"]:
        raise RuntimeError("Phase 8 manifest changed after preparation")
    controls_h5 = PROJECT_ROOT / manifest["controls_h5"]
    if sha256_file(controls_h5) != manifest["controls_h5_sha256"]:
        raise RuntimeError("Phase 8 controls HDF5 changed after preparation")

    phase7_rankings = _load_phase7_selected_rankings(source["phase7_rankings"])
    phase7_metric_values = _phase7_metric_map(source["phase7_metrics"])
    checkpoint_paths = sorted((interim / artifacts["checkpoints"]).glob("*_controls.csv"))
    metrics, control_rankings, max_phase7_error = build_metrics(
        phase7_rankings, checkpoint_paths, phase7_metric_values
    )
    if max_phase7_error > 1e-12:
        raise RuntimeError(f"Reused Phase 7 ranking AP mismatch: {max_phase7_error}")
    units = build_calibrated_units(metrics)
    statistics = contrast_statistics(units)
    fallback = fallback_sensitivity(units)
    pathways = pathway_summary(units)

    phase7b_units = {
        (row["pathway"], int(row["draw_id"])): float(row["calibrated_did"])
        for row in read_csv(source["phase7b_calibrated_units"])
        if row["scenario"] == "mean_shift"
    }
    replication_errors = [
        abs(float(row["C_total"]) - phase7b_units[(row["pathway"], int(row["draw_id"]))])
        for row in units
    ]
    max_replication_error = max(replication_errors)
    if max_replication_error > float(config["evaluation"]["phase7b_total_replication_tolerance"]):
        raise RuntimeError(f"Phase 7B C_total unit replication failed: {max_replication_error}")

    alpha = float(config["evaluation"]["family_alpha"])
    pattern = _scientific_pattern(statistics, alpha)
    means = _means(units)
    stats_by_name = {row["contrast"]: row for row in statistics}
    technical_qc = {
        "status": "PASS",
        "existing_phase_input_hashes_unchanged": True,
        "vanilla_ranking_reused": True,
        "true_genept_ranking_reused": True,
        "new_vanilla_mcm": 0,
        "new_true_genept_mcm": 0,
        "permuted_completed_mcm": qc["mcm_count_permuted"],
        "random_completed_mcm": qc["mcm_count_random"],
        "failed_mcm_calls": qc["failed_mcm_calls"],
        "partial_run": qc["partial_run"],
        "control_unit_count": manifest["control_unit_count"],
        "same_control_matrix_three_states": qc["control_qc"]["same_control_paths_for_three_states"],
        "permuted_fixed_point_max": qc["control_qc"]["permuted_fixed_point_max"],
        "permuted_vector_multiset_preserved": qc["control_qc"]["permuted_row_multiset_preserved_all"],
        "permuted_sorted_row_norm_max_abs_difference": qc["control_qc"]["permuted_sorted_row_norm_max_abs_difference"],
        "random_dimension": qc["control_qc"]["random_dimension"],
        "random_corresponding_row_norm_max_abs_difference": qc["control_qc"]["random_corresponding_row_norm_max_abs_difference"],
        "nonfinite_count": qc["control_qc"]["nonfinite_count"],
        "ground_truth_read_during_control_generation": qc["ground_truth_read_during_control_generation"],
        "ground_truth_read_during_masking": qc["ground_truth_read_during_masking"],
        "reused_phase7_ap_max_error": max_phase7_error,
        "phase7b_C_total_unit_max_error": max_replication_error,
    }
    if not all([
        technical_qc["permuted_fixed_point_max"] == 0,
        technical_qc["permuted_vector_multiset_preserved"] is True,
        technical_qc["random_dimension"] == 1536,
        technical_qc["nonfinite_count"] == 0,
        technical_qc["ground_truth_read_during_control_generation"] is False,
        technical_qc["ground_truth_read_during_masking"] is False,
    ]):
        raise RuntimeError("One or more Phase 8 final technical gates failed")

    summary = {
        "status": "COMPLETED_WITH_SINGLE_REALIZATION_WARNING",
        "technical_qc": technical_qc,
        "calibrated_ap_means": means,
        "primary_contrasts": stats_by_name,
        "scientific_pattern": pattern,
        "fallback_sensitivity": fallback,
        "pathway_summary": pathways,
        "interpretation_limits": config["interpretation_limits"],
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }

    write_csv_atomic(control_rankings, processed / artifacts["control_rankings"])
    write_csv_atomic(metrics, processed / artifacts["metrics"])
    write_csv_atomic(units, processed / artifacts["calibrated_units"])
    write_csv_atomic(statistics, processed / artifacts["statistics"])
    write_csv_atomic(fallback, processed / artifacts["fallback_sensitivity"])
    write_csv_atomic(pathways, processed / artifacts["pathway_summary"])
    write_json_atomic(summary, processed / artifacts["summary"])

    table = "\n".join(
        f"| {row['contrast']} | {row['mean_contrast']:.6f} | {row['median_contrast']:.6f} "
        f"| {row['wilcoxon_bonferroni_4_p_value']:.4g} "
        f"| {row['rank_biserial_positive_is_left']:.3f} |"
        for row in statistics
    )
    no_fallback = {row["contrast"]: row for row in fallback if row["fallback_stratum"] == "NO_FALLBACK"}
    yes_fallback = {row["contrast"]: row for row in fallback if row["fallback_stratum"] == "FALLBACK"}
    results = f"""# Phase 8 mean-shift mechanism decomposition 결과

Status: **{summary['status']}**
Technical QC: **PASS** — Vanilla/TRUE 새 MCM 0, PERMUTED {qc['mcm_count_permuted']:,}, RANDOM {qc['mcm_count_random']:,}, failed 0

## Mean calibrated AP

- Vanilla: {means['vanilla']:.6f}
- TRUE GenePT: {means['true_genept']:.6f}
- PERMUTED GenePT: {means['permuted_genept']:.6f}
- RANDOM projection: {means['random_projection']:.6f}

## Primary contrasts

| Contrast | Mean | Median | Bonferroni p | Rank-biserial |
|---|---:|---:|---:|---:|
{table}

Phase 7B mean-shift `TRUE−VANILLA` unit replication maximum error는 {max_replication_error:.3g}다.

## Fallback sensitivity

- NO_FALLBACK n={no_fallback['C_total']['n']}: mapping {no_fallback['C_mapping']['mean_contrast']:.6f}, geometry {no_fallback['C_geometry']['mean_contrast']:.6f}, projection {no_fallback['C_projection']['mean_contrast']:.6f}, total {no_fallback['C_total']['mean_contrast']:.6f}
- FALLBACK n={yes_fallback['C_total']['n']}: mapping {yes_fallback['C_mapping']['mean_contrast']:.6f}, geometry {yes_fallback['C_geometry']['mean_contrast']:.6f}, projection {yes_fallback['C_projection']['mean_contrast']:.6f}, total {yes_fallback['C_total']['mean_contrast']:.6f}

## Interpretation

사전 정의 statistical-support rule에 따른 pattern은 **{pattern}**이다. 이는 frozen mean-shift
synthetic benchmark의 한 control realization/pathway×draw 결과다. Correct correspondence,
GenePT geometry 또는 generic projection의 인과적 기여를 입증하지 않으며 biological
superiority, causality, general validity 또는 semantics proven을 주장하지 않는다.
"""
    write_text_atomic(results, processed / artifacts["results"])
    print(
        f"PHASE8_EVAL status={summary['status']} qc=PASS pattern={pattern} "
        f"C_mapping={stats_by_name['C_mapping']['mean_contrast']:.6f} "
        f"C_total={stats_by_name['C_total']['mean_contrast']:.6f}"
    )
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path,
                        default=PROJECT_ROOT / "config/phase8_mean_shift_mechanism.yaml")
    args = parser.parse_args()
    evaluate(args.config)
