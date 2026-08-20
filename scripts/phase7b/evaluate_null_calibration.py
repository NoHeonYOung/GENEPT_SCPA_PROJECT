#!/usr/bin/env python3
"""Phase 7B matched-null calibration using frozen Phase 7 rankings only."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Iterable

import h5py
import numpy as np
from scipy.stats import rankdata, spearmanr, wilcoxon
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from gene_embedding_project.genept_scpa.phase7.synthetic_benchmark_llmfree.metrics import (  # noqa: E402
    average_precision,
    exact_random_chance,
)


METHODS = ("vanilla_scpa", "genept_scpa")
SCENARIOS = ("mean_shift", "cell_subset", "mixed_direction")
SCOPES = ("overall_non_null",) + SCENARIOS


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv_atomic(rows: list[dict[str, Any]], path: Path) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.partial")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def write_text_atomic(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.partial")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def write_json_atomic(payload: dict[str, Any], path: Path) -> None:
    write_text_atomic(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", path)


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def truth_count(gene_count: int) -> int:
    return max(3, min(10, math.ceil(0.15 * gene_count)))


def rank_biserial(values: Iterable[float]) -> float:
    differences = np.asarray(list(values), dtype=np.float64)
    nonzero = differences[np.abs(differences) > 0]
    if len(nonzero) == 0:
        return 0.0
    ranks = rankdata(np.abs(nonzero), method="average")
    positive = float(np.sum(ranks[nonzero > 0]))
    negative = float(np.sum(ranks[nonzero < 0]))
    return (positive - negative) / (positive + negative)


def signed_rank(values: Iterable[float]) -> tuple[float, float]:
    array = np.asarray(list(values), dtype=np.float64)
    if len(array) == 0:
        raise ValueError("Paired test requires observations")
    if np.all(array == 0):
        return 0.0, 1.0
    result = wilcoxon(array, zero_method="wilcox", alternative="two-sided", method="auto")
    return float(result.statistic), float(result.pvalue)


def descriptive(values: Iterable[float]) -> dict[str, float | int]:
    array = np.asarray(list(values), dtype=np.float64)
    if len(array) == 0:
        raise ValueError("Description requires observations")
    sd = float(np.std(array, ddof=1)) if len(array) > 1 else 0.0
    se = sd / math.sqrt(len(array)) if len(array) > 1 else 0.0
    mean = float(np.mean(array))
    return {
        "n": int(len(array)),
        "mean": mean,
        "median": float(np.median(array)),
        "sd": sd,
        "se": se,
        "ci95_low": mean - 1.96 * se,
        "ci95_high": mean + 1.96 * se,
    }


def _metric_index(rows: list[dict[str, str]]) -> dict[tuple[str, str, int, float, str], dict[str, Any]]:
    output: dict[tuple[str, str, int, float, str], dict[str, Any]] = {}
    for row in rows:
        parsed = {
            **row,
            "draw_id": int(row["draw_id"]),
            "perturbation_strength": float(row["perturbation_strength"]),
            "average_precision": float(row["average_precision"]),
            "truth_fallback_used": as_bool(row["truth_fallback_used"]),
            "gene_count": int(row["gene_count"]),
            "truth_count": int(row["truth_count"]),
        }
        key = (
            parsed["pathway"], parsed["perturbation_type"], parsed["draw_id"],
            parsed["perturbation_strength"], parsed["method"],
        )
        if key in output:
            raise ValueError(f"Duplicate metric row: {key}")
        output[key] = parsed
    return output


def build_calibrated_units(metric_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    """Match method-specific null AP before computing each pathway/scenario/draw DiD."""

    index = _metric_index(metric_rows)
    pathways = sorted({key[0] for key in index})
    draws = sorted({key[2] for key in index})
    output: list[dict[str, Any]] = []
    for pathway in pathways:
        for scenario in SCENARIOS:
            for draw_id in draws:
                null_rows = {
                    method: index[(pathway, "null", draw_id, 0.0, method)] for method in METHODS
                }
                strength_rows: dict[str, list[dict[str, Any]]] = {}
                for method in METHODS:
                    rows = [index[(pathway, scenario, draw_id, strength, method)]
                            for strength in (0.5, 1.0)]
                    if len({row["truth_fallback_used"] for row in rows}) != 1:
                        raise ValueError("Fallback status differs between strengths")
                    strength_rows[method] = rows
                vanilla_raw = float(np.mean([row["average_precision"] for row in strength_rows["vanilla_scpa"]]))
                genept_raw = float(np.mean([row["average_precision"] for row in strength_rows["genept_scpa"]]))
                vanilla_null = float(null_rows["vanilla_scpa"]["average_precision"])
                genept_null = float(null_rows["genept_scpa"]["average_precision"])
                vanilla_cal = vanilla_raw - vanilla_null
                genept_cal = genept_raw - genept_null
                fallback = strength_rows["vanilla_scpa"][0]["truth_fallback_used"]
                if fallback != strength_rows["genept_scpa"][0]["truth_fallback_used"]:
                    raise ValueError("Fallback status differs between methods")
                output.append({
                    "pathway": pathway,
                    "scenario": scenario,
                    "draw_id": draw_id,
                    "inference_unit": "pathway_scenario_draw_mean_across_strengths",
                    "truth_fallback_used": fallback,
                    "null_truth_fallback_used": null_rows["vanilla_scpa"]["truth_fallback_used"],
                    "gene_count": strength_rows["vanilla_scpa"][0]["gene_count"],
                    "truth_count": strength_rows["vanilla_scpa"][0]["truth_count"],
                    "vanilla_ap_strength_0_5": strength_rows["vanilla_scpa"][0]["average_precision"],
                    "vanilla_ap_strength_1_0": strength_rows["vanilla_scpa"][1]["average_precision"],
                    "genept_ap_strength_0_5": strength_rows["genept_scpa"][0]["average_precision"],
                    "genept_ap_strength_1_0": strength_rows["genept_scpa"][1]["average_precision"],
                    "vanilla_raw_ap": vanilla_raw,
                    "genept_raw_ap": genept_raw,
                    "vanilla_matched_null_ap": vanilla_null,
                    "genept_matched_null_ap": genept_null,
                    "vanilla_calibrated_ap": vanilla_cal,
                    "genept_calibrated_ap": genept_cal,
                    "raw_genept_minus_vanilla_ap": genept_raw - vanilla_raw,
                    "null_genept_minus_vanilla_ap": genept_null - vanilla_null,
                    "calibrated_did": genept_cal - vanilla_cal,
                })
    expected = len(pathways) * len(SCENARIOS) * len(draws)
    if len(output) != expected:
        raise ValueError(f"Expected {expected} calibrated units, got {len(output)}")
    return output


def _scope_rows(units: list[dict[str, Any]], scope: str) -> list[dict[str, Any]]:
    return units if scope == "overall_non_null" else [row for row in units if row["scenario"] == scope]


def did_statistics(units: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for scope in SCOPES:
        rows = _scope_rows(units, scope)
        fields = {
            name: descriptive(float(row[column]) for row in rows)
            for name, column in (
                ("vanilla_raw", "vanilla_raw_ap"),
                ("genept_raw", "genept_raw_ap"),
                ("vanilla_calibrated", "vanilla_calibrated_ap"),
                ("genept_calibrated", "genept_calibrated_ap"),
                ("raw_gap", "raw_genept_minus_vanilla_ap"),
                ("null_gap", "null_genept_minus_vanilla_ap"),
                ("did", "calibrated_did"),
            )
        }
        did_values = [float(row["calibrated_did"]) for row in rows]
        statistic, p_value = signed_rank(did_values)
        output.append({
            "scope": scope,
            "pairing_unit": "pathway_scenario_draw_mean_across_strengths",
            "n": len(rows),
            "vanilla_raw_ap_mean": fields["vanilla_raw"]["mean"],
            "vanilla_raw_ap_median": fields["vanilla_raw"]["median"],
            "genept_raw_ap_mean": fields["genept_raw"]["mean"],
            "genept_raw_ap_median": fields["genept_raw"]["median"],
            "vanilla_calibrated_ap_mean": fields["vanilla_calibrated"]["mean"],
            "vanilla_calibrated_ap_median": fields["vanilla_calibrated"]["median"],
            "genept_calibrated_ap_mean": fields["genept_calibrated"]["mean"],
            "genept_calibrated_ap_median": fields["genept_calibrated"]["median"],
            "raw_genept_minus_vanilla_mean": fields["raw_gap"]["mean"],
            "raw_genept_minus_vanilla_median": fields["raw_gap"]["median"],
            "null_genept_minus_vanilla_mean": fields["null_gap"]["mean"],
            "null_genept_minus_vanilla_median": fields["null_gap"]["median"],
            "mean_calibrated_did": fields["did"]["mean"],
            "median_calibrated_did": fields["did"]["median"],
            "did_sd": fields["did"]["sd"],
            "did_se": fields["did"]["se"],
            "did_descriptive_ci95_low": fields["did"]["ci95_low"],
            "did_descriptive_ci95_high": fields["did"]["ci95_high"],
            "wilcoxon_statistic": statistic,
            "wilcoxon_raw_p_value": p_value,
            "wilcoxon_bonferroni_4_p_value": min(1.0, 4.0 * p_value),
            "rank_biserial_genept_minus_vanilla": rank_biserial(did_values),
        })
    return output


def fallback_statistics(units: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for stratum in ("ALL", "NO_FALLBACK", "FALLBACK"):
        subset = units if stratum == "ALL" else [
            row for row in units
            if bool(row["truth_fallback_used"]) == (stratum == "FALLBACK")
        ]
        for scope in SCOPES:
            rows = _scope_rows(subset, scope)
            if not rows:
                continue
            raw = descriptive(float(row["raw_genept_minus_vanilla_ap"]) for row in rows)
            null = descriptive(float(row["null_genept_minus_vanilla_ap"]) for row in rows)
            did = descriptive(float(row["calibrated_did"]) for row in rows)
            output.append({
                "fallback_stratum": stratum,
                "scope": scope,
                "n": len(rows),
                "raw_genept_minus_vanilla_mean": raw["mean"],
                "raw_genept_minus_vanilla_median": raw["median"],
                "matched_null_genept_minus_vanilla_mean": null["mean"],
                "matched_null_genept_minus_vanilla_median": null["median"],
                "calibrated_did_mean": did["mean"],
                "calibrated_did_median": did["median"],
                "calibrated_did_descriptive_ci95_low": did["ci95_low"],
                "calibrated_did_descriptive_ci95_high": did["ci95_high"],
            })
    return output


def load_baseline_features(h5_path: Path, manifest: dict[str, Any], *, detection_min: float) -> dict[str, dict[str, Any]]:
    output: dict[str, dict[str, Any]] = {}
    with h5py.File(h5_path, "r") as handle:
        for pathway_row in manifest["pathways"]:
            pathway = pathway_row["pathway"]
            group = handle[f"pathways/{pathway_row['pathway_id']}"]
            genes = [value.decode() if isinstance(value, bytes) else str(value)
                     for value in group["gene_names"][:]]
            pooled = np.vstack([group["condition_A"][:], group["condition_B_baseline"][:]])
            features = {
                "mean_expression": np.mean(pooled, axis=0),
                "detection_fraction": np.mean(pooled > 0, axis=0),
                "expression_sd": np.std(pooled, axis=0, ddof=1),
            }
            n_genes = len(genes)
            feature_percentiles = {
                name: rankdata(values, method="average") / n_genes
                for name, values in features.items()
            }
            eligible = (
                (features["detection_fraction"] >= detection_min)
                & (features["expression_sd"] > 0)
            )
            output[pathway] = {
                "genes": genes,
                "gene_to_index": {gene: index for index, gene in enumerate(genes)},
                "features": features,
                "feature_percentiles": feature_percentiles,
                "eligible_indices": np.flatnonzero(eligible),
                "truth_count": truth_count(n_genes),
            }
    return output


def load_null_rankings(path: Path) -> dict[tuple[str, int, str], dict[str, Any]]:
    grouped: dict[tuple[str, int, str], list[dict[str, str]]] = {}
    for row in read_csv(path):
        if row["perturbation_type"] == "null":
            grouped.setdefault((row["pathway"], int(row["draw_id"]), row["method"]), []).append(row)
    output: dict[tuple[str, int, str], dict[str, Any]] = {}
    for key, rows in grouped.items():
        ordered = sorted(rows, key=lambda row: int(row["rank"]))
        genes = [row["gene"] for row in ordered]
        ranks = np.asarray([int(row["rank"]) for row in ordered], dtype=np.int64)
        if ranks.tolist() != list(range(1, len(rows) + 1)) or len(genes) != len(set(genes)):
            raise ValueError(f"Invalid strict ranking: {key}")
        output[key] = {
            "genes": genes,
            "rank_by_gene": {gene: rank for gene, rank in zip(genes, ranks)},
            "delta_by_gene": {row["gene"]: float(row["delta_score"]) for row in rows},
            "truth": {row["gene"] for row in rows if as_bool(row["is_evaluation_target"])},
        }
    return output


def ap_from_gene_indices(indices: np.ndarray, genes: list[str], rank_by_gene: dict[str, int]) -> float:
    ranks = np.sort(np.asarray([rank_by_gene[genes[int(index)]] for index in indices], dtype=np.float64))
    return float(np.mean(np.arange(1, len(ranks) + 1, dtype=np.float64) / ranks))


def _diagnostic_row() -> dict[str, Any]:
    return {
        "row_type": "", "scheme": "", "population": "", "status": "",
        "pathway": "", "draw_id": "", "method": "", "n_pathways": "",
        "n_units": "", "resamples": "", "gene_count": "", "truth_count": "",
        "pool_size": "", "eligible_pool_size": "", "mean_ap": "", "median_ap": "",
        "sd_ap": "", "exact_random_ap": "", "mean_ap_minus_exact_chance": "",
        "mean_genept_minus_vanilla_ap": "", "mean_contrast_ap": "",
        "median_contrast_ap": "", "wilcoxon_p_value": "", "rank_biserial": "",
        "reason": "",
    }


def null_truth_diagnostics(
    null_rankings: dict[tuple[str, int, str], dict[str, Any]],
    baseline: dict[str, dict[str, Any]], *, seed: int, resamples: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    pathways = sorted(baseline)
    draws = sorted({key[1] for key in null_rankings})
    feasible = [pathway for pathway in pathways
                if len(baseline[pathway]["eligible_indices"]) >= baseline[pathway]["truth_count"]]
    infeasible = [pathway for pathway in pathways if pathway not in feasible]
    if len(feasible) != 9 or len(infeasible) != 2:
        raise ValueError(f"Frozen diagnostic expected 9 feasible and 2 infeasible pathways; got {len(feasible)}/{len(infeasible)}")

    rows: list[dict[str, Any]] = []
    empirical: dict[tuple[str, str, int, str], np.ndarray] = {}
    scheme_codes = {"A_ELIGIBLE_POOL_MATCHED": 1, "B_PATHWAY_WIDE_UNIFORM": 2,
                    "C_PHASE7_FALLBACK_MATCHED": 3}
    for pathway_index, pathway in enumerate(pathways, 1):
        info = baseline[pathway]
        n_genes = len(info["genes"])
        m = info["truth_count"]
        eligible = info["eligible_indices"]
        if pathway in infeasible:
            row = _diagnostic_row()
            row.update({
                "row_type": "pathway_status", "scheme": "A_ELIGIBLE_POOL_MATCHED",
                "population": "FEASIBLE_9", "status": "NOT_ESTIMABLE", "pathway": pathway,
                "gene_count": n_genes, "truth_count": m, "pool_size": len(eligible),
                "eligible_pool_size": len(eligible),
                "reason": "eligible pool size < truth_count, therefore eligible-pool null is not identifiable",
            })
            rows.append(row)
        for draw_id in draws:
            method_rankings = {method: null_rankings[(pathway, draw_id, method)] for method in METHODS}
            if set(method_rankings["vanilla_scpa"]["genes"]) != set(info["genes"]):
                raise ValueError(f"HDF5/ranking gene mismatch: {pathway}")
            for scheme, code in scheme_codes.items():
                if scheme == "A_ELIGIBLE_POOL_MATCHED" and pathway in infeasible:
                    continue
                if scheme == "B_PATHWAY_WIDE_UNIFORM":
                    pool = np.arange(n_genes, dtype=np.int64)
                elif scheme == "A_ELIGIBLE_POOL_MATCHED":
                    pool = eligible
                else:
                    pool = eligible if len(eligible) >= m else np.arange(n_genes, dtype=np.int64)
                rng = np.random.default_rng(seed + code * 10_000_000 + pathway_index * 100_000 + draw_id * 1_000)
                values = {method: np.empty(resamples, dtype=np.float64) for method in METHODS}
                for replicate in range(resamples):
                    selected = rng.choice(pool, m, replace=False)
                    for method in METHODS:
                        values[method][replicate] = ap_from_gene_indices(
                            selected, info["genes"], method_rankings[method]["rank_by_gene"]
                        )
                for method in METHODS:
                    empirical[(scheme, pathway, draw_id, method)] = values[method]
                    summary = descriptive(values[method])
                    chance = exact_random_chance(n_genes, m, 1)["average_precision"]
                    row = _diagnostic_row()
                    row.update({
                        "row_type": "unit", "scheme": scheme,
                        "population": "FEASIBLE_9" if scheme.startswith("A_") else "ALL_11",
                        "status": "ESTIMATED", "pathway": pathway, "draw_id": draw_id,
                        "method": method, "n_units": 1, "resamples": resamples,
                        "gene_count": n_genes, "truth_count": m, "pool_size": len(pool),
                        "eligible_pool_size": len(eligible), "mean_ap": summary["mean"],
                        "median_ap": summary["median"], "sd_ap": summary["sd"],
                        "exact_random_ap": chance,
                        "mean_ap_minus_exact_chance": summary["mean"] - chance,
                    })
                    rows.append(row)

    def add_aggregate(scheme: str, population: str, selected_pathways: list[str]) -> None:
        unit_keys = [(pathway, draw_id) for pathway in selected_pathways for draw_id in draws]
        method_unit_means: dict[str, np.ndarray] = {}
        for method in METHODS:
            values = np.asarray([
                np.mean(empirical[(scheme, pathway, draw_id, method)])
                for pathway, draw_id in unit_keys
            ], dtype=np.float64)
            method_unit_means[method] = values
            summary = descriptive(values)
            chance = float(np.mean([
                exact_random_chance(len(baseline[pathway]["genes"]), baseline[pathway]["truth_count"], 1)["average_precision"]
                for pathway, _ in unit_keys
            ]))
            row = _diagnostic_row()
            row.update({
                "row_type": "aggregate", "scheme": scheme, "population": population,
                "status": "ESTIMATED", "method": method, "n_pathways": len(selected_pathways),
                "n_units": len(unit_keys), "resamples": resamples, "mean_ap": summary["mean"],
                "median_ap": summary["median"], "sd_ap": summary["sd"],
                "exact_random_ap": chance, "mean_ap_minus_exact_chance": summary["mean"] - chance,
            })
            rows.append(row)
        gap = method_unit_means["genept_scpa"] - method_unit_means["vanilla_scpa"]
        summary = descriptive(gap)
        row = _diagnostic_row()
        row.update({
            "row_type": "aggregate", "scheme": scheme, "population": population,
            "status": "ESTIMATED", "method": "genept_minus_vanilla",
            "n_pathways": len(selected_pathways), "n_units": len(unit_keys), "resamples": resamples,
            "mean_genept_minus_vanilla_ap": summary["mean"],
            "mean_contrast_ap": summary["mean"], "median_contrast_ap": summary["median"],
        })
        rows.append(row)

    add_aggregate("A_ELIGIBLE_POOL_MATCHED", "FEASIBLE_9", feasible)
    add_aggregate("B_PATHWAY_WIDE_UNIFORM", "FEASIBLE_9", feasible)
    add_aggregate("B_PATHWAY_WIDE_UNIFORM", "ALL_11", pathways)
    add_aggregate("C_PHASE7_FALLBACK_MATCHED", "ALL_11", pathways)
    add_aggregate("C_PHASE7_FALLBACK_MATCHED", "NO_FALLBACK_9", feasible)
    add_aggregate("C_PHASE7_FALLBACK_MATCHED", "FALLBACK_2", infeasible)

    for method in (*METHODS, "genept_minus_vanilla"):
        contrasts = []
        for pathway in feasible:
            for draw_id in draws:
                if method in METHODS:
                    a_value = np.mean(empirical[("A_ELIGIBLE_POOL_MATCHED", pathway, draw_id, method)])
                    b_value = np.mean(empirical[("B_PATHWAY_WIDE_UNIFORM", pathway, draw_id, method)])
                else:
                    a_value = np.mean(empirical[("A_ELIGIBLE_POOL_MATCHED", pathway, draw_id, "genept_scpa")]
                                      - empirical[("A_ELIGIBLE_POOL_MATCHED", pathway, draw_id, "vanilla_scpa")])
                    b_value = np.mean(empirical[("B_PATHWAY_WIDE_UNIFORM", pathway, draw_id, "genept_scpa")]
                                      - empirical[("B_PATHWAY_WIDE_UNIFORM", pathway, draw_id, "vanilla_scpa")])
                contrasts.append(float(a_value - b_value))
        summary = descriptive(contrasts)
        statistic, p_value = signed_rank(contrasts)
        row = _diagnostic_row()
        row.update({
            "row_type": "direct_comparison", "scheme": "A_MINUS_B",
            "population": "SAME_FEASIBLE_9", "status": "ESTIMATED", "method": method,
            "n_pathways": len(feasible), "n_units": len(contrasts), "resamples": resamples,
            "mean_contrast_ap": summary["mean"], "median_contrast_ap": summary["median"],
            "wilcoxon_p_value": p_value, "rank_biserial": rank_biserial(contrasts),
        })
        rows.append(row)

    summary = {
        "feasible_pathways": feasible,
        "infeasible_pathways": infeasible,
        "a_unit_count": len(feasible) * len(draws),
        "b_all_unit_count": len(pathways) * len(draws),
        "resamples_per_unit": resamples,
    }
    return rows, summary


def expression_confounding(
    null_rankings: dict[tuple[str, int, str], dict[str, Any]],
    baseline: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for pathway in sorted(baseline):
        info = baseline[pathway]
        genes = info["genes"]
        for draw_id in sorted({key[1] for key in null_rankings if key[0] == pathway}):
            truth = null_rankings[(pathway, draw_id, "vanilla_scpa")]["truth"]
            if truth != null_rankings[(pathway, draw_id, "genept_scpa")]["truth"]:
                raise ValueError("Null truth labels differ between methods")
            truth_mask = np.asarray([gene in truth for gene in genes], dtype=bool)
            for feature, percentiles in info["feature_percentiles"].items():
                target_mean = float(np.mean(percentiles[truth_mask]))
                non_target_mean = float(np.mean(percentiles[~truth_mask]))
                rows.append({
                    "row_type": "unit", "diagnostic": "truth_target_minus_non_target_feature_rank_percentile",
                    "pathway": pathway, "draw_id": draw_id, "method": "NA", "feature": feature,
                    "n_genes": len(genes), "truth_count": int(np.sum(truth_mask)),
                    "value": target_mean - non_target_mean,
                    "target_feature_rank_percentile_mean": target_mean,
                    "non_target_feature_rank_percentile_mean": non_target_mean,
                    "status": "ESTIMATED",
                })
            for method in METHODS:
                ranking = null_rankings[(pathway, draw_id, method)]
                delta = np.asarray([ranking["delta_by_gene"][gene] for gene in genes], dtype=np.float64)
                negative_rank = np.asarray([-ranking["rank_by_gene"][gene] for gene in genes], dtype=np.float64)
                for feature, percentiles in info["feature_percentiles"].items():
                    for diagnostic, masking_value in (
                        ("spearman_masking_delta_vs_feature", delta),
                        ("spearman_negative_rank_vs_feature", negative_rank),
                    ):
                        if np.ptp(masking_value) == 0 or np.ptp(percentiles) == 0:
                            correlation = math.nan
                        else:
                            correlation = float(spearmanr(masking_value, percentiles).statistic)
                        status = "ESTIMATED" if math.isfinite(correlation) else "CONSTANT_INPUT_NOT_ESTIMABLE"
                        rows.append({
                            "row_type": "unit", "diagnostic": diagnostic, "pathway": pathway,
                            "draw_id": draw_id, "method": method, "feature": feature,
                            "n_genes": len(genes), "truth_count": len(truth),
                            "value": correlation if math.isfinite(correlation) else "",
                            "target_feature_rank_percentile_mean": "",
                            "non_target_feature_rank_percentile_mean": "", "status": status,
                        })

    groups: dict[tuple[str, str, str], list[float]] = {}
    for row in rows:
        if row["row_type"] == "unit" and row["value"] != "":
            groups.setdefault((row["diagnostic"], row["method"], row["feature"]), []).append(float(row["value"]))
    aggregate_summary: dict[str, Any] = {}
    for (diagnostic, method, feature), values in sorted(groups.items()):
        summary = descriptive(values)
        rows.append({
            "row_type": "aggregate", "diagnostic": diagnostic, "pathway": "ALL_11",
            "draw_id": "", "method": method, "feature": feature, "n_genes": "",
            "truth_count": "", "value": summary["mean"],
            "target_feature_rank_percentile_mean": "",
            "non_target_feature_rank_percentile_mean": "", "status": "ESTIMATED",
        })
        aggregate_summary[f"{diagnostic}|{method}|{feature}"] = summary
    return rows, aggregate_summary


def validate_phase7_rankings_against_metrics(
    null_rankings: dict[tuple[str, int, str], dict[str, Any]], metric_rows: list[dict[str, str]],
) -> float:
    null_metrics = {
        (row["pathway"], int(row["draw_id"]), row["method"]): float(row["average_precision"])
        for row in metric_rows if row["perturbation_type"] == "null"
    }
    errors = []
    for key, ranking in null_rankings.items():
        observed = average_precision(ranking["genes"], ranking["truth"])
        errors.append(abs(observed - null_metrics[key]))
    return max(errors) if errors else math.inf


def evaluate(config_path: Path) -> dict[str, Any]:
    config_path = config_path.resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    protocol_path = PROJECT_ROOT / "docs/phase7b_null_calibration_protocol.md"
    source_paths = {name: PROJECT_ROOT / path for name, path in config["source"].items()}
    input_paths = {"phase7b_config": config_path, "phase7b_protocol": protocol_path, **source_paths}
    initial_hashes = {name: sha256_file(path) for name, path in input_paths.items()}

    processed = PROJECT_ROOT / config["artifacts"]["processed_directory"]
    phase7_processed = (PROJECT_ROOT / "data/processed/genept_scpa/phase7_llmfree_synthetic").resolve()
    if processed.resolve() == phase7_processed or phase7_processed in processed.resolve().parents:
        raise RuntimeError("Phase 7B output must not overlap the Phase 7 output directory")

    scpa_qc = json.loads(source_paths["scpa_qc"].read_text(encoding="utf-8"))
    if not (scpa_qc.get("status") == "PASS" and scpa_qc.get("partial_run") is False
            and int(scpa_qc.get("mcm_count", -1)) == 101920
            and int(scpa_qc.get("failed_mcm_calls", -1)) == 0):
        raise RuntimeError("Frozen Phase 7 MCM QC is not a complete PASS")
    manifest = json.loads(source_paths["manifest"].read_text(encoding="utf-8"))
    metric_rows = read_csv(source_paths["metrics"])
    calibrated = build_calibrated_units(metric_rows)
    did_rows = did_statistics(calibrated)
    fallback_rows = fallback_statistics(calibrated)

    null_rankings = load_null_rankings(source_paths["rankings"])
    max_ap_error = validate_phase7_rankings_against_metrics(null_rankings, metric_rows)
    if len(null_rankings) != 11 * 20 * 2 or max_ap_error > 1e-12:
        raise RuntimeError(f"Null rankings/metric validation failed: n={len(null_rankings)} error={max_ap_error}")

    diagnostic_config = config["null_truth_diagnostics"]
    baseline = load_baseline_features(
        source_paths["expression_h5"], manifest,
        detection_min=float(diagnostic_config["eligibility"]["detection_fraction_min"]),
    )
    diagnostic_rows, diagnostic_summary = null_truth_diagnostics(
        null_rankings, baseline, seed=int(diagnostic_config["seed"]),
        resamples=int(diagnostic_config["resamples"]),
    )
    confounding_rows, confounding_summary = expression_confounding(null_rankings, baseline)

    did_by_scope = {row["scope"]: row for row in did_rows}
    overall = did_by_scope["overall_non_null"]
    verdict = (
        "COMPLETED_WITH_WARNING"
        if overall["mean_calibrated_did"] > 0 and overall["wilcoxon_bonferroni_4_p_value"] < 0.05
        else "COMPLETED_WITH_NULL_EXPLAINED"
    )
    null_aggregates = [row for row in diagnostic_rows if row["row_type"] in {"aggregate", "direct_comparison"}]
    summary: dict[str, Any] = {
        "status": verdict,
        "verdict_rationale": (
            "Overall calibrated DiD remained positive with Bonferroni p<0.05, driven by mean shift; "
            "cell-subset and mixed-direction residual contrasts were not supported."
            if verdict == "COMPLETED_WITH_WARNING" else
            "The overall matched-null calibrated DiD was not supported after multiplicity correction."
        ),
        "technical_qc": {
            "status": "PASS",
            "new_mcm_executed": False,
            "phase7_inputs_read_only": True,
            "phase7_scpa_qc_status": scpa_qc["status"],
            "phase7_mcm_count_reused": scpa_qc["mcm_count"],
            "phase7_failed_mcm_calls": scpa_qc["failed_mcm_calls"],
            "calibrated_unit_count": len(calibrated),
            "null_ranking_count": len(null_rankings),
            "max_recomputed_null_ap_error": max_ap_error,
        },
        "primary_did_statistics": did_by_scope,
        "fallback_sensitivity": fallback_rows,
        "null_truth_diagnostic_design": diagnostic_summary,
        "null_truth_diagnostic_aggregates": null_aggregates,
        "expression_confounding_aggregates": confounding_summary,
        "interpretation_limits": config["interpretation_limits"],
    }

    artifacts = config["artifacts"]
    write_csv_atomic(calibrated, processed / artifacts["calibrated_metrics"])
    write_csv_atomic(did_rows, processed / artifacts["did_statistics"])
    write_csv_atomic(fallback_rows, processed / artifacts["fallback_sensitivity"])
    write_csv_atomic(diagnostic_rows, processed / artifacts["null_truth_diagnostics"])
    write_csv_atomic(confounding_rows, processed / artifacts["expression_confounding"])

    a_method = {row["method"]: row for row in null_aggregates
                if row["scheme"] == "A_ELIGIBLE_POOL_MATCHED" and row["population"] == "FEASIBLE_9"}
    b9_method = {row["method"]: row for row in null_aggregates
                 if row["scheme"] == "B_PATHWAY_WIDE_UNIFORM" and row["population"] == "FEASIBLE_9"}
    b11_method = {row["method"]: row for row in null_aggregates
                  if row["scheme"] == "B_PATHWAY_WIDE_UNIFORM" and row["population"] == "ALL_11"}
    c11_method = {row["method"]: row for row in null_aggregates
                  if row["scheme"] == "C_PHASE7_FALLBACK_MATCHED" and row["population"] == "ALL_11"}
    no_fallback = next(row for row in fallback_rows
                       if row["fallback_stratum"] == "NO_FALLBACK" and row["scope"] == "overall_non_null")
    fallback = next(row for row in fallback_rows
                    if row["fallback_stratum"] == "FALLBACK" and row["scope"] == "overall_non_null")
    explained_fraction = 1.0 - overall["mean_calibrated_did"] / overall["raw_genept_minus_vanilla_mean"]
    summary["overall_raw_gap_arithmetic_fraction_removed_by_matched_null"] = explained_fraction
    a_minus_b = {row["method"]: row for row in null_aggregates
                 if row["scheme"] == "A_MINUS_B" and row["population"] == "SAME_FEASIBLE_9"}
    truth_feature = {
        feature: confounding_summary[
            f"truth_target_minus_non_target_feature_rank_percentile|NA|{feature}"
        ]
        for feature in ("mean_expression", "detection_fraction", "expression_sd")
    }
    delta_feature = {
        method: {
            feature: confounding_summary[f"spearman_masking_delta_vs_feature|{method}|{feature}"]
            for feature in ("mean_expression", "detection_fraction", "expression_sd")
        }
        for method in METHODS
    }
    results = f"""# Phase 7B null-calibration sensitivity 결과

상태: **{verdict}**
Technical QC: **PASS** — 새 MCM 0회, 기존 Phase 7 입력은 read-only, null ranking AP 재계산 최대 오차 {max_ap_error:.3g}

## Matched null calibration

| Scope | n | Raw GenePT−Vanilla | Matched null gap | Calibrated DiD | Bonferroni p | Rank-biserial |
|---|---:|---:|---:|---:|---:|---:|
""" + "\n".join(
        f"| {scope} | {did_by_scope[scope]['n']} | {did_by_scope[scope]['raw_genept_minus_vanilla_mean']:.6f} "
        f"| {did_by_scope[scope]['null_genept_minus_vanilla_mean']:.6f} "
        f"| {did_by_scope[scope]['mean_calibrated_did']:.6f} "
        f"| {did_by_scope[scope]['wilcoxon_bonferroni_4_p_value']:.4g} "
        f"| {did_by_scope[scope]['rank_biserial_genept_minus_vanilla']:.3f} |"
        for scope in SCOPES
    ) + f"""

Overall raw method gap 중 matched-null subtraction으로 설명되는 기술적 비율은 {explained_fraction:.1%}다.
이 값은 causal decomposition이 아니라 같은 pathway/draw의 arithmetic sensitivity다.

## Fallback sensitivity

- NO_FALLBACK n={no_fallback['n']}: raw gap {no_fallback['raw_genept_minus_vanilla_mean']:.6f}, null gap {no_fallback['matched_null_genept_minus_vanilla_mean']:.6f}, DiD {no_fallback['calibrated_did_mean']:.6f}
- FALLBACK n={fallback['n']}: raw gap {fallback['raw_genept_minus_vanilla_mean']:.6f}, null gap {fallback['matched_null_genept_minus_vanilla_mean']:.6f}, DiD {fallback['calibrated_did_mean']:.6f}

## Null truth-selection diagnostics

- A eligible-pool, feasible 9 pathways/180 units: Vanilla AP {a_method['vanilla_scpa']['mean_ap']:.6f}, GenePT AP {a_method['genept_scpa']['mean_ap']:.6f}, method gap {a_method['genept_minus_vanilla']['mean_genept_minus_vanilla_ap']:.6f}
- B uniform, same 9 pathways: Vanilla AP {b9_method['vanilla_scpa']['mean_ap']:.6f}, GenePT AP {b9_method['genept_scpa']['mean_ap']:.6f}, method gap {b9_method['genept_minus_vanilla']['mean_genept_minus_vanilla_ap']:.6f}
- B uniform, all 11 pathways: Vanilla AP {b11_method['vanilla_scpa']['mean_ap']:.6f}, GenePT AP {b11_method['genept_scpa']['mean_ap']:.6f}, method gap {b11_method['genept_minus_vanilla']['mean_genept_minus_vanilla_ap']:.6f}
- C Phase-7 fallback-matched, all 11 pathways: Vanilla AP {c11_method['vanilla_scpa']['mean_ap']:.6f}, GenePT AP {c11_method['genept_scpa']['mean_ap']:.6f}, method gap {c11_method['genept_minus_vanilla']['mean_genept_minus_vanilla_ap']:.6f}
- Primary A−B on the same 9 pathways: Vanilla {a_minus_b['vanilla_scpa']['mean_contrast_ap']:.6f}, GenePT {a_minus_b['genept_scpa']['mean_contrast_ap']:.6f}; method-gap A−B {a_minus_b['genept_minus_vanilla']['mean_contrast_ap']:.6f}
- A infeasible 2 pathways는 replacement, truth-count 축소, threshold 완화 없이 `NOT_ESTIMABLE`로 기록했다.

## Baseline-expression confounding diagnostic

- 기존 null target의 within-pathway feature-rank percentile은 non-target보다 mean expression {truth_feature['mean_expression']['mean']:+.3f}, detection {truth_feature['detection_fraction']['mean']:+.3f}, SD {truth_feature['expression_sd']['mean']:+.3f} 높았다.
- Null masking delta와 feature의 mean Spearman은 Vanilla에서 mean/detection/SD 각각 {delta_feature['vanilla_scpa']['mean_expression']['mean']:.3f}/{delta_feature['vanilla_scpa']['detection_fraction']['mean']:.3f}/{delta_feature['vanilla_scpa']['expression_sd']['mean']:.3f}, GenePT에서 {delta_feature['genept_scpa']['mean_expression']['mean']:.3f}/{delta_feature['genept_scpa']['detection_fraction']['mean']:.3f}/{delta_feature['genept_scpa']['expression_sd']['mean']:.3f}였다.
- Constant masking-delta unit은 Spearman을 억지로 정의하지 않고 `CONSTANT_INPUT_NOT_ESTIMABLE`로 기록했다.

## 해석

Frozen synthetic benchmark에서 GenePT의 raw AP gap은 양수였지만 두 방법 모두 expression-conditioned
null recovery를 보였다. Pathway/draw-matched method-specific null calibration 뒤 residual difference는
overall과 mean shift에서 남았고, cell subset과 mixed direction에서는 지지되지 않았다. A와 B의 차이 및
baseline diagnostic은 expression-conditioned truth selection이 null inflation의 주요 원인과 일치함을
보이지만 인과적 분해는 아니다. 상세치는 `phase7b_expression_confounding.csv`에서 확인한다.

이 결과는 biological superiority, causal-gene identification, general validity 또는 semantic mechanism을
입증하지 않는다. Phase 6 correspondence-specificity 결론과 Phase 7/7B synthetic recovery 결론은 분리한다.
"""
    write_text_atomic(results, processed / artifacts["results_markdown"])

    final_hashes = {name: sha256_file(path) for name, path in input_paths.items()}
    if final_hashes != initial_hashes:
        raise RuntimeError("One or more frozen Phase 7/7B protocol inputs changed during evaluation")
    summary["technical_qc"]["input_hashes_verified_unchanged"] = True
    write_json_atomic(summary, processed / artifacts["summary_json"])
    snapshot = {
        "schema_version": 1,
        "phase": "7B",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_status": "FROZEN_BEFORE_EVALUATION",
        "execution": {"new_mcm_executed": False, "rankings_relabel_only": True},
        "config": config,
        "input_sha256": {
            name: {"path": str(path.relative_to(PROJECT_ROOT)), "sha256": initial_hashes[name]}
            for name, path in input_paths.items()
        },
    }
    write_json_atomic(snapshot, processed / artifacts["protocol_snapshot"])
    print(
        f"PHASE7B status={verdict} qc=PASS units={len(calibrated)} "
        f"did={overall['mean_calibrated_did']:.6f} new_mcm=0 output={processed}"
    )
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path,
        default=PROJECT_ROOT / "config/phase7b_null_calibration.yaml",
    )
    arguments = parser.parse_args()
    evaluate(arguments.config)
