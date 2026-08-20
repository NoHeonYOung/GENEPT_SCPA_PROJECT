#!/usr/bin/env python3
"""Join masking checkpoints to truth, evaluate, aggregate, and report Phase 7."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Iterable

import numpy as np
from scipy.stats import rankdata, wilcoxon
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from gene_embedding_project.genept_scpa.io import write_json_atomic  # noqa: E402
from gene_embedding_project.genept_scpa.phase7.synthetic_benchmark_llmfree.metrics import (  # noqa: E402
    average_precision, exact_random_chance, ndcg_at_k, recall_at_k,
)


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
        writer.writeheader(); writer.writerows(rows)
    os.replace(temporary, path)


def as_bool(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes"}


def mean_sd(values: Iterable[float]) -> tuple[float, float, int]:
    array = np.asarray(list(values), dtype=np.float64)
    return float(np.mean(array)), float(np.std(array, ddof=1)) if len(array) > 1 else 0.0, len(array)


def rank_biserial(differences: np.ndarray) -> float:
    nonzero = differences[np.abs(differences) > 0]
    if len(nonzero) == 0:
        return 0.0
    ranks = rankdata(np.abs(nonzero), method="average")
    positive = float(np.sum(ranks[nonzero > 0]))
    negative = float(np.sum(ranks[nonzero < 0]))
    return (positive - negative) / (positive + negative)


def paired_test(rows: list[dict[str, Any]], scenario: str) -> dict[str, Any]:
    selected = [row for row in rows if row["perturbation_type"] != "null"
                and (scenario == "ALL_NON_NULL" or row["perturbation_type"] == scenario)]
    # Strengths share the same truth draw and are averaged before inference to
    # avoid counting the 0.5x/1.0x pair as independent replicates.
    grouped: dict[tuple[Any, ...], dict[str, list[float]]] = {}
    for row in selected:
        key = (row["pathway"], row["perturbation_type"], row["draw_id"])
        grouped.setdefault(key, {}).setdefault(row["method"], []).append(float(row["average_precision"]))
    differences = []
    for methods in grouped.values():
        if set(methods) != {"vanilla_scpa", "genept_scpa"}:
            raise ValueError("Incomplete paired method results")
        differences.append(np.mean(methods["genept_scpa"]) - np.mean(methods["vanilla_scpa"]))
    values = np.asarray(differences, dtype=np.float64)
    if len(values) == 0:
        raise ValueError("No non-null AP pairs")
    if np.all(values == 0):
        statistic, p_value = 0.0, 1.0
    else:
        test = wilcoxon(values, zero_method="wilcox", alternative="two-sided", method="auto")
        statistic, p_value = float(test.statistic), float(test.pvalue)
    return {
        "scope": scenario, "pairing_unit": "pathway_scenario_draw_mean_across_strengths",
        "pair_count": len(values), "mean_genept_minus_vanilla_ap": float(np.mean(values)),
        "wilcoxon_statistic": statistic, "raw_p_value": p_value,
        "rank_biserial_effect_genept_minus_vanilla": rank_biserial(values),
    }


def aggregate(metrics: list[dict[str, Any]]) -> list[dict[str, Any]]:
    keys = ["recall_at_3", "recall_at_5", "recall_at_10", "average_precision",
            "ndcg_at_3", "ndcg_at_5", "ndcg_at_10"]
    groups: list[tuple[str, str, str, list[dict[str, Any]]]] = []
    for method in ("vanilla_scpa", "genept_scpa"):
        method_rows = [row for row in metrics if row["method"] == method]
        groups.append(("all_non_null", "ALL_NON_NULL", "ALL", [row for row in method_rows if row["perturbation_type"] != "null"]))
        for fallback in (False, True):
            groups.append((
                "truth_fallback", "FALLBACK" if fallback else "NO_FALLBACK", "ALL",
                [row for row in method_rows
                 if row["perturbation_type"] != "null" and row["truth_fallback_used"] is fallback],
            ))
        for scenario in ("null", "mean_shift", "cell_subset", "mixed_direction"):
            scenario_rows = [row for row in method_rows if row["perturbation_type"] == scenario]
            groups.append(("scenario", scenario, "ALL", scenario_rows))
            for strength in sorted({row["perturbation_strength"] for row in scenario_rows}):
                groups.append(("scenario_strength", scenario, str(strength),
                               [row for row in scenario_rows if row["perturbation_strength"] == strength]))
    output = []
    for scope, scenario, strength, rows in groups:
        if not rows:
            continue
        base: dict[str, Any] = {"scope": scope, "scenario": scenario,
                                "strength": strength, "method": rows[0]["method"], "n": len(rows)}
        for metric in keys:
            mean, sd, _ = mean_sd(float(row[metric]) for row in rows)
            base[f"{metric}_mean"] = mean; base[f"{metric}_sd"] = sd
        output.append(base)
    return output


def evaluate(config_path: Path, *, allow_partial: bool = False) -> dict[str, Any]:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    artifacts = config["artifacts"]
    interim = PROJECT_ROOT / artifacts["interim_directory"]
    processed = PROJECT_ROOT / artifacts["processed_directory"]
    manifest = json.loads((processed / artifacts["manifest"]).read_text(encoding="utf-8"))
    truth_rows = read_csv(processed / artifacts["ground_truth"])
    truth_by_experiment: dict[str, set[str]] = {}
    universe_by_experiment: dict[str, set[str]] = {}
    truth_meta: dict[str, dict[str, str]] = {}
    fallback_by_experiment: dict[str, bool] = {}
    for row in truth_rows:
        experiment_id = row["experiment_id"]
        universe_by_experiment.setdefault(experiment_id, set()).add(row["gene"])
        truth_meta[experiment_id] = row
        fallback_by_experiment[experiment_id] = as_bool(row["truth_fallback_used"])
        if as_bool(row["is_evaluation_target"]):
            truth_by_experiment.setdefault(experiment_id, set()).add(row["gene"])
    expected_ids = {row["experiment_id"] for row in manifest["experiments"]}
    if not allow_partial and (
        int(manifest.get("experiment_count", 0)) != 1540
        or int(manifest.get("draw_count", 0)) < 20
    ):
        raise RuntimeError("Production evaluation requires the full frozen 1,540-experiment manifest")
    checkpoint_dir = interim / artifacts["checkpoints"]
    checkpoint_paths = sorted(checkpoint_dir.glob("*_masking.csv"))
    found_ids = {path.name.removesuffix("_masking.csv") for path in checkpoint_paths}
    missing = expected_ids - found_ids
    if missing and not allow_partial:
        raise RuntimeError(f"Missing {len(missing)} Phase 7 checkpoints; resume masking first")

    ranking_rows: list[dict[str, Any]] = []
    metrics: list[dict[str, Any]] = []
    k_values = [int(value) for value in config["evaluation"]["recall_k"]]
    for checkpoint in checkpoint_paths:
        rows = read_csv(checkpoint)
        if not rows:
            raise ValueError(f"Empty checkpoint: {checkpoint}")
        experiment_id = rows[0]["experiment_id"]
        if experiment_id not in expected_ids:
            raise ValueError(f"Unknown checkpoint experiment: {experiment_id}")
        genes = {row["gene"] for row in rows}
        if genes != universe_by_experiment[experiment_id]:
            raise ValueError(f"Checkpoint gene universe mismatch: {experiment_id}")
        truth = truth_by_experiment.get(experiment_id, set())
        if not truth:
            raise ValueError(f"Missing evaluation targets: {experiment_id}")
        meta = truth_meta[experiment_id]
        for method in ("vanilla", "genept"):
            delta_field = f"{method}_delta_score"
            ordered_rows = sorted(rows, key=lambda row: (-float(row[delta_field]), row["gene"]))
            ranking = [row["gene"] for row in ordered_rows]
            for rank, row in enumerate(ordered_rows, 1):
                ranking_rows.append({
                    "experiment_id": experiment_id, "draw_id": int(meta["draw_id"]),
                    "pathway": meta["pathway"], "perturbation_type": meta["perturbation_type"],
                    "perturbation_strength": float(meta["perturbation_strength"]),
                    "perturbation_seed": int(meta["perturbation_seed"]),
                    "truth_fallback_used": fallback_by_experiment[experiment_id],
                    "method": f"{method}_scpa", "gene": row["gene"], "rank": rank,
                    "delta_score": float(row[delta_field]),
                    "is_evaluation_target": row["gene"] in truth,
                    "is_ground_truth_perturbed": (
                        row["gene"] in truth and meta["perturbation_type"] != "null"
                    ),
                })
            metric_row: dict[str, Any] = {
                "experiment_id": experiment_id, "draw_id": int(meta["draw_id"]),
                "pathway": meta["pathway"], "source_database": meta["source_database"],
                "perturbation_type": meta["perturbation_type"],
                "perturbation_strength": float(meta["perturbation_strength"]),
                "perturbation_seed": int(meta["perturbation_seed"]),
                "method": f"{method}_scpa", "gene_count": len(ranking),
                "truth_fallback_used": fallback_by_experiment[experiment_id],
                "truth_count": len(truth), "average_precision": average_precision(ranking, truth),
            }
            for k in k_values:
                metric_row[f"recall_at_{k}"] = recall_at_k(ranking, truth, k)
                metric_row[f"ndcg_at_{k}"] = ndcg_at_k(ranking, truth, k)
                chance = exact_random_chance(len(ranking), len(truth), k)
                metric_row[f"chance_recall_at_{k}"] = chance["recall"]
                metric_row[f"chance_ndcg_at_{k}"] = chance["ndcg"]
            metric_row["chance_average_precision"] = exact_random_chance(
                len(ranking), len(truth), k_values[0]
            )["average_precision"]
            metrics.append(metric_row)
    if not metrics:
        raise ValueError("No completed masking checkpoints")

    aggregates = aggregate(metrics)
    overall_test = paired_test(metrics, "ALL_NON_NULL")
    overall_test["bonferroni_p_value"] = ""
    tests = [overall_test]
    required_scenarios = ("mean_shift", "cell_subset", "mixed_direction")
    present_scenarios = {row["perturbation_type"] for row in metrics}
    if not allow_partial and not set(required_scenarios) <= present_scenarios:
        raise RuntimeError("Production evaluation is missing one or more non-null scenarios")
    scenario_tests = [paired_test(metrics, scenario) for scenario in required_scenarios
                      if scenario in present_scenarios]
    for test in scenario_tests:
        test["bonferroni_p_value"] = min(1.0, 3 * test["raw_p_value"])
    tests.extend(scenario_tests)

    null_checks = []
    null_warning = False
    for method in ("vanilla_scpa", "genept_scpa"):
        rows = [row for row in metrics if row["method"] == method and row["perturbation_type"] == "null"]
        if not rows:
            if not allow_partial:
                raise RuntimeError("Production evaluation requires null checkpoints")
            null_warning = True
            null_checks.append({"method": method, "n": 0,
                                "mean_ap_minus_exact_chance": "",
                                "standard_error": "", "warning": True,
                                "reason": "null_not_evaluated_in_partial_run"})
            continue
        differences = np.asarray([row["average_precision"] - row["chance_average_precision"] for row in rows])
        mean_difference = float(np.mean(differences))
        se = float(np.std(differences, ddof=1) / math.sqrt(len(differences))) if len(differences) > 1 else 0.0
        warning = abs(mean_difference) > 2 * se if se > 0 else mean_difference != 0
        null_warning |= warning
        null_checks.append({"method": method, "n": len(rows),
                            "mean_ap_minus_exact_chance": mean_difference,
                            "standard_error": se, "warning": warning})

    write_csv_atomic(ranking_rows, processed / artifacts["rankings_csv"])
    write_csv_atomic(metrics, processed / artifacts["metrics_csv"])
    write_csv_atomic(aggregates, processed / artifacts["aggregate_csv"])
    write_csv_atomic(tests, processed / artifacts["statistics_csv"])
    payload = {"status": "PASS_WITH_NULL_WARNING" if null_warning else "PASS",
               "partial": bool(allow_partial or missing),
               "completed_experiments": len({row['experiment_id'] for row in metrics}),
               "metrics": metrics, "aggregates": aggregates, "paired_tests": tests,
               "null_checks": null_checks,
               "interpretation_limits": config["interpretation"]}
    write_json_atomic(payload, processed / artifacts["metrics_json"])

    nonnull = [row for row in metrics if row["perturbation_type"] != "null"]
    method_means = {method: float(np.mean([row["average_precision"] for row in nonnull
                                          if row["method"] == method]))
                    for method in ("vanilla_scpa", "genept_scpa")}
    winner = max(method_means, key=method_means.get)
    overall = tests[0]
    significant = overall["raw_p_value"] < 0.05
    fallback_experiment_count = len({row["experiment_id"] for row in metrics
                                     if row["truth_fallback_used"]})
    lead = (
        f"비-null synthetic benchmark에서 평균 AP는 {winner}가 더 높았고 "
        f"(Vanilla={method_means['vanilla_scpa']:.4f}, GenePT={method_means['genept_scpa']:.4f}), "
        f"paired Wilcoxon은 {'유의했다' if significant else '유의하지 않았다'} "
        f"(raw p={overall['raw_p_value']:.4g}, rank-biserial={overall['rank_biserial_effect_genept_minus_vanilla']:.3f})."
    )
    summary = f"""# Phase 7 LLM-free synthetic benchmark 결과

{lead}

- 상태: **{payload['status']}**
- 완료 experiment: {payload['completed_experiments']} / {len(expected_ids)}
- 방법: Vanilla SCPA zero-mask vs GenePT non-L2 exact subtraction mask
- 순위 신호: `-log10(raw_p_full) - -log10(raw_p_masked)` 내림차순
- 추론 단위: pathway × scenario × draw에서 두 strength AP를 평균한 paired unit
- Null 경고: **{null_warning}** — 경고가 있으면 비-null 해석보다 방법론 점검이 우선
- Truth-pool fallback experiment: {fallback_experiment_count}; FALLBACK/NO_FALLBACK 집계를 함께 확인

## 해석 제한

- 이 결과는 Naïve CD4 0h 한 집단과 동결된 synthetic perturbation에서의 injected-signal recovery만 뜻한다.
- GenePT의 biological superiority, 일반적 우월성, 인과성을 주장하지 않는다.
- Phase 6은 correct gene↔embedding correspondence의 특이성을 검정한 별도 질문이며 결론을 합치지 않는다.
- GenePT 열세가 나오더라도 pathway K genes와 1,536D projection의 geometry 차이가 confound일 수 있다.
"""
    (processed / artifacts["summary"]).write_text(summary, encoding="utf-8")
    print(f"PHASE7_EVAL status={payload['status']} metrics={len(metrics)} null_warning={null_warning}")
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path,
                        default=PROJECT_ROOT / "config/phase7_llmfree_synthetic.yaml")
    parser.add_argument("--allow-partial", action="store_true",
                        help="Smoke/debug only; production report must omit this flag")
    args = parser.parse_args()
    evaluate(args.config, allow_partial=args.allow_partial)
