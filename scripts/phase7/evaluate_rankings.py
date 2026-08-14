#!/usr/bin/env python3
"""Evaluate Phase 7 rankings; this is the only runner allowed to read truth labels."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from scripts.phase4.run_pathway_comparison import write_csv_atomic  # noqa: E402
from gene_embedding_project.genept_scpa.phase7.evaluation import (  # noqa: E402
    evaluate_ranking, prompt_order_spearman,
)
from gene_embedding_project.genept_scpa.phase7.ranking import aggregate_rankings  # noqa: E402


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _truth_tables(
    truth_path: Path,
) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    truth_by_experiment: dict[str, set[str]] = {}
    universe_by_experiment: dict[str, set[str]] = {}
    for row in _read_csv(truth_path):
        experiment_id = row["experiment_id"]
        universe_by_experiment.setdefault(experiment_id, set()).add(row["gene"])
        if row["is_ground_truth_perturbed"].strip().lower() in {"true", "1"}:
            truth_by_experiment.setdefault(experiment_id, set()).add(row["gene"])
    return truth_by_experiment, universe_by_experiment


def _metric_or_null(
    ranking: list[str], truth: set[str]
) -> dict[str, Any]:
    if truth:
        return evaluate_ranking(ranking, truth)
    return {
        "truth_k": 0, "recall_at_truth_k": "", "average_precision": "",
        "ndcg_at_n": "", "ndcg_at_truth_k": "",
    }


def evaluate_directory(
    response_dir: Path,
    request_dir: Path,
    mapping_dir: Path,
    truth_path: Path,
    output_path: Path,
    *,
    allow_mock_toy: bool = False,
) -> list[dict[str, Any]]:
    truth_by_experiment, universe_by_experiment = _truth_tables(truth_path)

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    backends: dict[tuple[str, str], str] = {}
    for artifact_path in sorted(response_dir.glob("*.json")):
        artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
        if artifact.get("scientific_evaluation_allowed") is False and not allow_mock_toy:
            raise RuntimeError("Mock rankings are forbidden from scientific evaluation")
        request = json.loads((request_dir / artifact_path.name).read_text(encoding="utf-8"))
        key = (request["experiment_id"], request["prompt_condition"])
        grouped.setdefault(key, []).append(artifact["response"])
        backends[key] = artifact["backend"]

    metric_rows: list[dict[str, Any]] = []
    for (experiment_id, prompt_condition), responses in sorted(grouped.items()):
        mapping = json.loads((mapping_dir / f"{experiment_id}.json").read_text(encoding="utf-8"))
        candidate_to_gene = mapping["candidate_to_gene"]
        aggregate = aggregate_rankings(responses, tie_seed=20290814)
        candidate_order = [row["candidate_id"] for row in aggregate]
        gene_order = [candidate_to_gene[candidate] for candidate in candidate_order]
        if set(gene_order) != universe_by_experiment.get(experiment_id, set()):
            raise ValueError(f"LLM candidate universe mismatch: {experiment_id}")
        truth = truth_by_experiment.get(experiment_id, set())
        row: dict[str, Any] = {
            "experiment_id": experiment_id,
            "prompt_condition": prompt_condition,
            "backend": backends[(experiment_id, prompt_condition)],
            "run_count": len(responses),
            "scientific_evaluation": False if allow_mock_toy else True,
        }
        row.update(_metric_or_null(gene_order, truth))
        stability = prompt_order_spearman(responses)
        row["prompt_order_mean_spearman"] = stability["mean_spearman"]
        row["prompt_order_min_spearman"] = stability["min_spearman"]
        metric_rows.append(row)
    if not metric_rows:
        raise ValueError("No rankings were evaluated")
    write_csv_atomic(metric_rows, list(metric_rows[0]), output_path)
    return metric_rows


def evaluate_scpa_directory(
    checkpoint_dir: Path,
    truth_path: Path,
    output_path: Path,
) -> list[dict[str, Any]]:
    """Evaluate completed R checkpoints; this function alone joins them to truth."""

    truth_by_experiment, universe_by_experiment = _truth_tables(truth_path)
    metric_rows: list[dict[str, Any]] = []
    for checkpoint in sorted(checkpoint_dir.glob("*_masking.csv")):
        rows = _read_csv(checkpoint)
        if not rows:
            raise ValueError(f"Empty SCPA checkpoint: {checkpoint}")
        experiment_ids = {row["experiment_id"] for row in rows}
        if len(experiment_ids) != 1:
            raise ValueError(f"Mixed experiment IDs in checkpoint: {checkpoint}")
        experiment_id = next(iter(experiment_ids))
        genes = {row["gene"] for row in rows}
        if genes != universe_by_experiment.get(experiment_id, set()):
            raise ValueError(f"SCPA gene universe mismatch: {experiment_id}")
        truth = truth_by_experiment.get(experiment_id, set())
        for method in ("vanilla", "genept"):
            rank_field = f"{method}_signed_rank"
            ordered = sorted(rows, key=lambda row: (float(row[rank_field]), row["gene"]))
            gene_order = [row["gene"] for row in ordered]
            result: dict[str, Any] = {
                "experiment_id": experiment_id,
                "method": f"{method}_scpa",
                "backend": "multicross::mcm",
                "ranking": "descending_signed_masking_delta",
            }
            result.update(_metric_or_null(gene_order, truth))
            metric_rows.append(result)
    if not metric_rows:
        raise ValueError("No SCPA masking checkpoints were evaluated")
    write_csv_atomic(metric_rows, list(metric_rows[0]), output_path)
    return metric_rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("llm", "scpa"), default="llm")
    parser.add_argument("--response-dir", type=Path)
    parser.add_argument("--request-dir", type=Path)
    parser.add_argument("--mapping-dir", type=Path)
    parser.add_argument("--scpa-checkpoint-dir", type=Path)
    parser.add_argument("--ground-truth", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-mock-toy", action="store_true")
    args = parser.parse_args()
    if args.mode == "llm":
        if any(value is None for value in (args.response_dir, args.request_dir, args.mapping_dir)):
            parser.error("LLM mode requires --response-dir, --request-dir and --mapping-dir")
        rows = evaluate_directory(
            args.response_dir, args.request_dir, args.mapping_dir,
            args.ground_truth, args.output, allow_mock_toy=args.allow_mock_toy,
        )
    else:
        if args.scpa_checkpoint_dir is None:
            parser.error("SCPA mode requires --scpa-checkpoint-dir")
        rows = evaluate_scpa_directory(
            args.scpa_checkpoint_dir, args.ground_truth, args.output,
        )
    print(
        f"PHASE7_EVALUATION status=PASS mode={args.mode} rows={len(rows)} "
        f"mock_toy={args.allow_mock_toy}"
    )
