"""Frozen perturbations in normalized log1p space.

The null scenario samples an uninjected evaluation-target set. This is not a
claim that those genes were perturbed; it makes the requested chance-level
negative-control metrics mathematically defined.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Sequence

import numpy as np


SCENARIOS = ("null", "mean_shift", "cell_subset", "mixed_direction")


@dataclass(frozen=True)
class PerturbationResult:
    condition_a: np.ndarray
    condition_b: np.ndarray
    truth_rows: tuple[dict[str, Any], ...]
    target_cell_indices: tuple[int, ...]


def ground_truth_gene_count(pathway_gene_count: int) -> int:
    if pathway_gene_count < 1:
        raise ValueError("Pathway must contain genes")
    return max(3, min(10, math.ceil(0.15 * pathway_gene_count)))


def pooled_gene_statistics(a: np.ndarray, b: np.ndarray) -> dict[str, np.ndarray]:
    pooled = np.vstack([a, b])
    return {
        "sd": np.std(pooled, axis=0, ddof=1),
        "median": np.median(pooled, axis=0),
        "detection_fraction": np.mean(pooled > 0, axis=0),
        "condition_b_detection_fraction": np.mean(b > 0, axis=0),
    }


def _validate(a: np.ndarray, b: np.ndarray, genes: Sequence[str]):
    xa = np.asarray(a, dtype=np.float64)
    xb = np.asarray(b, dtype=np.float64)
    names = tuple(str(gene) for gene in genes)
    if xa.ndim != 2 or xb.ndim != 2 or xa.shape[1] != xb.shape[1]:
        raise ValueError("Condition matrices must align as cells x genes")
    if xa.shape[1] != len(names) or len(names) != len(set(names)):
        raise ValueError("Gene axis must be unique and aligned")
    if xa.shape[0] < 2 or xb.shape[0] < 2:
        raise ValueError("Each condition needs at least two cells")
    if np.any(~np.isfinite(xa)) or np.any(~np.isfinite(xb)) or np.any(xa < 0) or np.any(xb < 0):
        raise ValueError("Expression must be finite and non-negative")
    return xa, xb, names


def _select_targets(
    statistics: dict[str, np.ndarray], scenario: str, count: int, seed: int,
    detection_fraction_min: float, negative_detection_fraction_min: float,
    negative_median_min: float,
) -> tuple[np.ndarray, dict[int, str], dict[int, str]]:
    eligible = np.flatnonzero(
        (statistics["detection_fraction"] >= detection_fraction_min)
        & (statistics["sd"] > 0)
    )
    rng = np.random.default_rng(seed)
    if scenario != "mixed_direction":
        if len(eligible) >= count:
            selected_pool = eligible
            rule = "general_truth_eligible"
        else:
            selected_pool = np.arange(len(statistics["sd"]), dtype=np.int64)
            rule = "all_pathway_genes_fallback"
        chosen = np.asarray(rng.choice(selected_pool, count, replace=False), dtype=np.int64)
        return (
            chosen,
            {int(index): "positive" for index in chosen},
            {int(index): rule for index in chosen},
        )

    negative_count = count // 2
    negative_pool = np.flatnonzero(
        (statistics["detection_fraction"] >= negative_detection_fraction_min)
        & (statistics["median"] > negative_median_min)
        & (statistics["sd"] > 0)
        & (statistics["condition_b_detection_fraction"] > 0)
    )
    if len(negative_pool) >= negative_count:
        selected_negative_pool = negative_pool
        negative_rule = "strict_negative_eligible"
    else:
        # The frozen 11-pathway universe contains pathways with too few genes
        # meeting the preferred median/detection rule.  Before any MCM result
        # is observed, fall back to the general nonzero-variance truth pool;
        # clipping is retained and logged per gene.
        selected_negative_pool = np.flatnonzero(
            statistics["condition_b_detection_fraction"] > 0
        )
        negative_count = min(negative_count, len(selected_negative_pool))
        negative_rule = "condition_b_detected_fallback"
    if negative_count < 1:
        raise ValueError("Mixed direction requires at least one gene detected in condition B")
    negative = np.asarray(
        rng.choice(selected_negative_pool, negative_count, replace=False), dtype=np.int64
    )
    negative_set = set(int(value) for value in negative)
    positive_count = count - negative_count
    positive_pool = np.asarray(
        [value for value in eligible if int(value) not in negative_set], dtype=np.int64
    )
    positive_rule = "general_truth_eligible"
    if len(positive_pool) < positive_count:
        positive_pool = np.asarray(
            [value for value in range(len(statistics["sd"])) if value not in negative_set],
            dtype=np.int64,
        )
        positive_rule = "all_pathway_genes_fallback"
    positive = np.asarray(rng.choice(positive_pool, positive_count, replace=False), dtype=np.int64)
    chosen = np.concatenate([positive, negative])
    directions = {int(index): "positive" for index in positive}
    directions.update({int(index): "negative" for index in negative})
    selection_rules = {int(index): positive_rule for index in positive}
    selection_rules.update({int(index): negative_rule for index in negative})
    return chosen, directions, selection_rules


def inject_perturbation(
    condition_a: np.ndarray, condition_b: np.ndarray, genes: Sequence[str], *,
    scenario: str, alpha: float, seed: int, cell_subset_fraction: float = 0.30,
    detection_fraction_min: float = 0.10,
    negative_detection_fraction_min: float = 0.50,
    negative_median_min: float = 0.0, scale_floor: float = 0.10,
) -> PerturbationResult:
    """Return copied matrices; only B is changed and it is never renormalized."""

    a, b, names = _validate(condition_a, condition_b, genes)
    if scenario not in SCENARIOS:
        raise ValueError(f"Unknown scenario: {scenario}")
    if scenario == "null" and alpha != 0:
        raise ValueError("Null scenario requires alpha=0")
    if scenario != "null" and alpha <= 0:
        raise ValueError("Non-null scenarios require alpha>0")
    if scale_floor <= 0 or not 0 < cell_subset_fraction <= 1:
        raise ValueError("Invalid scale floor or cell-subset fraction")

    statistics = pooled_gene_statistics(a, b)
    target_scenario = "mean_shift" if scenario == "null" else scenario
    targets, directions, selection_rules = _select_targets(
        statistics, target_scenario, ground_truth_gene_count(len(names)), seed,
        detection_fraction_min, negative_detection_fraction_min, negative_median_min,
    )
    rng = np.random.default_rng(seed + 1)
    if scenario == "cell_subset":
        cell_count = int(round(cell_subset_fraction * b.shape[0]))
        cells = np.sort(rng.choice(b.shape[0], cell_count, replace=False))
    elif scenario == "null":
        cells = np.asarray([], dtype=np.int64)
    else:
        cells = np.arange(b.shape[0], dtype=np.int64)

    output = b.copy()
    rows: list[dict[str, Any]] = []
    for raw_index in targets:
        index = int(raw_index)
        scale = max(float(statistics["sd"][index]), float(scale_floor))
        signed_delta = float(alpha * scale)
        direction = "none" if scenario == "null" else directions[index]
        clipped = 0
        applied = 0.0
        if scenario != "null":
            before = output[cells, index].copy()
            if direction == "negative":
                unclipped = before - signed_delta
                output[cells, index] = np.maximum(unclipped, 0.0)
                clipped = int(np.sum(unclipped < 0))
                applied = -signed_delta
            else:
                output[cells, index] = before + signed_delta
                applied = signed_delta
        rows.append({
            "gene": names[index], "gene_index": index,
            "is_evaluation_target": True,
            "is_ground_truth_perturbed": scenario != "null",
            "perturbation_direction": direction,
            "direction_selection_rule": selection_rules[index],
            "perturbation_strength": float(alpha),
            "target_cell_fraction": float(len(cells) / b.shape[0]),
            "target_cell_count": int(len(cells)),
            "pooled_baseline_sd": float(statistics["sd"][index]),
            "applied_log_delta": applied, "clipped_cell_count": clipped,
        })
    rows.sort(key=lambda row: row["gene_index"])
    return PerturbationResult(a.copy(), output, tuple(rows), tuple(int(i) for i in cells))
