"""Controlled Phase 7 perturbations in normalized log1p expression space."""

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
    ground_truth_rows: tuple[dict[str, Any], ...]
    perturbed_cell_indices: tuple[int, ...]


def ground_truth_gene_count(pathway_gene_count: int) -> int:
    if pathway_gene_count < 1:
        raise ValueError("Pathway must contain genes")
    return max(3, min(10, math.ceil(0.15 * pathway_gene_count)))


def _validated_arrays(
    condition_a: np.ndarray, condition_b: np.ndarray, genes: Sequence[str]
) -> tuple[np.ndarray, np.ndarray, tuple[str, ...]]:
    a = np.asarray(condition_a, dtype=np.float64)
    b = np.asarray(condition_b, dtype=np.float64)
    gene_names = tuple(str(gene) for gene in genes)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1] != b.shape[1]:
        raise ValueError("Condition matrices must be aligned cells x genes")
    if a.shape[1] != len(gene_names) or len(gene_names) != len(set(gene_names)):
        raise ValueError("Gene axis is not unique and aligned")
    if a.shape[0] < 2 or b.shape[0] < 2:
        raise ValueError("Each condition needs at least two cells")
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)) or np.any(a < 0) or np.any(b < 0):
        raise ValueError("Normalized expression must be finite and non-negative")
    return a, b, gene_names


def pooled_gene_statistics(
    condition_a: np.ndarray, condition_b: np.ndarray
) -> dict[str, np.ndarray]:
    pooled = np.vstack([condition_a, condition_b])
    return {
        "sd": np.std(pooled, axis=0, ddof=1),
        "median": np.median(pooled, axis=0),
        "detection_fraction": np.mean(pooled > 0, axis=0),
    }


def _select_truth(
    a: np.ndarray,
    b: np.ndarray,
    scenario: str,
    *,
    seed: int,
    detection_fraction_min: float,
    negative_detection_fraction_min: float,
    negative_median_min: float,
) -> tuple[np.ndarray, dict[int, str], dict[str, np.ndarray]]:
    statistics = pooled_gene_statistics(a, b)
    count = ground_truth_gene_count(a.shape[1])
    general = np.flatnonzero(
        (statistics["detection_fraction"] >= detection_fraction_min)
        & (statistics["sd"] > 0)
    )
    if len(general) < count:
        raise ValueError(f"Only {len(general)} truth-eligible genes; need {count}")
    generator = np.random.default_rng(seed)
    directions: dict[int, str] = {}
    if scenario != "mixed_direction":
        selected = generator.choice(general, size=count, replace=False)
        directions = {int(index): "positive" for index in selected}
        return np.asarray(selected, dtype=np.int64), directions, statistics

    negative_count = count // 2
    negative_pool = np.flatnonzero(
        (statistics["detection_fraction"] >= negative_detection_fraction_min)
        & (statistics["median"] > negative_median_min)
        & (statistics["sd"] > 0)
    )
    if len(negative_pool) < negative_count:
        raise ValueError(
            f"Only {len(negative_pool)} negative-eligible genes; need {negative_count}"
        )
    negative = np.asarray(
        generator.choice(negative_pool, size=negative_count, replace=False), dtype=np.int64
    )
    remaining = np.asarray([index for index in general if index not in set(negative)], dtype=np.int64)
    positive_count = count - negative_count
    if len(remaining) < positive_count:
        raise ValueError("Insufficient positive genes after negative selection")
    positive = np.asarray(
        generator.choice(remaining, size=positive_count, replace=False), dtype=np.int64
    )
    selected = np.concatenate([positive, negative])
    directions.update({int(index): "positive" for index in positive})
    directions.update({int(index): "negative" for index in negative})
    return selected, directions, statistics


def inject_perturbation(
    condition_a: np.ndarray,
    condition_b: np.ndarray,
    genes: Sequence[str],
    *,
    scenario: str,
    alpha: float,
    seed: int,
    cell_subset_fraction: float = 0.30,
    detection_fraction_min: float = 0.10,
    negative_detection_fraction_min: float = 0.50,
    negative_median_min: float = 0.0,
    scale_floor: float = 0.10,
) -> PerturbationResult:
    """Inject a frozen scenario into B without mutating or renormalizing inputs."""

    a, b, gene_names = _validated_arrays(condition_a, condition_b, genes)
    if scenario not in SCENARIOS:
        raise ValueError(f"Unknown perturbation scenario: {scenario}")
    if scenario == "null":
        if alpha != 0:
            raise ValueError("Null perturbation must use alpha=0")
        return PerturbationResult(a.copy(), b.copy(), tuple(), tuple())
    if alpha <= 0 or scale_floor <= 0:
        raise ValueError("Non-null alpha and scale floor must be positive")
    if not 0 < cell_subset_fraction <= 1:
        raise ValueError("Cell subset fraction must be in (0, 1]")

    truth, directions, statistics = _select_truth(
        a, b, scenario, seed=seed,
        detection_fraction_min=detection_fraction_min,
        negative_detection_fraction_min=negative_detection_fraction_min,
        negative_median_min=negative_median_min,
    )
    output = b.copy()
    generator = np.random.default_rng(seed + 1)
    if scenario == "cell_subset":
        cell_count = int(round(cell_subset_fraction * b.shape[0]))
        cells = np.sort(generator.choice(b.shape[0], size=cell_count, replace=False))
    else:
        cells = np.arange(b.shape[0], dtype=np.int64)

    rows: list[dict[str, Any]] = []
    for index in truth:
        gene_index = int(index)
        scale = max(float(statistics["sd"][gene_index]), float(scale_floor))
        signed_delta = float(alpha * scale)
        before = output[cells, gene_index].copy()
        if directions[gene_index] == "negative":
            after_unclipped = before - signed_delta
            after = np.maximum(after_unclipped, 0.0)
            clipped = int(np.sum(after_unclipped < 0))
            applied = -signed_delta
        else:
            after = before + signed_delta
            clipped = 0
            applied = signed_delta
        output[cells, gene_index] = after
        rows.append({
            "gene": gene_names[gene_index],
            "gene_index": gene_index,
            "is_ground_truth_perturbed": True,
            "perturbation_direction": directions[gene_index],
            "perturbation_strength": float(alpha),
            "target_cell_fraction": float(len(cells) / b.shape[0]),
            "target_cell_count": int(len(cells)),
            "pooled_baseline_sd": float(statistics["sd"][gene_index]),
            "applied_log_delta": applied,
            "clipped_cell_count": clipped,
        })
    rows.sort(key=lambda row: row["gene_index"])
    return PerturbationResult(
        a.copy(), output, tuple(rows), tuple(int(index) for index in cells)
    )
