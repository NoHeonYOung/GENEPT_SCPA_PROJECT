"""Direct NumPy ranking metrics frozen for Phase 7."""

from __future__ import annotations

import itertools
import math
from typing import Any, Mapping, Sequence

import numpy as np


def _validate_ranking(ranking: Sequence[str]) -> tuple[str, ...]:
    ordered = tuple(str(item) for item in ranking)
    if not ordered or len(ordered) != len(set(ordered)):
        raise ValueError("Ranking must be non-empty and unique")
    return ordered


def recall_at_k(ranking: Sequence[str], truth: set[str], k: int) -> float:
    ordered = _validate_ranking(ranking)
    if not truth or not 1 <= k <= len(ordered) or not truth <= set(ordered):
        raise ValueError("Recall@K requires non-empty aligned truth and valid K")
    return float(len(set(ordered[:k]) & truth) / len(truth))


def average_precision(ranking: Sequence[str], truth: set[str]) -> float:
    ordered = _validate_ranking(ranking)
    if not truth or not truth <= set(ordered):
        raise ValueError("Average Precision requires non-empty aligned truth")
    hits = 0
    total = 0.0
    for rank, item in enumerate(ordered, start=1):
        if item in truth:
            hits += 1
            total += hits / rank
    return float(total / len(truth))


def ndcg_at_k(ranking: Sequence[str], truth: set[str], k: int) -> float:
    ordered = _validate_ranking(ranking)
    if not truth or not 1 <= k <= len(ordered) or not truth <= set(ordered):
        raise ValueError("NDCG requires non-empty aligned truth and valid K")
    dcg = sum(
        (1.0 if item in truth else 0.0) / math.log2(rank + 1)
        for rank, item in enumerate(ordered[:k], start=1)
    )
    ideal_hits = min(len(truth), k)
    idcg = sum(1.0 / math.log2(rank + 1) for rank in range(1, ideal_hits + 1))
    return float(dcg / idcg)


def evaluate_ranking(ranking: Sequence[str], truth: set[str]) -> dict[str, float | int]:
    ordered = _validate_ranking(ranking)
    if not truth:
        raise ValueError("Null experiments do not have ground-truth recovery metrics")
    k = len(truth)
    return {
        "truth_k": k,
        "recall_at_truth_k": recall_at_k(ordered, truth, k),
        "average_precision": average_precision(ordered, truth),
        "ndcg_at_n": ndcg_at_k(ordered, truth, len(ordered)),
        "ndcg_at_truth_k": ndcg_at_k(ordered, truth, k),
    }


def _spearman_strict(a: Mapping[str, int], b: Mapping[str, int]) -> float:
    if set(a) != set(b) or len(a) < 2:
        raise ValueError("Spearman rankings must align and contain at least two candidates")
    names = sorted(a)
    values_a = np.asarray([a[name] for name in names], dtype=np.float64)
    values_b = np.asarray([b[name] for name in names], dtype=np.float64)
    return float(np.corrcoef(values_a, values_b)[0, 1])


def prompt_order_spearman(responses: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    if len(responses) < 2:
        raise ValueError("Prompt-order stability requires at least two runs")
    rank_maps = [
        {row["candidate_id"]: int(row["rank"]) for row in response["ranking"]}
        for response in responses
    ]
    correlations = [
        _spearman_strict(rank_maps[left], rank_maps[right])
        for left, right in itertools.combinations(range(len(rank_maps)), 2)
    ]
    return {
        "pair_count": len(correlations),
        "pairwise_spearman": correlations,
        "mean_spearman": float(np.mean(correlations)),
        "min_spearman": float(np.min(correlations)),
    }
