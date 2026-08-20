"""Hand-testable ranking metrics for synthetic ground-truth recovery.

Higher values mean better recovery in this synthetic setup only. They are not
evidence that a method is biologically superior or causal.
"""

from __future__ import annotations

import math
from typing import Collection, Optional, Sequence


def _validate(ranking: Sequence[str], truth: Collection[str], k: Optional[int] = None):
    ordered = tuple(str(item) for item in ranking)
    relevant = frozenset(str(item) for item in truth)
    if not ordered or len(ordered) != len(set(ordered)):
        raise ValueError("Ranking must be non-empty and unique")
    if not relevant or not relevant <= set(ordered):
        raise ValueError("Truth must be non-empty and contained in the ranking")
    if k is not None and not 1 <= k <= len(ordered):
        raise ValueError("K must be between 1 and ranking length")
    return ordered, relevant


def recall_at_k(ranking: Sequence[str], truth: Collection[str], k: int) -> float:
    ordered, relevant = _validate(ranking, truth, k)
    return len(set(ordered[:k]) & relevant) / len(relevant)


def average_precision(ranking: Sequence[str], truth: Collection[str]) -> float:
    ordered, relevant = _validate(ranking, truth)
    hits = 0
    precision_sum = 0.0
    for rank, gene in enumerate(ordered, start=1):
        if gene in relevant:
            hits += 1
            precision_sum += hits / rank
    return precision_sum / len(relevant)


def ndcg_at_k(ranking: Sequence[str], truth: Collection[str], k: int) -> float:
    ordered, relevant = _validate(ranking, truth, k)
    dcg = sum(
        (1.0 if gene in relevant else 0.0) / math.log2(rank + 1)
        for rank, gene in enumerate(ordered[:k], start=1)
    )
    ideal_hits = min(k, len(relevant))
    idcg = sum(1.0 / math.log2(rank + 1) for rank in range(1, ideal_hits + 1))
    return dcg / idcg


def exact_random_chance(n_genes: int, n_truth: int, k: int) -> dict[str, float]:
    """Exact expectations under a uniformly random complete ranking."""

    if not 1 <= n_truth <= n_genes or not 1 <= k <= n_genes:
        raise ValueError("Invalid universe, truth count, or K")
    harmonic = sum(1.0 / rank for rank in range(1, n_genes + 1))
    expected_ap = (
        harmonic + (n_truth - 1) / (n_genes - 1) * (n_genes - harmonic)
    ) / n_genes if n_genes > 1 else 1.0
    discounts = [1.0 / math.log2(rank + 1) for rank in range(1, k + 1)]
    ideal = sum(discounts[: min(k, n_truth)])
    return {
        "recall": min(k, n_genes) / n_genes,
        "average_precision": expected_ap,
        "ndcg": (n_truth / n_genes) * sum(discounts) / ideal,
    }
