"""Masking algebra and repeated strict-ranking aggregation for Phase 7."""

from __future__ import annotations

import hashlib
import math
from typing import Any, Callable, Mapping, Sequence

import numpy as np

from .schemas import validate_llm_response


def score_from_raw_p(raw_p: float, *, clip: float = 1e-300) -> float:
    value = float(raw_p)
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError("Raw p-value must be finite and in [0, 1]")
    if not 0 < clip <= 1:
        raise ValueError("Raw p-value clip must be in (0, 1]")
    return -math.log10(max(value, clip))


def vanilla_zero_mask(expression: np.ndarray, gene_index: int) -> np.ndarray:
    matrix = np.asarray(expression, dtype=np.float64)
    if matrix.ndim != 2 or not 0 <= gene_index < matrix.shape[1]:
        raise ValueError("Invalid Vanilla mask input")
    masked = matrix.copy()
    masked[:, gene_index] = 0.0
    return masked


def genept_projection(expression: np.ndarray, embeddings: np.ndarray) -> np.ndarray:
    matrix = np.asarray(expression, dtype=np.float64)
    embedding = np.asarray(embeddings, dtype=np.float64)
    if matrix.ndim != 2 or embedding.ndim != 2 or matrix.shape[1] != embedding.shape[0]:
        raise ValueError("Expression columns must align to embedding rows")
    result = matrix @ embedding
    if np.any(~np.isfinite(result)):
        raise ValueError("Non-finite GenePT projection")
    return result


def genept_subtraction_mask(
    full_projection: np.ndarray,
    expression_gene: np.ndarray,
    embedding_gene: np.ndarray,
) -> np.ndarray:
    full = np.asarray(full_projection, dtype=np.float64)
    values = np.asarray(expression_gene, dtype=np.float64).reshape(-1)
    vector = np.asarray(embedding_gene, dtype=np.float64).reshape(-1)
    if full.shape != (values.size, vector.size):
        raise ValueError("GenePT subtraction axes are not aligned")
    return full - np.outer(values, vector)


def _average_descending_ranks(values: Sequence[float]) -> list[float]:
    """Match R rank(-x, ties.method='average') for finite exact ties."""

    numeric = [float(value) for value in values]
    if any(not math.isfinite(value) for value in numeric):
        raise ValueError("Ranking values must be finite")
    order = sorted(range(len(numeric)), key=lambda index: -numeric[index])
    ranks = [0.0] * len(numeric)
    cursor = 0
    while cursor < len(order):
        end = cursor + 1
        while end < len(order) and numeric[order[end]] == numeric[order[cursor]]:
            end += 1
        average_rank = ((cursor + 1) + end) / 2.0
        for position in range(cursor, end):
            ranks[order[position]] = average_rank
        cursor = end
    return ranks


def compute_masking_rows(
    condition_a: np.ndarray,
    condition_b: np.ndarray,
    genes: Sequence[str],
    embeddings: np.ndarray,
    raw_p_function: Callable[[np.ndarray, np.ndarray], float],
    *,
    raw_p_clip: float = 1e-300,
) -> list[dict[str, Any]]:
    """Apply the Phase 5 definitions using an injected raw-p implementation."""

    a = np.asarray(condition_a, dtype=np.float64)
    b = np.asarray(condition_b, dtype=np.float64)
    ep = np.asarray(embeddings, dtype=np.float64)
    gene_names = tuple(str(gene) for gene in genes)
    if a.shape[1] != len(gene_names) or b.shape[1] != len(gene_names) or ep.shape[0] != len(gene_names):
        raise ValueError("Masking gene axes are not aligned")
    za = genept_projection(a, ep)
    zb = genept_projection(b, ep)
    vanilla_full_p = float(raw_p_function(a, b))
    genept_full_p = float(raw_p_function(za, zb))
    vanilla_full_score = score_from_raw_p(vanilla_full_p, clip=raw_p_clip)
    genept_full_score = score_from_raw_p(genept_full_p, clip=raw_p_clip)
    rows: list[dict[str, Any]] = []
    for index, gene in enumerate(gene_names):
        vanilla_p = float(raw_p_function(
            vanilla_zero_mask(a, index), vanilla_zero_mask(b, index)
        ))
        genept_p = float(raw_p_function(
            genept_subtraction_mask(za, a[:, index], ep[index]),
            genept_subtraction_mask(zb, b[:, index], ep[index]),
        ))
        rows.append({
            "gene": gene,
            "gene_index": index,
            "vanilla_raw_p_full": vanilla_full_p,
            "vanilla_raw_p_masked": vanilla_p,
            "vanilla_delta_score": vanilla_full_score - score_from_raw_p(vanilla_p, clip=raw_p_clip),
            "genept_raw_p_full": genept_full_p,
            "genept_raw_p_masked": genept_p,
            "genept_delta_score": genept_full_score - score_from_raw_p(genept_p, clip=raw_p_clip),
        })
    for method in ("vanilla", "genept"):
        signed = _average_descending_ranks(
            [row[f"{method}_delta_score"] for row in rows]
        )
        absolute = _average_descending_ranks(
            [abs(row[f"{method}_delta_score"]) for row in rows]
        )
        for index, row in enumerate(rows):
            row[f"{method}_signed_rank"] = signed[index]
            row[f"{method}_absolute_rank"] = absolute[index]
    return rows


def _tie_hash(seed: int, candidate_id: str) -> str:
    return hashlib.sha256(f"{seed}:{candidate_id}".encode("utf-8")).hexdigest()


def aggregate_rankings(
    responses: Sequence[Mapping[str, Any]], *, tie_seed: int
) -> list[dict[str, Any]]:
    if not responses:
        raise ValueError("At least one response is required")
    for response in responses:
        validate_llm_response(response)
    first_ids = {row["candidate_id"] for row in responses[0]["ranking"]}
    if any({row["candidate_id"] for row in response["ranking"]} != first_ids for response in responses):
        raise ValueError("Repeated responses have different candidate sets")
    ranks = {
        candidate: [
            next(row["rank"] for row in response["ranking"] if row["candidate_id"] == candidate)
            for response in responses
        ]
        for candidate in first_ids
    }
    ordered = sorted(
        first_ids,
        key=lambda candidate: (float(np.mean(ranks[candidate])), _tie_hash(tie_seed, candidate)),
    )
    return [
        {
            "candidate_id": candidate,
            "mean_rank": float(np.mean(ranks[candidate])),
            "aggregate_rank": index + 1,
            "individual_ranks": list(ranks[candidate]),
        }
        for index, candidate in enumerate(ordered)
    ]
