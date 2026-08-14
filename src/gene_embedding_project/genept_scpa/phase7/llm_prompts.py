"""Leakage-controlled Phase 7 LLM request construction."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np

from .pathway_selection import sanitize_gene_description
from .schemas import LLM_INPUT_SCHEMA_VERSION, validate_llm_request


PROMPT_CONDITIONS = ("stats_only", "true_description", "shuffled_description")


@dataclass(frozen=True)
class LLMRequestBundle:
    request: dict[str, Any]
    candidate_to_gene: dict[str, str]
    description_source_candidate: dict[str, str]
    description_changed_fraction: float


def expression_summary(values: np.ndarray) -> dict[str, float]:
    array = np.asarray(values, dtype=np.float64).reshape(-1)
    if array.size == 0 or np.any(~np.isfinite(array)) or np.any(array < 0):
        raise ValueError("Expression summary requires finite non-negative values")
    return {
        "mean": float(np.mean(array)),
        "sd": float(np.std(array, ddof=1)) if array.size > 1 else 0.0,
        "median": float(np.median(array)),
        "q10": float(np.quantile(array, 0.10)),
        "q25": float(np.quantile(array, 0.25)),
        "q75": float(np.quantile(array, 0.75)),
        "q90": float(np.quantile(array, 0.90)),
        "nonzero_fraction": float(np.mean(array > 0)),
    }


def _deranged_source_indices(size: int, seed: int) -> np.ndarray:
    if size < 2:
        raise ValueError("Description derangement requires at least two candidates")
    generator = np.random.default_rng(seed)
    offset = int(generator.integers(1, size))
    return np.roll(np.arange(size, dtype=np.int64), offset)


def build_llm_request(
    *,
    experiment_id: str,
    run_id: str,
    pathway: str,
    source_database: str,
    genes: Sequence[str],
    condition_a: np.ndarray,
    condition_b: np.ndarray,
    descriptions: Mapping[str, str],
    prompt_condition: str,
    candidate_order_seed: int,
    description_shuffle_seed: int,
    backend: str,
) -> LLMRequestBundle:
    """Build an LLM-visible request with no gene symbol or ground-truth field."""

    if prompt_condition not in PROMPT_CONDITIONS:
        raise ValueError(f"Unknown prompt condition: {prompt_condition}")
    gene_names = tuple(str(gene) for gene in genes)
    a = np.asarray(condition_a, dtype=np.float64)
    b = np.asarray(condition_b, dtype=np.float64)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1] != len(gene_names) or b.shape[1] != len(gene_names):
        raise ValueError("LLM expression matrices and gene axis are not aligned")
    width = max(3, len(str(len(gene_names))))
    candidate_ids = tuple(f"C{index + 1:0{width}d}" for index in range(len(gene_names)))
    candidate_to_gene = dict(zip(candidate_ids, gene_names))
    sanitized = [sanitize_gene_description(gene, descriptions.get(gene, "")) for gene in gene_names]
    if prompt_condition != "stats_only" and not all(sanitized):
        raise ValueError("Every description condition candidate needs a usable description")

    if prompt_condition == "shuffled_description":
        source_indices = _deranged_source_indices(len(gene_names), description_shuffle_seed)
    else:
        source_indices = np.arange(len(gene_names), dtype=np.int64)
    source_mapping = {
        candidate_ids[index]: candidate_ids[int(source_indices[index])]
        for index in range(len(gene_names))
    }
    changed = float(np.mean(source_indices != np.arange(len(gene_names))))
    if prompt_condition == "shuffled_description" and changed <= 0.9:
        raise RuntimeError("Description shuffle changed fraction must exceed 0.9")

    candidates: list[dict[str, Any]] = []
    for index, candidate_id in enumerate(candidate_ids):
        record: dict[str, Any] = {
            "candidate_id": candidate_id,
            "condition_a": expression_summary(a[:, index]),
            "condition_b": expression_summary(b[:, index]),
        }
        if prompt_condition != "stats_only":
            record["description"] = sanitized[int(source_indices[index])]
        candidates.append(record)
    order = np.random.default_rng(candidate_order_seed).permutation(len(candidates))
    visible_candidates = [candidates[int(index)] for index in order]
    request = {
        "schema_version": LLM_INPUT_SCHEMA_VERSION,
        "run_id": str(run_id),
        "experiment_id": str(experiment_id),
        "backend": str(backend),
        "pathway": {
            "name": str(pathway),
            "source_database": str(source_database),
            "candidate_count": len(candidates),
        },
        "comparison": {
            "condition_a": "A", "condition_b": "B",
            "cells_a": int(a.shape[0]), "cells_b": int(b.shape[0]),
        },
        "prompt_condition": prompt_condition,
        "candidate_order_seed": int(candidate_order_seed),
        "candidates": visible_candidates,
        "output_instruction": (
            "Rank every opaque candidate from most to least likely to drive the "
            "difference between A and B. Return only the required JSON schema."
        ),
    }
    validate_llm_request(request)
    return LLMRequestBundle(request, candidate_to_gene, source_mapping, changed)
