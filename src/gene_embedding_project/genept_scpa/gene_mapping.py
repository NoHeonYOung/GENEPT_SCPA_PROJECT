"""Strict GenePT artifact loading and gene-symbol matching helpers."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import pickle
from typing import Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class GeneMatch:
    """One dataset feature's match to the official GenePT lookup."""

    dataset_index: int
    dataset_gene: str
    match_type: str
    embedding_key: str | None


def load_official_genept_embeddings(
    path: str | Path,
    *,
    expected_dimension: int = 1536,
) -> dict[str, np.ndarray]:
    """Load the checksum-validated official pickle and validate its schema.

    Pickle is intentionally accepted only for the pinned official Zenodo artifact.
    The acquisition script verifies the archive checksum before extracting it.
    """

    with Path(path).open("rb") as handle:
        payload = pickle.load(handle)  # noqa: S301 - pinned official artifact format
    if not isinstance(payload, Mapping) or not payload:
        raise ValueError("GenePT embedding artifact must be a non-empty dictionary")

    embeddings: dict[str, np.ndarray] = {}
    for key, value in payload.items():
        if not isinstance(key, str) or not key:
            raise ValueError("GenePT embedding keys must be non-empty strings")
        vector = np.asarray(value, dtype=np.float32).reshape(-1)
        if vector.shape != (expected_dimension,):
            raise ValueError(
                f"Embedding {key!r} has shape {vector.shape}; "
                f"expected ({expected_dimension},)"
            )
        if not np.isfinite(vector).all():
            raise ValueError(f"Embedding {key!r} contains NaN or Inf")
        embeddings[key] = vector
    return embeddings


def load_primary_gene_keys(path: str | Path) -> set[str]:
    """Return primary NCBI gene keys used to distinguish official alias keys."""

    with Path(path).open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, Mapping):
        raise ValueError("NCBI_summary_of_genes.json must contain a JSON object")
    keys = set(payload)
    if not keys or any(not isinstance(key, str) or not key for key in keys):
        raise ValueError("Primary NCBI gene keys must be non-empty strings")
    return keys


def classify_gene_matches(
    dataset_genes: Sequence[str],
    embedding_keys: set[str],
    primary_gene_keys: set[str],
) -> list[GeneMatch]:
    """Match exact artifact keys only; never case-fold or fuzzy-match symbols.

    The official GenePT methods add HGNC aliases to the lookup. An exact artifact
    key absent from the primary NCBI-summary keys is therefore recorded as an
    ``official_alias`` match rather than silently treated as a primary symbol.
    """

    matches: list[GeneMatch] = []
    for index, gene in enumerate(dataset_genes):
        if not isinstance(gene, str) or not gene:
            raise ValueError(f"Dataset gene at index {index} is empty or not a string")
        if gene not in embedding_keys:
            matches.append(GeneMatch(index, gene, "unmatched", None))
        elif gene in primary_gene_keys:
            matches.append(GeneMatch(index, gene, "exact", gene))
        else:
            matches.append(GeneMatch(index, gene, "official_alias", gene))
    return matches


def build_aligned_embedding_matrix(
    matches: Sequence[GeneMatch],
    embeddings: Mapping[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    """Return matched dataset-column indices and the aligned embedding matrix."""

    matched = [match for match in matches if match.embedding_key is not None]
    if not matched:
        raise ValueError("No dataset genes match the official GenePT artifact")
    indices = np.asarray([match.dataset_index for match in matched], dtype=np.int64)
    matrix = np.stack(
        [embeddings[match.embedding_key] for match in matched], axis=0
    ).astype(np.float32, copy=False)
    return indices, matrix


def mapping_counts(matches: Sequence[GeneMatch]) -> dict[str, int]:
    """Summarize exact, official-alias, unmatched, and duplicate mappings."""

    dataset_genes = [match.dataset_gene for match in matches]
    duplicate_dataset_genes = len(dataset_genes) - len(set(dataset_genes))
    return {
        "dataset_genes": len(matches),
        "exact_matches": sum(match.match_type == "exact" for match in matches),
        "alias_matches": sum(
            match.match_type == "official_alias" for match in matches
        ),
        "unmatched_dataset_genes": sum(
            match.match_type == "unmatched" for match in matches
        ),
        "dataset_duplicate_genes": duplicate_dataset_genes,
        "embedding_duplicate_keys": 0,
        "duplicate_mapping_count": duplicate_dataset_genes,
    }
