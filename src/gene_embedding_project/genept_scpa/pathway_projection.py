"""Pathway-specific GenePT projection helpers for the Phase 4 paired analysis."""

from __future__ import annotations

from dataclasses import dataclass
import csv
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from scipy import sparse
from scipy.stats import kendalltau, spearmanr


@dataclass(frozen=True)
class PathwayDefinition:
    """One pathway and its ordered, unique gene symbols."""

    name: str
    source_database: str
    genes: tuple[str, ...]


@dataclass(frozen=True)
class PairedPathway:
    """Auditable gene-set counts and the frozen primary paired genes."""

    definition: PathwayDefinition
    shared_genes: tuple[str, ...]
    genept_mappable_genes: tuple[str, ...]
    paired_genes: tuple[str, ...]


def pathway_source(name: str) -> str:
    for source in ("HALLMARK", "KEGG", "REACTOME"):
        if name.startswith(f"{source}_"):
            return source
    return "OTHER"


def read_wide_pathway_csv(path: str | Path) -> list[PathwayDefinition]:
    """Read the headerless wide SCPA collection without changing gene order."""

    pathways: list[PathwayDefinition] = []
    names: set[str] = set()
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        for row_number, row in enumerate(csv.reader(handle), start=1):
            if not row or not row[0].strip():
                continue
            name = row[0].strip()
            if name in names:
                raise ValueError(f"Duplicate pathway name at row {row_number}: {name}")
            names.add(name)
            genes = tuple(dict.fromkeys(value.strip() for value in row[1:] if value.strip()))
            if not genes:
                raise ValueError(f"Pathway has no genes at row {row_number}: {name}")
            pathways.append(PathwayDefinition(name, pathway_source(name), genes))
    if not pathways:
        raise ValueError(f"No pathways found in {path}")
    return pathways


def build_paired_pathways(
    pathways: Sequence[PathwayDefinition],
    cd4_genes: Sequence[str],
    cd8_genes: Sequence[str],
    embedding_keys: set[str],
) -> list[PairedPathway]:
    """Freeze each primary gene set as pathway ∩ CD4 ∩ CD8 ∩ GenePT."""

    cd4 = set(cd4_genes)
    cd8 = set(cd8_genes)
    output: list[PairedPathway] = []
    for pathway in pathways:
        shared = tuple(sorted(gene for gene in pathway.genes if gene in cd4 and gene in cd8))
        mappable = tuple(sorted(gene for gene in pathway.genes if gene in embedding_keys))
        paired = tuple(gene for gene in shared if gene in embedding_keys)
        output.append(PairedPathway(pathway, shared, mappable, paired))
    return output


def filter_eligible_pathways(
    pathways: Sequence[PairedPathway], *, min_genes: int, max_genes: int
) -> list[PairedPathway]:
    if min_genes < 1 or max_genes < min_genes:
        raise ValueError("Invalid pathway gene thresholds")
    return [
        pathway
        for pathway in pathways
        if min_genes <= len(pathway.paired_genes) <= max_genes
    ]


def project_pathway(
    expression: np.ndarray | sparse.spmatrix,
    expression_genes: Sequence[str],
    embeddings: np.ndarray,
    embedding_keys: Sequence[str],
    *,
    l2_normalize: bool = False,
) -> np.ndarray:
    """Compute X_P @ E_P with strict row/column gene-order assertions."""

    expression_genes = tuple(expression_genes)
    embedding_keys = tuple(embedding_keys)
    if expression_genes != embedding_keys:
        raise ValueError("Expression columns and embedding rows have different gene order")
    if len(expression_genes) != len(set(expression_genes)):
        raise ValueError("Pathway gene order contains duplicate symbols")
    if expression.shape[1] != len(expression_genes):
        raise ValueError("Expression column count does not match pathway gene order")
    embedding_array = np.asarray(embeddings, dtype=np.float64)
    if embedding_array.ndim != 2 or embedding_array.shape[0] != len(embedding_keys):
        raise ValueError("Embedding rows do not match pathway gene order")
    if not np.isfinite(embedding_array).all():
        raise ValueError("Embedding matrix contains NaN or Inf")
    projected = np.asarray(expression @ embedding_array, dtype=np.float64)
    if not np.isfinite(projected).all():
        raise ValueError("Pathway projection contains NaN or Inf")
    if l2_normalize:
        norms = np.linalg.norm(projected, axis=1, keepdims=True)
        nonzero = norms[:, 0] > 0
        projected[nonzero] /= norms[nonzero]
    return projected


def expression_mass_coverage(
    full_log_expression: sparse.spmatrix,
    retained_indices: np.ndarray,
    reference_indices: np.ndarray,
) -> float:
    """Fraction of pathway log-expression mass retained by the paired gene set."""

    numerator = float(full_log_expression[:, retained_indices].sum())
    denominator = float(full_log_expression[:, reference_indices].sum())
    if denominator == 0:
        return 1.0 if numerator == 0 else 0.0
    return numerator / denominator


def scpa_bonferroni_qvalues(raw_p: Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
    """Reproduce SCPA 1.6.2: Bonferroni then sqrt(-log10(adjPval))."""

    values = np.asarray(raw_p, dtype=np.float64)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("raw_p must be a non-empty vector")
    if np.any(~np.isfinite(values)) or np.any(values < 0) or np.any(values > 1):
        raise ValueError("raw p-values must be finite and within [0, 1]")
    adjusted = np.minimum(values * values.size, 1.0)
    with np.errstate(divide="ignore"):
        qvalues = np.sqrt(-np.log10(adjusted))
    return adjusted, qvalues


def rank_descending(values: Sequence[float], names: Sequence[str]) -> np.ndarray:
    """Return deterministic 1-based ranks; larger values rank first."""

    numeric = np.asarray(values, dtype=np.float64)
    labels = np.asarray(names, dtype=str)
    if numeric.shape != labels.shape or np.any(~np.isfinite(numeric)):
        raise ValueError("Rank inputs must be aligned and finite")
    order = np.lexsort((labels, -numeric))
    ranks = np.empty(order.size, dtype=np.int64)
    ranks[order] = np.arange(1, order.size + 1)
    return ranks


def average_rank_descending(values: Sequence[float]) -> np.ndarray:
    """Average 1-based descending ranks, matching R ties.method='average'."""

    numeric = np.asarray(values, dtype=np.float64)
    if numeric.ndim != 1 or numeric.size == 0 or np.any(~np.isfinite(numeric)):
        raise ValueError("Rank values must be a non-empty finite vector")
    order = np.argsort(-numeric, kind="stable")
    ranks = np.empty(numeric.size, dtype=np.float64)
    start = 0
    while start < order.size:
        end = start + 1
        while end < order.size and numeric[order[end]] == numeric[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + 1 + end) / 2.0
        start = end
    return ranks


def significance_state(vanilla_adjusted_p: float, genept_adjusted_p: float) -> str:
    """Classify paired primary branches at the frozen adjusted-p threshold."""

    vanilla = float(vanilla_adjusted_p) < 0.05
    genept = float(genept_adjusted_p) < 0.05
    if vanilla and genept:
        return "Both significant"
    if vanilla:
        return "Vanilla-only significant"
    if genept:
        return "GenePT-only significant"
    return "Neither significant"


def ranking_agreement(
    vanilla_ranks: Sequence[int], genept_ranks: Sequence[int]
) -> dict[str, float | int]:
    """Agreement metrics only; these values are not accuracy measurements."""

    vanilla = np.asarray(vanilla_ranks, dtype=np.int64)
    genept = np.asarray(genept_ranks, dtype=np.int64)
    if vanilla.shape != genept.shape or vanilla.ndim != 1 or vanilla.size == 0:
        raise ValueError("Rank vectors must be non-empty and aligned")
    spearman = float(spearmanr(vanilla, genept).statistic)
    kendall = float(kendalltau(vanilla, genept).statistic)
    result: dict[str, float | int] = {"spearman": spearman, "kendall": kendall}
    for n in (10, 20):
        actual_n = min(n, vanilla.size)
        vanilla_top = set(np.flatnonzero(vanilla <= actual_n).tolist())
        genept_top = set(np.flatnonzero(genept <= actual_n).tolist())
        overlap = len(vanilla_top & genept_top)
        union = len(vanilla_top | genept_top)
        result[f"top{n}_n"] = actual_n
        result[f"top{n}_overlap"] = overlap
        result[f"top{n}_jaccard"] = float(overlap / union) if union else 1.0
    return result
