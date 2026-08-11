"""Published GenePT-w sparse projection with explicit numerical QC."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy import sparse


@dataclass
class ProjectionDiagnostics:
    """Per-cell values required by the Phase 2 QC gate."""

    library_sizes: np.ndarray
    matched_raw_mass: np.ndarray
    expression_coverage: np.ndarray
    pre_l2_norms: np.ndarray
    post_l2_norms: np.ndarray
    finite_values: int
    zero_vectors: int


def _validate_projection_inputs(
    counts: sparse.spmatrix,
    matched_gene_indices: np.ndarray,
    embedding_matrix: np.ndarray,
) -> sparse.csr_matrix:
    if not sparse.issparse(counts):
        raise TypeError("GenePT-w counts input must remain a scipy sparse matrix")
    matrix = counts.tocsr()
    if matrix.ndim != 2 or matrix.shape[0] == 0 or matrix.shape[1] == 0:
        raise ValueError("Counts must be a non-empty cells x genes matrix")
    if matrix.data.size and (
        not np.isfinite(matrix.data).all() or np.any(matrix.data < 0)
    ):
        raise ValueError("Counts must be finite and non-negative")
    indices = np.asarray(matched_gene_indices, dtype=np.int64)
    if indices.ndim != 1 or indices.size == 0:
        raise ValueError("matched_gene_indices must be a non-empty vector")
    if np.any(indices < 0) or np.any(indices >= matrix.shape[1]):
        raise ValueError("matched_gene_indices contains an out-of-range column")
    embeddings = np.asarray(embedding_matrix)
    if embeddings.ndim != 2 or embeddings.shape[0] != indices.size:
        raise ValueError("Embedding rows must align one-to-one with matched genes")
    if not np.isfinite(embeddings).all():
        raise ValueError("Embedding matrix contains NaN or Inf")
    return matrix


def _normalize_log1p_block(
    block: sparse.csr_matrix,
    library_sizes: np.ndarray,
    normalization_target: float,
) -> sparse.csr_matrix:
    normalized = block.astype(np.float64, copy=True)
    scale = np.zeros_like(library_sizes, dtype=np.float64)
    nonzero = library_sizes > 0
    scale[nonzero] = normalization_target / library_sizes[nonzero]
    normalized = sparse.diags(scale, format="csr") @ normalized
    np.log1p(normalized.data, out=normalized.data)
    return normalized.tocsr()


def normalize_log1p_sparse(
    counts: sparse.spmatrix,
    *,
    normalization_target: float = 10_000.0,
) -> sparse.csr_matrix:
    """Normalize each sparse cell over all input genes, then apply log1p."""

    if not sparse.issparse(counts):
        raise TypeError("Counts input must remain a scipy sparse matrix")
    matrix = counts.tocsr()
    if matrix.ndim != 2 or matrix.shape[0] == 0 or matrix.shape[1] == 0:
        raise ValueError("Counts must be a non-empty cells x genes matrix")
    if matrix.data.size and (
        not np.isfinite(matrix.data).all() or np.any(matrix.data < 0)
    ):
        raise ValueError("Counts must be finite and non-negative")
    if normalization_target <= 0:
        raise ValueError("normalization_target must be positive")
    library_sizes = np.asarray(matrix.sum(axis=1)).reshape(-1).astype(np.float64)
    return _normalize_log1p_block(matrix, library_sizes, normalization_target)


def project_genept_w_sparse(
    counts: sparse.spmatrix,
    matched_gene_indices: np.ndarray,
    embedding_matrix: np.ndarray,
    *,
    normalization_target: float = 10_000.0,
    batch_size: int = 512,
    output_path: str | Path | None = None,
) -> tuple[np.ndarray, ProjectionDiagnostics]:
    """Compute GenePT-w without densifying the cell-by-gene expression matrix.

    Order follows the published method: normalize raw counts over *all* dataset
    genes, apply log1p, align to official GenePT keys (unmatched genes therefore
    have zero projection contribution), take the expression-weighted average,
    and finally L2-normalize each non-zero cell vector. The official notebook's
    division by total dataset gene count is retained; it is a scalar cancelled by
    final L2 normalization but is reflected in the reported pre-L2 norms.
    """

    matrix = _validate_projection_inputs(
        counts, matched_gene_indices, embedding_matrix
    )
    if normalization_target <= 0:
        raise ValueError("normalization_target must be positive")
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")

    matched_gene_indices = np.asarray(matched_gene_indices, dtype=np.int64)
    embedding_matrix = np.asarray(embedding_matrix, dtype=np.float64)
    cells, dataset_gene_count = matrix.shape
    dimension = embedding_matrix.shape[1]
    library_sizes = np.asarray(matrix.sum(axis=1)).reshape(-1).astype(np.float64)
    matched_raw_mass = np.asarray(
        matrix[:, matched_gene_indices].sum(axis=1)
    ).reshape(-1).astype(np.float64)
    coverage = np.divide(
        matched_raw_mass,
        library_sizes,
        out=np.zeros_like(matched_raw_mass),
        where=library_sizes > 0,
    )

    if output_path is None:
        output: np.ndarray = np.empty((cells, dimension), dtype=np.float32)
    else:
        output = np.lib.format.open_memmap(
            Path(output_path), mode="w+", dtype=np.float32, shape=(cells, dimension)
        )
    pre_norms = np.zeros(cells, dtype=np.float64)
    post_norms = np.zeros(cells, dtype=np.float64)
    finite_values = 0

    for start in range(0, cells, batch_size):
        stop = min(start + batch_size, cells)
        normalized = _normalize_log1p_block(
            matrix[start:stop],
            library_sizes[start:stop],
            normalization_target,
        )
        projected = (
            normalized[:, matched_gene_indices] @ embedding_matrix
        ) / float(dataset_gene_count)
        projected = np.asarray(projected, dtype=np.float64)
        norms = np.linalg.norm(projected, axis=1)
        nonzero = norms > 0
        projected[nonzero] /= norms[nonzero, None]
        pre_norms[start:stop] = norms
        post_norms[start:stop] = np.linalg.norm(projected, axis=1)
        finite_values += int(np.isfinite(projected).sum())
        output[start:stop] = projected.astype(np.float32)

    if isinstance(output, np.memmap):
        output.flush()
    diagnostics = ProjectionDiagnostics(
        library_sizes=library_sizes,
        matched_raw_mass=matched_raw_mass,
        expression_coverage=coverage,
        pre_l2_norms=pre_norms,
        post_l2_norms=post_norms,
        finite_values=finite_values,
        zero_vectors=int(np.sum(pre_norms == 0)),
    )
    return output, diagnostics


def project_genept_w_direct(
    counts: sparse.spmatrix,
    matched_gene_indices: Sequence[int],
    embedding_matrix: np.ndarray,
    *,
    normalization_target: float = 10_000.0,
) -> np.ndarray:
    """Small explicit-loop reference implementation for correctness checks."""

    matrix = _validate_projection_inputs(
        counts, np.asarray(matched_gene_indices), embedding_matrix
    )
    matched_lookup = {
        int(dataset_index): np.asarray(embedding_matrix[row], dtype=np.float64)
        for row, dataset_index in enumerate(matched_gene_indices)
    }
    output = np.zeros((matrix.shape[0], embedding_matrix.shape[1]), dtype=np.float64)
    for cell_index in range(matrix.shape[0]):
        row = matrix.getrow(cell_index)
        library_size = float(row.data.sum())
        if library_size == 0:
            continue
        scale = normalization_target / library_size
        for gene_index, raw_count in zip(row.indices, row.data):
            vector = matched_lookup.get(int(gene_index))
            if vector is not None:
                weight = np.log1p(float(raw_count) * scale)
                output[cell_index] += weight * vector
        output[cell_index] /= float(matrix.shape[1])
        norm = np.linalg.norm(output[cell_index])
        if norm > 0:
            output[cell_index] /= norm
    return output


def numeric_summary(values: np.ndarray) -> dict[str, float]:
    """Return the frozen six-number summary used in Phase 2 QC."""

    array = np.asarray(values, dtype=np.float64)
    return {
        "min": float(np.min(array)),
        "q1": float(np.quantile(array, 0.25)),
        "median": float(np.median(array)),
        "mean": float(np.mean(array)),
        "q3": float(np.quantile(array, 0.75)),
        "max": float(np.max(array)),
    }
