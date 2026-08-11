#!/usr/bin/env python3
"""Build and QC published-workflow GenePT-w for a configured Seurat dataset."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable

import numpy as np
from scipy import io as scipy_io
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

DATASET_SETTINGS = {
    "naive_cd4": {
        "phase": 2,
        "phase_label": "Phase 2 - Published GenePT-w reproduction",
        "expected_cells": 14_894,
        "expected_features": 17_856,
        "output_directory": "data/processed/genept_scpa/phase2",
        "mapping_filename": "genept_gene_mapping.csv",
        "qc": "data/interim/genept_scpa/phase2_genept_w_qc.json",
        "summary": "data/interim/genept_scpa/phase2_genept_w_summary.md",
    },
    "naive_cd8": {
        "phase": 3,
        "phase_label": "Phase 3 preparation - naïve CD8 GenePT-w reproduction",
        "expected_cells": 7_041,
        "expected_features": 17_942,
        "output_directory": "data/processed/genept_scpa/phase3",
        "mapping_filename": "naive_cd8_genept_gene_mapping.csv",
        "qc": "data/interim/genept_scpa/phase3_naive_cd8_genept_w_qc.json",
        "summary": "data/interim/genept_scpa/phase3_naive_cd8_genept_w_summary.md",
    },
}

from gene_embedding_project.genept_scpa.config import load_config  # noqa: E402
from gene_embedding_project.genept_scpa.gene_mapping import (  # noqa: E402
    GeneMatch,
    build_aligned_embedding_matrix,
    classify_gene_matches,
    load_official_genept_embeddings,
    load_primary_gene_keys,
    mapping_counts,
)
from gene_embedding_project.genept_scpa.genept_projection import (  # noqa: E402
    numeric_summary,
    project_genept_w_direct,
    project_genept_w_sparse,
)
from gene_embedding_project.genept_scpa.io import (  # noqa: E402
    sha256_file,
    write_json_atomic,
)


def read_lines(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8") as handle:
        return [line.rstrip("\n\r") for line in handle]


def write_lines_atomic(lines: Iterable[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            for line in lines:
                handle.write(str(line))
                handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def copy_atomic(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as target, source.open("rb") as handle:
            temporary_name = target.name
            shutil.copyfileobj(handle, target)
            target.flush()
            os.fsync(target.fileno())
        os.replace(temporary_name, destination)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def write_mapping_atomic(matches: list[GeneMatch], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            writer = csv.writer(handle)
            writer.writerow(
                ["dataset_index", "dataset_gene", "match_type", "embedding_key"]
            )
            for match in matches:
                writer.writerow(
                    [
                        match.dataset_index,
                        match.dataset_gene,
                        match.match_type,
                        match.embedding_key or "",
                    ]
                )
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, destination)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def write_markdown_atomic(lines: list[str], destination: Path) -> None:
    write_lines_atomic(lines, destination)


def require_export(dataset: str, force: bool) -> Path:
    export_dir = (
        PROJECT_ROOT / "data/interim/genept_scpa/phase2_export" / dataset
    )
    manifest = export_dir / f"{dataset}_export_manifest.json"
    if force or not manifest.exists():
        subprocess.run(
            [
                "Rscript",
                str(PROJECT_ROOT / "scripts/data/export_seurat_for_genept.R"),
                "--dataset",
                dataset,
            ],
            check=True,
            cwd=PROJECT_ROOT,
        )
    if not manifest.exists():
        raise FileNotFoundError(f"Seurat export manifest was not created: {manifest}")
    return manifest


def load_and_validate_export(manifest_path: Path) -> tuple[dict, sparse.csr_matrix, list[str], list[str]]:
    with manifest_path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    if manifest.get("assay") != "RNA" or manifest.get("layer") != "counts":
        raise ValueError("GenePT-w requires the Seurat RNA/counts layer")
    if manifest.get("matrix_orientation") != "genes_by_cells":
        raise ValueError("Unexpected Seurat export matrix orientation")
    files = {key: Path(value) for key, value in manifest["files"].items()}
    for name, path in files.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing Seurat export file {name}: {path}")
        expected_sha = manifest["sha256"][name]
        if sha256_file(path) != expected_sha:
            raise ValueError(f"Seurat export checksum mismatch: {path}")
    genes = read_lines(files["genes"])
    cell_ids = read_lines(files["cell_ids"])
    if len(genes) != int(manifest["genes"]) or len(cell_ids) != int(manifest["cells"]):
        raise ValueError("Seurat export ID lengths do not match manifest dimensions")
    if len(set(genes)) != len(genes) or len(set(cell_ids)) != len(cell_ids):
        raise ValueError("Seurat export contains duplicate gene or cell IDs")
    genes_by_cells = scipy_io.mmread(files["counts"])
    if not sparse.issparse(genes_by_cells):
        genes_by_cells = sparse.coo_matrix(genes_by_cells)
    if genes_by_cells.shape != (len(genes), len(cell_ids)):
        raise ValueError("Matrix Market dimensions do not match exported IDs")
    counts = genes_by_cells.transpose().tocsr()
    return manifest, counts, genes, cell_ids


def validate_metadata(path: Path, cell_ids: list[str]) -> None:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not {"cell_id", "Hour", "Cell_Type"}.issubset(reader.fieldnames):
            raise ValueError("Exported metadata lacks cell_id, Hour, or Cell_Type")
        observed = [row["cell_id"] for row in reader]
    if observed != cell_ids:
        raise ValueError("Exported metadata cell IDs/order do not match counts")


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def run_synthetic_correctness_check() -> dict[str, float]:
    """Run the hand-calculable GenePT-w gate before any full-data projection."""

    toy_counts = sparse.csr_matrix(
        np.array([[1.0, 1.0, 0.0], [0.0, 2.0, 2.0]], dtype=np.float64)
    )
    toy_indices = np.array([0, 1, 2], dtype=np.int64)
    toy_embeddings = np.array(
        [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]], dtype=np.float32
    )
    weight = np.log1p(5_000.0)
    expected = np.array([[weight, weight], [weight, 2.0 * weight]]) / 3.0
    expected /= np.linalg.norm(expected, axis=1, keepdims=True)
    optimized, _ = project_genept_w_sparse(
        toy_counts,
        toy_indices,
        toy_embeddings,
        normalization_target=10_000.0,
        batch_size=2,
    )
    direct = project_genept_w_direct(
        toy_counts,
        toy_indices,
        toy_embeddings,
        normalization_target=10_000.0,
    )
    repeated, _ = project_genept_w_sparse(
        toy_counts,
        toy_indices,
        toy_embeddings,
        normalization_target=10_000.0,
        batch_size=1,
    )
    result = {
        "hand_calculation_max_abs_error": float(
            np.max(np.abs(optimized.astype(np.float64) - expected))
        ),
        "optimized_vs_direct_max_abs_error": float(
            np.max(np.abs(optimized.astype(np.float64) - direct))
        ),
        "batch_determinism_max_abs_error": float(
            np.max(np.abs(optimized - repeated))
        ),
    }
    if (
        result["hand_calculation_max_abs_error"] >= 1e-6
        or result["optimized_vs_direct_max_abs_error"] >= 1e-6
        or result["batch_determinism_max_abs_error"] != 0.0
    ):
        raise RuntimeError(f"Synthetic GenePT-w correctness gate failed: {result}")
    return result


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=("naive_cd4", "naive_cd8"), default="naive_cd4")
    parser.add_argument("--force-export", action="store_true")
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument(
        "--embedding-provenance",
        type=Path,
        default=PROJECT_ROOT
        / "data/interim/genept_scpa/phase2_genept_embedding_provenance.json",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    settings = DATASET_SETTINGS[args.dataset]
    config = load_config(PROJECT_ROOT / "config/genept_scpa.yaml")
    config.require_phase(int(settings["phase"]))
    cd8_dataset_qc_pass = args.dataset != "naive_cd8"
    if args.dataset == "naive_cd4" and config.active_phase != 2:
        raise RuntimeError(
            "Phase 2 is already PASS; refusing to overwrite completed naïve CD4 outputs"
        )
    if args.dataset == "naive_cd8":
        if config.active_phase != 3:
            raise RuntimeError("Phase 3 must be active to generate naïve CD8 GenePT-w")
        if config.values["phase2"]["status"] != "passed":
            raise RuntimeError("Phase 2 must be PASS before generating naïve CD8 GenePT-w")
        cd8_dataset_qc = (
            PROJECT_ROOT / "data/interim/genept_scpa/naive_cd8_dataset_qc.json"
        )
        with cd8_dataset_qc.open("r", encoding="utf-8") as handle:
            cd8_dataset_qc_pass = json.load(handle)["gate"]["status"] == "PASS"
        if not cd8_dataset_qc_pass:
            raise RuntimeError("Naïve CD8 acquisition/QC gate must be PASS")
    phase1b_qc_path = PROJECT_ROOT / "data/interim/genept_scpa/phase1b_scpa_qc.json"
    with phase1b_qc_path.open("r", encoding="utf-8") as handle:
        phase1b_qc = json.load(handle)
    if phase1b_qc["gate"]["status"] != "PASS":
        raise RuntimeError("Phase 1B gate must be PASS before Phase 2")

    if not args.embedding_provenance.is_file():
        raise FileNotFoundError(
            "Official embedding provenance is missing. Run "
            "scripts/genept/prepare_genept_embeddings.py --download first."
        )
    with args.embedding_provenance.open("r", encoding="utf-8") as handle:
        provenance = json.load(handle)
    embedding_path = Path(provenance["primary_embedding"]["path"])
    summary_path = Path(provenance["primary_gene_summaries"]["path"])
    if sha256_file(embedding_path) != provenance["primary_embedding"]["sha256"]:
        raise ValueError("Official GenePT embedding SHA-256 mismatch")
    if sha256_file(summary_path) != provenance["primary_gene_summaries"]["sha256"]:
        raise ValueError("Official NCBI summary SHA-256 mismatch")

    manifest_path = require_export(args.dataset, args.force_export)
    manifest, counts, genes, cell_ids = load_and_validate_export(manifest_path)
    metadata_source = Path(manifest["files"]["metadata"])
    validate_metadata(metadata_source, cell_ids)

    embeddings = load_official_genept_embeddings(embedding_path, expected_dimension=1536)
    primary_keys = load_primary_gene_keys(summary_path)
    matches = classify_gene_matches(genes, set(embeddings), primary_keys)
    match_summary = mapping_counts(matches)
    matched_indices, embedding_matrix = build_aligned_embedding_matrix(matches, embeddings)

    synthetic_check = run_synthetic_correctness_check()
    crosscheck_cells = min(3, counts.shape[0])
    optimized_small, _ = project_genept_w_sparse(
        counts[:crosscheck_cells],
        matched_indices,
        embedding_matrix,
        normalization_target=10_000.0,
        batch_size=crosscheck_cells,
    )
    direct_small = project_genept_w_direct(
        counts[:crosscheck_cells],
        matched_indices,
        embedding_matrix,
        normalization_target=10_000.0,
    )
    direct_max_abs_error = float(
        np.max(np.abs(optimized_small.astype(np.float64) - direct_small))
    )
    repeated_small, _ = project_genept_w_sparse(
        counts[:crosscheck_cells],
        matched_indices,
        embedding_matrix,
        normalization_target=10_000.0,
        batch_size=1,
    )
    deterministic_max_abs_error = float(
        np.max(np.abs(optimized_small - repeated_small))
    )
    if direct_max_abs_error >= 1e-6 or deterministic_max_abs_error != 0.0:
        raise RuntimeError(
            "Actual-cell GenePT-w preflight gate failed: "
            f"direct_error={direct_max_abs_error}, "
            f"determinism_error={deterministic_max_abs_error}"
        )

    output_dir = PROJECT_ROOT / str(settings["output_directory"])
    output_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = output_dir / f"{args.dataset}_genept_w.npy"
    temporary_matrix = output_dir / f".{args.dataset}_genept_w.tmp.npy"
    if temporary_matrix.exists():
        temporary_matrix.unlink()
    output, diagnostics = project_genept_w_sparse(
        counts,
        matched_indices,
        embedding_matrix,
        normalization_target=10_000.0,
        batch_size=args.batch_size,
        output_path=temporary_matrix,
    )
    del output
    os.replace(temporary_matrix, matrix_path)

    cell_ids_output = output_dir / f"{args.dataset}_genept_w_cell_ids.txt"
    metadata_output = output_dir / f"{args.dataset}_genept_w_metadata.csv"
    mapping_output = output_dir / str(settings["mapping_filename"])
    write_lines_atomic(cell_ids, cell_ids_output)
    copy_atomic(metadata_source, metadata_output)
    write_mapping_atomic(matches, mapping_output)

    expected_values = counts.shape[0] * embedding_matrix.shape[1]
    nonzero_norms = diagnostics.post_l2_norms[diagnostics.pre_l2_norms > 0]
    l2_tolerance = 1e-5
    l2_pass = bool(
        nonzero_norms.size > 0
        and np.all(np.abs(nonzero_norms - 1.0) < l2_tolerance)
    )
    criteria = {
        "phase1b_pass": True,
        "phase2_pass_for_cd8": (
            args.dataset != "naive_cd8"
            or config.values["phase2"]["status"] == "passed"
        ),
        "cd8_dataset_qc_pass": cd8_dataset_qc_pass,
        "official_embedding_dimension": embedding_matrix.shape[1] == 1536,
        "rna_counts_source": manifest["assay"] == "RNA" and manifest["layer"] == "counts",
        "source_cell_count_preserved": counts.shape[0] == settings["expected_cells"],
        "source_feature_count_preserved": counts.shape[1] == settings["expected_features"],
        "cell_ids_unique_and_complete": len(set(cell_ids)) == counts.shape[0],
        "dataset_gene_ids_unique": match_summary["dataset_duplicate_genes"] == 0,
        "duplicate_mapping_absent": match_summary["duplicate_mapping_count"] == 0,
        "finite_output": diagnostics.finite_values == expected_values,
        "unexpected_zero_vectors_absent": diagnostics.zero_vectors == 0,
        "post_l2_norms": l2_pass,
        "optimized_vs_direct": direct_max_abs_error < 1e-6,
        "deterministic_reproduction": deterministic_max_abs_error == 0.0,
        "synthetic_correctness": all(
            error < 1e-6 for error in synthetic_check.values()
        ),
        "source_object_unmodified": manifest.get("source_object_modified") is False,
        "no_downstream_analysis": True,
    }
    failed_checks = [name for name, passed in criteria.items() if not passed]
    gate_status = "READY_FOR_GPT_REVIEW" if not failed_checks else "NEEDS_REVIEW"

    cd4_coverage_comparison = None
    if args.dataset == "naive_cd8":
        cd4_qc_path = PROJECT_ROOT / "data/interim/genept_scpa/phase2_genept_w_qc.json"
        with cd4_qc_path.open("r", encoding="utf-8") as handle:
            cd4_qc = json.load(handle)
        cd4_median_coverage = float(
            cd4_qc["gene_matching"]["expression_coverage_summary"]["median"]
        )
        cd8_median_coverage = float(np.median(diagnostics.expression_coverage))
        cd4_coverage_comparison = {
            "cd4_median_expression_coverage": cd4_median_coverage,
            "cd8_median_expression_coverage": cd8_median_coverage,
            "cd8_minus_cd4_median_coverage": (
                cd8_median_coverage - cd4_median_coverage
            ),
            "exclusion_threshold_applied": False,
        }

    qc = {
        "phase": settings["phase_label"],
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": {
            "accession": "GSE212270",
            "source_file": manifest["source_file"],
            "cells": counts.shape[0],
            "features": counts.shape[1],
        },
        "genept_embedding": {
            "source": provenance["source"],
            "model": "text-embedding-ada-002",
            "dimension": embedding_matrix.shape[1],
            "gene_count": len(embeddings),
            "sha256": provenance["primary_embedding"]["sha256"],
            "file": str(embedding_path),
        },
        "expression_input": {
            "assay": "RNA",
            "layer": "counts",
            "source_dimensions": {"genes": counts.shape[1], "cells": counts.shape[0]},
            "sparse": True,
            "matrix_market_manifest": str(manifest_path),
        },
        "preprocessing": {
            "normalization_target": 10_000,
            "log_transform": "log1p",
            "gene_filtering_order": "normalize_all_dataset_genes_then_project_exact_official_lookup_matches",
            "aggregation": "log_normalized_expression_weighted_sum_divided_by_total_dataset_gene_count",
            "l2_normalization": "rowwise_unit_norm_for_nonzero_vectors",
            "unmatched_gene_handling": "included_in_library_size_normalization_then_zero_projection_contribution",
            "alias_policy": "exact artifact-key match only; official HGNC alias keys distinguished using NCBI primary keys",
        },
        "gene_matching": {
            **match_summary,
            "embedding_genes": len(embeddings),
            "matched_dataset_genes": int(matched_indices.size),
            "unmatched": match_summary["unmatched_dataset_genes"],
            "duplicates": match_summary["duplicate_mapping_count"],
            "expression_coverage_summary": numeric_summary(diagnostics.expression_coverage),
            "mapping_file": str(mapping_output),
            "comparison_with_naive_cd4": cd4_coverage_comparison,
        },
        "output": {
            "matrix_file": str(matrix_path),
            "matrix_sha256": sha256_file(matrix_path),
            "cell_ids_file": str(cell_ids_output),
            "metadata_file": str(metadata_output),
            "cells": counts.shape[0],
            "dimensions": embedding_matrix.shape[1],
            "dtype": "float32",
            "finite_values": diagnostics.finite_values,
            "expected_values": expected_values,
            "zero_vectors": diagnostics.zero_vectors,
            "pre_l2_norm_summary": numeric_summary(diagnostics.pre_l2_norms),
            "post_l2_norm_summary": numeric_summary(diagnostics.post_l2_norms),
            "l2_norm_summary": numeric_summary(diagnostics.post_l2_norms),
            "l2_tolerance": l2_tolerance,
        },
        "correctness": {
            "synthetic_preflight": synthetic_check,
            "optimized_vs_direct_cells": crosscheck_cells,
            "optimized_vs_direct_max_abs_error": direct_max_abs_error,
            "deterministic_max_abs_error": deterministic_max_abs_error,
        },
        "compatibility": {
            "seurat_export": manifest.get("compatibility", {}),
            "official_notebook_adaptation": (
                "The notebook reads preprocessed AnnData.X and performs X@E/number_of_dataset_genes. "
                "This pipeline makes the paper's preceding counts->normalize_total(10000)->log1p "
                "and final rowwise L2 steps explicit while preserving the notebook aggregation."
            ),
        },
        "reproducibility": {
            "Python_version": platform.python_version(),
            "package_versions": {
                "numpy": package_version("numpy"),
                "scipy": package_version("scipy"),
                "PyYAML": package_version("PyYAML"),
            },
            "code_path": str(Path(__file__).resolve()),
            "batch_size": args.batch_size,
        },
        "scope": {
            "SCPA_run": False,
            "CD4_vs_CD8_comparison": False,
            "classifier_training": False,
            "downstream_metrics": False,
        },
        "gate": {
            "status": gate_status,
            "failed_checks": failed_checks,
            "warnings": [],
            "criteria": criteria,
        },
    }
    qc_path = PROJECT_ROOT / str(settings["qc"])
    write_json_atomic(qc, qc_path)
    summary_path = PROJECT_ROOT / str(settings["summary"])
    write_markdown_atomic(
        [
            f"# {settings['phase_label']} summary",
            "",
            f"Gate status: `{gate_status}`",
            "",
            "## Official source and primary embedding",
            "",
            "- Final paper: Chen & Zou, Nature Biomedical Engineering 9, 483–493 (2025).",
            "- Official code: `yiqunchen/GenePT`, including `aorta_data_analysis.ipynb`.",
            "- Artifact: Zenodo DOI `10.5281/zenodo.10833191`.",
            "- Primary model: `text-embedding-ada-002`, 1,536 dimensions.",
            "- OpenAI API calls: none; the author-provided precomputed artifact was used.",
            "",
            "## Input and published preprocessing",
            "",
            f"- Source: `{manifest['source_file']}` (`RNA/counts`, {counts.shape[0]} cells x {counts.shape[1]} genes).",
            "- Raw counts -> cell-wise total-count normalization to 10,000 -> log1p.",
            "- Dataset genes are then aligned by exact official artifact key; no fuzzy or case conversion.",
            "- Expression-weighted GenePT aggregation uses the official notebook denominator (all dataset genes), followed by row-wise unit L2 normalization.",
            "- Unmatched genes participate in library-size normalization and then contribute a zero embedding vector.",
            "",
            "## Gene matching and expression coverage",
            "",
            f"- Embedding keys: {len(embeddings)}",
            f"- Exact primary-symbol matches: {match_summary['exact_matches']}",
            f"- Official HGNC-alias-key matches: {match_summary['alias_matches']}",
            f"- Unmatched dataset genes: {match_summary['unmatched_dataset_genes']}",
            f"- Duplicate mappings: {match_summary['duplicate_mapping_count']}",
            f"- Raw-count mass coverage summary: `{numeric_summary(diagnostics.expression_coverage)}`",
            (
                f"- CD8-minus-CD4 median coverage: "
                f"{cd4_coverage_comparison['cd8_minus_cd4_median_coverage']:.6g} "
                "(recorded as a potential confounder; no exclusion threshold applied)."
                if cd4_coverage_comparison is not None
                else "- Cross-dataset coverage comparison: deferred to naïve CD8 generation."
            ),
            "",
            "## Output and correctness",
            "",
            f"- Matrix: `{matrix_path}` ({counts.shape[0]} x {embedding_matrix.shape[1]}, float32).",
            f"- Finite values: {diagnostics.finite_values}/{expected_values}",
            f"- Zero vectors: {diagnostics.zero_vectors}",
            f"- Post-L2 norm summary: `{numeric_summary(diagnostics.post_l2_norms)}`",
            f"- Optimized-vs-direct maximum absolute error: {direct_max_abs_error:.3g}",
            f"- Deterministic repeat maximum absolute error: {deterministic_max_abs_error:.3g}",
            "",
            "## Compatibility and scope",
            "",
            "The official notebook loads an already-preprocessed AnnData `.X`; this pipeline makes the paper's explicit normalization/log1p and final L2 steps reproducible from Seurat raw counts. The notebook's lookup/aggregation rule is retained.",
            "",
            "No SCPA, CD4-vs-CD8 comparison, classifier, separability metric, or pathway/gene-level interpretation was run by this builder.",
        ],
        summary_path,
    )
    print(
        "GENEPT_W_SUMMARY "
        f"status={gate_status} cells={counts.shape[0]} dimension={embedding_matrix.shape[1]} "
        f"matched_genes={matched_indices.size} qc_json={qc_path.resolve()}"
    )
    return 0 if gate_status == "READY_FOR_GPT_REVIEW" else 1


if __name__ == "__main__":
    raise SystemExit(main())
