#!/usr/bin/env python3
"""Run the Phase 3 CD4 0h vs CD8 0h SCPA-core benchmark."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Sequence

import numpy as np
from scipy import io as scipy_io
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from scripts.genept.build_genept_w import (  # noqa: E402
    load_and_validate_export,
    require_export,
    write_lines_atomic,
    write_markdown_atomic,
)
from gene_embedding_project.genept_scpa.config import load_config  # noqa: E402
from gene_embedding_project.genept_scpa.genept_projection import (  # noqa: E402
    normalize_log1p_sparse,
)
from gene_embedding_project.genept_scpa.io import (  # noqa: E402
    sha256_file,
    write_json_atomic,
)


def canonical_hour(value: object) -> str:
    text = str(value).strip().lower().replace(" ", "")
    if text.endswith(".0"):
        text = text[:-2]
    if text.endswith("hours"):
        text = text[:-5]
    elif text.endswith("hour"):
        text = text[:-4]
    elif text.endswith("hrs"):
        text = text[:-3]
    elif text.endswith("hr"):
        text = text[:-2]
    elif text.endswith("h"):
        text = text[:-1]
    return f"{text}h"


def read_metadata(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not {"cell_id", "Hour", "Cell_Type"}.issubset(
            reader.fieldnames
        ):
            raise ValueError(f"Metadata schema is incomplete: {path}")
        rows = list(reader)
    ids = [row["cell_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError(f"Metadata contains duplicate cell IDs: {path}")
    return rows


def select_canonical_cells(
    metadata: Sequence[dict[str, str]],
    *,
    hour: str,
    sample_size: int,
    seed: int,
) -> tuple[list[str], list[str]]:
    candidates = sorted(
        row["cell_id"]
        for row in metadata
        if canonical_hour(row["Hour"]) == canonical_hour(hour)
    )
    if len(candidates) < sample_size:
        raise ValueError(
            f"Only {len(candidates)} cells are available for {hour}; "
            f"need {sample_size}"
        )
    generator = np.random.default_rng(seed)
    selected_indices = generator.choice(
        len(candidates), size=sample_size, replace=False
    )
    selected = [candidates[int(index)] for index in selected_indices]
    return candidates, selected


def aligned_shared_gene_indices(
    cd4_genes: Sequence[str], cd8_genes: Sequence[str]
) -> tuple[list[str], np.ndarray, np.ndarray]:
    if len(cd4_genes) != len(set(cd4_genes)) or len(cd8_genes) != len(set(cd8_genes)):
        raise ValueError("Original-expression gene IDs must be unique")
    shared = sorted(set(cd4_genes).intersection(cd8_genes))
    if not shared:
        raise ValueError("CD4 and CD8 have no shared gene symbols")
    cd4_lookup = {gene: index for index, gene in enumerate(cd4_genes)}
    cd8_lookup = {gene: index for index, gene in enumerate(cd8_genes)}
    return (
        shared,
        np.asarray([cd4_lookup[gene] for gene in shared], dtype=np.int64),
        np.asarray([cd8_lookup[gene] for gene in shared], dtype=np.int64),
    )


def rows_for_ids(all_ids: Sequence[str], selected_ids: Sequence[str]) -> np.ndarray:
    lookup = {cell_id: index for index, cell_id in enumerate(all_ids)}
    missing = [cell_id for cell_id in selected_ids if cell_id not in lookup]
    if missing:
        raise ValueError(f"Selected cell IDs are missing from a representation: {missing[:3]}")
    return np.asarray([lookup[cell_id] for cell_id in selected_ids], dtype=np.int64)


def write_matrix_market_atomic(matrix: np.ndarray | sparse.spmatrix, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            scipy_io.mmwrite(handle, sparse.csr_matrix(matrix))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--sample-size", type=int, default=500)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_config(PROJECT_ROOT / "config/genept_scpa.yaml")
    config.require_phase(3)
    if config.active_phase != 3 or config.values["phase3"]["status"] != "in_progress":
        raise RuntimeError("Phase 3 must be active and in_progress")
    frozen = config.values["phase3"]["primary_comparison"]
    if args.seed != int(frozen["seed"]) or args.sample_size != int(frozen["sample_size"]):
        raise RuntimeError("Phase 3 seed and sample size are protocol-frozen")

    cd4_qc_path = PROJECT_ROOT / "data/interim/genept_scpa/phase2_genept_w_qc.json"
    cd8_qc_path = PROJECT_ROOT / "data/interim/genept_scpa/phase3_naive_cd8_genept_w_qc.json"
    cd4_qc = load_json(cd4_qc_path)
    cd8_qc = load_json(cd8_qc_path)
    if config.values["phase2"]["status"] != "passed":
        raise RuntimeError("Phase 2 status must be passed")
    if cd4_qc["gate"]["status"] != "READY_FOR_GPT_REVIEW":
        raise RuntimeError("The reviewed Phase 2 QC artifact is not intact")
    if cd8_qc["gate"]["status"] != "READY_FOR_GPT_REVIEW":
        raise RuntimeError("Naïve CD8 GenePT-w must be READY_FOR_GPT_REVIEW")

    cd4_matrix_path = Path(cd4_qc["output"]["matrix_file"])
    cd8_matrix_path = Path(cd8_qc["output"]["matrix_file"])
    cd4_genept = np.load(cd4_matrix_path, mmap_mode="r")
    cd8_genept = np.load(cd8_matrix_path, mmap_mode="r")
    cd4_genept_ids = Path(cd4_qc["output"]["cell_ids_file"]).read_text(
        encoding="utf-8"
    ).splitlines()
    cd8_genept_ids = Path(cd8_qc["output"]["cell_ids_file"]).read_text(
        encoding="utf-8"
    ).splitlines()
    cd4_metadata = read_metadata(Path(cd4_qc["output"]["metadata_file"]))
    cd8_metadata = read_metadata(Path(cd8_qc["output"]["metadata_file"]))
    if cd4_genept.shape != (len(cd4_genept_ids), 1_536):
        raise ValueError("CD4 GenePT matrix shape does not match its cell IDs")
    if cd8_genept.shape != (len(cd8_genept_ids), 1_536):
        raise ValueError("CD8 GenePT matrix shape does not match its cell IDs")
    if [row["cell_id"] for row in cd4_metadata] != cd4_genept_ids:
        raise ValueError("CD4 GenePT and metadata cell order differ")
    if [row["cell_id"] for row in cd8_metadata] != cd8_genept_ids:
        raise ValueError("CD8 GenePT and metadata cell order differ")

    cd4_all_0h, cd4_selected = select_canonical_cells(
        cd4_metadata, hour="0h", sample_size=args.sample_size, seed=args.seed
    )
    cd8_all_0h, cd8_selected = select_canonical_cells(
        cd8_metadata, hour="0h", sample_size=args.sample_size, seed=args.seed
    )
    sampling_dir = PROJECT_ROOT / "data/interim/genept_scpa/phase3_sampling"
    sampling_files = {
        "cd4_all_0h": sampling_dir / "cd4_0h_all_cells.txt",
        "cd8_all_0h": sampling_dir / "cd8_0h_all_cells.txt",
        "cd4_canonical": sampling_dir / "cd4_0h_cells.txt",
        "cd8_canonical": sampling_dir / "cd8_0h_cells.txt",
    }
    write_lines_atomic(cd4_all_0h, sampling_files["cd4_all_0h"])
    write_lines_atomic(cd8_all_0h, sampling_files["cd8_all_0h"])
    write_lines_atomic(cd4_selected, sampling_files["cd4_canonical"])
    write_lines_atomic(cd8_selected, sampling_files["cd8_canonical"])

    cd4_genept_rows = rows_for_ids(cd4_genept_ids, cd4_selected)
    cd8_genept_rows = rows_for_ids(cd8_genept_ids, cd8_selected)
    selected_cd4_genept = np.asarray(cd4_genept[cd4_genept_rows], dtype=np.float64)
    selected_cd8_genept = np.asarray(cd8_genept[cd8_genept_rows], dtype=np.float64)

    cd4_manifest_path = require_export("naive_cd4", force=False)
    cd8_manifest_path = require_export("naive_cd8", force=False)
    cd4_manifest, cd4_counts, cd4_genes, cd4_count_ids = load_and_validate_export(
        cd4_manifest_path
    )
    cd8_manifest, cd8_counts, cd8_genes, cd8_count_ids = load_and_validate_export(
        cd8_manifest_path
    )
    shared_genes, cd4_gene_columns, cd8_gene_columns = aligned_shared_gene_indices(
        cd4_genes, cd8_genes
    )
    cd4_count_rows = rows_for_ids(cd4_count_ids, cd4_selected)
    cd8_count_rows = rows_for_ids(cd8_count_ids, cd8_selected)
    cd4_log = normalize_log1p_sparse(cd4_counts[cd4_count_rows])[:, cd4_gene_columns]
    cd8_log = normalize_log1p_sparse(cd8_counts[cd8_count_rows])[:, cd8_gene_columns]

    adapter_input_dir = PROJECT_ROOT / "data/interim/genept_scpa/phase3_adapter_inputs"
    adapter_inputs = {
        "genept_cd4": adapter_input_dir / "genept_w_cd4_0h_cells_by_features.mtx",
        "genept_cd8": adapter_input_dir / "genept_w_cd8_0h_cells_by_features.mtx",
        "expression_cd4": adapter_input_dir / "expression_cd4_0h_cells_by_genes.mtx",
        "expression_cd8": adapter_input_dir / "expression_cd8_0h_cells_by_genes.mtx",
        "shared_genes": adapter_input_dir / "original_expression_shared_genes.txt",
    }
    write_matrix_market_atomic(selected_cd4_genept, adapter_inputs["genept_cd4"])
    write_matrix_market_atomic(selected_cd8_genept, adapter_inputs["genept_cd8"])
    write_matrix_market_atomic(cd4_log, adapter_inputs["expression_cd4"])
    write_matrix_market_atomic(cd8_log, adapter_inputs["expression_cd8"])
    write_lines_atomic(shared_genes, adapter_inputs["shared_genes"])

    processed_dir = PROJECT_ROOT / "data/processed/genept_scpa/phase3"
    core_result_path = processed_dir / "phase3_scpa_core_results.json"
    subprocess.run(
        [
            "Rscript",
            str(PROJECT_ROOT / "scripts/scpa/run_phase3_scpa_core.R"),
            "--genept-cd4",
            str(adapter_inputs["genept_cd4"]),
            "--genept-cd8",
            str(adapter_inputs["genept_cd8"]),
            "--expression-cd4",
            str(adapter_inputs["expression_cd4"]),
            "--expression-cd8",
            str(adapter_inputs["expression_cd8"]),
            "--output",
            str(core_result_path),
        ],
        check=True,
        cwd=PROJECT_ROOT,
    )
    core_results = load_json(core_result_path)

    expected_cd4_0h = 4_428
    expected_cd8_0h = 1_048
    no_core_warnings = all(
        not analysis["warnings"] for analysis in core_results["analyses"].values()
    )
    criteria = {
        "phase2_pass": True,
        "cd8_genept_ready_for_review": True,
        "cd4_0h_cell_count": len(cd4_all_0h) == expected_cd4_0h,
        "cd8_0h_cell_count": len(cd8_all_0h) == expected_cd8_0h,
        "canonical_sample_size": (
            len(cd4_selected) == args.sample_size == len(cd8_selected)
        ),
        "same_cells_across_representations": True,
        "genept_dimensions": (
            selected_cd4_genept.shape == (args.sample_size, 1_536)
            and selected_cd8_genept.shape == (args.sample_size, 1_536)
        ),
        "selected_genept_values_finite": bool(
            np.isfinite(selected_cd4_genept).all()
            and np.isfinite(selected_cd8_genept).all()
        ),
        "original_expression_exact_shared_genes": len(shared_genes) == 17_085,
        "scpa_core_toy_test": core_results["toy_test"]["passed"] is True,
        "scpa_core_runtime_warnings_absent": no_core_warnings,
        "source_objects_unmodified": (
            cd4_manifest.get("source_object_modified") is False
            and cd8_manifest.get("source_object_modified") is False
        ),
        "classifier_not_run": True,
        "pathway_gene_interpretation_not_run": True,
    }
    failed_checks = [name for name, passed in criteria.items() if not passed]
    status = "READY_FOR_GPT_REVIEW" if not failed_checks else "NEEDS_REVIEW"
    cd4_coverage = cd4_qc["gene_matching"]["expression_coverage_summary"]
    cd8_coverage = cd8_qc["gene_matching"]["expression_coverage_summary"]
    qc = {
        "phase": "Phase 3 - Primary CD4 0h vs CD8 0h GenePT-w benchmark",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "datasets": {
            "cd4": {"accession": "GSE212270", "cells": 14_894, "features": 17_856},
            "cd8": {"accession": "GSE212270", "cells": 7_041, "features": 17_942},
        },
        "genept": {
            "embedding_model": "text-embedding-ada-002",
            "dimension": 1_536,
            "artifact_sha256": cd4_qc["genept_embedding"]["sha256"],
        },
        "cd4_genept": {
            "cells": cd4_qc["output"]["cells"],
            "matched_genes": cd4_qc["gene_matching"]["matched_dataset_genes"],
            "expression_coverage": cd4_coverage,
        },
        "cd8_genept": {
            "cells": cd8_qc["output"]["cells"],
            "matched_genes": cd8_qc["gene_matching"]["matched_dataset_genes"],
            "expression_coverage": cd8_coverage,
        },
        "coverage_comparison": {
            "cd4_median": cd4_coverage["median"],
            "cd8_median": cd8_coverage["median"],
            "cd8_minus_cd4_median": cd8_coverage["median"] - cd4_coverage["median"],
            "exclusion_threshold_applied": False,
        },
        "primary_comparison": {
            "cd4_definition": "naive_cd4 AND Hour=0",
            "cd8_definition": "naive_cd8 AND Hour=0",
            "full_cell_counts": {"cd4": len(cd4_all_0h), "cd8": len(cd8_all_0h)},
            "sampled_cell_counts": {"cd4": len(cd4_selected), "cd8": len(cd8_selected)},
            "seed": args.seed,
            "sampling_method": "numpy Generator(PCG64), without replacement, from sorted 0h IDs",
            "sample_size_basis": "SCPA default and frozen Phase 1 convention: 500 cells/group",
            "canonical_cell_id_files": {
                key: str(path) for key, path in sampling_files.items()
            },
        },
        "scpa_core": {
            "implementation_source": "SCPA 1.6.2 single_comparison calls multicross::mcm; adapter calls the same multicross::mcm directly",
            "adapter": "SCPA-core multivariate framework adaptation, not standard pathway analysis",
            "input_orientation": "cells_by_features",
            "feature_handling": "no pathway min/max filter; use all 1536 GenePT dimensions or all 17085 aligned shared genes within each separate global test",
            "sampling": "explicit canonical 500/group before adapter; no hidden random_cells call",
            "toy_test": core_results["toy_test"],
            "warnings": [
                warning
                for analysis in core_results["analyses"].values()
                for warning in analysis["warnings"]
            ],
        },
        "analyses": core_results["analyses"],
        "original_expression_reference": {
            "source": "sparse RNA/counts",
            "preprocessing": "per-cell total-count normalization to 10000 over all genes, then log1p, then exact shared-symbol alignment",
            "shared_gene_count": len(shared_genes),
            "feature_order": "lexicographically sorted exact shared gene symbols",
            "gene_file": str(adapter_inputs["shared_genes"]),
            "role": "confirm population difference in original space; not rank representation quality",
        },
        "adapter_inputs": {
            key: {"path": str(path), "sha256": sha256_file(path)}
            for key, path in adapter_inputs.items()
        },
        "scope": {
            "classifier_run": False,
            "pathway_gene_interpretation_run": False,
            "timepoint_extension_run": False,
            "robustness_analysis_run": False,
        },
        "gate": {"status": status, "failed_checks": failed_checks, "warnings": [], "criteria": criteria},
    }
    qc_path = PROJECT_ROOT / "data/interim/genept_scpa/phase3_cd4_cd8_qc.json"
    write_json_atomic(qc, qc_path)
    summary_path = PROJECT_ROOT / "data/interim/genept_scpa/phase3_cd4_cd8_summary.md"
    write_markdown_atomic(
        [
            "# Phase 3 CD4 0h vs CD8 0h summary",
            "",
            f"Gate status: `{status}`",
            "",
            "## Research question",
            "",
            "Does GenePT-w preserve a detectable multivariate difference between naïve CD4 and naïve CD8 cells, and can the SCPA core framework detect it?",
            "",
            "## Why 0h",
            "",
            "The primary comparison uses 0h first to reduce activation/time confounding and isolate the cell-type distinction.",
            "",
            "## GenePT-w and coverage",
            "",
            f"- CD4: {cd4_qc['output']['cells']} cells; {cd4_qc['gene_matching']['matched_dataset_genes']} matched genes; median coverage {cd4_coverage['median']:.6g}.",
            f"- CD8: {cd8_qc['output']['cells']} cells; {cd8_qc['gene_matching']['matched_dataset_genes']} matched genes; median coverage {cd8_coverage['median']:.6g}.",
            f"- CD8-minus-CD4 median coverage: {cd8_coverage['median'] - cd4_coverage['median']:.6g}; recorded as a potential confounder without an exclusion threshold.",
            "",
            "## Canonical cohort",
            "",
            f"Actual 0h counts were CD4={len(cd4_all_0h)} and CD8={len(cd8_all_0h)}. A fixed seed ({args.seed}) selected {args.sample_size}/group, and the same IDs were used for both representations.",
            "",
            "## SCPA-core adaptation",
            "",
            "SCPA 1.6.2 calls `multicross::mcm()` on cells-by-genes pathway matrices. This adapter calls that same core function on aligned cells-by-features matrices without pretending the 1,536 GenePT dimensions are genes or a pathway.",
            "",
            "## Original-expression reference",
            "",
            f"RNA/counts were normalized to 10,000 over each dataset's full gene set, log1p transformed, then aligned by {len(shared_genes)} exact shared symbols. This reference only checks that a population difference exists in original space.",
            "",
            "## Interpretation limits",
            "",
            "The two representations have different dimensions and geometry. Their raw MCM p/q values are not comparable representation-quality scores. This phase does not establish that GenePT is better, measure classifier accuracy, or identify important pathways/genes.",
            "",
            "## Next planned work",
            "",
            "Phase 4 will separately freeze separability metrics and True/Permuted/Random controls. The core pathway/gene-level interpretation question is reserved for its own later methodological decision gate.",
        ],
        summary_path,
    )
    print(
        "PHASE3_SUMMARY "
        f"status={status} cd4_0h={len(cd4_all_0h)} cd8_0h={len(cd8_all_0h)} "
        f"sampled={args.sample_size} qc_json={qc_path.resolve()}"
    )
    return 0 if status == "READY_FOR_GPT_REVIEW" else 1


if __name__ == "__main__":
    raise SystemExit(main())
