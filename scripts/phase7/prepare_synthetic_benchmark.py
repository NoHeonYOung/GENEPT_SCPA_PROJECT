#!/usr/bin/env python3
"""Prepare the frozen Phase 7 benchmark after its production gate is unlocked.

The current frozen configuration deliberately refuses production preparation.
This module is implemented now for audit and toy testing, but must not be run on
GSE212270 until a later explicit approval changes the execution gate.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Sequence

import h5py
import numpy as np
from scipy import io as scipy_io
from scipy import sparse
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
import sys
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from scripts.phase3.run_cd4_cd8_benchmark import canonical_hour, rows_for_ids  # noqa: E402
from scripts.phase4.run_pathway_comparison import read_json, write_csv_atomic  # noqa: E402
from scripts.phase4.run_timecourse_validation import read_metadata  # noqa: E402
from gene_embedding_project.genept_scpa.config import load_config  # noqa: E402
from gene_embedding_project.genept_scpa.gene_mapping import load_official_genept_embeddings  # noqa: E402
from gene_embedding_project.genept_scpa.genept_projection import normalize_log1p_sparse  # noqa: E402
from gene_embedding_project.genept_scpa.io import sha256_file, write_json_atomic  # noqa: E402
from gene_embedding_project.genept_scpa.phase7.cohort import split_pseudo_conditions  # noqa: E402
from gene_embedding_project.genept_scpa.phase7.pathway_selection import select_phase7_pathways  # noqa: E402
from gene_embedding_project.genept_scpa.phase7.schemas import (  # noqa: E402
    SYNTHETIC_SCHEMA_VERSION, validate_synthetic_manifest,
)
from gene_embedding_project.genept_scpa.phase7.synthetic_perturbation import inject_perturbation  # noqa: E402


def load_phase7_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict) or config.get("phase") != 7:
        raise ValueError("Invalid Phase 7 configuration")
    return config


def load_portable_counts_export(
    manifest_path: Path, configured_files: dict[str, str],
    *, repository_root: Path = PROJECT_ROOT,
) -> tuple[dict[str, Any], sparse.csr_matrix, list[str], list[str]]:
    """Load exact repo-relative transfer paths and verify the frozen export hashes."""

    manifest = read_json(manifest_path)
    if manifest.get("assay") != "RNA" or manifest.get("layer") != "counts":
        raise ValueError("Phase 7 requires the RNA/counts export")
    if manifest.get("matrix_orientation") != "genes_by_cells":
        raise ValueError("Phase 7 export must remain genes_by_cells")
    required = {"counts", "genes", "cell_ids", "metadata"}
    if set(configured_files) != required:
        raise ValueError("Phase 7 configured export file map is incomplete")
    paths = {key: repository_root / configured_files[key] for key in required}
    for key, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"Exact Phase 7 transfer path is missing: {path}")
        if sha256_file(path) != manifest["sha256"][key]:
            raise ValueError(f"Phase 7 transferred export checksum mismatch: {key}")
    genes = paths["genes"].read_text(encoding="utf-8").splitlines()
    cell_ids = paths["cell_ids"].read_text(encoding="utf-8").splitlines()
    if len(genes) != int(manifest["genes"]) or len(cell_ids) != int(manifest["cells"]):
        raise ValueError("Phase 7 transferred export dimensions differ from manifest")
    if len(genes) != len(set(genes)) or len(cell_ids) != len(set(cell_ids)):
        raise ValueError("Phase 7 transferred gene/cell IDs must remain unique")
    genes_by_cells = scipy_io.mmread(paths["counts"])
    if not sparse.issparse(genes_by_cells):
        genes_by_cells = sparse.coo_matrix(genes_by_cells)
    if genes_by_cells.shape != (len(genes), len(cell_ids)):
        raise ValueError("Phase 7 transferred counts axes differ from ID files")
    return manifest, genes_by_cells.transpose().tocsr(), genes, cell_ids


def _experiment_specs(config: dict[str, Any]) -> list[tuple[str, float]]:
    specifications: list[tuple[str, float]] = []
    for scenario in config["perturbations"]["scenarios"]:
        strengths = (
            config["perturbations"]["null_strengths_sd"]
            if scenario == "null" else config["perturbations"]["strengths_sd"]
        )
        specifications.extend((str(scenario), float(alpha)) for alpha in strengths)
    return specifications


def _write_h5_atomic(
    path: Path,
    experiment_payloads: Sequence[dict[str, Any]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    string_dtype = h5py.string_dtype(encoding="utf-8")
    try:
        with tempfile.NamedTemporaryFile(
            dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", delete=False
        ) as temporary:
            temporary_name = temporary.name
        with h5py.File(temporary_name, "w") as handle:
            handle.attrs["schema_version"] = SYNTHETIC_SCHEMA_VERSION
            handle.attrs["orientation"] = "cells_by_genes"
            for item in experiment_payloads:
                group = handle.require_group(f"experiments/{item['experiment_id']}")
                group.create_dataset(
                    "gene_names", data=np.asarray(item["genes"], dtype=object), dtype=string_dtype
                )
                group.create_dataset("embeddings", data=item["embeddings"], compression="gzip", shuffle=True)
                for condition in ("A", "B"):
                    condition_group = group.require_group(f"condition_{condition}")
                    condition_group.create_dataset(
                        "expression", data=item[f"condition_{condition.lower()}"],
                        compression="gzip", shuffle=True,
                    )
                    condition_group.create_dataset(
                        "cell_ids", data=np.asarray(item[f"condition_{condition.lower()}_ids"], dtype=object),
                        dtype=string_dtype,
                    )
        os.replace(temporary_name, path)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def prepare(config_path: Path) -> dict[str, Any]:
    protocol = load_config(PROJECT_ROOT / "config/genept_scpa.yaml")
    protocol.require_phase(7)
    config = load_phase7_config(config_path)
    if not config["execution_gate"]["production_synthetic_generation_allowed"]:
        raise RuntimeError(
            "Phase 7 production synthetic generation remains locked; toy smoke only"
        )

    source = config["source"]
    export_path = PROJECT_ROOT / source["counts_export"]
    export_manifest, counts, dataset_genes, all_cell_ids = load_portable_counts_export(
        export_path, source["counts_export_files"]
    )
    metadata = read_metadata(PROJECT_ROOT / source["counts_export_files"]["metadata"])
    cd4_zero_ids = [row["cell_id"] for row in metadata if canonical_hour(row["Hour"]) == "0h"]
    split = split_pseudo_conditions(
        cd4_zero_ids,
        cells_per_condition=int(config["cohort"]["cells_per_pseudo_condition"]),
        seed=int(config["cohort"]["split_seed"]),
    )
    selected_rows = rows_for_ids(all_cell_ids, split.all_ids)
    normalized = normalize_log1p_sparse(
        counts[selected_rows],
        normalization_target=float(config["normalization"]["full_transcriptome_total"]),
    )
    n = len(split.condition_a_ids)
    condition_a_all = normalized[:n]
    condition_b_all = normalized[n:]

    phase4_manifest_path = PROJECT_ROOT / source["phase4_manifest"]
    phase4_manifest = read_json(phase4_manifest_path)
    with (PROJECT_ROOT / source["descriptions"]).open("r", encoding="utf-8") as handle:
        descriptions = json.load(handle)
    selected_pathways, selection_audit = select_phase7_pathways(
        phase4_manifest, descriptions,
        size_bins=config["pathways"]["size_bins"],
        source_bin_quota=config["pathways"]["source_bin_quota"],
        seed=int(config["pathways"]["selection_seed"]),
    )
    if len(selected_pathways) != int(config["pathways"]["expected_selected_pathways"]):
        raise RuntimeError("Frozen Phase 7 pathway count changed")
    embeddings_all = load_official_genept_embeddings(PROJECT_ROOT / source["embeddings"])
    gene_index = {gene: index for index, gene in enumerate(dataset_genes)}

    cell_rows = [
        {
            "split_id": "split_001", "source_cell_id": cell_id,
            "pseudo_condition": condition, "within_condition_index": index,
            "sampling_seed": split.seed,
        }
        for condition, cells in (("A", split.condition_a_ids), ("B", split.condition_b_ids))
        for index, cell_id in enumerate(cells, start=1)
    ]
    specifications = _experiment_specs(config)
    experiments: list[dict[str, Any]] = []
    ground_truth_rows: list[dict[str, Any]] = []
    h5_payloads: list[dict[str, Any]] = []
    truth_base = int(config["ground_truth"]["selection_seed_base"])
    for pathway_index, pathway in enumerate(selected_pathways, start=1):
        genes = list(pathway["analysis_genes"])
        columns = np.asarray([gene_index[gene] for gene in genes], dtype=np.int64)
        baseline_a = condition_a_all[:, columns].toarray().astype(np.float64)
        baseline_b = condition_b_all[:, columns].toarray().astype(np.float64)
        embedding = np.stack([embeddings_all[gene] for gene in genes]).astype(np.float64)
        scenario_seed_index = {
            scenario: index
            for index, scenario in enumerate(config["perturbations"]["scenarios"], start=1)
        }
        for scenario, alpha in specifications:
            # Strength is deliberately excluded: the same pathway/scenario truth
            # set is reused at every alpha, as required by the frozen protocol.
            truth_seed = truth_base + pathway_index * 1000 + scenario_seed_index[scenario]
            result = inject_perturbation(
                baseline_a, baseline_b, genes, scenario=scenario, alpha=alpha,
                seed=truth_seed,
                cell_subset_fraction=float(config["perturbations"]["cell_subset_fraction"]),
                detection_fraction_min=float(config["ground_truth"]["detection_fraction_min"]),
                negative_detection_fraction_min=float(config["ground_truth"]["negative_detection_fraction_min"]),
                negative_median_min=float(config["ground_truth"]["negative_median_min_exclusive"]),
                scale_floor=float(config["ground_truth"]["scale_floor"]),
            )
            experiment_id = (
                f"p{pathway_index:02d}_{pathway['pathway'].lower()}_"
                f"{scenario}_a{int(round(alpha * 100)):03d}"
            )
            experiments.append({
                "experiment_id": experiment_id,
                "split_id": "split_001",
                "pathway": pathway["pathway"],
                "source_database": pathway["source_database"],
                "analysis_gene_count": len(genes),
                "analysis_genes": genes,
                "perturbation_type": scenario,
                "perturbation_strength": alpha,
                "random_seed": truth_seed,
                "condition_labels": {"A": "pseudo_A", "B": "pseudo_B"},
                "injection_space": "normalized_log1p_expression",
                "post_injection_renormalization": False,
                "expression_h5_group": f"/experiments/{experiment_id}",
            })
            truth_lookup = {row["gene"]: row for row in result.ground_truth_rows}
            for gene_position, gene in enumerate(genes):
                base = {
                    "experiment_id": experiment_id, "pathway": pathway["pathway"],
                    "source_database": pathway["source_database"], "gene": gene,
                    "candidate_id": f"C{gene_position + 1:03d}", "gene_index": gene_position,
                    "is_ground_truth_perturbed": False, "perturbation_direction": "none",
                    "perturbation_strength": alpha, "target_cell_fraction": 0.0,
                    "target_cell_count": 0, "pooled_baseline_sd": "",
                    "applied_log_delta": 0.0, "clipped_cell_count": 0,
                    "ground_truth_seed": truth_seed, "description_available": True,
                }
                if gene in truth_lookup:
                    base.update(truth_lookup[gene])
                ground_truth_rows.append(base)
            h5_payloads.append({
                "experiment_id": experiment_id, "genes": genes, "embeddings": embedding,
                "condition_a": result.condition_a, "condition_b": result.condition_b,
                "condition_a_ids": split.condition_a_ids, "condition_b_ids": split.condition_b_ids,
            })

    processed = PROJECT_ROOT / config["artifacts"]["processed_directory"]
    interim = PROJECT_ROOT / config["artifacts"]["interim_directory"]
    processed.mkdir(parents=True, exist_ok=True)
    interim.mkdir(parents=True, exist_ok=True)
    cells_path = processed / config["artifacts"]["cell_assignments"]
    truth_path = processed / config["artifacts"]["ground_truth"]
    h5_path = interim / config["artifacts"]["expression_h5"]
    write_csv_atomic(cell_rows, list(cell_rows[0]), cells_path)
    write_csv_atomic(ground_truth_rows, list(ground_truth_rows[0]), truth_path)
    _write_h5_atomic(h5_path, h5_payloads)
    manifest = {
        "schema_version": SYNTHETIC_SCHEMA_VERSION,
        "phase": 7,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "backend_scope": "method_inputs_only",
        "ground_truth_embedded": False,
        "source": {"dataset": "GSE212270", "population": "naive_cd4_0h"},
        "cell_split": {
            "split_id": "split_001", "random_seed": split.seed,
            "condition_a_count": len(split.condition_a_ids),
            "condition_b_count": len(split.condition_b_ids),
            "assignment_file": str(cells_path), "assignment_sha256": sha256_file(cells_path),
        },
        "normalization": config["normalization"],
        "pathway_selection": selection_audit,
        "pathways": selected_pathways,
        "experiments": experiments,
        "expression_h5": str(h5_path),
        "expression_h5_sha256": sha256_file(h5_path),
        "ground_truth_file": str(truth_path),
        "ground_truth_sha256": sha256_file(truth_path),
        "phase4_manifest_sha256": sha256_file(phase4_manifest_path),
        "method_runners_may_read_ground_truth": False,
        "execution_gate": config["execution_gate"],
    }
    validate_synthetic_manifest(manifest)
    manifest_path = processed / config["artifacts"]["manifest"]
    write_json_atomic(manifest, manifest_path)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path,
        default=PROJECT_ROOT / "config/phase7_gpt_oss_synthetic.yaml",
    )
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    prepare(arguments.config)
