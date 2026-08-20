#!/usr/bin/env python3
"""Prepare the frozen Phase 7 LLM-free synthetic benchmark (no SCPA run)."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
from typing import Any, Optional

import h5py
import numpy as np
from scipy import io as scipy_io
from scipy import sparse
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from scripts.phase3.run_cd4_cd8_benchmark import canonical_hour, rows_for_ids  # noqa: E402
from scripts.phase4.run_pathway_comparison import read_json  # noqa: E402
from scripts.phase4.run_timecourse_validation import read_metadata  # noqa: E402
from gene_embedding_project.genept_scpa.gene_mapping import load_official_genept_embeddings  # noqa: E402
from gene_embedding_project.genept_scpa.genept_projection import normalize_log1p_sparse  # noqa: E402
from gene_embedding_project.genept_scpa.io import sha256_file, write_json_atomic  # noqa: E402
from gene_embedding_project.genept_scpa.phase7.synthetic_benchmark_llmfree.cohort import (  # noqa: E402
    split_pseudo_conditions,
)
from gene_embedding_project.genept_scpa.phase7.synthetic_benchmark_llmfree.perturbation import (  # noqa: E402
    inject_perturbation,
)


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("phase") != 7:
        raise ValueError("Invalid Phase 7 LLM-free configuration")
    if not config["execution"]["llm_backend_forbidden"]:
        raise ValueError("Phase 7 replacement must forbid LLM backends")
    if int(config["ground_truth"]["draw_count"]) < 20:
        raise ValueError("Frozen protocol requires at least 20 draws")
    return config


def write_csv_atomic(rows: list[dict[str, Any]], path: Path) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.partial")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def load_counts_export(config: dict[str, Any]):
    source = config["source"]
    manifest = read_json(PROJECT_ROOT / source["counts_export"])
    if (manifest.get("assay"), manifest.get("layer"), manifest.get("matrix_orientation")) != (
        "RNA", "counts", "genes_by_cells"
    ):
        raise ValueError("Phase 7 requires the RNA counts genes-by-cells export")
    paths = {key: PROJECT_ROOT / value for key, value in source["counts_export_files"].items()}
    for key, path in paths.items():
        if not path.is_file() or sha256_file(path) != manifest["sha256"][key]:
            raise ValueError(f"Missing or checksum-mismatched Phase 2 export: {key}")
    genes = paths["genes"].read_text(encoding="utf-8").splitlines()
    cell_ids = paths["cell_ids"].read_text(encoding="utf-8").splitlines()
    matrix = scipy_io.mmread(paths["counts"])
    if not sparse.issparse(matrix):
        matrix = sparse.coo_matrix(matrix)
    if matrix.shape != (len(genes), len(cell_ids)):
        raise ValueError("Counts axes do not match gene/cell identifiers")
    return manifest, matrix.transpose().tocsr(), genes, cell_ids, paths["metadata"]


def frozen_pathways(config: dict[str, Any], phase4_manifest: dict[str, Any]):
    lookup = {row["pathway"]: row for row in phase4_manifest["pathways"]}
    names = list(config["pathways"]["names"])
    if len(names) != int(config["pathways"]["expected_count"]) or len(names) != len(set(names)):
        raise ValueError("Frozen pathway names are not exactly 11 unique pathways")
    selected = []
    for index, name in enumerate(names, start=1):
        if name not in lookup:
            raise ValueError(f"Frozen pathway missing from Phase 4 manifest: {name}")
        source = str(lookup[name]["source_database"])
        genes = list(lookup[name]["paired_genes"])
        selected.append({
            "pathway_id": f"p{index:02d}", "pathway_index": index,
            "pathway": name, "source_database": source,
            "analysis_genes": genes, "analysis_gene_count": len(genes),
        })
    counts = {source: sum(row["source_database"] == source for row in selected)
              for source in ("KEGG", "REACTOME", "HALLMARK")}
    if counts != {"KEGG": 6, "REACTOME": 5, "HALLMARK": 0}:
        raise ValueError(f"Frozen 6/5/0 pathway composition changed: {counts}")
    return selected


def experiment_specs(config: dict[str, Any]):
    for scenario in config["perturbations"]["scenarios"]:
        strengths = (config["perturbations"]["null_strengths_sd"] if scenario == "null"
                     else config["perturbations"]["strengths_sd"])
        for strength in strengths:
            yield str(scenario), float(strength)


def array_input_hash(prefix: str, *arrays: np.ndarray) -> str:
    digest = hashlib.sha256(prefix.encode("utf-8"))
    for array in arrays:
        contiguous = np.ascontiguousarray(array)
        digest.update(str(contiguous.shape).encode("ascii"))
        digest.update(str(contiguous.dtype).encode("ascii"))
        digest.update(contiguous.tobytes(order="C"))
    return digest.hexdigest()


def prepare(config_path: Path) -> dict[str, Any]:
    config = load_config(config_path)
    _, counts, dataset_genes, all_cell_ids, metadata_path = load_counts_export(config)
    metadata = read_metadata(metadata_path)
    zero_hour_ids = [row["cell_id"] for row in metadata if canonical_hour(row["Hour"]) == "0h"]
    split = split_pseudo_conditions(
        zero_hour_ids,
        cells_per_condition=int(config["cohort"]["cells_per_pseudo_condition"]),
        seed=int(config["cohort"]["split_seed"]),
    )
    selected_rows = rows_for_ids(all_cell_ids, split.all_ids)
    normalized = normalize_log1p_sparse(
        counts[selected_rows],
        normalization_target=float(config["normalization"]["full_transcriptome_total"]),
    )
    n = len(split.condition_a_ids)
    all_a, all_b = normalized[:n], normalized[n:]

    phase4_path = PROJECT_ROOT / config["source"]["phase4_manifest"]
    pathways = frozen_pathways(config, read_json(phase4_path))
    embeddings = load_official_genept_embeddings(PROJECT_ROOT / config["source"]["embeddings"])
    dataset_lookup = {gene: index for index, gene in enumerate(dataset_genes)}
    for pathway in pathways:
        missing = [gene for gene in pathway["analysis_genes"]
                   if gene not in dataset_lookup or gene not in embeddings]
        if missing:
            raise ValueError(f"Frozen pathway has unmapped genes: {pathway['pathway']} {missing[:3]}")

    artifact_cfg = config["artifacts"]
    interim = PROJECT_ROOT / artifact_cfg["interim_directory"]
    processed = PROJECT_ROOT / artifact_cfg["processed_directory"]
    interim.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)
    h5_path = interim / artifact_cfg["expression_h5"]
    truth_path = processed / artifact_cfg["ground_truth"]
    cells_path = processed / artifact_cfg["cell_assignments"]
    temporary_h5: Optional[str] = None
    truth_rows: list[dict[str, Any]] = []
    experiments: list[dict[str, Any]] = []
    scenario_index = {name: index for index, name in enumerate(config["perturbations"]["scenarios"], 1)}
    base_seed = int(config["ground_truth"]["perturbation_seed_base"])
    specs = list(experiment_specs(config))
    total_experiments = len(pathways) * int(config["ground_truth"]["draw_count"]) * len(specs)
    completed = 0

    try:
        with tempfile.NamedTemporaryFile(dir=interim, prefix=".phase7.", suffix=".h5", delete=False) as handle:
            temporary_h5 = handle.name
        string_dtype = h5py.string_dtype(encoding="utf-8")
        with h5py.File(temporary_h5, "w") as h5:
            h5.attrs["schema_version"] = 2
            h5.attrs["orientation"] = "cells_by_genes"
            h5.attrs["post_injection_renormalization"] = False
            for pathway in pathways:
                genes = pathway["analysis_genes"]
                columns = np.asarray([dataset_lookup[gene] for gene in genes], dtype=np.int64)
                baseline_a = all_a[:, columns].toarray().astype(np.float64)
                baseline_b = all_b[:, columns].toarray().astype(np.float64)
                embedding = np.stack([embeddings[gene] for gene in genes]).astype(np.float64)
                pathway_input_hash = array_input_hash("\0".join(genes), baseline_a, embedding)
                pg = h5.require_group(f"pathways/{pathway['pathway_id']}")
                pg.create_dataset("gene_names", data=np.asarray(genes, dtype=object), dtype=string_dtype)
                pg.create_dataset("embeddings", data=embedding, compression="gzip", shuffle=True)
                pg.create_dataset("condition_A", data=baseline_a, compression="gzip", shuffle=True)
                pg.create_dataset("condition_B_baseline", data=baseline_b, compression="gzip", shuffle=True)

                for draw_id in range(1, int(config["ground_truth"]["draw_count"]) + 1):
                    for scenario, strength in specs:
                        seed = (base_seed + 100000 * pathway["pathway_index"]
                                + 1000 * scenario_index[scenario] + draw_id)
                        result = inject_perturbation(
                            baseline_a, baseline_b, genes, scenario=scenario, alpha=strength,
                            seed=seed,
                            cell_subset_fraction=float(config["perturbations"]["cell_subset_fraction"]),
                            detection_fraction_min=float(config["ground_truth"]["detection_fraction_min"]),
                            negative_detection_fraction_min=float(config["ground_truth"]["negative_detection_fraction_min"]),
                            negative_median_min=float(config["ground_truth"]["negative_median_min_exclusive"]),
                            scale_floor=float(config["ground_truth"]["scale_floor"]),
                        )
                        fallback_count = sum(
                            "fallback" in str(row["direction_selection_rule"])
                            for row in result.truth_rows
                        )
                        input_hash = array_input_hash(pathway_input_hash, result.condition_b)
                        alpha_code = int(round(strength * 100))
                        experiment_id = (
                            f"d{draw_id:02d}_{pathway['pathway_id']}_{scenario}_a{alpha_code:03d}"
                        )
                        eg = h5.require_group(f"experiments/{experiment_id}")
                        eg.create_dataset("condition_B", data=result.condition_b,
                                          compression="gzip", shuffle=True)
                        experiments.append({
                            "experiment_id": experiment_id, "draw_id": draw_id,
                            "pathway_id": pathway["pathway_id"], "pathway": pathway["pathway"],
                            "source_database": pathway["source_database"],
                            "analysis_gene_count": len(genes), "perturbation_type": scenario,
                            "perturbation_strength": strength, "perturbation_seed": seed,
                            "truth_fallback_used": fallback_count > 0,
                            "truth_fallback_gene_count": fallback_count,
                            "experiment_input_sha256": input_hash,
                            "condition_a_h5": f"/pathways/{pathway['pathway_id']}/condition_A",
                            "condition_b_h5": f"/experiments/{experiment_id}/condition_B",
                            "gene_names_h5": f"/pathways/{pathway['pathway_id']}/gene_names",
                            "embeddings_h5": f"/pathways/{pathway['pathway_id']}/embeddings",
                            "post_injection_renormalization": False,
                        })
                        target_lookup = {row["gene"]: row for row in result.truth_rows}
                        for gene_index, gene in enumerate(genes):
                            row = {
                                "experiment_id": experiment_id, "draw_id": draw_id,
                                "pathway": pathway["pathway"], "source_database": pathway["source_database"],
                                "perturbation_type": scenario, "perturbation_strength": strength,
                                "perturbation_seed": seed, "gene": gene, "gene_index": gene_index,
                                "experiment_input_sha256": input_hash,
                                "truth_fallback_used": fallback_count > 0,
                                "truth_fallback_gene_count": fallback_count,
                                "is_evaluation_target": False, "is_ground_truth_perturbed": False,
                                "perturbation_direction": "none",
                                "direction_selection_rule": "not_selected",
                                "target_cell_fraction": 0.0,
                                "target_cell_count": 0, "pooled_baseline_sd": "",
                                "applied_log_delta": 0.0, "clipped_cell_count": 0,
                            }
                            if gene in target_lookup:
                                row.update(target_lookup[gene])
                            truth_rows.append(row)
                        completed += 1
                        if completed % 50 == 0 or completed == total_experiments:
                            print(f"[Phase 7 prepare] {completed}/{total_experiments} experiments written", flush=True)
        os.replace(temporary_h5, h5_path)
        temporary_h5 = None
    finally:
        if temporary_h5 and os.path.exists(temporary_h5):
            os.unlink(temporary_h5)

    cell_rows = [
        {"source_cell_id": cell, "pseudo_condition": condition,
         "within_condition_index": index, "cohort_seed": split.seed}
        for condition, ids in (("A", split.condition_a_ids), ("B", split.condition_b_ids))
        for index, cell in enumerate(ids, 1)
    ]
    write_csv_atomic(cell_rows, cells_path)
    write_csv_atomic(truth_rows, truth_path)
    expected_mcm = sum(2 * (int(row["analysis_gene_count"]) + 1) for row in experiments)
    fallback_rows = [row for row in truth_rows if "fallback" in str(row["direction_selection_rule"])]
    fallback_experiments = sorted({row["experiment_id"] for row in fallback_rows})
    manifest = {
        "schema_version": 2, "phase": 7,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "llm_free_vanilla_vs_genept_ground_truth_recovery",
        "source": {"dataset": "GSE212270", "population": "naive_cd4_0h"},
        "cell_split": {"cohort_seed": split.seed, "condition_a_count": n,
                       "condition_b_count": n, "assignment_file": str(cells_path),
                       "assignment_sha256": sha256_file(cells_path)},
        "normalization": config["normalization"], "pathways": pathways,
        "draw_count": int(config["ground_truth"]["draw_count"]),
        "experiments": experiments, "experiment_count": len(experiments),
        "expected_mcm_count": expected_mcm,
        "expression_h5": str(h5_path), "expression_h5_sha256": sha256_file(h5_path),
        "ground_truth_file": str(truth_path), "ground_truth_sha256": sha256_file(truth_path),
        "phase4_manifest_sha256": sha256_file(phase4_path),
        "method_runner_reads_ground_truth": False,
        "masking_implementation": config["scpa_masking"],
        "selection_fallback_audit": {
            "fallback_gene_rows": len(fallback_rows),
            "fallback_experiment_count": len(fallback_experiments),
            "fallback_experiment_ids": fallback_experiments,
            "scientific_results_seen_before_freeze": False,
        },
        "llm_backend_present": False,
    }
    manifest_path = processed / artifact_cfg["manifest"]
    write_json_atomic(manifest, manifest_path)
    print(f"PHASE7_PREP status=PASS experiments={len(experiments)} expected_mcm={expected_mcm}")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path,
                        default=PROJECT_ROOT / "config/phase7_llmfree_synthetic.yaml")
    args = parser.parse_args()
    prepare(args.config)
