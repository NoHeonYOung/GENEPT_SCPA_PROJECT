#!/usr/bin/env python3
"""Prepare truth-blind Phase 8 control matrices and execution manifest."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any

import h5py
import numpy as np
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))


CONTROL_PROTOCOL = "phase8_deranged_and_norm_matched_random_v1"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_array(array: np.ndarray) -> str:
    value = np.ascontiguousarray(array, dtype="<f8")
    return hashlib.sha256(value.tobytes(order="C")).hexdigest()


def write_json_atomic(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.partial")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def make_derangement(size: int, seed: int) -> np.ndarray:
    if size < 2:
        raise ValueError("A derangement requires at least two rows")
    rng = np.random.default_rng(seed)
    identity = np.arange(size, dtype=np.int64)
    for _ in range(10_000):
        permutation = rng.permutation(size)
        if not np.any(permutation == identity):
            return permutation
    raise RuntimeError(f"Could not generate a derangement for K={size}")


def make_random_projection(true_embeddings: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    random = rng.standard_normal(true_embeddings.shape)
    generated_norms = np.linalg.norm(random, axis=1)
    if np.any(generated_norms == 0):
        raise RuntimeError("Generated a zero-norm random direction")
    true_norms = np.linalg.norm(true_embeddings, axis=1)
    return random / generated_norms[:, None] * true_norms[:, None]


def _experiment_hash(experiment: dict[str, Any], permuted_sha: str, random_sha: str) -> str:
    fields = {
        "control_protocol": CONTROL_PROTOCOL,
        "phase7_experiment_input_sha256": experiment["experiment_input_sha256"],
        "permuted_sha256": permuted_sha,
        "random_sha256": random_sha,
    }
    encoded = json.dumps(fields, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def prepare(config_path: Path) -> dict[str, Any]:
    config_path = config_path.resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    phase8_protocol = PROJECT_ROOT / "docs/phase8_mean_shift_mechanism_protocol.md"
    source = {name: PROJECT_ROOT / value for name, value in config["source"].items()}
    phase7_qc = json.loads(source["phase7_scpa_qc"].read_text(encoding="utf-8"))
    if not (phase7_qc.get("status") == "PASS" and phase7_qc.get("partial_run") is False
            and int(phase7_qc.get("mcm_count", -1)) == 101920
            and int(phase7_qc.get("failed_mcm_calls", -1)) == 0):
        raise RuntimeError("Phase 7 production MCM is not a complete PASS")
    phase7_manifest = json.loads(source["phase7_manifest"].read_text(encoding="utf-8"))
    artifacts = config["artifacts"]
    interim = PROJECT_ROOT / artifacts["interim_directory"]
    processed = PROJECT_ROOT / artifacts["processed_directory"]
    phase7_processed = (PROJECT_ROOT / "data/processed/genept_scpa/phase7_llmfree_synthetic").resolve()
    phase7b_processed = (PROJECT_ROOT / "data/processed/genept_scpa/phase7b_null_calibration").resolve()
    if processed.resolve() in {phase7_processed, phase7b_processed}:
        raise RuntimeError("Phase 8 output overlaps an existing phase")
    interim.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)

    selected = [
        experiment for experiment in phase7_manifest["experiments"]
        if experiment["perturbation_type"] == "null"
        or experiment["perturbation_type"] == "mean_shift"
    ]
    pathways = phase7_manifest["pathways"]
    sum_k = sum(int(row["analysis_gene_count"]) for row in pathways)
    sum_k_plus_one = sum(int(row["analysis_gene_count"]) + 1 for row in pathways)
    expected_per_control = 20 * 3 * sum_k_plus_one
    frozen = config["workload"]
    if (len(selected) != int(config["scope_filter"]["experiment_count"])
            or sum_k != int(config["scope_filter"]["pathway_gene_count_sum"])
            or expected_per_control != int(frozen["expected_mcm_per_control"])):
        raise RuntimeError(
            f"Unexpected workload: experiments={len(selected)} sumK={sum_k} "
            f"per_control={expected_per_control}"
        )

    controls_path = interim / artifacts["controls_h5"]
    temporary_h5 = controls_path.with_name(f".{controls_path.name}.{os.getpid()}.partial")
    control_units: list[dict[str, Any]] = []
    perm_seed_base = int(config["representations"]["permuted_genept"]["seed_base"])
    random_seed_base = int(config["representations"]["random_projection"]["seed_base"])
    with h5py.File(source["phase7_expression_h5"], "r") as source_h5, h5py.File(temporary_h5, "w") as output_h5:
        output_h5.attrs["control_protocol"] = CONTROL_PROTOCOL
        output_h5.attrs["ground_truth_labels_parsed"] = False
        for pathway in pathways:
            pathway_id = pathway["pathway_id"]
            pathway_index = int(pathway["pathway_index"])
            true_embeddings = np.asarray(
                source_h5[f"pathways/{pathway_id}/embeddings"][:], dtype=np.float64
            )
            gene_count, dimension = true_embeddings.shape
            if gene_count != int(pathway["analysis_gene_count"]) or dimension != 1536:
                raise ValueError(f"Unexpected TRUE embedding shape for {pathway_id}: {true_embeddings.shape}")
            true_norms = np.linalg.norm(true_embeddings, axis=1)
            for draw_id in range(1, 21):
                permuted_seed = perm_seed_base + 100_000 * pathway_index + draw_id
                random_seed = random_seed_base + 100_000 * pathway_index + draw_id
                permutation = make_derangement(gene_count, permuted_seed)
                permuted = true_embeddings[permutation].copy()
                random = make_random_projection(true_embeddings, random_seed)
                fixed_points = int(np.sum(permutation == np.arange(gene_count)))
                perm_norm_diff = float(np.max(np.abs(
                    np.sort(np.linalg.norm(permuted, axis=1)) - np.sort(true_norms)
                )))
                random_norm_diff = float(np.max(np.abs(np.linalg.norm(random, axis=1) - true_norms)))
                if fixed_points != 0 or perm_norm_diff != 0.0:
                    raise RuntimeError(f"Invalid derangement control for {pathway_id}/draw {draw_id}")
                if random_norm_diff > float(config["representations"]["random_projection"]["row_norm_tolerance"]):
                    raise RuntimeError(f"Random row norm mismatch for {pathway_id}/draw {draw_id}")
                group = output_h5.require_group(f"controls/{pathway_id}/d{draw_id:02d}")
                group.create_dataset("permuted", data=permuted, compression="gzip", compression_opts=1)
                group.create_dataset("random", data=random, compression="gzip", compression_opts=1)
                group.create_dataset("permutation_zero_based", data=permutation)
                control_units.append({
                    "pathway_id": pathway_id,
                    "pathway_index": pathway_index,
                    "pathway": pathway["pathway"],
                    "draw_id": draw_id,
                    "gene_count": gene_count,
                    "dimension": dimension,
                    "permuted_seed": permuted_seed,
                    "random_seed": random_seed,
                    "permuted_h5": f"/controls/{pathway_id}/d{draw_id:02d}/permuted",
                    "random_h5": f"/controls/{pathway_id}/d{draw_id:02d}/random",
                    "permuted_sha256": sha256_array(permuted),
                    "random_sha256": sha256_array(random),
                    "permutation_fixed_point_count": fixed_points,
                    "permuted_vector_multiset_preserved": bool(np.array_equal(
                        np.sort(permutation), np.arange(gene_count)
                    )),
                    "permuted_sorted_row_norm_max_abs_difference": perm_norm_diff,
                    "random_corresponding_row_norm_max_abs_difference": random_norm_diff,
                    "nonfinite_count": int(np.sum(~np.isfinite(permuted)) + np.sum(~np.isfinite(random))),
                })
    os.replace(temporary_h5, controls_path)
    control_by_unit = {(row["pathway_id"], row["draw_id"]): row for row in control_units}

    execution_experiments = []
    for experiment in selected:
        control = control_by_unit[(experiment["pathway_id"], int(experiment["draw_id"]))]
        # Deliberately omit truth genes and truth-fallback fields. The masking runner
        # receives expression paths and truth-blind control paths only.
        execution_experiments.append({
            "experiment_id": experiment["experiment_id"],
            "draw_id": int(experiment["draw_id"]),
            "pathway_id": experiment["pathway_id"],
            "pathway_index": control["pathway_index"],
            "pathway": experiment["pathway"],
            "source_database": experiment["source_database"],
            "analysis_gene_count": int(experiment["analysis_gene_count"]),
            "perturbation_type": experiment["perturbation_type"],
            "perturbation_strength": float(experiment["perturbation_strength"]),
            "perturbation_seed": int(experiment["perturbation_seed"]),
            "condition_a_h5": experiment["condition_a_h5"],
            "condition_b_h5": experiment["condition_b_h5"],
            "gene_names_h5": experiment["gene_names_h5"],
            "permuted_h5": control["permuted_h5"],
            "random_h5": control["random_h5"],
            "permuted_seed": control["permuted_seed"],
            "random_seed": control["random_seed"],
            "permuted_sha256": control["permuted_sha256"],
            "random_sha256": control["random_sha256"],
            "phase8_experiment_input_sha256": _experiment_hash(
                experiment, control["permuted_sha256"], control["random_sha256"]
            ),
        })

    manifest = {
        "schema_version": 1,
        "phase": 8,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "mean_shift_mechanism_decomposition",
        "control_protocol": CONTROL_PROTOCOL,
        "source_phase7_manifest": str(source["phase7_manifest"].relative_to(PROJECT_ROOT)),
        "source_expression_h5": str(source["phase7_expression_h5"].relative_to(PROJECT_ROOT)),
        "controls_h5": str(controls_path.relative_to(PROJECT_ROOT)),
        "controls_h5_sha256": sha256_file(controls_path),
        "pathways": pathways,
        "control_units": control_units,
        "control_unit_count": len(control_units),
        "experiments": execution_experiments,
        "experiment_count": len(execution_experiments),
        "sum_pathway_gene_count": sum_k,
        "sum_pathway_gene_count_plus_one": sum_k_plus_one,
        "expected_mcm_per_control": expected_per_control,
        "expected_total_mcm": 2 * expected_per_control,
        "masking": config["masking"],
        "anti_leakage": {
            "ground_truth_file_read": False,
            "ground_truth_labels_parsed": False,
            "truth_fields_in_execution_manifest": False,
        },
    }
    manifest_path = processed / artifacts["manifest"]
    write_json_atomic(manifest, manifest_path)

    snapshot_paths = {
        "phase8_config": config_path,
        "phase8_protocol": phase8_protocol,
        **source,
    }
    snapshot = {
        "schema_version": 1,
        "phase": 8,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "note": "Files containing evaluation labels were byte-hashed only; labels were not parsed during control generation.",
        "input_sha256": {
            name: {"path": str(path.relative_to(PROJECT_ROOT)), "sha256": sha256_file(path)}
            for name, path in snapshot_paths.items()
        },
        "generated": {
            "manifest": {"path": str(manifest_path.relative_to(PROJECT_ROOT)), "sha256": sha256_file(manifest_path)},
            "controls_h5": {"path": str(controls_path.relative_to(PROJECT_ROOT)), "sha256": sha256_file(controls_path)},
        },
    }
    write_json_atomic(snapshot, processed / artifacts["input_snapshot"])
    print(
        f"PHASE8_PREP status=PASS experiments={len(execution_experiments)} controls={len(control_units)} "
        f"mcm_per_control={expected_per_control} total_mcm={2 * expected_per_control}"
    )
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path,
                        default=PROJECT_ROOT / "config/phase8_mean_shift_mechanism.yaml")
    arguments = parser.parse_args()
    prepare(arguments.config)
