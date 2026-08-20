#!/usr/bin/env python3
"""Run Phase 6 True/Permuted/Random semantic-specificity controls."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import time
from typing import Any, Iterable, Sequence

import h5py
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from scripts.genept.build_genept_w import load_and_validate_export, require_export, write_lines_atomic  # noqa: E402
from scripts.phase3.run_cd4_cd8_benchmark import rows_for_ids  # noqa: E402
from scripts.phase4.run_pathway_comparison import read_json, write_csv_atomic, write_markdown_atomic  # noqa: E402
from scripts.phase4.run_timecourse_validation import (  # noqa: E402
    PRIMARY_COMPARISONS,
    choose_cells,
    create_source_audit,
    prepare_canonical_sampling,
    prepare_core_inputs,
    read_metadata,
)
from scripts.phase5.run_gene_contribution import (  # noqa: E402
    PHASE4_DIR,
    OUTPUT_DIR as PHASE5_DIR,
    phase4b_gate,
    representative_targets,
    select_targets,
)
from gene_embedding_project.genept_scpa.config import load_config  # noqa: E402
from gene_embedding_project.genept_scpa.genept_projection import normalize_log1p_sparse  # noqa: E402
from gene_embedding_project.genept_scpa.io import sha256_file, write_json_atomic  # noqa: E402


OUTPUT_DIR = PROJECT_ROOT / "data/processed/genept_scpa/phase6_semantic_controls"
CHECKPOINT_DIR = PROJECT_ROOT / "data/interim/genept_scpa/phase6_semantic_controls_checkpoint"
PREP_DIR = PROJECT_ROOT / "data/interim/genept_scpa/phase6_semantic_controls_prepared"
PATHWAY_REPS = 100
GENE_REPS = 20
RESAMPLING_REPS = 10
TARGET_COUNT = 30
GENE_TARGET_COUNT = 6
PATHWAY_COUNT = 123
RAW_P_CLIP = 1e-300
TOTAL_PATHWAY_CONTROL_MCM = TARGET_COUNT * 2 * PATHWAY_REPS
TOTAL_GENE_MASK_CONTROL_MCM = 11960
TOTAL_GENE_CONTROL_BASELINE_MCM = GENE_TARGET_COUNT * 2 * GENE_REPS
TOTAL_ROBUSTNESS_MCM = RESAMPLING_REPS * TARGET_COUNT * 3
PHASE5_SECONDS_PER_MCM = 95.2 * 60 / 2270


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--smoke-test", action="store_true")
    parser.add_argument("--cores", type=int, default=min(6, os.cpu_count() or 1))
    parser.add_argument(
        "--resume-from", choices=("pathway", "gene", "robustness"), default="pathway",
        help="Resume production at this stage; earlier stage QC/checkpoints must already be complete.",
    )
    parser.add_argument(
        "--progress-seconds", type=int, default=30,
        help="Seconds between aggregate task/MCM/ETA progress reports (minimum 5).",
    )
    return parser.parse_args(argv)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def phase5_gate() -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], dict[str, str]]:
    qc_path = PHASE5_DIR / "phase5_gene_contribution_qc.json"
    manifest_path = PHASE5_DIR / "phase5_gene_contribution_manifest.json"
    all_path = PHASE5_DIR / "phase5_gene_masking_all_results.csv"
    summary_path = PHASE5_DIR / "phase5_pathway_summary.csv"
    qc = read_json(qc_path)
    if qc.get("gate", {}).get("status") != "READY_FOR_GPT_REVIEW":
        raise RuntimeError("Phase 5 is not READY_FOR_GPT_REVIEW")
    if qc["gate"].get("failed_checks") or qc["gate"].get("warnings") or not all(qc["gate"]["criteria"].values()):
        raise RuntimeError("Phase 5 has failed checks or warnings")
    rows: list[dict[str, Any]] = []
    text_fields = {"comparison", "group_a", "group_b", "pathway", "detection_state", "gene"}
    for source in read_csv(all_path):
        rows.append({
            key: (value if key in text_fields else parse_value(value))
            for key, value in source.items()
        })
    summary = [{key: value for key, value in row.items()} for row in read_csv(summary_path)]
    hashes = {name: sha256_file(path) for name, path in {
        "qc": qc_path, "manifest": manifest_path, "all_results": all_path, "summary": summary_path,
        "targets": PHASE5_DIR / "phase5_target_pathways.csv",
    }.items()}
    return qc, rows, summary, hashes


def frozen_targets() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    _, phase4_rows, _ = phase4b_gate()
    targets = select_targets(phase4_rows)
    representatives = representative_targets(targets)
    lookup = {(row["comparison"], row["pathway"]): row for row in targets}
    gene_targets = [dict(lookup[(row["comparison"], row["pathway"])]) for row in representatives]
    expected = [
        ("cd4_0h_vs_12h", "KEGG_NITROGEN_METABOLISM"),
        ("cd4_0h_vs_12h", "REACTOME_INOSITOL_PHOSPHATE_METABOLISM"),
        ("cd4_12h_vs_24h", "KEGG_DRUG_METABOLISM_CYTOCHROME_P450"),
        ("cd4_12h_vs_24h", "REACTOME_PHOSPHOLIPID_METABOLISM"),
        ("cd4_0h_vs_24h", "REACTOME_SYNTHESIS_OF_VERY_LONG_CHAIN_FATTY_ACYL_COAS"),
        ("cd4_0h_vs_24h", "REACTOME_FATTY_ACYL_COA_BIOSYNTHESIS"),
    ]
    if [(row["comparison"], row["pathway"]) for row in gene_targets] != expected:
        raise RuntimeError("Frozen Phase 6 gene-control representative set changed")
    if sum(int(row["n_paired_genes"]) for row in gene_targets) != 299:
        raise RuntimeError("Expected 299 representative gene instances")
    return targets, gene_targets


def target_csv_rows(targets: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    fields = ("comparison", "group_a", "group_b", "pathway", "detection_state", "n_paired_genes",
              "vanilla_raw_p", "vanilla_adjusted_p", "genept_raw_p", "genept_adjusted_p")
    return [{field: row[field] for field in fields} for row in targets]


def build_manifest(base: dict[str, Any], targets: list[dict[str, Any]], gene_targets: list[dict[str, Any]],
                   phase5_hashes: dict[str, str], cores: int) -> dict[str, Any]:
    base["phase"] = "Phase 6 semantic-specific True/Permuted/Random controls"
    base["phase6"] = {
        "question": "Are correct gene-to-GenePT correspondences unusual relative to correspondence-destroyed and semantic-free controls?",
        "pathway_targets": targets,
        "gene_targets": gene_targets,
        "eligible_pathway_count": PATHWAY_COUNT,
        "raw_p_clip": RAW_P_CLIP,
        "representations": {
            "true": "exact official text-embedding-ada-002 rows assigned to their genes",
            "permuted": "within-pathway row-assignment permutation; vector multiset unchanged; changed fraction >0.9",
            "random": "independent 1536D Gaussian direction normalized and scaled to the corresponding True-row norm",
        },
        "same_control_for_both_groups_full_and_masks": True,
        "l2_normalization": False,
        "seeds": {
            "base": 20260812, "permuted": 20260812, "random": 20270812,
            "formula": "control_base + 1000 * one_based_target_index + replicate_id; robustness adds 200000",
        },
        "replicates": {"pathway_per_control": PATHWAY_REPS, "gene_per_control": GENE_REPS, "resampling": RESAMPLING_REPS},
        "resampling": {
            "sampling_seed_formula": "20280812 + 1000 * replicate_id + hour_index",
            "without_replacement": True, "cells_per_timepoint": 500,
            "same_cells_true_permuted_random_within_replicate": True,
        },
        "runtime": {
            "cores": cores, "phase5_seconds_per_mcm": PHASE5_SECONDS_PER_MCM,
            "single_core_estimate_hours": (TOTAL_PATHWAY_CONTROL_MCM + TOTAL_GENE_MASK_CONTROL_MCM + TOTAL_GENE_CONTROL_BASELINE_MCM + TOTAL_ROBUSTNESS_MCM) * PHASE5_SECONDS_PER_MCM / 3600,
            "production_threshold_hours": 24,
        },
        "expected_mcm": {
            "pathway_controls": TOTAL_PATHWAY_CONTROL_MCM,
            "gene_masks": TOTAL_GENE_MASK_CONTROL_MCM,
            "gene_control_baselines": TOTAL_GENE_CONTROL_BASELINE_MCM,
            "robustness": TOTAL_ROBUSTNESS_MCM,
            "total_including_control_baselines": TOTAL_PATHWAY_CONTROL_MCM + TOTAL_GENE_MASK_CONTROL_MCM + TOTAL_GENE_CONTROL_BASELINE_MCM + TOTAL_ROBUSTNESS_MCM,
        },
        "phase5_artifact_sha256": phase5_hashes,
        "scope": {"phase7_run": False, "cd8_generalization_run": False, "classifier_run": False,
                  "genept_l2_gene_masking_run": False, "external_network_run": False, "new_pathway_database_run": False},
    }
    return base


def prepare_resampling(h5_path: Path, manifest: dict[str, Any]) -> None:
    interim = PROJECT_ROOT / "data/interim/genept_scpa"
    metadata = read_metadata(interim / "phase2_export/naive_cd4/naive_cd4_metadata.csv")
    _, counts, genes, cell_ids = load_and_validate_export(require_export("naive_cd4", force=False))
    gene_index = {gene: index for index, gene in enumerate(genes)}
    global_genes = list(manifest["global_gene_order"])
    columns = np.asarray([gene_index[gene] for gene in global_genes], dtype=np.int64)
    cell_dir = OUTPUT_DIR / "robustness" / "cells"
    cell_dir.mkdir(parents=True, exist_ok=True)
    records: dict[str, Any] = {}
    with h5py.File(h5_path, "a") as handle:
        for replicate_id in range(1, RESAMPLING_REPS + 1):
            records[str(replicate_id)] = {}
            for hour_index, hour in enumerate(("0h", "12h", "24h"), start=1):
                seed = 20280812 + 1000 * replicate_id + hour_index
                selected = choose_cells(metadata, hour, sample_size=500, seed=seed)
                destination = cell_dir / f"rep_{replicate_id:02d}_cd4_{hour}_cells.txt"
                write_lines_atomic(selected, destination)
                rows = rows_for_ids(cell_ids, selected)
                normalized = normalize_log1p_sparse(counts[rows], normalization_target=10_000.0)
                matrix = normalized[:, columns].toarray().astype(np.float64)
                handle.create_dataset(
                    f"resampling/rep_{replicate_id:02d}/expression/cd4_{hour}",
                    data=matrix, compression="gzip", shuffle=True,
                )
                records[str(replicate_id)][f"cd4_{hour}"] = {
                    "seed": seed, "cell_count": 500, "cell_file": str(destination), "sha256": sha256_file(destination),
                }
    manifest["phase6"]["resampling"]["samples"] = records


def _csv_data_row_count(path: Path) -> int:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def _expected_work(manifest: dict[str, Any], mode: str, max_targets: int, max_reps: int,
                   max_genes: int) -> tuple[int, int]:
    if mode == "pathway":
        targets = manifest["phase6"]["pathway_targets"][:max_targets or None]
        reps = min(PATHWAY_REPS, max_reps) if max_reps else PATHWAY_REPS
        return len(targets) * 2 * reps, len(targets) * 2 * reps
    if mode == "gene":
        targets = manifest["phase6"]["gene_targets"][:max_targets or None]
        reps = min(GENE_REPS, max_reps) if max_reps else GENE_REPS
        gene_counts = [min(int(target["n_paired_genes"]), max_genes) if max_genes else int(target["n_paired_genes"]) for target in targets]
        return len(targets) * 2 * reps, sum((count + 1) * 2 * reps for count in gene_counts)
    if mode == "robustness":
        targets = manifest["phase6"]["pathway_targets"][:max_targets or None]
        reps = min(RESAMPLING_REPS, max_reps) if max_reps else RESAMPLING_REPS
        return len(targets) * 3 * reps, len(targets) * 3 * reps
    return 0, 0


def _progress_snapshot(checkpoint_dir: Path, mode: str) -> tuple[int, int, list[str]]:
    checkpoints = sorted(checkpoint_dir.rglob("*.csv")) if checkpoint_dir.exists() else []
    completed_tasks = len(checkpoints)
    completed_mcm = sum((_csv_data_row_count(path) + 1) if mode == "gene" else 1 for path in checkpoints)
    active: list[str] = []
    for marker in sorted(checkpoint_dir.rglob("*.progress.json")) if checkpoint_dir.exists() else []:
        checkpoint = Path(str(marker).removesuffix(".progress.json"))
        if checkpoint.exists():
            continue
        try:
            item = read_json(marker)
        except (OSError, ValueError, json.JSONDecodeError):
            continue
        completed_mcm += int(item["completed_mcm"])
        active.append(
            f"T{int(item['target_index']):02d}/{item['control']}/R{int(item['replicate_id']):02d}"
            f" gene={item['current_gene']} ({int(item['completed_mcm'])}/{int(item['total_mcm'])} MCM)"
        )
    return completed_tasks, completed_mcm, active


def _duration(seconds: float | None) -> str:
    if seconds is None or not math.isfinite(seconds):
        return "calculating"
    seconds = max(0, int(seconds))
    hours, remainder = divmod(seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def run_r(h5_path: Path, manifest_path: Path, checkpoint_dir: Path, qc_path: Path, *,
          mode: str, cores: int, max_targets: int = 0, max_reps: int = 0, max_genes: int = 0,
          progress_seconds: int = 30) -> dict[str, Any]:
    command = [
        "Rscript", str(PROJECT_ROOT / "scripts/scpa/run_phase6_semantic_controls_core.R"),
        "--input-h5", str(h5_path), "--manifest", str(manifest_path),
        "--checkpoint-dir", str(checkpoint_dir), "--output-json", str(qc_path),
        "--mode", mode, "--cores", str(cores),
    ]
    for flag, value in (("--max-targets", max_targets), ("--max-reps", max_reps), ("--max-genes", max_genes)):
        if value:
            command += [flag, str(value)]
    started = time.monotonic()
    print(f"[Phase 6] R core mode={mode} cores={cores} START", flush=True)
    if mode == "preflight":
        subprocess.run(command, cwd=PROJECT_ROOT, check=True)
    else:
        manifest = read_json(manifest_path)
        total_tasks, total_mcm = _expected_work(manifest, mode, max_targets, max_reps, max_genes)
        initial_tasks, initial_mcm, _ = _progress_snapshot(checkpoint_dir, mode)
        print(
            f"[Phase 6 progress] mode={mode} RESUME completed_tasks={initial_tasks}/{total_tasks} "
            f"completed_mcm={initial_mcm}/{total_mcm} remaining_mcm={max(total_mcm - initial_mcm, 0)}",
            flush=True,
        )
        process = subprocess.Popen(command, cwd=PROJECT_ROOT)
        interval = max(5, int(progress_seconds))
        try:
            while process.poll() is None:
                time.sleep(interval)
                completed_tasks, completed_mcm, active = _progress_snapshot(checkpoint_dir, mode)
                elapsed = time.monotonic() - started
                new_mcm = max(completed_mcm - initial_mcm, 0)
                rate = new_mcm / elapsed if elapsed > 0 else 0.0
                eta = (total_mcm - completed_mcm) / rate if rate > 0 else None
                action = " | ".join(active[:3]) if active else "checkpoint scan / worker startup"
                print(
                    f"[Phase 6 progress] mode={mode} tasks={completed_tasks}/{total_tasks} "
                    f"mcm={completed_mcm}/{total_mcm} remaining={max(total_mcm - completed_mcm, 0)} "
                    f"elapsed={_duration(elapsed)} ETA={_duration(eta)} active={len(active)} "
                    f"doing=[{action}]",
                    flush=True,
                )
            if process.returncode:
                raise subprocess.CalledProcessError(process.returncode, command)
        except KeyboardInterrupt:
            print(
                f"\n[Phase 6] Interrupted by user. Completed CSV checkpoints are safe. "
                f"Run the same command to resume mode={mode}.", flush=True,
            )
            process.send_signal(2)
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.terminate()
            raise
    print(f"[Phase 6] R core mode={mode} DONE elapsed={(time.monotonic() - started) / 60:.1f} min", flush=True)
    return read_json(qc_path)


def parse_value(value: str) -> Any:
    # R writes logical values as TRUE/FALSE and prior Python artifacts use
    # True/False. Keep lowercase representation labels such as "true" as text.
    if value in {"TRUE", "FALSE", "True", "False"}:
        return value in {"TRUE", "True"}
    if value in {"NA", "NaN", ""}:
        return None
    try:
        return float(value)
    except ValueError:
        return value


def collect_files(directory: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(directory.rglob("*.csv")):
        rows.extend({key: parse_value(value) for key, value in row.items()} for row in read_csv(path))
    return rows


def bh_adjust(values: Sequence[float]) -> list[float]:
    order = np.argsort(np.asarray(values, dtype=float))
    adjusted = np.empty(len(values), dtype=float)
    running = 1.0
    for reverse_rank, index in enumerate(order[::-1], start=1):
        rank = len(values) - reverse_rank + 1
        running = min(running, values[int(index)] * len(values) / rank)
        adjusted[int(index)] = min(running, 1.0)
    return adjusted.tolist()


def summarize_pathway_controls(rows: Sequence[dict[str, Any]], targets: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    for control in ("permuted", "random"):
        control_rows: list[dict[str, Any]] = []
        for target_index, target in enumerate(targets, start=1):
            subset = [row for row in rows if int(row["target_index"]) == target_index and row["control"] == control]
            if len(subset) != PATHWAY_REPS:
                raise RuntimeError(f"Incomplete Phase 6A target={target_index} control={control}: {len(subset)}")
            scores = np.asarray([row["score"] for row in subset], dtype=float)
            true_score = -math.log10(max(float(target["genept_raw_p"]), RAW_P_CLIP))
            greater_equal = int(np.sum(scores >= true_score))
            less_equal = int(np.sum(scores <= true_score))
            control_rows.append({
                "target_index": target_index, "comparison": target["comparison"], "pathway": target["pathway"],
                "detection_state": target["detection_state"], "control": control,
                "true_raw_p": target["genept_raw_p"], "true_score": true_score,
                "control_median_score": float(np.median(scores)),
                "true_minus_control_median_score": true_score - float(np.median(scores)),
                "true_percentile_among_control": float((1 + less_equal) / (len(scores) + 1)),
                "fraction_control_ge_true": float(greater_equal / len(scores)),
                "fraction_control_significant": float(np.mean([row["significant"] for row in subset])),
                "empirical_upper_p": float((1 + greater_equal) / (len(scores) + 1)),
                "empirical_lower_p": float((1 + less_equal) / (len(scores) + 1)),
            })
        upper_bh = bh_adjust([row["empirical_upper_p"] for row in control_rows])
        lower_bh = bh_adjust([row["empirical_lower_p"] for row in control_rows])
        for row, upper, lower in zip(control_rows, upper_bh, lower_bh):
            row["empirical_upper_p_bh30"] = upper
            row["empirical_lower_p_bh30"] = lower
        summaries.extend(control_rows)
    return sorted(summaries, key=lambda row: (row["target_index"], row["control"]))


def rank_desc(values: Sequence[float]) -> np.ndarray:
    values_array = np.asarray(values, dtype=float)
    order = np.argsort(-values_array, kind="mergesort")
    ranks = np.empty(len(values_array), dtype=float)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values_array[order[end]] == values_array[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + 1 + end) / 2
        start = end
    return ranks


def spearman(a: Sequence[float], b: Sequence[float]) -> float | None:
    ra, rb = rank_desc(a), rank_desc(b)
    if np.std(ra) == 0 or np.std(rb) == 0:
        return None
    return float(np.corrcoef(ra, rb)[0, 1])


def top_genes(rows: Sequence[dict[str, Any]], field: str, k: int) -> set[str]:
    return {row["gene"] for row in sorted(rows, key=lambda row: (-abs(float(row[field])), row["gene"]))[:min(k, len(rows))]}


def summarize_gene_controls(control_rows: list[dict[str, Any]], phase5_rows: Sequence[dict[str, Any]],
                            gene_targets: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    phase5_lookup = {(row["comparison"], row["pathway"], row["gene"]): row for row in phase5_rows}
    for row in control_rows:
        true = phase5_lookup[(row["comparison"], row["pathway"], row["gene"])]
        row["vanilla_delta_score"] = true["vanilla_delta_score"]
        row["true_genept_delta_score"] = true["genept_delta_score"]
        row["true_genept_detection_flip"] = true["genept_detection_flip"]
    for target_index, target in enumerate(gene_targets, start=1):
        true_subset = [row for row in phase5_rows if row["comparison"] == target["comparison"] and row["pathway"] == target["pathway"]]
        vanilla_d = [row["vanilla_delta_score"] for row in true_subset]
        true_d = [row["genept_delta_score"] for row in true_subset]
        true_metrics = {
            "true_signed_spearman_vs_vanilla": spearman(vanilla_d, true_d),
            "true_absolute_spearman_vs_vanilla": spearman([abs(x) for x in vanilla_d], [abs(x) for x in true_d]),
            "true_top5_overlap_vs_vanilla": len(top_genes(true_subset, "vanilla_delta_score", 5) & top_genes(true_subset, "genept_delta_score", 5)),
            "true_top10_overlap_vs_vanilla": len(top_genes(true_subset, "vanilla_delta_score", 10) & top_genes(true_subset, "genept_delta_score", 10)),
            "true_detection_flip_count": sum(bool(row["genept_detection_flip"]) for row in true_subset),
        }
        for control in ("permuted", "random"):
            replicate_metrics = []
            for replicate_id in range(1, GENE_REPS + 1):
                subset = [row for row in control_rows if int(row["target_index"]) == target_index and row["control"] == control and int(row["replicate_id"]) == replicate_id]
                if len(subset) != len(true_subset):
                    raise RuntimeError(f"Incomplete Phase 6B target={target_index} control={control} rep={replicate_id}")
                control_d = [row["control_delta_score"] for row in subset]
                replicate_metrics.append({
                    "signed_spearman_vs_vanilla": spearman(vanilla_d, control_d),
                    "absolute_spearman_vs_vanilla": spearman([abs(x) for x in vanilla_d], [abs(x) for x in control_d]),
                    "top5_overlap_vs_vanilla": len(top_genes(subset, "vanilla_delta_score", 5) & top_genes(subset, "control_delta_score", 5)),
                    "top10_overlap_vs_vanilla": len(top_genes(subset, "vanilla_delta_score", 10) & top_genes(subset, "control_delta_score", 10)),
                    "detection_flip_count": sum(bool(row["control_detection_flip"]) for row in subset),
                })
            summary = {"target_index": target_index, "comparison": target["comparison"], "pathway": target["pathway"],
                       "detection_state": target["detection_state"], "n_genes": len(true_subset), "control": control, **true_metrics}
            for metric in replicate_metrics[0]:
                values = [row[metric] for row in replicate_metrics if row[metric] is not None]
                summary[f"control_median_{metric}"] = float(np.median(values)) if values else None
                true_key = "true_" + metric
                if true_key in true_metrics and true_metrics[true_key] is not None and values:
                    summary[f"true_percentile_{metric}"] = float((1 + sum(value <= true_metrics[true_key] for value in values)) / (len(values) + 1))
            summaries.append(summary)
    return summaries


def augment_robustness(rows: list[dict[str, Any]], targets: Sequence[dict[str, Any]]) -> dict[str, Any]:
    canonical = np.asarray([-math.log10(max(float(row["genept_raw_p"]), RAW_P_CLIP)) for row in targets])
    sign_records: dict[str, list[bool]] = {"permuted": [], "random": []}
    rank_stability: list[float] = []
    for rep in range(1, RESAMPLING_REPS + 1):
        true_rows = sorted([row for row in rows if int(row["resample_id"]) == rep and row["representation"] == "true"], key=lambda row: int(row["target_index"]))
        true_scores = np.asarray([row["score"] for row in true_rows])
        rank_stability.append(float(spearman(canonical, true_scores)))
        true_lookup = {int(row["target_index"]): row for row in true_rows}
        for row in rows:
            if int(row["resample_id"]) != rep:
                continue
            row["canonical_true_score"] = canonical[int(row["target_index"]) - 1]
            row["resampled_true_score"] = true_lookup[int(row["target_index"])]["score"]
            if row["representation"] != "true":
                row["true_minus_control_score"] = row["resampled_true_score"] - row["score"]
                canonical_difference = canonical[int(row["target_index"]) - 1] - next(
                    summary["control_median_score"] for summary in _PATHWAY_SUMMARY_CACHE
                    if int(summary["target_index"]) == int(row["target_index"]) and summary["control"] == row["representation"]
                )
                sign_records[row["representation"]].append(np.sign(row["true_minus_control_score"]) == np.sign(canonical_difference))
            else:
                row["true_minus_control_score"] = 0.0
    return {
        "median_true_rank_stability": float(np.median(rank_stability)),
        "true_rank_stability_by_resample": rank_stability,
        "true_minus_permuted_sign_consistency": float(np.mean(sign_records["permuted"])),
        "true_minus_random_sign_consistency": float(np.mean(sign_records["random"])),
    }


_PATHWAY_SUMMARY_CACHE: list[dict[str, Any]] = []


def create_figures(pathway_summary: Sequence[dict[str, Any]], gene_rows: Sequence[dict[str, Any]],
                   gene_summary: Sequence[dict[str, Any]], robustness_rows: Sequence[dict[str, Any]]) -> list[Path]:
    cache = Path("/tmp/genept_scpa_plot_cache")
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache / "xdg"))
    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415

    figure_dir = OUTPUT_DIR / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    files: list[Path] = []

    def pathway_label(value: str, width: int = 32) -> str:
        return "\n".join(textwrap.wrap(value.replace("_", " "), width=width))

    def comparison_label(value: str) -> str:
        return value.replace("cd4_", "CD4 ").replace("h_vs_", "–").replace("h", " h")

    def save(name: str, fig: Any) -> None:
        path = figure_dir / name
        fig.savefig(path, dpi=180, bbox_inches="tight", pad_inches=0.15)
        plt.close(fig); files.append(path)

    fig, ax = plt.subplots(figsize=(8, 6), layout="constrained")
    for control, marker in (("permuted", "o"), ("random", "s")):
        subset = [row for row in pathway_summary if row["control"] == control]
        ax.scatter([row["control_median_score"] for row in subset], [row["true_score"] for row in subset], label=control, marker=marker, alpha=.75)
    limit = max(ax.get_xlim()[1], ax.get_ylim()[1]); ax.plot([0, limit], [0, limit], "--", color="grey")
    ax.set(xlabel="Control median pathway score", ylabel="True pathway score", title="True vs semantic controls"); ax.legend()
    save("01_true_vs_control_pathway_scores.png", fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
    for ax, control in zip(axes, ("permuted", "random")):
        values = [row["true_percentile_among_control"] for row in pathway_summary if row["control"] == control]
        ax.hist(values, bins=10, edgecolor="black"); ax.set(title=control, xlabel="True percentile", ylabel="Targets")
    save("02_true_percentile_distributions.png", fig)

    fig, ax = plt.subplots(figsize=(8, 9), layout="constrained")
    matrix = np.asarray([[next(row["true_minus_control_median_score"] for row in pathway_summary if int(row["target_index"]) == i and row["control"] == control) for control in ("permuted", "random")] for i in range(1, 31)])
    image = ax.imshow(matrix, aspect="auto", cmap="coolwarm")
    ax.set(xticks=(0, 1), xticklabels=("Permuted", "Random"), yticks=np.arange(30), yticklabels=np.arange(1, 31), ylabel="Target index", title="True minus control median score")
    fig.colorbar(image, ax=ax)
    save("03_control_difference_heatmap.png", fig)

    fig, axes = plt.subplots(2, 3, figsize=(18, 12), layout="constrained")
    for target_index, ax in enumerate(axes.flat, start=1):
        subset = [row for row in gene_rows if int(row["target_index"]) == target_index and row["control"] == "permuted" and int(row["replicate_id"]) == 1]
        ax.scatter([row["vanilla_delta_score"] for row in subset], [row["control_delta_score"] for row in subset], s=15)
        if subset: ax.set(title=f"{comparison_label(subset[0]['comparison'])}\n{pathway_label(subset[0]['pathway'], 34)}")
        ax.set(xlabel="Vanilla delta", ylabel="Permuted delta")
    save("04_gene_rank_control_comparison.png", fig)

    ordered_summary = sorted(gene_summary, key=lambda row: (int(row["target_index"]), row["control"]))
    fig, ax = plt.subplots(figsize=(10, 7), layout="constrained")
    y = np.arange(len(ordered_summary))
    ax.barh(y - .18, [row["true_top10_overlap_vs_vanilla"] for row in ordered_summary], .36, label="True")
    ax.barh(y + .18, [row["control_median_top10_overlap_vs_vanilla"] for row in ordered_summary], .36, label="Control median")
    ax.set(yticks=y, yticklabels=[f"Target {int(row['target_index']):02d} · {row['control'].title()}" for row in ordered_summary], xlabel="Top-10 overlap", title="Gene-ranking overlap with Vanilla")
    ax.invert_yaxis(); ax.legend()
    save("05_top_gene_overlap_controls.png", fig)

    fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
    values = {representation: [row["score"] for row in robustness_rows if row["representation"] == representation] for representation in ("true", "permuted", "random")}
    ax.boxplot([values[key] for key in values], labels=[key.title() for key in values]); ax.set(ylabel="Pathway score", title="Paired resampling robustness")
    save("06_resampling_robustness.png", fig)
    return files


def build_qc(preflight: dict[str, Any], cores_qc: dict[str, dict[str, Any]], pathway_rows: list[dict[str, Any]],
             pathway_summary: list[dict[str, Any]], gene_rows: list[dict[str, Any]], gene_summary: list[dict[str, Any]],
             robustness_rows: list[dict[str, Any]], robustness_metrics: dict[str, Any], manifest: dict[str, Any],
             phase5_hashes: dict[str, str], smoke_resume: bool, figures: Sequence[Path]) -> dict[str, Any]:
    current_hashes = {name: sha256_file(path) for name, path in {
        "qc": PHASE5_DIR / "phase5_gene_contribution_qc.json",
        "manifest": PHASE5_DIR / "phase5_gene_contribution_manifest.json",
        "all_results": PHASE5_DIR / "phase5_gene_masking_all_results.csv",
        "summary": PHASE5_DIR / "phase5_pathway_summary.csv",
        "targets": PHASE5_DIR / "phase5_target_pathways.csv",
    }.items()}
    criteria = {
        "phase5_pass": True,
        "exact_30_pathway_targets": len({(row["comparison"], row["pathway"]) for row in pathway_summary}) == 30,
        "true_baseline_reproduced": preflight["checks"]["true_baseline_reproduced"],
        "permuted_vector_multiset_preserved": preflight["checks"]["permuted_vector_multiset_preserved"],
        "permuted_mapping_changed_gt_0_9": min(row["mapping_changed_fraction"] for row in pathway_rows if row["control"] == "permuted") > .9,
        "random_dimension_1536": preflight["checks"]["random_dimension_1536"],
        "row_norms_preserved": max(row["norm_max_difference"] for row in pathway_rows + gene_rows + robustness_rows) <= 1e-9,
        "seed_determinism": preflight["checks"]["same_seed_deterministic"] and preflight["checks"]["different_seed_changes_control"],
        "same_cells_genes_order": preflight["checks"]["same_cells_and_gene_order"],
        "no_l2": preflight["checks"]["no_l2_normalization"],
        "pathway_control_complete_6000": len(pathway_rows) == TOTAL_PATHWAY_CONTROL_MCM,
        "gene_mask_control_complete_11960": len(gene_rows) == TOTAL_GENE_MASK_CONTROL_MCM,
        "gene_control_baseline_complete_240": cores_qc["gene"]["task_count"] == TOTAL_GENE_CONTROL_BASELINE_MCM,
        "resampling_complete_900": len(robustness_rows) == TOTAL_ROBUSTNESS_MCM,
        "failed_mcm_zero": all(item["failed_mcm_calls"] == 0 for item in cores_qc.values()),
        "runtime_warnings_absent": not any(item["warnings"] for item in cores_qc.values()),
        "checkpoint_resume_pass": smoke_resume and all(item["checkpoint_resume"] for item in cores_qc.values()),
        "historical_phase5_outputs_unmodified": current_hashes == phase5_hashes,
        "scope_respected": not any(manifest["phase6"]["scope"].values()),
    }
    failed = [name for name, passed in criteria.items() if not passed]
    status = "READY_FOR_GPT_REVIEW" if not failed else "NEEDS_REVIEW"
    return {
        "phase": "Phase 6 semantic-specific controls", "status": status,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "mcm_counts": manifest["phase6"]["expected_mcm"], "preflight": preflight,
        "core": cores_qc, "robustness": robustness_metrics,
        "figures": [str(path) for path in figures], "scope": manifest["phase6"]["scope"],
        "phase5_artifact_sha256_before": phase5_hashes, "phase5_artifact_sha256_after": current_hashes,
        "gate": {"status": status, "criteria": criteria, "failed_checks": failed, "warnings": []},
    }


def summary_lines(qc: dict[str, Any], pathway_summary: Sequence[dict[str, Any]], gene_summary: Sequence[dict[str, Any]]) -> list[str]:
    return [
        "# Phase 6 semantic-specific controls", "",
        f"Gate status: `{qc['gate']['status']}`", "",
        "## Frozen method", "",
        "True uses the correct official GenePT gene-row correspondence. Permuted preserves the exact within-pathway vector multiset and row norms while changing more than 90% of gene assignments. Random uses independent 1536-dimensional Gaussian directions scaled to each corresponding True-row norm. The same control realization is used for both groups and for full/masked calculations; no L2 normalization is used.", "",
        "## Completed scope", "",
        f"- Phase 6A pathway-control MCM: {qc['mcm_counts']['pathway_controls']}",
        f"- Phase 6B gene-mask control MCM: {qc['mcm_counts']['gene_masks']}",
        f"- Phase 6B control-baseline MCM: {qc['mcm_counts']['gene_control_baselines']}",
        f"- Paired resampling MCM: {qc['mcm_counts']['robustness']}",
        f"- Pathway target/control summaries: {len(pathway_summary)}",
        f"- Representative target/control summaries: {len(gene_summary)}", "",
        "## Interpretation limit", "",
        "These controls test whether results depend on the correct GenePT gene-to-vector correspondence relative to specified null constructions. They do not establish causal genes, biological correctness, predictive superiority, or generalization to CD8 or other datasets. Phase 7, CD8 generalization, classifiers, GenePT L2 gene masking, external networks and new pathway databases were not run.",
    ]


def run(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.cores < 1:
        raise ValueError("--cores must be positive")
    protocol = load_config(PROJECT_ROOT / "config/genept_scpa.yaml")
    protocol.require_phase(6)
    if protocol.active_phase != 6:
        raise RuntimeError("Phase 6 must be the active phase")
    phase5_qc, phase5_rows, _, phase5_hashes = phase5_gate()
    targets, gene_targets = frozen_targets()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for name in ("pathway_controls", "gene_controls", "robustness", "figures"):
        (OUTPUT_DIR / name).mkdir(parents=True, exist_ok=True)
    write_csv_atomic(target_csv_rows(targets), list(target_csv_rows(targets)[0]), OUTPUT_DIR / "phase6_control_targets.csv")
    write_csv_atomic(target_csv_rows(gene_targets), list(target_csv_rows(gene_targets)[0]), OUTPUT_DIR / "phase6_gene_control_targets.csv")

    estimate = (TOTAL_PATHWAY_CONTROL_MCM + TOTAL_GENE_MASK_CONTROL_MCM + TOTAL_GENE_CONTROL_BASELINE_MCM + TOTAL_ROBUSTNESS_MCM) * PHASE5_SECONDS_PER_MCM / 3600
    print(f"[Phase 6 runtime gate] MCM=19100 single-core estimate={estimate:.2f}h threshold=24h cores={args.cores} decision=RUN", flush=True)
    if estimate >= 24:
        print("READY_FOR_RUNTIME_REVIEW")
        return 2

    audit = create_source_audit()
    if audit["gate"]["status"] != "PASS":
        raise RuntimeError("Source audit failed")
    sampling = prepare_canonical_sampling(groups=("cd4_0h", "cd4_12h", "cd4_24h"))
    PREP_DIR.mkdir(parents=True, exist_ok=True)
    base, h5_path, _ = prepare_core_inputs(dict(protocol.values), PREP_DIR, sampling, PRIMARY_COMPARISONS, comparison_set_name="cd4_activation")
    manifest = build_manifest(base, targets, gene_targets, phase5_hashes, args.cores)
    print("[Phase 6] Preparing 10 paired resampling sets from all available CD4 cells...", flush=True)
    prepare_resampling(h5_path, manifest)
    manifest_path = PREP_DIR / "phase6_semantic_control_manifest.json"
    write_json_atomic(manifest, manifest_path)
    if args.prepare_only:
        print("PHASE6_PREP status=PASS runtime_gate=RUN production_run=false")
        return 0

    preflight_path = PREP_DIR / "phase6_preflight_qc.json"
    preflight = run_r(h5_path, manifest_path, PREP_DIR / "preflight_checkpoint", preflight_path, mode="preflight", cores=1)
    if preflight["status"] != "PASS":
        raise RuntimeError("Phase 6 preflight failed")

    with tempfile.TemporaryDirectory(prefix="genept_scpa_phase6_smoke_") as directory:
        smoke = Path(directory)
        smoke_ok = True
        for mode, max_genes in (("pathway", 0), ("gene", 2), ("robustness", 0)):
            first = run_r(h5_path, manifest_path, smoke / mode, smoke / f"{mode}_first.json", mode=mode, cores=1, max_targets=1, max_reps=1, max_genes=max_genes, progress_seconds=5)
            second = run_r(h5_path, manifest_path, smoke / mode, smoke / f"{mode}_second.json", mode=mode, cores=1, max_targets=1, max_reps=1, max_genes=max_genes, progress_seconds=5)
            smoke_ok &= first["reused_checkpoint_count"] == 0 and second["reused_checkpoint_count"] == first["task_count"]
        if not smoke_ok:
            raise RuntimeError("Phase 6 checkpoint/resume smoke failed")
    if args.smoke_test:
        print("PHASE6_SMOKE status=PASS preflight=PASS checkpoint_resume=PASS production_run=false")
        return 0

    core_qc: dict[str, dict[str, Any]] = {}
    production_modes = ("pathway", "gene", "robustness")
    start_index = production_modes.index(args.resume_from)
    for mode in production_modes[:start_index]:
        prior_qc_path = CHECKPOINT_DIR / f"phase6_{mode}_core_qc.json"
        if not prior_qc_path.is_file():
            raise RuntimeError(
                f"--resume-from {args.resume_from} requires completed earlier QC: {prior_qc_path}"
            )
        core_qc[mode] = read_json(prior_qc_path)
        if core_qc[mode].get("status") != "PASS":
            raise RuntimeError(f"Earlier Phase 6 mode is not PASS: {mode}")
        print(f"[Phase 6 production] mode={mode} already PASS; skipped by --resume-from {args.resume_from}.", flush=True)
    for mode in production_modes[start_index:]:
        print(f"[Phase 6 production] mode={mode} starting with atomic resume checkpoints...", flush=True)
        core_qc[mode] = run_r(
            h5_path, manifest_path, CHECKPOINT_DIR / mode, CHECKPOINT_DIR / f"phase6_{mode}_core_qc.json",
            mode=mode, cores=args.cores, progress_seconds=args.progress_seconds,
        )

    pathway_rows = collect_files(CHECKPOINT_DIR / "pathway")
    pathway_summary = summarize_pathway_controls(pathway_rows, targets)
    global _PATHWAY_SUMMARY_CACHE
    _PATHWAY_SUMMARY_CACHE = pathway_summary
    gene_rows = collect_files(CHECKPOINT_DIR / "gene")
    gene_summary = summarize_gene_controls(gene_rows, phase5_rows, gene_targets)
    robustness_rows = collect_files(CHECKPOINT_DIR / "robustness")
    robustness_metrics = augment_robustness(robustness_rows, targets)

    write_csv_atomic(pathway_rows, list(pathway_rows[0]), OUTPUT_DIR / "phase6_pathway_control_all_results.csv")
    write_csv_atomic(pathway_summary, list(pathway_summary[0]), OUTPUT_DIR / "phase6_pathway_control_summary.csv")
    write_csv_atomic(gene_rows, list(gene_rows[0]), OUTPUT_DIR / "phase6_gene_control_all_results.csv")
    write_csv_atomic(gene_summary, list(gene_summary[0]), OUTPUT_DIR / "phase6_gene_control_summary.csv")
    write_csv_atomic(robustness_rows, list(robustness_rows[0]), OUTPUT_DIR / "phase6_resampling_results.csv")
    write_json_atomic(manifest, OUTPUT_DIR / "phase6_semantic_control_manifest.json")
    figures = create_figures(pathway_summary, gene_rows, gene_summary, robustness_rows)
    qc = build_qc(preflight, core_qc, pathway_rows, pathway_summary, gene_rows, gene_summary,
                  robustness_rows, robustness_metrics, manifest, phase5_hashes, smoke_ok, figures)
    write_json_atomic(qc, OUTPUT_DIR / "phase6_semantic_control_qc.json")
    write_markdown_atomic(summary_lines(qc, pathway_summary, gene_summary), OUTPUT_DIR / "phase6_semantic_control_summary.md")
    print(f"PHASE6 status={qc['gate']['status']} total_mcm=19100 failed={len(qc['gate']['failed_checks'])}")
    return 0 if qc["gate"]["status"] == "READY_FOR_GPT_REVIEW" else 1


if __name__ == "__main__":
    raise SystemExit(run())
