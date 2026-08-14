#!/usr/bin/env python3
"""Build leakage-controlled Phase 7 LLM requests without reading ground truth."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import h5py
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from gene_embedding_project.genept_scpa.io import write_json_atomic  # noqa: E402
from gene_embedding_project.genept_scpa.phase7.gpt_oss_backend import (  # noqa: E402
    TransformersGPTOSSBackend,
)
from gene_embedding_project.genept_scpa.phase7.llm_backend import MockLLMBackend  # noqa: E402
from gene_embedding_project.genept_scpa.phase7.llm_prompts import (  # noqa: E402
    PROMPT_CONDITIONS, build_llm_request,
)


def build_requests(
    manifest_path: Path, config_path: Path, output_dir: Path,
    *, backend_name: str = MockLLMBackend.name,
) -> int:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    descriptions = json.loads((PROJECT_ROOT / config["source"]["descriptions"]).read_text(encoding="utf-8"))
    h5_path = Path(manifest["expression_h5"])
    output_dir.mkdir(parents=True, exist_ok=True)
    mapping_dir = output_dir.parent / "candidate_mappings"
    mapping_dir.mkdir(parents=True, exist_ok=True)
    order_base = int(config["llm"]["prompt_order_seed_base"])
    shuffle_base = int(config["llm"]["description_shuffle_seed_base"])
    repeat_count = int(config["llm"]["prompt_order_repeats"])
    count = 0
    with h5py.File(h5_path, "r") as handle:
        for experiment_index, experiment in enumerate(manifest["experiments"], start=1):
            group = handle[f"experiments/{experiment['experiment_id']}"]
            genes = [value.decode() if isinstance(value, bytes) else str(value) for value in group["gene_names"][:]]
            a = group["condition_A/expression"][:]
            b = group["condition_B/expression"][:]
            mapping_saved = False
            for condition_index, prompt_condition in enumerate(PROMPT_CONDITIONS, start=1):
                for repeat in range(1, repeat_count + 1):
                    run_id = f"{experiment['experiment_id']}__{prompt_condition}__order{repeat:02d}"
                    bundle = build_llm_request(
                        experiment_id=experiment["experiment_id"], run_id=run_id,
                        pathway=experiment["pathway"], source_database=experiment["source_database"],
                        genes=genes, condition_a=a, condition_b=b, descriptions=descriptions,
                        prompt_condition=prompt_condition,
                        candidate_order_seed=order_base + experiment_index * 100 + repeat,
                        description_shuffle_seed=shuffle_base + experiment_index * 100 + condition_index,
                        backend=backend_name,
                    )
                    write_json_atomic(bundle.request, output_dir / f"{run_id}.json")
                    if not mapping_saved:
                        write_json_atomic(
                            {
                                "experiment_id": experiment["experiment_id"],
                                "candidate_to_gene": bundle.candidate_to_gene,
                                "llm_visible": False,
                            },
                            mapping_dir / f"{experiment['experiment_id']}.json",
                        )
                        mapping_saved = True
                    count += 1
    return count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument(
        "--config", type=Path,
        default=PROJECT_ROOT / "config/phase7_gpt_oss_synthetic.yaml",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--backend",
        choices=(MockLLMBackend.name, TransformersGPTOSSBackend.name),
        default=MockLLMBackend.name,
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    total = build_requests(
        args.manifest, args.config, args.output_dir, backend_name=args.backend
    )
    print(f"PHASE7_LLM_REQUESTS status=PASS backend={args.backend} requests={total}")
