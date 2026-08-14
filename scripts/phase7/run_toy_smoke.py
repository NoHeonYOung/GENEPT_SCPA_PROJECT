#!/usr/bin/env python3
"""Run the cheap Phase 7 end-to-end toy pipeline with the mock backend only."""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from gene_embedding_project.genept_scpa.phase7.evaluation import (  # noqa: E402
    evaluate_ranking, prompt_order_spearman,
)
from gene_embedding_project.genept_scpa.phase7.llm_backend import MockLLMBackend  # noqa: E402
from gene_embedding_project.genept_scpa.phase7.llm_prompts import (  # noqa: E402
    PROMPT_CONDITIONS, build_llm_request,
)
from gene_embedding_project.genept_scpa.phase7.ranking import (  # noqa: E402
    aggregate_rankings, compute_masking_rows,
)
from gene_embedding_project.genept_scpa.phase7.synthetic_perturbation import (  # noqa: E402
    inject_perturbation,
)


def toy_raw_p(condition_a: np.ndarray, condition_b: np.ndarray) -> float:
    """Deterministic distance surrogate for wiring tests; this is not SCPA."""

    distance = float(np.linalg.norm(np.mean(condition_b, axis=0) - np.mean(condition_a, axis=0)))
    return float(math.exp(-distance))


def run_toy_smoke() -> dict[str, object]:
    generator = np.random.default_rng(20260814)
    genes = [f"G{index}" for index in range(1, 7)]
    a = generator.gamma(shape=2.0, scale=0.5, size=(20, len(genes)))
    b = generator.gamma(shape=2.0, scale=0.5, size=(20, len(genes)))
    perturbed = inject_perturbation(
        a, b, genes, scenario="mean_shift", alpha=1.0, seed=20261814
    )
    truth = {row["gene"] for row in perturbed.ground_truth_rows}
    embeddings = generator.normal(size=(len(genes), 8))
    masking = compute_masking_rows(
        perturbed.condition_a, perturbed.condition_b, genes, embeddings, toy_raw_p
    )
    vanilla_order = [
        row["gene"] for row in sorted(masking, key=lambda row: row["vanilla_signed_rank"])
    ]
    genept_order = [
        row["gene"] for row in sorted(masking, key=lambda row: row["genept_signed_rank"])
    ]
    descriptions = {
        gene: f"Gene Symbol {gene} toy metabolic description number {index}"
        for index, gene in enumerate(genes, start=1)
    }
    backend = MockLLMBackend()
    llm_results: dict[str, object] = {}
    for condition_index, prompt_condition in enumerate(PROMPT_CONDITIONS, start=1):
        responses = []
        mapping: dict[str, str] | None = None
        for repeat in range(1, 4):
            bundle = build_llm_request(
                experiment_id="toy_exp", run_id=f"toy__{prompt_condition}__{repeat}",
                pathway="TOY_PATHWAY", source_database="TOY", genes=genes,
                condition_a=perturbed.condition_a, condition_b=perturbed.condition_b,
                descriptions=descriptions, prompt_condition=prompt_condition,
                candidate_order_seed=20280814 + repeat,
                description_shuffle_seed=20270814 + condition_index,
                backend=backend.name,
            )
            mapping = bundle.candidate_to_gene
            responses.append(backend.rank(bundle.request))
        assert mapping is not None
        aggregate = aggregate_rankings(responses, tie_seed=20290814)
        gene_order = [mapping[row["candidate_id"]] for row in aggregate]
        llm_results[prompt_condition] = {
            "metrics": evaluate_ranking(gene_order, truth),
            "stability": prompt_order_spearman(responses),
            "individual_run_count": len(responses),
        }
    return {
        "status": "PASS",
        "backend": backend.name,
        "scientific_evaluation_allowed": backend.scientific_evaluation_allowed,
        "raw_p_backend": "toy_mean_distance_surrogate_not_scpa",
        "truth_gene_count": len(truth),
        "vanilla_metrics": evaluate_ranking(vanilla_order, truth),
        "genept_metrics": evaluate_ranking(genept_order, truth),
        "llm": llm_results,
        "real_llm_inference": False,
        "production_scpa": False,
        "real_dataset_used": False,
    }


if __name__ == "__main__":
    print(json.dumps(run_toy_smoke(), indent=2, ensure_ascii=False))
