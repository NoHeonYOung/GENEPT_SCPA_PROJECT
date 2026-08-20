import math
import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from gene_embedding_project.genept_scpa.phase7.synthetic_benchmark_llmfree.cohort import (
    split_pseudo_conditions,
)
from gene_embedding_project.genept_scpa.phase7.synthetic_benchmark_llmfree.metrics import (
    average_precision,
    exact_random_chance,
    ndcg_at_k,
    recall_at_k,
)
from gene_embedding_project.genept_scpa.phase7.synthetic_benchmark_llmfree.perturbation import (
    ground_truth_gene_count,
    inject_perturbation,
)
from scripts.phase7.evaluate_llmfree_benchmark import rank_biserial
from scripts.phase7.prepare_llmfree_benchmark import experiment_specs, frozen_pathways


ROOT = Path(__file__).resolve().parents[1]


def test_hand_calculated_ranking_metrics():
    ranking = ["A", "B", "C", "D", "E"]
    truth = {"A", "C"}
    assert recall_at_k(ranking, truth, 1) == pytest.approx(0.5)
    assert recall_at_k(ranking, truth, 3) == pytest.approx(1.0)
    assert average_precision(ranking, truth) == pytest.approx((1.0 + 2.0 / 3.0) / 2.0)
    expected_ndcg = (1.0 + 1.0 / math.log2(4)) / (1.0 + 1.0 / math.log2(3))
    assert ndcg_at_k(ranking, truth, 3) == pytest.approx(expected_ndcg)


def test_metrics_reject_invalid_truth_and_ranking():
    with pytest.raises(ValueError):
        average_precision(["A", "A"], {"A"})
    with pytest.raises(ValueError):
        recall_at_k(["A", "B"], set(), 1)
    with pytest.raises(ValueError):
        ndcg_at_k(["A", "B"], {"C"}, 1)


def test_exact_random_chance_matches_simple_cases():
    chance = exact_random_chance(4, 1, 2)
    assert chance["recall"] == pytest.approx(0.5)
    assert chance["average_precision"] == pytest.approx(
        (1 + 1 / 2 + 1 / 3 + 1 / 4) / 4
    )
    assert chance["ndcg"] == pytest.approx(
        0.25 * (1 + 1 / math.log2(3))
    )


def test_cohort_split_is_disjoint_deterministic_and_seeded():
    cells = [f"cell_{index:04d}" for index in range(1200)]
    first = split_pseudo_conditions(cells, cells_per_condition=500, seed=20260810)
    second = split_pseudo_conditions(reversed(cells), cells_per_condition=500, seed=20260810)
    assert first == second
    assert len(first.condition_a_ids) == len(first.condition_b_ids) == 500
    assert not set(first.condition_a_ids) & set(first.condition_b_ids)


def _toy_expression():
    rng = np.random.default_rng(7)
    a = rng.gamma(2.0, 1.0, size=(50, 20))
    b = rng.gamma(2.0, 1.0, size=(50, 20))
    return a, b, [f"G{index:02d}" for index in range(20)]


def test_null_uses_uninjected_evaluation_targets():
    a, b, genes = _toy_expression()
    result = inject_perturbation(a, b, genes, scenario="null", alpha=0.0, seed=20260901)
    np.testing.assert_array_equal(result.condition_a, a)
    np.testing.assert_array_equal(result.condition_b, b)
    assert len(result.truth_rows) == ground_truth_gene_count(len(genes))
    assert all(row["is_evaluation_target"] for row in result.truth_rows)
    assert not any(row["is_ground_truth_perturbed"] for row in result.truth_rows)


def test_truth_and_subset_structure_are_reused_across_strengths():
    a, b, genes = _toy_expression()
    weak = inject_perturbation(a, b, genes, scenario="cell_subset", alpha=0.5, seed=20262001)
    strong = inject_perturbation(a, b, genes, scenario="cell_subset", alpha=1.0, seed=20262001)
    assert [row["gene"] for row in weak.truth_rows] == [row["gene"] for row in strong.truth_rows]
    assert weak.target_cell_indices == strong.target_cell_indices
    assert len(weak.target_cell_indices) == 15
    np.testing.assert_array_equal(weak.condition_a, strong.condition_a)
    assert np.all(strong.condition_b >= 0)


def test_mixed_direction_logs_pre_result_negative_pool_fallback():
    a = np.zeros((50, 20)); b = np.zeros((50, 20))
    for gene_index in range(20):
        a[:10, gene_index] = 1.0 + gene_index / 100
        b[10:20, gene_index] = 1.0 + gene_index / 100
    result = inject_perturbation(
        a, b, [f"G{i:02d}" for i in range(20)],
        scenario="mixed_direction", alpha=0.5, seed=44,
    )
    negative = [row for row in result.truth_rows if row["perturbation_direction"] == "negative"]
    assert negative
    assert all(row["direction_selection_rule"] == "condition_b_detected_fallback"
               for row in negative)
    assert all(row["clipped_cell_count"] >= 0 for row in negative)


def test_frozen_pathways_and_workload_match_protocol():
    config = yaml.safe_load((ROOT / "config/phase7_llmfree_synthetic.yaml").read_text())
    manifest = json.loads((ROOT / config["source"]["phase4_manifest"]).read_text())
    pathways = frozen_pathways(config, manifest)
    assert len(pathways) == 11
    assert sum(row["source_database"] == "KEGG" for row in pathways) == 6
    assert sum(row["source_database"] == "REACTOME" for row in pathways) == 5
    assert len(list(experiment_specs(config))) == 7
    expected_experiments = len(pathways) * config["ground_truth"]["draw_count"] * 7
    expected_mcm = config["ground_truth"]["draw_count"] * 14 * sum(
        row["analysis_gene_count"] + 1 for row in pathways
    )
    assert expected_experiments == 1540
    assert expected_mcm == 101920


def test_rank_biserial_direction_is_genept_minus_vanilla():
    assert rank_biserial(np.asarray([1.0, 2.0, 3.0])) == pytest.approx(1.0)
    assert rank_biserial(np.asarray([-1.0, -2.0, -3.0])) == pytest.approx(-1.0)
    assert rank_biserial(np.asarray([0.0, 0.0])) == pytest.approx(0.0)
