import csv
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from gene_embedding_project.genept_scpa.phase7.evaluation import (
    average_precision,
    evaluate_ranking,
    ndcg_at_k,
    prompt_order_spearman,
    recall_at_k,
)
from gene_embedding_project.genept_scpa.phase7.ranking import (
    _average_descending_ranks,
    aggregate_rankings,
    genept_projection,
    genept_subtraction_mask,
    vanilla_zero_mask,
)
from scripts.phase7.evaluate_rankings import evaluate_scpa_directory


ROOT = Path(__file__).resolve().parents[1]


class Phase7RankingMetricTests(unittest.TestCase):
    def test_phase5_average_tie_rank_rule(self):
        self.assertEqual(_average_descending_ranks([3.0, 1.0, 3.0]), [1.5, 3.0, 1.5])

    def test_hand_calculated_metrics(self):
        ranking = ["A", "B", "C", "D"]
        truth = {"A", "C"}
        self.assertEqual(recall_at_k(ranking, truth, 2), 0.5)
        self.assertAlmostEqual(average_precision(ranking, truth), (1.0 + 2 / 3) / 2)
        expected_ndcg = (1.0 + 1 / math.log2(4)) / (1.0 + 1 / math.log2(3))
        self.assertAlmostEqual(ndcg_at_k(ranking, truth, 4), expected_ndcg)
        metrics = evaluate_ranking(ranking, truth)
        self.assertEqual(metrics["truth_k"], 2)

    def test_phase5_masking_algebra_equivalence(self):
        expression = np.array([[1.0, 2.0, 3.0], [0.5, 0.0, 1.5]])
        embeddings = np.array([[1.0, 0.0], [0.0, 2.0], [1.0, 1.0]])
        full = genept_projection(expression, embeddings)
        for gene_index in range(3):
            direct = vanilla_zero_mask(expression, gene_index) @ embeddings
            subtract = genept_subtraction_mask(
                full, expression[:, gene_index], embeddings[gene_index]
            )
            np.testing.assert_allclose(subtract, direct, atol=1e-12)
            self.assertTrue(np.all(vanilla_zero_mask(expression, gene_index)[:, gene_index] == 0))
        phase5 = (ROOT / "scripts/scpa/run_phase5_gene_masking_core.R").read_text()
        phase7 = (ROOT / "scripts/scpa/run_phase7_synthetic_masking_core.R").read_text()
        self.assertIn("xa_masked[, gene_index] <- 0", phase5)
        self.assertIn("xa_masked[, gene_index] <- 0", phase7)
        self.assertIn("za - tcrossprod(xa[, gene_index], ep[gene_index, ])", phase5)
        self.assertIn("za - tcrossprod(xa[, gene_index], ep[gene_index, ])", phase7)
        self.assertNotIn("row_l2", phase7)

    def test_aggregation_preserves_individual_ranks(self):
        responses = []
        for suffix, order in (("1", ["C001", "C002", "C003"]), ("2", ["C002", "C001", "C003"])):
            responses.append({
                "schema_version": "phase7.llm-ranking.v1", "pathway": "P",
                "run_id": suffix, "backend": "mock",
                "ranking": [{"candidate_id": candidate, "rank": rank} for rank, candidate in enumerate(order, start=1)],
            })
        aggregate = aggregate_rankings(responses, tie_seed=7)
        self.assertEqual({tuple(row["individual_ranks"]) for row in aggregate}, {(1, 2), (2, 1), (3, 3)})
        stability = prompt_order_spearman(responses)
        self.assertEqual(stability["pair_count"], 1)
        self.assertAlmostEqual(stability["mean_spearman"], 0.5)

    def test_scpa_checkpoint_evaluator_joins_truth_only_at_evaluation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            checkpoint_dir = root / "checkpoints"
            checkpoint_dir.mkdir()
            truth_path = root / "truth.csv"
            checkpoint_path = checkpoint_dir / "E1_masking.csv"
            output_path = root / "metrics.csv"
            with truth_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["experiment_id", "gene", "is_ground_truth_perturbed"],
                )
                writer.writeheader()
                writer.writerows([
                    {"experiment_id": "E1", "gene": "A", "is_ground_truth_perturbed": True},
                    {"experiment_id": "E1", "gene": "B", "is_ground_truth_perturbed": False},
                    {"experiment_id": "E1", "gene": "C", "is_ground_truth_perturbed": True},
                ])
            with checkpoint_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "experiment_id", "gene", "vanilla_signed_rank", "genept_signed_rank"
                    ],
                )
                writer.writeheader()
                writer.writerows([
                    {"experiment_id": "E1", "gene": "A", "vanilla_signed_rank": 1, "genept_signed_rank": 3},
                    {"experiment_id": "E1", "gene": "B", "vanilla_signed_rank": 3, "genept_signed_rank": 2},
                    {"experiment_id": "E1", "gene": "C", "vanilla_signed_rank": 2, "genept_signed_rank": 1},
                ])
            rows = evaluate_scpa_directory(checkpoint_dir, truth_path, output_path)
            self.assertEqual(len(rows), 2)
            by_method = {row["method"]: row for row in rows}
            self.assertEqual(by_method["vanilla_scpa"]["recall_at_truth_k"], 1.0)
            self.assertEqual(by_method["genept_scpa"]["recall_at_truth_k"], 0.5)
            self.assertTrue(output_path.is_file())


if __name__ == "__main__":
    unittest.main()
