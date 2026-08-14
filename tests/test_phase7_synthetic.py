import unittest

import numpy as np

from gene_embedding_project.genept_scpa.phase7.synthetic_perturbation import (
    ground_truth_gene_count,
    inject_perturbation,
)


class Phase7SyntheticTests(unittest.TestCase):
    def setUp(self):
        generator = np.random.default_rng(11)
        self.a = generator.gamma(2.0, 0.5, size=(20, 20))
        self.b = generator.gamma(2.0, 0.5, size=(20, 20))
        self.genes = [f"G{index}" for index in range(20)]

    def test_truth_count_formula(self):
        self.assertEqual(ground_truth_gene_count(15), 3)
        self.assertEqual(ground_truth_gene_count(40), 6)
        self.assertEqual(ground_truth_gene_count(299), 10)

    def test_mean_shift_is_deterministic_and_does_not_mutate_inputs(self):
        before_a = self.a.copy()
        before_b = self.b.copy()
        first = inject_perturbation(
            self.a, self.b, self.genes, scenario="mean_shift", alpha=0.5, seed=101
        )
        second = inject_perturbation(
            self.a, self.b, self.genes, scenario="mean_shift", alpha=0.5, seed=101
        )
        np.testing.assert_array_equal(first.condition_b, second.condition_b)
        np.testing.assert_array_equal(self.a, before_a)
        np.testing.assert_array_equal(self.b, before_b)
        self.assertEqual(len(first.ground_truth_rows), ground_truth_gene_count(20))
        self.assertTrue(all(row["perturbation_direction"] == "positive" for row in first.ground_truth_rows))

    def test_truth_set_is_reused_across_strengths_for_same_seed(self):
        weak = inject_perturbation(
            self.a, self.b, self.genes, scenario="mean_shift", alpha=0.5, seed=909
        )
        strong = inject_perturbation(
            self.a, self.b, self.genes, scenario="mean_shift", alpha=1.0, seed=909
        )
        weak_genes = [row["gene"] for row in weak.ground_truth_rows]
        strong_genes = [row["gene"] for row in strong.ground_truth_rows]
        self.assertEqual(weak_genes, strong_genes)
        weak_delta = {row["gene"]: row["applied_log_delta"] for row in weak.ground_truth_rows}
        strong_delta = {row["gene"]: row["applied_log_delta"] for row in strong.ground_truth_rows}
        for gene in weak_genes:
            self.assertAlmostEqual(strong_delta[gene], 2.0 * weak_delta[gene])

    def test_cell_subset_changes_exactly_seeded_thirty_percent(self):
        result = inject_perturbation(
            self.a, self.b, self.genes, scenario="cell_subset", alpha=1.0,
            seed=202, cell_subset_fraction=0.30,
        )
        self.assertEqual(len(result.perturbed_cell_indices), 6)
        self.assertTrue(all(row["target_cell_count"] == 6 for row in result.ground_truth_rows))
        changed_rows = set(np.flatnonzero(np.any(result.condition_b != self.b, axis=1)).tolist())
        self.assertEqual(changed_rows, set(result.perturbed_cell_indices))

    def test_mixed_direction_clips_at_zero_and_records_directions(self):
        result = inject_perturbation(
            self.a, self.b, self.genes, scenario="mixed_direction", alpha=10.0, seed=303
        )
        self.assertGreaterEqual(float(result.condition_b.min()), 0.0)
        directions = {row["perturbation_direction"] for row in result.ground_truth_rows}
        self.assertEqual(directions, {"positive", "negative"})
        self.assertGreater(sum(row["clipped_cell_count"] for row in result.ground_truth_rows), 0)

    def test_null_changes_nothing_and_has_empty_truth(self):
        result = inject_perturbation(
            self.a, self.b, self.genes, scenario="null", alpha=0.0, seed=404
        )
        np.testing.assert_array_equal(result.condition_a, self.a)
        np.testing.assert_array_equal(result.condition_b, self.b)
        self.assertEqual(result.ground_truth_rows, tuple())


if __name__ == "__main__":
    unittest.main()
