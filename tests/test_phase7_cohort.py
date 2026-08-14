import unittest

from gene_embedding_project.genept_scpa.phase7.cohort import split_pseudo_conditions


class Phase7CohortTests(unittest.TestCase):
    def test_split_is_deterministic_unique_and_disjoint(self):
        cells = [f"cell_{index:04d}" for index in range(1200)]
        first = split_pseudo_conditions(cells, cells_per_condition=500, seed=20260814)
        second = split_pseudo_conditions(reversed(cells), cells_per_condition=500, seed=20260814)
        self.assertEqual(first, second)
        self.assertEqual(len(first.condition_a_ids), 500)
        self.assertEqual(len(first.condition_b_ids), 500)
        self.assertFalse(set(first.condition_a_ids) & set(first.condition_b_ids))
        self.assertEqual(len(set(first.all_ids)), 1000)

    def test_split_rejects_duplicate_or_insufficient_source_ids(self):
        with self.assertRaises(ValueError):
            split_pseudo_conditions(["A", "A", "B", "C"], cells_per_condition=2, seed=1)
        with self.assertRaises(ValueError):
            split_pseudo_conditions(["A", "B", "C"], cells_per_condition=2, seed=1)


if __name__ == "__main__":
    unittest.main()
