import json
from pathlib import Path
import unittest

import numpy as np

from scripts.phase6.run_semantic_controls import (
    GENE_REPS,
    PATHWAY_REPS,
    RESAMPLING_REPS,
    TOTAL_GENE_CONTROL_BASELINE_MCM,
    TOTAL_GENE_MASK_CONTROL_MCM,
    TOTAL_PATHWAY_CONTROL_MCM,
    TOTAL_ROBUSTNESS_MCM,
    bh_adjust,
    frozen_targets,
    parse_value,
)


ROOT = Path(__file__).resolve().parents[1]


class Phase6SemanticControlTests(unittest.TestCase):
    def test_frozen_target_counts_and_representatives(self):
        targets, gene_targets = frozen_targets()
        self.assertEqual(len(targets), 30)
        self.assertEqual(len(gene_targets), 6)
        self.assertEqual(sum(int(row["n_paired_genes"]) for row in gene_targets), 299)

    def test_exact_replication_and_mcm_budget(self):
        self.assertEqual(PATHWAY_REPS, 100)
        self.assertEqual(GENE_REPS, 20)
        self.assertEqual(RESAMPLING_REPS, 10)
        self.assertEqual(TOTAL_PATHWAY_CONTROL_MCM, 6000)
        self.assertEqual(TOTAL_GENE_MASK_CONTROL_MCM, 11960)
        self.assertEqual(TOTAL_GENE_CONTROL_BASELINE_MCM, 240)
        self.assertEqual(TOTAL_ROBUSTNESS_MCM, 900)

    def test_bh_adjustment_is_valid_and_monotone_in_p_order(self):
        values = [0.01, 0.04, 0.03, 0.9]
        adjusted = bh_adjust(values)
        self.assertTrue(all(0 <= value <= 1 for value in adjusted))
        ordered = np.argsort(values)
        self.assertTrue(all(adjusted[int(ordered[i])] <= adjusted[int(ordered[i + 1])] for i in range(3)))
        json.dumps(adjusted)

    def test_lowercase_true_representation_is_not_parsed_as_boolean(self):
        self.assertEqual(parse_value("true"), "true")
        self.assertIs(parse_value("TRUE"), True)
        self.assertIs(parse_value("False"), False)

    def test_r_core_freezes_semantic_controls_and_no_l2(self):
        core = (ROOT / "scripts/scpa/run_phase6_semantic_controls_core.R").read_text()
        self.assertIn("ep[permutation, , drop = FALSE]", core)
        self.assertIn("stats::rnorm(p * ncol(ep))", core)
        self.assertIn("mean(permutation != seq_len(p)) <= 0.9", core)
        self.assertIn("atomic_csv", core)
        self.assertNotIn("row_l2", core)


if __name__ == "__main__":
    unittest.main()
