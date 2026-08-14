import unittest

from scripts.phase7.run_toy_smoke import run_toy_smoke


class Phase7ToySmokeTests(unittest.TestCase):
    def test_end_to_end_toy_smoke_uses_no_real_model_or_scpa(self):
        result = run_toy_smoke()
        self.assertEqual(result["status"], "PASS")
        self.assertFalse(result["scientific_evaluation_allowed"])
        self.assertFalse(result["real_llm_inference"])
        self.assertFalse(result["production_scpa"])
        self.assertFalse(result["real_dataset_used"])
        self.assertEqual(result["truth_gene_count"], 3)
        self.assertEqual(set(result["llm"]), {"stats_only", "true_description", "shuffled_description"})
        self.assertTrue(all(value["individual_run_count"] == 3 for value in result["llm"].values()))


if __name__ == "__main__":
    unittest.main()
