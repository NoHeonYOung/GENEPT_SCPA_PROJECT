import json
from pathlib import Path
import unittest

from scripts.phase5.run_gene_contribution import (
    add_ranks,
    build_pathway_summary,
    phase4b_gate,
    representative_targets,
    select_targets,
)


ROOT = Path(__file__).resolve().parents[1]


class Phase5GeneContributionTests(unittest.TestCase):
    def test_phase4b_gate_and_frozen_target_counts(self):
        qc, rows, _ = phase4b_gate()
        self.assertEqual(qc["gate"]["status"], "READY_FOR_GPT_REVIEW")
        targets = select_targets(rows)
        self.assertEqual(len(targets), 30)
        self.assertEqual(
            {comparison: sum(t["comparison"] == comparison for t in targets) for comparison in {
                "cd4_0h_vs_12h", "cd4_12h_vs_24h", "cd4_0h_vs_24h"
            }},
            {"cd4_0h_vs_12h": 11, "cd4_12h_vs_24h": 9, "cd4_0h_vs_24h": 10},
        )

    def test_representative_selection_is_deterministic_and_has_six_panels(self):
        _, rows, _ = phase4b_gate()
        targets = select_targets(rows)
        first = representative_targets(targets)
        self.assertEqual(first, representative_targets(list(reversed(targets))))
        self.assertEqual(len(first), 6)

    def test_average_tied_gene_ranks_and_summary_are_deterministic(self):
        rows = []
        for index, gene in enumerate(("A", "B", "C", "D")):
            rows.append({
                "comparison": "cd4_0h_vs_12h", "pathway": "P", "gene": gene,
                "detection_state": "Vanilla-only significant",
                "vanilla_delta_score": (1.0, 0.0, 0.0, -1.0)[index],
                "genept_delta_score": (0.5, 0.0, 0.0, -0.5)[index],
                "vanilla_detection_flip": False, "genept_detection_flip": False,
            })
        ranked = add_ranks(rows)
        self.assertEqual(ranked[1]["vanilla_supporting_rank"], 2.5)
        self.assertEqual(ranked[2]["vanilla_supporting_rank"], 2.5)
        self.assertEqual(ranked, add_ranks(rows))
        summary = build_pathway_summary(ranked)
        self.assertEqual(len(summary), 1)
        json.dumps({"rows": ranked, "summary": summary})

    def test_r_core_freezes_paired_masking_and_checkpoint(self):
        core = (ROOT / "scripts/scpa/run_phase5_gene_masking_core.R").read_text()
        self.assertIn("xa_masked[, gene_index] <- 0", core)
        self.assertIn("za - tcrossprod(xa[, gene_index], ep[gene_index, ])", core)
        self.assertIn("CHECKPOINT SAVED", core)
        self.assertIn("RESUME checkpoint reused", core)
        self.assertNotIn("row_l2", core)


if __name__ == "__main__":
    unittest.main()
