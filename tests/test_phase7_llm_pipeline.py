import json
from pathlib import Path
import unittest

import numpy as np

from gene_embedding_project.genept_scpa.phase7.llm_backend import MockLLMBackend
from gene_embedding_project.genept_scpa.phase7.llm_prompts import build_llm_request
from gene_embedding_project.genept_scpa.phase7.pathway_selection import (
    sanitize_gene_description,
    select_phase7_pathways,
)
from gene_embedding_project.genept_scpa.phase7.schemas import (
    validate_llm_request,
    validate_llm_response,
)


ROOT = Path(__file__).resolve().parents[1]


class Phase7LLMPipelineTests(unittest.TestCase):
    def setUp(self):
        self.genes = ["G1", "G2", "G3", "G4"]
        self.a = np.arange(24, dtype=float).reshape(6, 4) / 10
        self.b = self.a + np.array([0.0, 0.2, 0.5, 0.1])
        self.descriptions = {
            gene: f"Gene Symbol {gene} enzyme description {index}"
            for index, gene in enumerate(self.genes)
        }

    def build(self, condition, order_seed=1, shuffle_seed=2):
        return build_llm_request(
            experiment_id="E1", run_id=f"E1_{condition}_{order_seed}", pathway="P1",
            source_database="KEGG", genes=self.genes, condition_a=self.a,
            condition_b=self.b, descriptions=self.descriptions,
            prompt_condition=condition, candidate_order_seed=order_seed,
            description_shuffle_seed=shuffle_seed, backend=MockLLMBackend.name,
        )

    def test_gene_symbol_prefix_is_removed_and_mapping_stays_internal(self):
        self.assertEqual(sanitize_gene_description("G1", "Gene Symbol G1 useful enzyme"), "useful enzyme")
        self.assertEqual(
            sanitize_gene_description("G1", "Gene Symbol G1 G1 regulates metabolism"),
            "this gene regulates metabolism",
        )
        bundle = self.build("true_description")
        encoded = json.dumps(bundle.request)
        self.assertNotIn('"gene"', encoded)
        self.assertNotIn("Gene Symbol", encoded)
        for gene in self.genes:
            self.assertNotIn(gene, encoded)
        self.assertEqual(set(bundle.candidate_to_gene.values()), set(self.genes))
        validate_llm_request(bundle.request)

    def test_stats_only_has_no_description(self):
        bundle = self.build("stats_only")
        self.assertTrue(all("description" not in row for row in bundle.request["candidates"]))

    def test_shuffled_descriptions_are_a_deranged_exact_multiset(self):
        true_bundle = self.build("true_description", order_seed=9)
        shuffled = self.build("shuffled_description", order_seed=9)
        true_descriptions = sorted(row["description"] for row in true_bundle.request["candidates"])
        shuffled_descriptions = sorted(row["description"] for row in shuffled.request["candidates"])
        self.assertEqual(true_descriptions, shuffled_descriptions)
        self.assertEqual(shuffled.description_changed_fraction, 1.0)
        self.assertTrue(all(key != value for key, value in shuffled.description_source_candidate.items()))

    def test_mock_backend_is_deterministic_and_schema_valid(self):
        backend = MockLLMBackend()
        request = self.build("stats_only").request
        first = backend.rank(request)
        self.assertEqual(first, backend.rank(request))
        self.assertFalse(backend.scientific_evaluation_allowed)
        validate_llm_response(first, expected_candidate_ids=[row["candidate_id"] for row in request["candidates"]])

    def test_real_frozen_pathway_selection_is_outcome_independent_and_deterministic(self):
        manifest = json.loads((ROOT / "data/processed/genept_scpa/phase4/pathway_projection_manifest.json").read_text())
        descriptions = json.loads((ROOT / "data/reference/genept_scpa/genept_ada002/NCBI_summary_of_genes.json").read_text())
        bins = {"small": [15, 25], "medium": [26, 40], "large": [41, 60]}
        quotas = {
            "KEGG": {"small": 2, "medium": 2, "large": 2},
            "REACTOME": {"small": 2, "medium": 2, "large": 1},
            "HALLMARK": {"small": 0, "medium": 0, "large": 0},
        }
        first, audit = select_phase7_pathways(manifest, descriptions, size_bins=bins, source_bin_quota=quotas, seed=20260815)
        second, _ = select_phase7_pathways(manifest, descriptions, size_bins=bins, source_bin_quota=quotas, seed=20260815)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 11)
        self.assertFalse(audit["selection_uses_phase4_or_phase5_scores"])
        self.assertEqual({row["source_database"] for row in first}, {"KEGG", "REACTOME"})


if __name__ == "__main__":
    unittest.main()
