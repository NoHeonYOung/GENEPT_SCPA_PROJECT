from pathlib import Path
import json
import tempfile
import unittest

import numpy as np
from scipy import sparse

from gene_embedding_project.genept_scpa.genept_projection import normalize_log1p_sparse
from gene_embedding_project.genept_scpa.pathway_projection import (
    PathwayDefinition,
    average_rank_descending,
    build_paired_pathways,
    filter_eligible_pathways,
    project_pathway,
    ranking_agreement,
    read_wide_pathway_csv,
    scpa_bonferroni_qvalues,
    significance_state,
)
from scripts.phase4.run_pathway_comparison import build_qc, create_figures


ROOT = Path(__file__).resolve().parents[1]


class Phase4PathwayProjectionTests(unittest.TestCase):
    def test_hand_calculated_projection_and_shape(self):
        expression = np.array([[1.0, 2.0, 0.0, 1.0], [0.0, 1.0, 3.0, 2.0]])
        embeddings = np.array(
            [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [-1.0, 2.0]]
        )
        genes = ["A", "B", "C", "D"]
        expected = expression @ embeddings
        observed = project_pathway(expression, genes, embeddings, genes)
        self.assertEqual(observed.shape, (2, 2))
        np.testing.assert_array_equal(observed, expected)

    def test_projection_is_deterministic_and_does_not_mutate_inputs(self):
        expression = np.array([[1.0, 2.0], [3.0, 4.0]])
        embeddings = np.array([[1.0, 0.0], [0.0, 1.0]])
        before_expression = expression.copy()
        before_embeddings = embeddings.copy()
        first = project_pathway(expression, ["A", "B"], embeddings, ["A", "B"])
        second = project_pathway(expression, ["A", "B"], embeddings, ["A", "B"])
        np.testing.assert_array_equal(first, second)
        np.testing.assert_array_equal(expression, before_expression)
        np.testing.assert_array_equal(embeddings, before_embeddings)

    def test_gene_order_mismatch_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "gene order"):
            project_pathway(
                np.ones((3, 2)), ["A", "B"], np.eye(2), ["B", "A"]
            )

    def test_paired_gene_set_is_identical_by_construction(self):
        pathway = PathwayDefinition("KEGG_TOY", "KEGG", ("A", "B", "C", "D"))
        paired = build_paired_pathways(
            [pathway], ["A", "B", "C"], ["A", "B", "D"], {"A", "C", "D"}
        )[0]
        self.assertEqual(paired.shared_genes, ("A", "B"))
        self.assertEqual(paired.genept_mappable_genes, ("A", "C", "D"))
        self.assertEqual(paired.paired_genes, ("A",))
        self.assertEqual(filter_eligible_pathways([paired], min_genes=1, max_genes=10), [paired])

    def test_normalize_full_transcriptome_then_subset_differs_from_wrong_order(self):
        counts = sparse.csr_matrix([[1.0, 1.0, 8.0]])
        correct = normalize_log1p_sparse(counts)[:, :2].toarray()
        wrong = normalize_log1p_sparse(counts[:, :2]).toarray()
        self.assertFalse(np.allclose(correct, wrong))
        np.testing.assert_allclose(correct, np.log1p([[1000.0, 1000.0]]))

    def test_l2_sensitivity_is_optional_and_rowwise(self):
        expression = np.array([[2.0, 0.0], [0.0, 3.0]])
        embeddings = np.eye(2)
        non_l2 = project_pathway(expression, ["A", "B"], embeddings, ["A", "B"])
        l2 = project_pathway(
            expression, ["A", "B"], embeddings, ["A", "B"], l2_normalize=True
        )
        np.testing.assert_array_equal(non_l2, expression)
        np.testing.assert_allclose(np.linalg.norm(l2, axis=1), 1.0)

    def test_scpa_official_bonferroni_and_base10_qval(self):
        adjusted, qval = scpa_bonferroni_qvalues([0.001, 0.2])
        np.testing.assert_allclose(adjusted, [0.002, 0.4])
        np.testing.assert_allclose(qval, np.sqrt(-np.log10(adjusted)))

    def test_average_rank_does_not_break_qval_zero_ties(self):
        observed = average_rank_descending([3.0, 2.0, 0.0, 0.0, 0.0])
        np.testing.assert_array_equal(observed, [1.0, 2.0, 4.0, 4.0, 4.0])

    def test_significance_categories(self):
        self.assertEqual(significance_state(0.01, 0.01), "Both significant")
        self.assertEqual(significance_state(0.01, 1.0), "Vanilla-only significant")
        self.assertEqual(significance_state(1.0, 0.01), "GenePT-only significant")
        self.assertEqual(significance_state(1.0, 1.0), "Neither significant")

    def test_ranking_agreement_top_overlap(self):
        result = ranking_agreement([1, 2, 3, 4], [1, 3, 2, 4])
        self.assertAlmostEqual(result["spearman"], 0.8)
        self.assertEqual(result["top10_overlap"], 4)
        self.assertEqual(result["top20_jaccard"], 1.0)

    def test_official_pathway_collection_is_read_without_header(self):
        pathways = read_wide_pathway_csv(
            ROOT / "data/reference/genept_scpa/combined_metabolic_pathways.csv"
        )
        self.assertEqual(len(pathways), 243)
        self.assertEqual(pathways[0].name, "HALLMARK_BILE_ACID_METABOLISM")
        self.assertEqual(pathways[0].source_database, "HALLMARK")

    def test_all_predeclared_figures_render(self):
        rows = [
            {
                "pathway": f"PATH_{index}",
                "vanilla_rank": index,
                "genept_rank": 5 - index,
                "vanilla_qval": float(6 - index),
                "genept_qval": float(index),
            }
            for index in range(1, 5)
        ]
        with tempfile.TemporaryDirectory() as directory:
            files = create_figures(rows, Path(directory))
            self.assertEqual(len(files), 4)
            self.assertTrue(all(path.is_file() and path.stat().st_size > 0 for path in files))

    def test_qc_payload_is_standard_json_serializable(self):
        manifest = {
            "cohort": {"canonical_sampling_reused": True},
            "pathway_collection": {"eligible_full_count": 1},
            "paired_gene_policy": {"identical_between_branches": True},
            "preprocessing": {"pathway_renormalization": False},
            "scope": {},
        }
        rows = [{
            "n_primary_paired_genes": 15,
            "primary_gene_set_identical": True,
            "vanilla_raw_p": 0.1,
            "vanilla_adjusted_p": 0.1,
            "vanilla_qval": 1.0,
            "genept_raw_p": 0.2,
            "genept_adjusted_p": 0.2,
            "genept_qval": 0.8,
            "embedding_rank": 15,
            "projected_rank": 14,
        }]
        core_qc = {
            "warnings": [], "scpa_version": "1.6.2", "multicross_version": "2.1.0",
            "raw_p_source": "multicross::mcm result[[1]]",
            "multiple_testing": "bonferroni", "qval_formula": "sqrt(-log10(adjusted_p))",
            "log_base": 10, "effective_rank": [],
        }
        qc = build_qc(manifest, rows, core_qc, {"top10_overlap": 1}, [], run_l2=True)
        encoded = json.dumps(qc)
        self.assertIn('"valid_results": 1', encoded)


if __name__ == "__main__":
    unittest.main()
