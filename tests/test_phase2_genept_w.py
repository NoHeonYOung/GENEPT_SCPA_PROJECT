import json
from pathlib import Path
import pickle
import tempfile
import unittest

import numpy as np
from scipy import sparse

from scripts.genept.build_genept_w import (
    DATASET_SETTINGS,
    run_synthetic_correctness_check,
)
from gene_embedding_project.genept_scpa.gene_mapping import (
    build_aligned_embedding_matrix,
    classify_gene_matches,
    load_official_genept_embeddings,
    load_primary_gene_keys,
    mapping_counts,
)
from gene_embedding_project.genept_scpa.genept_projection import (
    project_genept_w_direct,
    project_genept_w_sparse,
)


ROOT = Path(__file__).resolve().parents[1]


class GenePTWCorrectnessTests(unittest.TestCase):
    def setUp(self):
        self.genes = ["G1", "G2", "G3"]
        self.embeddings = {
            "G1": np.array([1.0, 0.0], dtype=np.float32),
            "G2": np.array([0.0, 1.0], dtype=np.float32),
            "G3": np.array([1.0, 1.0], dtype=np.float32),
            "ALIAS3": np.array([1.0, 1.0], dtype=np.float32),
        }
        self.counts = sparse.csr_matrix(
            np.array([[1.0, 1.0, 0.0], [0.0, 2.0, 2.0]])
        )
        matches = classify_gene_matches(
            self.genes, set(self.embeddings), {"G1", "G2", "G3"}
        )
        self.indices, self.embedding_matrix = build_aligned_embedding_matrix(
            matches, self.embeddings
        )

    def test_hand_calculated_toy_example(self):
        observed, diagnostics = project_genept_w_sparse(
            self.counts,
            self.indices,
            self.embedding_matrix,
            normalization_target=10_000,
            batch_size=1,
        )
        weight_cell_1 = np.log1p(5_000.0)
        raw_cell_1 = np.array([weight_cell_1, weight_cell_1]) / 3.0
        expected_cell_1 = raw_cell_1 / np.linalg.norm(raw_cell_1)

        weight_cell_2 = np.log1p(5_000.0)
        raw_cell_2 = np.array([weight_cell_2, 2.0 * weight_cell_2]) / 3.0
        expected_cell_2 = raw_cell_2 / np.linalg.norm(raw_cell_2)
        expected = np.vstack([expected_cell_1, expected_cell_2])
        np.testing.assert_allclose(observed, expected, rtol=0, atol=1e-7)
        np.testing.assert_allclose(diagnostics.post_l2_norms, 1.0, atol=1e-12)

    def test_build_preflight_synthetic_gate(self):
        result = run_synthetic_correctness_check()
        self.assertLess(result["hand_calculation_max_abs_error"], 1e-6)
        self.assertLess(result["optimized_vs_direct_max_abs_error"], 1e-6)
        self.assertEqual(result["batch_determinism_max_abs_error"], 0.0)

    def test_sparse_matches_explicit_loop(self):
        optimized, _ = project_genept_w_sparse(
            self.counts, self.indices, self.embedding_matrix, batch_size=2
        )
        direct = project_genept_w_direct(
            self.counts, self.indices, self.embedding_matrix
        )
        np.testing.assert_allclose(optimized, direct, rtol=0, atol=1e-7)

    def test_deterministic_across_batch_sizes(self):
        one, _ = project_genept_w_sparse(
            self.counts, self.indices, self.embedding_matrix, batch_size=1
        )
        two, _ = project_genept_w_sparse(
            self.counts, self.indices, self.embedding_matrix, batch_size=2
        )
        np.testing.assert_array_equal(one, two)

    def test_unmatched_gene_is_kept_in_normalization_denominator(self):
        genes = ["G1", "UNMATCHED"]
        counts = sparse.csr_matrix([[1.0, 9.0]])
        matches = classify_gene_matches(genes, set(self.embeddings), {"G1"})
        indices, matrix = build_aligned_embedding_matrix(matches, self.embeddings)
        _, diagnostics = project_genept_w_sparse(counts, indices, matrix)
        self.assertAlmostEqual(diagnostics.library_sizes[0], 10.0)
        self.assertAlmostEqual(diagnostics.expression_coverage[0], 0.1)

    def test_gene_order_invariance(self):
        original, _ = project_genept_w_sparse(
            self.counts, self.indices, self.embedding_matrix
        )
        order = np.array([2, 0, 1])
        reordered_genes = [self.genes[index] for index in order]
        reordered_counts = self.counts[:, order]
        matches = classify_gene_matches(
            reordered_genes, set(self.embeddings), {"G1", "G2", "G3"}
        )
        indices, matrix = build_aligned_embedding_matrix(matches, self.embeddings)
        reordered, _ = project_genept_w_sparse(reordered_counts, indices, matrix)
        np.testing.assert_allclose(original, reordered, rtol=0, atol=1e-7)

    def test_zero_count_cell_is_reported(self):
        counts = sparse.vstack([self.counts, sparse.csr_matrix((1, 3))])
        observed, diagnostics = project_genept_w_sparse(
            counts, self.indices, self.embedding_matrix
        )
        self.assertEqual(diagnostics.zero_vectors, 1)
        np.testing.assert_array_equal(observed[-1], np.zeros(2, dtype=np.float32))

    def test_strict_matching_and_official_alias_classification(self):
        matches = classify_gene_matches(
            ["G1", "ALIAS3", "g2", "G1"],
            set(self.embeddings),
            {"G1", "G2", "G3"},
        )
        self.assertEqual(
            [match.match_type for match in matches],
            ["exact", "official_alias", "unmatched", "exact"],
        )
        counts = mapping_counts(matches)
        self.assertEqual(counts["exact_matches"], 2)
        self.assertEqual(counts["alias_matches"], 1)
        self.assertEqual(counts["unmatched_dataset_genes"], 1)
        self.assertEqual(counts["duplicate_mapping_count"], 1)

    def test_pinned_artifact_schema_loaders(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            embedding_path = root / "embedding.pickle"
            summary_path = root / "summaries.json"
            with embedding_path.open("wb") as handle:
                pickle.dump(
                    {
                        "G1": np.array([1.0, 0.0], dtype=np.float32),
                        "ALIAS3": np.array([0.0, 1.0], dtype=np.float32),
                    },
                    handle,
                )
            summary_path.write_text(json.dumps({"G1": "summary"}), encoding="utf-8")
            observed = load_official_genept_embeddings(
                embedding_path, expected_dimension=2
            )
            self.assertEqual(set(observed), {"G1", "ALIAS3"})
            self.assertEqual(load_primary_gene_keys(summary_path), {"G1"})

    def test_invalid_artifact_dimension_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.pickle"
            with path.open("wb") as handle:
                pickle.dump({"G1": [1.0, 2.0, 3.0]}, handle)
            with self.assertRaisesRegex(ValueError, "expected"):
                load_official_genept_embeddings(path, expected_dimension=2)

    def test_phase2_scripts_do_not_run_downstream_analysis(self):
        script = (ROOT / "scripts/genept/build_genept_w.py").read_text(
            encoding="utf-8"
        )
        forbidden_calls = ("compare_pathways(", "LogisticRegression(", "KNeighborsClassifier(")
        for call in forbidden_calls:
            self.assertNotIn(call, script)

    def test_acquisition_is_pinned_and_does_not_call_openai(self):
        script = (ROOT / "scripts/genept/prepare_genept_embeddings.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("GenePT_gene_embedding_ada_text.pickle", script)
        self.assertIn("3f6ce4317e3a0091978ae5cb8fbf05a3", script)
        self.assertIn("--continue-at", script)
        self.assertIn("archive_incomplete", script)
        self.assertNotIn("OpenAI(", script)
        self.assertNotIn("embeddings.create(", script)

    def test_cd8_reuses_builder_without_overwriting_cd4_outputs(self):
        cd4 = DATASET_SETTINGS["naive_cd4"]
        cd8 = DATASET_SETTINGS["naive_cd8"]
        self.assertEqual(cd4["output_directory"], "data/processed/genept_scpa/phase2")
        self.assertEqual(cd8["output_directory"], "data/processed/genept_scpa/phase3")
        self.assertNotEqual(cd4["mapping_filename"], cd8["mapping_filename"])
        self.assertNotEqual(cd4["qc"], cd8["qc"])


if __name__ == "__main__":
    unittest.main()
