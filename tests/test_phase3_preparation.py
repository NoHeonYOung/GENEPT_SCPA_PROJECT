from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy import io as scipy_io

from scripts.phase3.run_cd4_cd8_benchmark import (
    aligned_shared_gene_indices,
    canonical_hour,
    rows_for_ids,
    select_canonical_cells,
    write_matrix_market_atomic,
)


class Phase3PreparationTests(unittest.TestCase):
    def test_hour_canonicalization(self):
        self.assertEqual(canonical_hour("0"), "0h")
        self.assertEqual(canonical_hour("0.0"), "0h")
        self.assertEqual(canonical_hour("0 hours"), "0h")

    def test_canonical_sampling_is_deterministic(self):
        metadata = [
            {"cell_id": f"cell{i}", "Hour": "0", "Cell_Type": "Resting"}
            for i in range(20)
        ]
        all_a, selected_a = select_canonical_cells(
            metadata, hour="0h", sample_size=5, seed=17
        )
        all_b, selected_b = select_canonical_cells(
            list(reversed(metadata)), hour="0", sample_size=5, seed=17
        )
        self.assertEqual(all_a, all_b)
        self.assertEqual(selected_a, selected_b)
        self.assertEqual(len(set(selected_a)), 5)

    def test_shared_gene_alignment_uses_symbols_not_positions(self):
        shared, cd4, cd8 = aligned_shared_gene_indices(
            ["B", "A", "CD4_ONLY"], ["C", "B", "A"]
        )
        self.assertEqual(shared, ["A", "B"])
        self.assertEqual(cd4.tolist(), [1, 0])
        self.assertEqual(cd8.tolist(), [2, 1])

    def test_rows_follow_canonical_id_order(self):
        rows = rows_for_ids(["b", "a", "c"], ["c", "b"])
        self.assertEqual(rows.tolist(), [2, 0])

    def test_adapter_matrix_market_roundtrip(self):
        matrix = np.array([[1.0, 0.0], [2.0, 3.0]])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "adapter.mtx"
            write_matrix_market_atomic(matrix, path)
            observed = scipy_io.mmread(path).toarray()
        np.testing.assert_array_equal(observed, matrix)


if __name__ == "__main__":
    unittest.main()
