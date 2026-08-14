import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy import io as scipy_io
from scipy import sparse

from gene_embedding_project.genept_scpa.phase7.transfer import verify_transfer_manifest
from gene_embedding_project.genept_scpa.io import sha256_file
from scripts.phase7.prepare_synthetic_benchmark import load_portable_counts_export


class Phase7TransferIntegrityTests(unittest.TestCase):
    def manifest(self, root: Path, expected_size: int, expected_hash: str):
        payload = {
            "schema_version": "phase7.transfer.v1",
            "resources": [{
                "id": "tiny", "repository_relative_path": "data/tiny.txt",
                "required_status": "required", "size_bytes": expected_size,
                "sha256": expected_hash,
            }],
        }
        path = root / "manifest.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def test_exact_file_passes_without_modification(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "data/tiny.txt"
            path.parent.mkdir()
            path.write_text("phase7", encoding="utf-8")
            before = path.stat().st_mtime_ns
            digest = hashlib.sha256(b"phase7").hexdigest()
            report = verify_transfer_manifest(
                self.manifest(root, 6, digest), root
            )
            self.assertEqual(report["status"], "PASS")
            self.assertEqual(path.stat().st_mtime_ns, before)
            self.assertFalse(report["files_modified"])
            self.assertFalse(report["downloads_performed"])
            self.assertFalse(report["alternate_paths_inferred"])

    def test_missing_or_changed_required_file_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            digest = hashlib.sha256(b"expected").hexdigest()
            missing = verify_transfer_manifest(
                self.manifest(root, 8, digest), root
            )
            self.assertEqual(missing["status"], "FAIL")
            path = root / "data/tiny.txt"
            path.parent.mkdir()
            path.write_text("changed!", encoding="utf-8")
            changed = verify_transfer_manifest(
                self.manifest(root, 8, digest), root
            )
            self.assertEqual(changed["status"], "FAIL")

    def test_absolute_manifest_path_is_rejected_without_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.json"
            manifest.write_text(json.dumps({
                "schema_version": "phase7.transfer.v1",
                "resources": [{
                    "id": "bad", "repository_relative_path": "/tmp/alternate",
                    "required_status": "required", "size_bytes": 0, "sha256": "0" * 64,
                }],
            }), encoding="utf-8")
            with self.assertRaises(ValueError):
                verify_transfer_manifest(manifest, root)

    def test_portable_counts_loader_uses_only_configured_relative_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = root / "portable"
            data.mkdir()
            files = {
                "counts": "portable/counts.mtx",
                "genes": "portable/genes.txt",
                "cell_ids": "portable/cells.txt",
                "metadata": "portable/metadata.csv",
            }
            scipy_io.mmwrite(root / files["counts"], sparse.coo_matrix([[1, 0], [2, 3]]))
            (root / files["genes"]).write_text("G1\nG2\n", encoding="utf-8")
            (root / files["cell_ids"]).write_text("C1\nC2\n", encoding="utf-8")
            (root / files["metadata"]).write_text(
                "cell_id,Hour,Cell_Type\nC1,0,Naive CD4 T\nC2,0,Naive CD4 T\n",
                encoding="utf-8",
            )
            manifest = {
                "assay": "RNA", "layer": "counts",
                "matrix_orientation": "genes_by_cells", "genes": 2, "cells": 2,
                "files": {key: f"/obsolete/lab/path/{key}" for key in files},
                "sha256": {key: sha256_file(root / relative) for key, relative in files.items()},
            }
            manifest_path = data / "manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            _, counts, genes, cells = load_portable_counts_export(
                manifest_path, files, repository_root=root
            )
            self.assertEqual(genes, ["G1", "G2"])
            self.assertEqual(cells, ["C1", "C2"])
            np.testing.assert_array_equal(counts.toarray(), np.array([[1, 2], [0, 3]]))


if __name__ == "__main__":
    unittest.main()
