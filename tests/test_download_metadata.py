import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from gene_embedding_project.genept_scpa.io import (
    build_download_metadata,
    write_json_atomic,
)


class DownloadMetadataTests(unittest.TestCase):
    def test_download_metadata_schema_uses_observed_values(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            archive = root / "GSE212270_integrated_naive_cd4.rds.gz"
            with gzip.open(archive, "wb") as handle:
                handle.write(b"small synthetic fixture")

            metadata = build_download_metadata(
                archive,
                geo_accession="GSE212270",
                download_source="https://ftp.ncbi.nlm.nih.gov/example",
            )
            output = root / "interim" / "phase1_download_metadata.json"
            write_json_atomic(metadata, output)
            observed = json.loads(output.read_text(encoding="utf-8"))

            required = {
                "geo_accession",
                "filename",
                "download_source",
                "file_size_bytes",
                "sha256",
                "gzip_integrity",
                "recorded_at",
                "recorded_at_utc",
            }
            self.assertTrue(required.issubset(observed))
            self.assertEqual(observed["file_size_bytes"], archive.stat().st_size)
            self.assertEqual(len(observed["sha256"]), 64)
            self.assertTrue(observed["gzip_integrity"])

    def test_invalid_gzip_fails(self):
        with TemporaryDirectory() as directory:
            archive = Path(directory) / "invalid.rds.gz"
            archive.write_bytes(b"not gzip")
            with self.assertRaisesRegex(ValueError, "Gzip integrity check failed"):
                build_download_metadata(
                    archive,
                    geo_accession="GSE212270",
                    download_source="https://ftp.ncbi.nlm.nih.gov/example",
                )


if __name__ == "__main__":
    unittest.main()
