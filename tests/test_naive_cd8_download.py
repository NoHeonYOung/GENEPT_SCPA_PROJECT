from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class NaiveCD8DownloadScriptTests(unittest.TestCase):
    def test_script_is_pinned_and_does_not_overwrite_cd4(self):
        script = (ROOT / "scripts/data/download_naive_cd8.sh").read_text(encoding="utf-8")
        self.assertIn("ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212270/suppl/", script)
        self.assertIn("GSE212270_integrated_naive_cd8.rds.gz", script)
        self.assertIn("naive_cd8_download_metadata.json", script)
        self.assertNotIn("GSE212270_integrated_naive_cd4.rds.gz", script)

    def test_script_has_resume_and_partial_file_safeguards(self):
        script = (ROOT / "scripts/data/download_naive_cd8.sh").read_text(encoding="utf-8")
        self.assertIn("--continue-at -", script)
        self.assertIn('partial="${archive}.part"', script)
        self.assertIn("gzip -t", script)
        self.assertIn("download-metadata", script)


if __name__ == "__main__":
    unittest.main()
