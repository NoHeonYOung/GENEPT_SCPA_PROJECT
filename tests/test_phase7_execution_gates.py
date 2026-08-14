from pathlib import Path
import tempfile
import unittest

from scripts.phase7.prepare_synthetic_benchmark import prepare
from scripts.phase7.run_gpt_oss_inference import run as run_gpt_oss
from scripts.phase7.run_scpa_masking import run as run_scpa


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config/phase7_gpt_oss_synthetic.yaml"


class Phase7ExecutionGateTests(unittest.TestCase):
    def test_production_synthetic_generation_is_locked_before_data_access(self):
        with self.assertRaisesRegex(RuntimeError, "production synthetic generation remains locked"):
            prepare(CONFIG)

    def test_production_scpa_is_locked_before_subprocess(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(RuntimeError, "production SCPA remains locked"):
                run_scpa(
                    CONFIG, root / "missing.h5", root / "missing.json",
                    root / "checkpoints", root / "output.json",
                )

    def test_real_gpt_oss_is_locked_before_runtime_or_model_access(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(RuntimeError, "Real gpt-oss inference remains locked"):
                run_gpt_oss(
                    CONFIG, root / "requests", root / "traces", root / "invalid",
                    root / "runtime.json",
                )
            self.assertFalse((root / "runtime.json").exists())


if __name__ == "__main__":
    unittest.main()
