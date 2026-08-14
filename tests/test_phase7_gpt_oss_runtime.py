import unittest

from gene_embedding_project.genept_scpa.phase7.runtime import (
    PILOT_ONLY,
    SUPPORTED_PRIMARY,
    UNSUPPORTED_PRIMARY,
    classify_runtime,
)


REQUIREMENTS = {
    "min_compute_capability": [7, 5],
    "min_vram_gib_primary": 16,
    "min_vram_gib_pilot_class": 8,
    "min_system_ram_gib": 32,
    "min_free_disk_gib": 30,
    "required_packages": ["torch", "transformers", "accelerate", "kernels", "triton"],
    "min_triton_version": "3.4",
}


def snapshot(vram_gib=24, disk_gib=50, *, missing=(), cuda=True):
    packages = {
        "torch": "2.8.0", "transformers": "4.57.6", "accelerate": "1.10.0",
        "kernels": "0.11.0", "triton": "3.4.0", "openai-harmony": None,
    }
    for name in missing:
        packages[name] = None
    return {
        "gpu": {
            "name": "test GPU", "vram_total_bytes": vram_gib * 1024**3,
            "compute_capability": [8, 0],
        },
        "system_ram_bytes": 64 * 1024**3,
        "disk": {"free_bytes": disk_gib * 1024**3},
        "cuda": {"available": cuda, "runtime": "12.8"},
        "packages": packages,
    }


class Phase7GPTOSSRuntimeTests(unittest.TestCase):
    def test_supported_primary_requires_all_primary_resources(self):
        report = classify_runtime(snapshot(), REQUIREMENTS)
        self.assertEqual(report["status"], SUPPORTED_PRIMARY)

    def test_eight_gib_gpu_is_pilot_only_when_every_other_requirement_passes(self):
        report = classify_runtime(snapshot(vram_gib=8), REQUIREMENTS)
        self.assertEqual(report["status"], PILOT_ONLY)
        self.assertFalse(report["primary_supported"])

    def test_missing_packages_and_disk_fail_closed(self):
        report = classify_runtime(
            snapshot(vram_gib=8, disk_gib=4, missing=("accelerate", "kernels")),
            REQUIREMENTS,
        )
        self.assertEqual(report["status"], UNSUPPORTED_PRIMARY)
        self.assertTrue(any("free disk" in reason for reason in report["primary_failure_reasons"]))
        self.assertTrue(any("accelerate" in reason for reason in report["primary_failure_reasons"]))


if __name__ == "__main__":
    unittest.main()
