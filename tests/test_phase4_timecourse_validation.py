import json
from pathlib import Path
import tempfile
import unittest

from scripts.phase4.run_timecourse_validation import (
    ALL_COMPARISONS,
    COMPARISONS,
    GROUPS,
    PRIMARY_COMPARISONS,
    build_reporting,
    choose_cells,
    comparison_set,
    create_validation_figures,
    load_phase4a_historical_reference,
    method_counts,
    output_profile,
    parse_args,
    ties_are_valid,
)


ROOT = Path(__file__).resolve().parents[1]


class Phase4TimecourseValidationTests(unittest.TestCase):
    def test_primary_defaults_to_three_cd4_comparisons_and_all_nine_are_preserved(self):
        self.assertEqual(len(GROUPS), 6)
        self.assertEqual(len(COMPARISONS), 9)
        self.assertEqual(COMPARISONS, ALL_COMPARISONS)
        self.assertEqual(len(PRIMARY_COMPARISONS), 3)
        self.assertEqual(parse_args([]).comparison_set, "cd4_activation")
        self.assertEqual(comparison_set("cd4_activation"), PRIMARY_COMPARISONS)
        self.assertEqual(len(comparison_set("all_9")), 9)
        self.assertTrue(all(item["id"].startswith("cd4_") for item in PRIMARY_COMPARISONS))
        self.assertIn("cd4_vs_cd8_12h", {item["id"] for item in ALL_COMPARISONS})

    def test_timepoint_sampling_is_deterministic(self):
        metadata = [
            {"cell_id": f"cell_{index:03d}", "Hour": "12"}
            for index in range(30)
        ]
        first = choose_cells(metadata, "12h", sample_size=10, seed=20260810)
        second = choose_cells(list(reversed(metadata)), "12h", sample_size=10, seed=20260810)
        self.assertEqual(first, second)
        self.assertEqual(len(set(first)), 10)

    def test_reporting_counts_and_detection_states(self):
        rows = []
        for comparison in PRIMARY_COMPARISONS:
            for index in range(123):
                significant = index == 0
                rows.append({
                    "comparison": comparison["id"],
                    "pathway": f"PATH_{index:03d}",
                    "vanilla_raw_p": 0.0001 if significant else 1.0,
                    "vanilla_adjusted_p": 0.01 if significant else 1.0,
                    "vanilla_qval": 1.0 if significant else 0.0,
                    "vanilla_rank": 1.0 if significant else 62.5,
                    "genept_raw_p": 0.0001 if significant else 1.0,
                    "genept_adjusted_p": 0.01 if significant else 1.0,
                    "genept_qval": 1.0 if significant else 0.0,
                    "genept_rank": 1.0 if significant else 62.5,
                    "l2_raw_p": 1.0,
                    "l2_adjusted_p": 1.0,
                    "l2_qval": 0.0,
                    "l2_rank": 62.0,
                })
        overview, detection, timing = build_reporting(rows)
        self.assertEqual((overview, detection, timing), build_reporting(rows))
        self.assertEqual(len(overview), 3)
        self.assertEqual(len(detection), 3 * 123)
        self.assertEqual(len(timing), 123 * 3)
        self.assertEqual(overview[0]["both_significant"], 1)
        self.assertEqual(overview[0]["neither_significant"], 122)
        self.assertTrue(ties_are_valid(rows))
        self.assertEqual(method_counts(rows[:123], "vanilla")["vanilla_qval_zero_count"], 122)
        self.assertEqual(overview[0]["vanilla_n_pathways"], 123)
        self.assertAlmostEqual(overview[0]["vanilla_qval_floor_fraction"], 122 / 123)
        self.assertEqual(overview[0]["genept_non_l2_l2_significant_overlap"], 0)
        json.dumps({"overview": overview, "detection": detection, "timing": timing})

    def test_primary_output_names_are_frozen(self):
        profile = output_profile("cd4_activation")
        self.assertEqual(profile["all_results"], "phase4_cd4_activation_all_results.csv")
        self.assertEqual(profile["qc"], "phase4_cd4_activation_qc.json")
        self.assertEqual(profile["comparison_filenames"]["cd4_0h_vs_12h"], "cd4_0_vs_12.csv")

    def test_zero_hour_sampling_exactly_reuses_phase3(self):
        timecourse = ROOT / "data/interim/genept_scpa/phase4_timecourse_sampling"
        phase3 = ROOT / "data/interim/genept_scpa/phase3_sampling"
        for lineage in ("cd4", "cd8"):
            self.assertEqual(
                (timecourse / f"{lineage}_0h_cells.txt").read_bytes(),
                (phase3 / f"{lineage}_0h_cells.txt").read_bytes(),
            )

    def test_source_audit_reports_no_download(self):
        audit = json.loads(
            (ROOT / "data/interim/genept_scpa/phase4_timecourse_source_audit.json").read_text()
        )
        self.assertEqual(audit["source_status"], "ALREADY_PRESENT_IN_RDS")
        self.assertFalse(audit["download_performed"])
        self.assertEqual(audit["gate"]["status"], "PASS")

    def test_runner_contains_no_download_or_phase5_execution(self):
        script = (ROOT / "scripts/phase4/run_timecourse_validation.py").read_text()
        for forbidden in ("curl", "wget", "Phase5_gene_contribution", "leave_one_gene_out"):
            self.assertNotIn(forbidden, script)

    def test_r_core_uses_official_log10_and_average_ties(self):
        core = (ROOT / "scripts/scpa/run_phase4_timecourse_core.R").read_text()
        self.assertIn('ties.method = "average"', core)
        adapter = (ROOT / "scripts/scpa/scpa_core_adapter.R").read_text()
        self.assertIn("sqrt(-log10(p_value))", adapter)

    def test_validation_figures_render(self):
        rows = []
        overview = []
        detection = []
        for comparison in PRIMARY_COMPARISONS:
            overview.append({
                "comparison": comparison["id"],
                "vanilla_qval_zero_count": 1, "genept_qval_zero_count": 1,
                "l2_qval_zero_count": 1, "vanilla_adj_p_lt_0_05": 1,
                "genept_adj_p_lt_0_05": 1, "l2_adj_p_lt_0_05": 0,
            })
            for index in range(2):
                pathway = f"PATH_{index}"
                rows.append({
                    "comparison": comparison["id"], "pathway": pathway,
                    "vanilla_qval": float(index + 0.5),
                    "genept_qval": float(index + 1),
                    "vanilla_rank": float(index + 1),
                    "genept_rank": float(index + 1),
                })
                detection.append({
                    "comparison": comparison["id"], "pathway": pathway,
                    "detection_state": "Both significant" if index == 0 else "Neither significant",
                })
        with tempfile.TemporaryDirectory() as directory:
            files = create_validation_figures(
                overview,
                detection,
                rows,
                Path(directory),
                load_phase4a_historical_reference(),
            )
            self.assertEqual(len(files), 5)
            self.assertTrue(all(path.is_file() and path.stat().st_size > 0 for path in files))


if __name__ == "__main__":
    unittest.main()
