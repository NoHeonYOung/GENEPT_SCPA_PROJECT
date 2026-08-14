from pathlib import Path
import unittest

from gene_embedding_project.genept_scpa.config import ConfigError, load_config


ROOT = Path(__file__).resolve().parents[1]


class ProtocolConfigTests(unittest.TestCase):
    def setUp(self):
        self.config = load_config(ROOT / "config/genept_scpa.yaml")

    def test_phase_seven_is_active_after_phase_six_pass(self):
        self.assertEqual(self.config.active_phase, 7)
        self.assertEqual(self.config.max_phase_allowed, 7)
        self.assertEqual(self.config.values["phase0"]["status"], "passed")
        self.assertEqual(self.config.values["phase1"]["status"], "passed")
        self.assertEqual(
            self.config.values["phase1"]["dataset_gate_status"],
            "passed",
        )
        self.assertEqual(
            self.config.values["phase1"]["stage"],
            "phase1b_passed",
        )
        self.config.require_phase(1)
        self.config.require_phase(2)
        self.config.require_phase(3)
        self.config.require_phase(4)
        self.config.require_phase(5)
        self.config.require_phase(6)
        self.config.require_phase(7)
        with self.assertRaises(ConfigError):
            self.config.require_phase(8)

    def test_phase_one_uses_full_official_geo_object(self):
        dataset = self.config.values["phase1"]["dataset"]
        self.assertEqual(dataset["accession"], "GSE212270")
        self.assertFalse(dataset["reduced"])
        self.assertEqual(dataset["required_hours"], [0, 12, 24])
        self.assertEqual(
            dataset["filename"],
            "GSE212270_integrated_naive_cd4.rds.gz",
        )
        self.assertTrue(dataset["url"].startswith("https://ftp.ncbi.nlm.nih.gov/geo/"))

    def test_download_script_uses_only_the_pinned_geo_dataset(self):
        script = (ROOT / "scripts/data/download_phase1_data.sh").read_text(
            encoding="utf-8"
        )
        self.assertIn("ftp.ncbi.nlm.nih.gov/geo/", script)
        self.assertNotIn("drive.google.com", script)
        self.assertNotIn("raw.githubusercontent.com", script)

    def test_phase_one_protocol_defaults_are_frozen(self):
        scpa = self.config.values["phase1"]["scpa"]
        self.assertEqual(scpa["seed"], 20260810)
        self.assertEqual(scpa["downsample"], 500)
        self.assertEqual(scpa["min_genes"], 15)
        self.assertEqual(scpa["max_genes"], 500)

    def test_phase_one_b_protocol_is_frozen(self):
        phase1b = self.config.values["phase1"]["phase1b"]
        self.assertEqual(phase1b["metadata_column"], "Hour")
        self.assertEqual(phase1b["time_values"], [0, 12, 24])
        self.assertEqual(phase1b["expression"]["assay"], "RNA")
        self.assertEqual(phase1b["expression"]["layer_or_slot"], "data")
        self.assertEqual(phase1b["pathways"]["input_pathway_count"], 243)
        self.assertEqual(
            phase1b["pathways"]["sha256"],
            "6bc5977da3fa60f86d5ffb59fc938740bf418fa4d976182a314d65479eb8b744",
        )
        self.assertEqual(
            set(phase1b["analyses"]),
            {"global", "0_vs_12", "12_vs_24", "0_vs_24", "reference"},
        )
        self.assertEqual(phase1b["qval_convention"], "larger_is_stronger")
        reference = phase1b["analyses"]["reference"]
        self.assertEqual(reference["population_1"], {"Cell_Type": "Resting", "Hour": 0})
        self.assertEqual(reference["population_2"], {"Cell_Type": "Activated", "Hour": 24})
        self.assertFalse(reference["parameter_tuning"])
        self.assertEqual(phase1b["visualization"]["heatmap_top_n"], 30)
        self.assertEqual(
            phase1b["visualization"]["comparison_scope"],
            "qualitative_only",
        )

    def test_phase_two_published_genept_w_protocol_is_frozen(self):
        phase2 = self.config.values["phase2"]
        self.assertEqual(phase2["status"], "passed")
        self.assertEqual(phase2["dataset"], "naive_cd4")
        self.assertEqual(phase2["embedding"]["model"], "text-embedding-ada-002")
        self.assertEqual(phase2["embedding"]["dimension"], 1536)
        self.assertFalse(phase2["embedding"]["generate_via_openai_api"])
        self.assertEqual(phase2["expression"], {
            "assay": "RNA",
            "layer": "counts",
            "sparse": True,
            "integrated_assay_allowed": False,
            "scale_data_allowed": False,
        })
        self.assertEqual(phase2["preprocessing"]["normalization_target"], 10000)
        self.assertEqual(phase2["preprocessing"]["log_transform"], "log1p")
        self.assertEqual(phase2["preprocessing"]["final_normalization"], "rowwise_unit_l2")

    def test_phase_three_primary_protocol_is_frozen(self):
        phase3 = self.config.values["phase3"]
        self.assertEqual(phase3["status"], "passed")
        self.assertEqual(phase3["genept"]["model"], "text-embedding-ada-002")
        self.assertTrue(phase3["genept"]["reuse_phase2_artifact"])
        self.assertEqual(phase3["primary_comparison"]["hour"], 0)
        self.assertEqual(phase3["primary_comparison"]["sample_size"], 500)
        self.assertEqual(phase3["original_expression"]["expected_shared_genes"], 17085)
        self.assertEqual(phase3["scpa_core"]["implementation"], "multicross::mcm")
        self.assertFalse(phase3["scpa_core"]["standard_pathway_analysis"])

    def test_phase_four_pathway_protocol_is_frozen(self):
        phase4 = self.config.values["phase4"]
        self.assertEqual(phase4["status"], "passed")
        self.assertEqual(phase4["stage"], "cd4_activation_passed")
        self.assertEqual(phase4["primary_comparison"]["seed"], 20260810)
        self.assertEqual(phase4["primary_comparison"]["sample_size"], 500)
        self.assertTrue(phase4["primary_comparison"]["canonical_sampling_reused"])
        self.assertEqual(phase4["expression"]["normalization_target"], 10000)
        self.assertFalse(phase4["expression"]["pathway_renormalization"])
        self.assertEqual(phase4["pathways"]["min_genes"], 15)
        self.assertEqual(phase4["pathways"]["max_genes"], 500)
        self.assertEqual(phase4["pathways"]["preflight_input_count"], 243)
        self.assertEqual(phase4["pathways"]["preflight_eligible_paired_count"], 123)
        self.assertTrue(phase4["paired_gene_policy"]["identical_between_branches"])
        self.assertEqual(phase4["genept_projection"]["formula"], "X_P @ E_P")
        self.assertFalse(phase4["genept_projection"]["primary_l2_normalization"])
        extension = phase4["validation_extension"]
        self.assertEqual(extension["source_status"], "ALREADY_PRESENT_IN_RDS")
        self.assertEqual(extension["comparison_count"], 3)
        self.assertEqual(extension["default_comparison_set"], "cd4_activation")
        self.assertTrue(extension["preserved_all_9_capability"])
        self.assertFalse(extension["cd8_current_production"])
        self.assertEqual(extension["pathway_count"], 123)
        self.assertEqual(extension["qval_log_base"], 10)
        self.assertEqual(extension["tie_ranking"], "average")

    def test_phase_five_gene_masking_protocol_is_frozen(self):
        phase5 = self.config.values["phase5"]
        self.assertEqual(phase5["status"], "passed")
        self.assertEqual(phase5["stage"], "passed_after_gpt_review")
        self.assertEqual(phase5["expected_target_count"], 30)
        self.assertEqual(phase5["expected_by_comparison"], {
            "cd4_0h_vs_12h": 11,
            "cd4_12h_vs_24h": 9,
            "cd4_0h_vs_24h": 10,
        })
        self.assertTrue(phase5["masking"]["same_gene_between_branches"])
        self.assertTrue(phase5["masking"]["genept_l2_deferred"])
        self.assertEqual(float(phase5["contribution_metric"]["raw_p_clip"]), 1e-300)

    def test_phase_six_semantic_control_protocol_is_frozen(self):
        phase6 = self.config.values["phase6"]
        self.assertEqual(phase6["status"], "passed")
        self.assertEqual(phase6["stage"], "passed_after_scientific_review")
        self.assertEqual(phase6["targets"]["pathway_control_count"], 30)
        self.assertEqual(phase6["targets"]["gene_control_representative_count"], 6)
        self.assertFalse(phase6["representations"]["l2_normalization"])
        self.assertTrue(phase6["representations"]["same_realization_for_both_groups_full_and_masks"])
        self.assertEqual(phase6["replicates"], {
            "pathway_per_control": 100,
            "gene_per_control": 20,
            "paired_resampling": 10,
        })
        self.assertEqual(phase6["expected_mcm"]["total_including_control_baselines"], 19100)

    def test_phase_seven_synthetic_protocol_is_frozen_and_mock_only(self):
        phase7 = self.config.values["phase7"]
        self.assertEqual(phase7["status"], "in_progress")
        self.assertEqual(phase7["source_population"], "naive_cd4_0h")
        self.assertEqual(phase7["pseudo_condition_cells"], 500)
        self.assertEqual(phase7["pathway_universe"], "frozen_phase4_paired_pathways")
        self.assertFalse(phase7["genept_primary_l2"])
        self.assertEqual(phase7["llm_backend_current"], "mock_only")
        self.assertEqual(
            phase7["real_backend_adapter"],
            "transformers_gpt_oss_mxfp4_v1_not_activated",
        )


if __name__ == "__main__":
    unittest.main()
