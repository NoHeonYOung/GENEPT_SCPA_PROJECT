import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from scripts.phase8.evaluate_mean_shift_mechanism import (
    build_calibrated_units,
    contrast_statistics,
)
from scripts.phase8.prepare_mean_shift_controls import (
    make_derangement,
    make_random_projection,
)


ROOT = Path(__file__).resolve().parents[1]


def test_derangement_is_deterministic_fixed_point_free_and_complete():
    first = make_derangement(17, 12345)
    second = make_derangement(17, 12345)
    np.testing.assert_array_equal(first, second)
    assert not np.any(first == np.arange(17))
    np.testing.assert_array_equal(np.sort(first), np.arange(17))


def test_random_projection_is_deterministic_1536d_and_norm_matched():
    rng = np.random.default_rng(9)
    true = rng.normal(size=(12, 1536))
    first = make_random_projection(true, 77)
    second = make_random_projection(true, 77)
    assert first.shape == (12, 1536)
    np.testing.assert_array_equal(first, second)
    np.testing.assert_allclose(np.linalg.norm(first, axis=1), np.linalg.norm(true, axis=1), atol=1e-10)
    assert np.isfinite(first).all()


def test_frozen_manifest_workload_is_exact():
    config = yaml.safe_load((ROOT / "config/phase8_mean_shift_mechanism.yaml").read_text())
    manifest = json.loads((ROOT / config["source"]["phase7_manifest"]).read_text())
    sum_k = sum(row["analysis_gene_count"] for row in manifest["pathways"])
    sum_k_plus_one = sum(row["analysis_gene_count"] + 1 for row in manifest["pathways"])
    selected = [row for row in manifest["experiments"]
                if row["perturbation_type"] in {"null", "mean_shift"}]
    assert len(manifest["pathways"]) == 11
    assert sum_k == 353
    assert sum_k_plus_one == 364
    assert len(selected) == 660
    assert 20 * 3 * sum_k_plus_one == config["workload"]["expected_mcm_per_control"] == 21840
    assert config["workload"]["expected_total_new_mcm"] == 43680


def test_calibration_and_contrast_algebra_on_full_toy_unit_grid():
    metrics = []
    base = {
        "vanilla": (0.20, 0.40, 0.60),
        "true_genept": (0.30, 0.70, 0.90),
        "permuted_genept": (0.25, 0.55, 0.75),
        "random_projection": (0.22, 0.48, 0.68),
    }
    for pathway_index in range(11):
        pathway = f"P{pathway_index:02d}"
        for draw in range(1, 21):
            for representation, values in base.items():
                for scenario, strength, ap in (
                    ("null", 0.0, values[0]),
                    ("mean_shift", 0.5, values[1]),
                    ("mean_shift", 1.0, values[2]),
                ):
                    metrics.append({
                        "pathway": pathway, "draw_id": draw,
                        "perturbation_type": scenario, "perturbation_strength": strength,
                        "truth_fallback_used": False, "representation": representation,
                        "average_precision": ap,
                    })
    units = build_calibrated_units(metrics)
    assert len(units) == 220
    first = units[0]
    assert first["vanilla_calibrated_ap"] == pytest.approx(0.30)
    assert first["true_genept_calibrated_ap"] == pytest.approx(0.50)
    assert first["permuted_genept_calibrated_ap"] == pytest.approx(0.40)
    assert first["random_projection_calibrated_ap"] == pytest.approx(0.36)
    assert first["C_mapping"] == pytest.approx(0.10)
    assert first["C_geometry"] == pytest.approx(0.04)
    assert first["C_projection"] == pytest.approx(0.06)
    assert first["C_total"] == pytest.approx(0.20)
    statistics = {row["contrast"]: row for row in contrast_statistics(units)}
    assert statistics["C_total"]["mean_contrast"] == pytest.approx(0.20)
    assert statistics["C_total"]["rank_biserial_positive_is_left"] == pytest.approx(1.0)


def test_phase8_paths_are_isolated_and_old_methods_have_zero_new_mcm():
    config = yaml.safe_load((ROOT / "config/phase8_mean_shift_mechanism.yaml").read_text())
    output = (ROOT / config["artifacts"]["processed_directory"]).resolve()
    assert output != (ROOT / "data/processed/genept_scpa/phase7_llmfree_synthetic").resolve()
    assert output != (ROOT / "data/processed/genept_scpa/phase7b_null_calibration").resolve()
    assert config["representations"]["vanilla"]["new_mcm"] == 0
    assert config["representations"]["true_genept"]["new_mcm"] == 0
    assert config["control_pairing"]["truth_labels_used_during_generation"] is False
