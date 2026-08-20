from pathlib import Path

import json
import numpy as np
import pytest
import yaml

from scripts.phase7b.evaluate_null_calibration import (
    ap_from_gene_indices,
    build_calibrated_units,
    load_baseline_features,
    rank_biserial,
)


ROOT = Path(__file__).resolve().parents[1]


def _metric(pathway, scenario, draw, strength, method, ap, fallback=False):
    return {
        "experiment_id": f"{pathway}_{scenario}_{draw}_{strength}",
        "draw_id": str(draw),
        "pathway": pathway,
        "source_database": "TEST",
        "perturbation_type": scenario,
        "perturbation_strength": str(strength),
        "perturbation_seed": "1",
        "method": method,
        "gene_count": "20",
        "truth_fallback_used": str(fallback),
        "truth_count": "3",
        "average_precision": str(ap),
    }


def test_method_specific_null_is_matched_before_did():
    rows = []
    for method, null_ap in (("vanilla_scpa", 0.20), ("genept_scpa", 0.30)):
        rows.append(_metric("P", "null", 1, 0.0, method, null_ap))
    for scenario in ("mean_shift", "cell_subset", "mixed_direction"):
        for method, values in (("vanilla_scpa", (0.40, 0.60)),
                               ("genept_scpa", (0.70, 0.90))):
            rows.extend(_metric("P", scenario, 1, strength, method, ap)
                        for strength, ap in zip((0.5, 1.0), values))
    units = build_calibrated_units(rows)
    assert len(units) == 3
    assert units[0]["vanilla_raw_ap"] == pytest.approx(0.50)
    assert units[0]["genept_raw_ap"] == pytest.approx(0.80)
    assert units[0]["raw_genept_minus_vanilla_ap"] == pytest.approx(0.30)
    assert units[0]["null_genept_minus_vanilla_ap"] == pytest.approx(0.10)
    assert units[0]["calibrated_did"] == pytest.approx(0.20)


def test_rank_biserial_and_ap_from_fixed_ranking_are_directional():
    assert rank_biserial([1, 2, 3]) == pytest.approx(1.0)
    assert rank_biserial([-1, -2, -3]) == pytest.approx(-1.0)
    genes = ["A", "B", "C", "D"]
    ranks = {"A": 1, "B": 2, "C": 3, "D": 4}
    assert ap_from_gene_indices(np.asarray([0, 2]), genes, ranks) == pytest.approx((1 + 2 / 3) / 2)


def test_frozen_A_has_exactly_nine_feasible_and_two_not_estimable_pathways():
    config = yaml.safe_load((ROOT / "config/phase7b_null_calibration.yaml").read_text())
    manifest = json.loads((ROOT / config["source"]["manifest"]).read_text())
    baseline = load_baseline_features(
        ROOT / config["source"]["expression_h5"], manifest,
        detection_min=config["null_truth_diagnostics"]["eligibility"]["detection_fraction_min"],
    )
    feasible = [name for name, info in baseline.items()
                if len(info["eligible_indices"]) >= info["truth_count"]]
    infeasible = [name for name in baseline if name not in feasible]
    assert len(feasible) == 9
    assert set(infeasible) == set(config["null_truth_diagnostics"]["infeasible_pathways"])


def test_phase7b_output_is_isolated_and_new_mcm_is_forbidden():
    config = yaml.safe_load((ROOT / "config/phase7b_null_calibration.yaml").read_text())
    output = (ROOT / config["artifacts"]["processed_directory"]).resolve()
    phase7 = (ROOT / "data/processed/genept_scpa/phase7_llmfree_synthetic").resolve()
    assert output != phase7
    assert phase7 not in output.parents
    assert config["execution"]["new_mcm_forbidden"] is True
    assert config["execution"]["reuse_existing_rankings_only"] is True
    assert config["null_truth_diagnostics"]["A_ELIGIBLE_POOL_MATCHED"]["insufficient_pool_status"] == "NOT_ESTIMABLE"
