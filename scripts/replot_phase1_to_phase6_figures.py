#!/usr/bin/env python3
"""Regenerate Phase 4-6 figures from existing result tables only.

Phase 1 uses ``scripts/scpa/replot_phase1b_figures.R``. Phases 2 and 3 do not
currently have dedicated scientific-result figures in the repository.
"""

from __future__ import annotations

import csv
from pathlib import Path
import sys
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from scripts.phase4.run_pathway_comparison import create_figures  # noqa: E402
from scripts.phase4.run_timecourse_validation import (  # noqa: E402
    create_validation_figures,
    load_phase4a_historical_reference,
)
from scripts.phase5.run_gene_contribution import (  # noqa: E402
    create_figures as create_phase5_figures,
    representative_targets,
)
from scripts.phase6.run_semantic_controls import (  # noqa: E402
    create_figures as create_phase6_figures,
)


def parse_value(value: str) -> Any:
    if value in {"TRUE", "True"}:
        return True
    if value in {"FALSE", "False"}:
        return False
    if value in {"", "NA", "NaN"}:
        return None
    try:
        return float(value)
    except ValueError:
        return value


def read_typed_csv(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [
            {key: parse_value(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]


def main() -> int:
    processed = PROJECT_ROOT / "data/processed/genept_scpa"

    phase4a = processed / "phase4"
    phase4a_rows = read_typed_csv(phase4a / "vanilla_vs_genept_pathway_comparison.csv")
    phase4a_files = create_figures(phase4a_rows, phase4a / "figures")
    print(f"[Replot] Phase 4A: {len(phase4a_files)} figures", flush=True)

    phase4b = processed / "phase4_cd4_activation"
    phase4b_files = create_validation_figures(
        read_typed_csv(phase4b / "phase4_cd4_activation_overview.csv"),
        read_typed_csv(phase4b / "phase4_cd4_activation_detection_states.csv"),
        read_typed_csv(phase4b / "phase4_cd4_activation_all_results.csv"),
        phase4b / "figures",
        load_phase4a_historical_reference(),
    )
    print(f"[Replot] Phase 4B: {len(phase4b_files)} figures (L2 omitted)", flush=True)

    phase5 = processed / "phase5_gene_contribution"
    targets = read_typed_csv(phase5 / "phase5_target_pathways.csv")
    phase5_files = create_phase5_figures(
        targets,
        read_typed_csv(phase5 / "phase5_gene_masking_all_results.csv"),
        read_typed_csv(phase5 / "phase5_pathway_summary.csv"),
        representative_targets(targets),
        phase5 / "figures",
    )
    print(f"[Replot] Phase 5: {len(phase5_files)} figures", flush=True)

    phase6 = processed / "phase6_semantic_controls"
    phase6_files = create_phase6_figures(
        read_typed_csv(phase6 / "phase6_pathway_control_summary.csv"),
        read_typed_csv(phase6 / "phase6_gene_control_all_results.csv"),
        read_typed_csv(phase6 / "phase6_gene_control_summary.csv"),
        read_typed_csv(phase6 / "phase6_resampling_results.csv"),
    )
    print(f"[Replot] Phase 6: {len(phase6_files)} figures", flush=True)
    print("REPLOT_PHASE1_TO_PHASE6 status=PASS analysis_rerun=false", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
