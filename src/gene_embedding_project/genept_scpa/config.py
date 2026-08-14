"""Load and validate the frozen GenePT × SCPA experiment configuration."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping
from urllib.parse import urlparse

import yaml


class ConfigError(ValueError):
    """Raised when the protocol configuration is incomplete or inconsistent."""


@dataclass(frozen=True)
class ProtocolConfig:
    """Validated protocol settings with explicit phase gating."""

    path: Path
    values: Mapping[str, Any]

    @property
    def active_phase(self) -> int:
        return int(self.values["project"]["active_phase"])

    @property
    def max_phase_allowed(self) -> int:
        return int(self.values["project"]["max_phase_allowed"])

    def require_phase(self, phase: int) -> None:
        """Refuse work that has not been unlocked in the frozen protocol."""
        if phase > self.max_phase_allowed:
            raise ConfigError(
                f"Phase {phase} is locked; max_phase_allowed="
                f"{self.max_phase_allowed}. Update the decision log before unlocking it."
            )


def _require(mapping: Mapping[str, Any], key: str, context: str) -> Any:
    if key not in mapping:
        raise ConfigError(f"Missing required key: {context}.{key}")
    return mapping[key]


def load_config(path: str | Path = "config/genept_scpa.yaml") -> ProtocolConfig:
    """Read YAML, validate protocol invariants, and return a phase-aware config."""
    config_path = Path(path)
    with config_path.open("r", encoding="utf-8") as handle:
        values = yaml.safe_load(handle)

    if not isinstance(values, Mapping):
        raise ConfigError("The configuration root must be a mapping")

    schema_version = _require(values, "schema_version", "root")
    if schema_version != 1:
        raise ConfigError(f"Unsupported schema_version: {schema_version!r}")

    project = _require(values, "project", "root")
    if not isinstance(project, Mapping):
        raise ConfigError("project must be a mapping")

    active_phase = int(_require(project, "active_phase", "project"))
    max_phase_allowed = int(
        _require(project, "max_phase_allowed", "project")
    )
    if not 0 <= active_phase <= 10:
        raise ConfigError("project.active_phase must be between 0 and 10")
    if not active_phase <= max_phase_allowed <= 10:
        raise ConfigError(
            "project.max_phase_allowed must be >= active_phase and <= 10"
        )

    phase1 = _require(values, "phase1", "root")
    dataset = _require(phase1, "dataset", "phase1")
    for key in ("accession", "filename", "url", "reduced"):
        _require(dataset, key, "phase1.dataset")
    if dataset["accession"] != "GSE212270":
        raise ConfigError("Phase 1 dataset must remain pinned to GEO GSE212270")
    if dataset["reduced"] is not False:
        raise ConfigError("Phase 1 must use the non-reduced dataset")
    expected_filename = "GSE212270_integrated_naive_cd4.rds.gz"
    if dataset["filename"] != expected_filename:
        raise ConfigError(f"Phase 1 filename must remain pinned to {expected_filename}")
    dataset_url = urlparse(str(dataset["url"]))
    if dataset_url.scheme != "https" or dataset_url.hostname != "ftp.ncbi.nlm.nih.gov":
        raise ConfigError("Phase 1 download URL must use the official NCBI GEO HTTPS host")
    if not dataset_url.path.endswith(f"/GSE212270/suppl/{expected_filename}"):
        raise ConfigError("Phase 1 download URL does not match the pinned GEO file")

    phase1_status = _require(phase1, "status", "phase1")
    if active_phase == 1 and phase1_status != "in_progress":
        raise ConfigError("Active Phase 1 must remain in_progress until its PASS gate")
    if active_phase >= 2 and phase1_status != "passed":
        raise ConfigError("Phase 1 must be passed before Phase 2 can be active")

    scpa = _require(phase1, "scpa", "phase1")
    for key in ("seed", "downsample", "min_genes", "max_genes"):
        _require(scpa, key, "phase1.scpa")

    if active_phase >= 2:
        phase2 = _require(values, "phase2", "root")
        phase2_status = _require(phase2, "status", "phase2")
        expected_phase2_status = "in_progress" if active_phase == 2 else "passed"
        if phase2_status != expected_phase2_status:
            raise ConfigError(
                f"Phase 2 must be {expected_phase2_status} when Phase {active_phase} is active"
            )
        embedding = _require(phase2, "embedding", "phase2")
        if embedding.get("model") != "text-embedding-ada-002":
            raise ConfigError("Phase 2 primary model must remain text-embedding-ada-002")
        if int(embedding.get("dimension", 0)) != 1536:
            raise ConfigError("Phase 2 primary embedding dimension must remain 1536")
        expression = _require(phase2, "expression", "phase2")
        if expression.get("assay") != "RNA" or expression.get("layer") != "counts":
            raise ConfigError("Phase 2 must use the RNA/counts expression source")

    if active_phase == 3:
        phase3 = _require(values, "phase3", "root")
        if _require(phase3, "status", "phase3") != "in_progress":
            raise ConfigError("Active Phase 3 must be in_progress")
        genept = _require(phase3, "genept", "phase3")
        if genept.get("model") != "text-embedding-ada-002":
            raise ConfigError("Phase 3 must reuse the Phase 2 ada-002 artifact")
        if int(genept.get("dimension", 0)) != 1536:
            raise ConfigError("Phase 3 GenePT-w dimension must remain 1536")
        if genept.get("reuse_phase2_artifact") is not True:
            raise ConfigError("Phase 3 must reuse the Phase 2 GenePT artifact")
        comparison = _require(phase3, "primary_comparison", "phase3")
        if int(comparison.get("hour", -1)) != 0:
            raise ConfigError("Phase 3 primary comparison must remain 0h")
        if int(comparison.get("sample_size", 0)) != 500:
            raise ConfigError("Phase 3 sample size must remain 500/group")
        scpa_core = _require(phase3, "scpa_core", "phase3")
        if scpa_core.get("implementation") != "multicross::mcm":
            raise ConfigError("Phase 3 must use the core function called by SCPA")
        if scpa_core.get("standard_pathway_analysis") is not False:
            raise ConfigError("Phase 3 is not standard SCPA pathway analysis")

    if active_phase >= 4:
        phase3 = _require(values, "phase3", "root")
        if _require(phase3, "status", "phase3") != "passed":
            raise ConfigError("Phase 3 must be passed before Phase 4 can be active")
        phase4 = _require(values, "phase4", "root")
        if active_phase == 4 and _require(phase4, "status", "phase4") != "in_progress":
            raise ConfigError("Active Phase 4 must remain in_progress until review")
        if active_phase >= 5 and _require(phase4, "status", "phase4") != "passed":
            raise ConfigError("Phase 4 must be passed before Phase 5 can be active")
        comparison = _require(phase4, "primary_comparison", "phase4")
        if int(comparison.get("seed", -1)) != 20260810:
            raise ConfigError("Phase 4 must reuse the frozen Phase 3 seed")
        if int(comparison.get("sample_size", -1)) != 500:
            raise ConfigError("Phase 4 must reuse 500 canonical cells per group")
        if comparison.get("canonical_sampling_reused") is not True:
            raise ConfigError("Phase 4 must reuse Phase 3 canonical sampling")
        projection = _require(phase4, "genept_projection", "phase4")
        if projection.get("primary_l2_normalization") is not False:
            raise ConfigError("Phase 4 primary projection must remain non-L2")
        if projection.get("formula") != "X_P @ E_P":
            raise ConfigError("Phase 4 projection formula must remain X_P @ E_P")
        extension = _require(phase4, "validation_extension", "phase4")
        if int(extension.get("comparison_count", 0)) != 3:
            raise ConfigError("Phase 4B primary validation must contain three comparisons")
        if extension.get("default_comparison_set") != "cd4_activation":
            raise ConfigError("Phase 4B must default to the CD4 activation comparison set")
        if extension.get("comparisons") != [
            "cd4_0h_vs_12h", "cd4_12h_vs_24h", "cd4_0h_vs_24h"
        ]:
            raise ConfigError("Phase 4B comparison identities have changed")
        if extension.get("groups") != ["cd4_0h", "cd4_12h", "cd4_24h"]:
            raise ConfigError("Phase 4B primary groups must remain CD4 0h/12h/24h")
        if extension.get("preserved_all_9_capability") is not True:
            raise ConfigError("The previous nine-comparison capability must be preserved")
        if int(extension.get("pathway_count", 0)) != 123:
            raise ConfigError("Phase 4 validation must reuse the 123-pathway universe")
        if extension.get("qval_log_base") != 10:
            raise ConfigError("Official SCPA 1.6.2 qval uses base-10 logarithms")
        if extension.get("tie_ranking") != "average":
            raise ConfigError("Phase 4 validation must use average tie ranks")
        scpa_core = _require(phase3, "scpa_core", "phase3")
        if scpa_core.get("implementation") != "multicross::mcm":
            raise ConfigError("Phase 3 must use the core function called by SCPA")
        if scpa_core.get("standard_pathway_analysis") is not False:
            raise ConfigError("Phase 3 is not standard SCPA pathway analysis")

    if active_phase >= 5:
        phase5 = _require(values, "phase5", "root")
        expected_status = "in_progress" if active_phase == 5 else "passed"
        if _require(phase5, "status", "phase5") != expected_status:
            raise ConfigError(
                f"Phase 5 must be {expected_status} when Phase {active_phase} is active"
            )
        if int(phase5.get("expected_target_count", 0)) != 30:
            raise ConfigError("Phase 5 must freeze exactly 30 discordant targets")
        if phase5.get("expected_by_comparison") != {
            "cd4_0h_vs_12h": 11,
            "cd4_12h_vs_24h": 9,
            "cd4_0h_vs_24h": 10,
        }:
            raise ConfigError("Phase 5 target counts by comparison have changed")
        masking = _require(phase5, "masking", "phase5")
        if masking.get("same_gene_between_branches") is not True:
            raise ConfigError("Phase 5 must mask the same gene in both branches")
        if masking.get("genept_l2_deferred") is not True:
            raise ConfigError("GenePT L2 gene masking must remain deferred")
        metric = _require(phase5, "contribution_metric", "phase5")
        if float(metric.get("raw_p_clip", 0)) != 1e-300:
            raise ConfigError("Phase 5 raw-p clipping must remain 1e-300")

    if active_phase >= 6:
        phase6 = _require(values, "phase6", "root")
        expected_status = "in_progress" if active_phase == 6 else "passed"
        if _require(phase6, "status", "phase6") != expected_status:
            raise ConfigError(
                f"Phase 6 must be {expected_status} when Phase {active_phase} is active"
            )
        if int(phase6.get("expected_mcm", {}).get(
            "total_including_control_baselines", 0
        )) != 19100:
            raise ConfigError("Phase 6 expected MCM total must remain 19100")

    if active_phase >= 7:
        phase7 = _require(values, "phase7", "root")
        if _require(phase7, "status", "phase7") != "in_progress":
            raise ConfigError("Active Phase 7 must remain in_progress until review")
        if phase7.get("source_population") != "naive_cd4_0h":
            raise ConfigError("Phase 7 source population must remain naive CD4 0h")
        if int(phase7.get("pseudo_condition_cells", 0)) != 500:
            raise ConfigError("Phase 7 must use 500 cells per pseudo-condition")
        if phase7.get("pathway_universe") != "frozen_phase4_paired_pathways":
            raise ConfigError("Phase 7 must reuse the frozen Phase 4 pathway universe")
        if phase7.get("genept_primary_l2") is not False:
            raise ConfigError("Phase 7 primary GenePT projection must remain non-L2")
        if phase7.get("llm_backend_current") != "mock_only":
            raise ConfigError("Phase 7 scaffolding must remain mock-only")

    return ProtocolConfig(path=config_path.resolve(), values=values)
