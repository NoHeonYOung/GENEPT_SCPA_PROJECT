"""Strict, dependency-free Phase 7 artifact schema validation."""

from __future__ import annotations

import re
from typing import Any, Mapping, Sequence


LLM_INPUT_SCHEMA_VERSION = "phase7.llm-input.v1"
LLM_RANKING_SCHEMA_VERSION = "phase7.llm-ranking.v1"
SYNTHETIC_SCHEMA_VERSION = "phase7.synthetic.v1"
CANDIDATE_PATTERN = re.compile(r"^C[0-9]{3,}$")
SUMMARY_FIELDS = (
    "mean", "sd", "median", "q10", "q25", "q75", "q90", "nonzero_fraction"
)
FORBIDDEN_LLM_KEYS = {
    "gene", "gene_symbol", "ground_truth", "is_ground_truth_perturbed",
    "perturbation_type", "perturbation_scenario", "perturbation_strength",
    "selected_perturbed_genes", "p_value", "fold_change", "fold_change_rank",
}


def _require_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a mapping")
    return value


def _forbidden_keys(value: Any, location: str = "root") -> list[str]:
    found: list[str] = []
    if isinstance(value, Mapping):
        for key, child in value.items():
            if str(key).lower() in FORBIDDEN_LLM_KEYS:
                found.append(f"{location}.{key}")
            found.extend(_forbidden_keys(child, f"{location}.{key}"))
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        for index, child in enumerate(value):
            found.extend(_forbidden_keys(child, f"{location}[{index}]"))
    return found


def validate_llm_request(payload: Mapping[str, Any]) -> None:
    request = _require_mapping(payload, "LLM request")
    if request.get("schema_version") != LLM_INPUT_SCHEMA_VERSION:
        raise ValueError("Unexpected LLM input schema version")
    for key in ("run_id", "experiment_id", "backend", "prompt_condition"):
        if not isinstance(request.get(key), str) or not request[key]:
            raise ValueError(f"LLM request field {key} must be a non-empty string")
    if request["prompt_condition"] not in {
        "stats_only", "true_description", "shuffled_description"
    }:
        raise ValueError("Unknown prompt condition")
    pathway = _require_mapping(request.get("pathway"), "pathway")
    if not isinstance(pathway.get("name"), str) or not pathway["name"]:
        raise ValueError("Pathway name is required")
    candidates = request.get("candidates")
    if not isinstance(candidates, list) or not candidates:
        raise ValueError("LLM request candidates must be a non-empty list")
    expected = int(pathway.get("candidate_count", -1))
    if expected != len(candidates):
        raise ValueError("Pathway candidate_count does not match candidates")
    identifiers: list[str] = []
    for candidate in candidates:
        item = _require_mapping(candidate, "candidate")
        identifier = item.get("candidate_id")
        if not isinstance(identifier, str) or not CANDIDATE_PATTERN.fullmatch(identifier):
            raise ValueError(f"Invalid opaque candidate ID: {identifier!r}")
        identifiers.append(identifier)
        for condition in ("condition_a", "condition_b"):
            summary = _require_mapping(item.get(condition), condition)
            if set(summary) != set(SUMMARY_FIELDS):
                raise ValueError(f"{condition} summary fields are not frozen")
            if not all(isinstance(summary[field], (int, float)) for field in SUMMARY_FIELDS):
                raise ValueError(f"{condition} summaries must be numeric")
        has_description = "description" in item
        if request["prompt_condition"] == "stats_only" and has_description:
            raise ValueError("stats_only candidates must not contain descriptions")
        if request["prompt_condition"] != "stats_only":
            if not isinstance(item.get("description"), str) or not item["description"].strip():
                raise ValueError("Description condition requires non-empty descriptions")
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("Duplicate candidate IDs")
    forbidden = _forbidden_keys(request)
    if forbidden:
        raise ValueError(f"LLM-visible request contains forbidden keys: {forbidden}")


def validate_llm_response(
    payload: Mapping[str, Any], *, expected_candidate_ids: Sequence[str] | None = None
) -> None:
    response = _require_mapping(payload, "LLM response")
    if response.get("schema_version") != LLM_RANKING_SCHEMA_VERSION:
        raise ValueError("Unexpected LLM ranking schema version")
    for key in ("pathway", "run_id", "backend"):
        if not isinstance(response.get(key), str) or not response[key]:
            raise ValueError(f"LLM response field {key} must be a non-empty string")
    ranking = response.get("ranking")
    if not isinstance(ranking, list) or not ranking:
        raise ValueError("LLM ranking must be non-empty")
    ids: list[str] = []
    ranks: list[int] = []
    for row in ranking:
        item = _require_mapping(row, "ranking row")
        if set(item) != {"candidate_id", "rank"}:
            raise ValueError("Ranking rows may contain only candidate_id and rank")
        candidate = item["candidate_id"]
        rank = item["rank"]
        if not isinstance(candidate, str) or not CANDIDATE_PATTERN.fullmatch(candidate):
            raise ValueError("Invalid ranked candidate ID")
        if not isinstance(rank, int) or isinstance(rank, bool):
            raise ValueError("Ranks must be integers")
        ids.append(candidate)
        ranks.append(rank)
    if len(ids) != len(set(ids)) or sorted(ranks) != list(range(1, len(ids) + 1)):
        raise ValueError("Ranking must contain unique candidates and ranks 1..N")
    if expected_candidate_ids is not None and set(ids) != set(expected_candidate_ids):
        raise ValueError("Response candidate set differs from the request")


def validate_synthetic_manifest(payload: Mapping[str, Any]) -> None:
    manifest = _require_mapping(payload, "synthetic manifest")
    if manifest.get("schema_version") != SYNTHETIC_SCHEMA_VERSION:
        raise ValueError("Unexpected synthetic manifest schema version")
    if manifest.get("ground_truth_embedded") is not False:
        raise ValueError("Ground truth must remain separate from the method manifest")
    experiments = manifest.get("experiments")
    if not isinstance(experiments, list) or not experiments:
        raise ValueError("Synthetic manifest requires experiments")
    for experiment in experiments:
        item = _require_mapping(experiment, "experiment")
        forbidden = {"ground_truth_genes", "perturbed_genes", "truth_indices"} & set(item)
        if forbidden:
            raise ValueError("Method manifest must not embed ground-truth labels")
