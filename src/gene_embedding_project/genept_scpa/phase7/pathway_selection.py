"""Outcome-independent Phase 7 pathway selection and description handling."""

from __future__ import annotations

import re
from typing import Any, Mapping, Sequence

import numpy as np


def sanitize_gene_description(gene: str, description: str) -> str:
    """Remove the official label and redact later primary-symbol mentions."""

    if not isinstance(description, str):
        return ""
    pattern = re.compile(
        rf"^\s*Gene\s+Symbol\s+{re.escape(gene)}\b[\s:;,-]*",
        flags=re.IGNORECASE,
    )
    sanitized = pattern.sub("", description, count=1).strip()
    symbol_pattern = re.compile(
        rf"(?<![A-Za-z0-9]){re.escape(gene)}(?![A-Za-z0-9])",
        flags=re.IGNORECASE,
    )
    return symbol_pattern.sub("this gene", sanitized).strip()


def description_is_usable(gene: str, descriptions: Mapping[str, str]) -> bool:
    return gene in descriptions and bool(sanitize_gene_description(gene, descriptions[gene]))


def _size_bin(count: int, bins: Mapping[str, Sequence[int]]) -> str | None:
    for name, bounds in bins.items():
        if len(bounds) != 2:
            raise ValueError(f"Invalid size bin {name}")
        if int(bounds[0]) <= count <= int(bounds[1]):
            return name
    return None


def select_phase7_pathways(
    phase4_manifest: Mapping[str, Any],
    descriptions: Mapping[str, str],
    *,
    size_bins: Mapping[str, Sequence[int]],
    source_bin_quota: Mapping[str, Mapping[str, int]],
    seed: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Select only by frozen pathway metadata, description coverage and seed."""

    pathways = phase4_manifest.get("pathways")
    if not isinstance(pathways, list) or not pathways:
        raise ValueError("Frozen Phase 4 manifest has no pathways")
    eligible: dict[tuple[str, str], list[dict[str, Any]]] = {}
    excluded = {"size": 0, "description": 0, "unrequested_source": 0}
    for pathway in pathways:
        genes = tuple(pathway["paired_genes"])
        bin_name = _size_bin(len(genes), size_bins)
        if bin_name is None:
            excluded["size"] += 1
            continue
        if not all(description_is_usable(gene, descriptions) for gene in genes):
            excluded["description"] += 1
            continue
        source = str(pathway["source_database"])
        if source not in source_bin_quota:
            excluded["unrequested_source"] += 1
            continue
        record = {
            "pathway": str(pathway["pathway"]),
            "source_database": source,
            "size_bin": bin_name,
            "analysis_genes": list(genes),
            "analysis_gene_count": len(genes),
            "original_pathway_genes": list(pathway.get("original_pathway_genes", genes)),
        }
        eligible.setdefault((source, bin_name), []).append(record)

    generator = np.random.default_rng(seed)
    selected: list[dict[str, Any]] = []
    counts: dict[str, int] = {}
    for source, bin_quotas in source_bin_quota.items():
        for bin_name, requested_value in bin_quotas.items():
            requested = int(requested_value)
            candidates = sorted(
                eligible.get((source, bin_name), []), key=lambda row: row["pathway"]
            )
            counts[f"{source}:{bin_name}:eligible"] = len(candidates)
            counts[f"{source}:{bin_name}:selected"] = requested
            if requested > len(candidates):
                raise ValueError(
                    f"Pathway quota {source}/{bin_name} requests {requested}; "
                    f"only {len(candidates)} eligible"
                )
            if requested:
                indices = generator.choice(len(candidates), size=requested, replace=False)
                selected.extend(candidates[int(index)] for index in indices)
    selected.sort(key=lambda row: (row["source_database"], row["size_bin"], row["pathway"]))
    if len({row["pathway"] for row in selected}) != len(selected):
        raise RuntimeError("Duplicate selected pathways")
    audit = {
        "selection_seed": int(seed),
        "phase4_pathway_count": len(pathways),
        "selected_pathway_count": len(selected),
        "selection_uses_phase4_or_phase5_scores": False,
        "counts": counts,
        "excluded": excluded,
    }
    return selected, audit
