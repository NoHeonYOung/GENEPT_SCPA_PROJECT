"""Deterministic, disjoint pseudo-condition sampling for Phase 7."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class PseudoConditionSplit:
    seed: int
    condition_a_ids: tuple[str, ...]
    condition_b_ids: tuple[str, ...]

    @property
    def all_ids(self) -> tuple[str, ...]:
        return self.condition_a_ids + self.condition_b_ids


def split_pseudo_conditions(
    source_cell_ids: Iterable[str], *, cells_per_condition: int = 500, seed: int
) -> PseudoConditionSplit:
    """Select 2N unique cells once and split them into disjoint A/B groups."""

    if cells_per_condition < 2:
        raise ValueError("Each pseudo-condition needs at least two cells")
    cells = sorted(str(cell_id) for cell_id in source_cell_ids)
    if any(not cell_id for cell_id in cells):
        raise ValueError("Source cell IDs must be non-empty")
    if len(cells) != len(set(cells)):
        raise ValueError("Source cell IDs must be unique")
    total = cells_per_condition * 2
    if len(cells) < total:
        raise ValueError(f"Need {total} source cells; found {len(cells)}")
    generator = np.random.default_rng(seed)
    chosen = generator.choice(len(cells), size=total, replace=False)
    selected = np.asarray([cells[int(index)] for index in chosen], dtype=object)
    selected = selected[generator.permutation(total)]
    condition_a = tuple(str(value) for value in selected[:cells_per_condition])
    condition_b = tuple(str(value) for value in selected[cells_per_condition:])
    if set(condition_a) & set(condition_b):
        raise RuntimeError("Pseudo-condition split unexpectedly overlaps")
    return PseudoConditionSplit(int(seed), condition_a, condition_b)
