"""Deterministic, disjoint pseudo-condition sampling."""

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
    """Select 2N cells once, then reproducibly split them into disjoint A/B."""

    cells = sorted(str(cell_id) for cell_id in source_cell_ids)
    if cells_per_condition < 2:
        raise ValueError("Each pseudo-condition needs at least two cells")
    if any(not cell for cell in cells) or len(cells) != len(set(cells)):
        raise ValueError("Source cell IDs must be non-empty and unique")
    total = 2 * cells_per_condition
    if len(cells) < total:
        raise ValueError(f"Need {total} source cells; found {len(cells)}")
    rng = np.random.default_rng(seed)
    selected = np.asarray(cells, dtype=object)[rng.choice(len(cells), total, replace=False)]
    selected = selected[rng.permutation(total)]
    a = tuple(str(value) for value in selected[:cells_per_condition])
    b = tuple(str(value) for value in selected[cells_per_condition:])
    if set(a) & set(b):
        raise RuntimeError("Pseudo-condition split overlaps")
    return PseudoConditionSplit(int(seed), a, b)
