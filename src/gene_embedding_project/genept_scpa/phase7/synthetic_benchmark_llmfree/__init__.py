"""Vanilla versus non-L2 GenePT synthetic recovery benchmark.

Metrics quantify recovery of injected genes in this frozen synthetic setup.
They do not establish biological superiority, causality, or transfer to real
unknown signals, and they answer a different question from Phase 6.
"""

from .cohort import PseudoConditionSplit, split_pseudo_conditions
from .metrics import average_precision, ndcg_at_k, recall_at_k
from .perturbation import PerturbationResult, inject_perturbation

__all__ = [
    "PseudoConditionSplit", "PerturbationResult", "average_precision",
    "inject_perturbation", "ndcg_at_k", "recall_at_k", "split_pseudo_conditions",
]
