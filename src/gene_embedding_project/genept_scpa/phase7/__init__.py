"""LLM-free Phase 7 synthetic ground-truth benchmark."""

from .synthetic_benchmark_llmfree.metrics import average_precision, ndcg_at_k, recall_at_k

__all__ = ["average_precision", "ndcg_at_k", "recall_at_k"]
