"""Phase 7 LLM backend interface and deterministic mock implementation."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Mapping

from .schemas import (
    LLM_RANKING_SCHEMA_VERSION,
    validate_llm_request,
    validate_llm_response,
)


class LLMBackend(ABC):
    """Interface shared by mock and any future explicitly approved backend."""

    name: str
    scientific_evaluation_allowed: bool = False

    @abstractmethod
    def rank(self, request: Mapping[str, Any]) -> dict[str, Any]:
        raise NotImplementedError


class MockLLMBackend(LLMBackend):
    """Deterministic schema test backend; never a scientific method result."""

    name = "mock_expression_difference_v1"
    scientific_evaluation_allowed = False

    def rank(self, request: Mapping[str, Any]) -> dict[str, Any]:
        validate_llm_request(request)
        if request["backend"] != self.name:
            raise ValueError("Request backend does not match MockLLMBackend")
        ordered = sorted(
            request["candidates"],
            key=lambda item: (
                -abs(float(item["condition_b"]["mean"]) - float(item["condition_a"]["mean"])),
                item["candidate_id"],
            ),
        )
        response = {
            "schema_version": LLM_RANKING_SCHEMA_VERSION,
            "pathway": request["pathway"]["name"],
            "run_id": request["run_id"],
            "backend": self.name,
            "ranking": [
                {"candidate_id": item["candidate_id"], "rank": index + 1}
                for index, item in enumerate(ordered)
            ],
        }
        validate_llm_response(
            response,
            expected_candidate_ids=[item["candidate_id"] for item in request["candidates"]],
        )
        return response
