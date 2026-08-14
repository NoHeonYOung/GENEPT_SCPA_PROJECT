"""Lazy, fail-closed Transformers adapter for the frozen gpt-oss protocol."""

from __future__ import annotations

import hashlib
import importlib
import importlib.metadata
import json
import math
from pathlib import Path
import re
import time
from typing import Any, Mapping, Sequence

from gene_embedding_project.genept_scpa.io import write_json_atomic

from .llm_backend import LLMBackend
from .runtime import SUPPORTED_PRIMARY, build_runtime_report
from .schemas import LLM_RANKING_SCHEMA_VERSION, validate_llm_request, validate_llm_response


class InvalidLLMOutputError(RuntimeError):
    """Raised by the compatibility rank() method after the frozen retries fail."""


def _canonical_hash(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def deterministic_decoding_config(inference: Mapping[str, Any]) -> dict[str, Any]:
    decoding = inference["decoding"]
    if decoding.get("do_sample") is not False or int(decoding.get("num_beams", 0)) != 1:
        raise ValueError("Primary gpt-oss decoding must remain deterministic greedy")
    if decoding.get("temperature") is not None or decoding.get("top_p") is not None:
        raise ValueError("Primary gpt-oss decoding must not set temperature or top_p")
    return {
        "do_sample": False,
        "num_beams": 1,
        "max_new_tokens": int(decoding["max_new_tokens"]),
        "use_cache": bool(decoding["use_cache"]),
    }


def candidate_order_hash(request: Mapping[str, Any]) -> str:
    return _canonical_hash([row["candidate_id"] for row in request["candidates"]])


def prompt_token_hash(token_ids: Sequence[Sequence[int]]) -> str:
    return _canonical_hash([[int(token) for token in row] for row in token_ids])


def exact_ranking_developer_instruction(request: Mapping[str, Any], backend_name: str) -> str:
    candidate_ids = [row["candidate_id"] for row in request["candidates"]]
    schema = {
        "type": "object",
        "additionalProperties": False,
        "required": ["schema_version", "pathway", "run_id", "backend", "ranking"],
        "properties": {
            "schema_version": {"const": LLM_RANKING_SCHEMA_VERSION},
            "pathway": {"const": request["pathway"]["name"]},
            "run_id": {"const": request["run_id"]},
            "backend": {"const": backend_name},
            "ranking": {
                "type": "array",
                "minItems": len(candidate_ids),
                "maxItems": len(candidate_ids),
                "items": {
                    "type": "object",
                    "additionalProperties": False,
                    "required": ["candidate_id", "rank"],
                    "properties": {
                        "candidate_id": {"enum": candidate_ids},
                        "rank": {"type": "integer", "minimum": 1, "maximum": len(candidate_ids)},
                    },
                },
            },
        },
    }
    return (
        "Rank every candidate from most to least likely to drive the difference between A and B. "
        "Return JSON only: no Markdown, commentary, code fence, or additional key. "
        "The response must conform exactly to this Phase 7 schema: "
        f"{json.dumps(schema, ensure_ascii=False, separators=(',', ':'))}. "
        f"Use every candidate ID exactly once with ranks 1..{len(candidate_ids)}. "
        f"The allowed candidate IDs are {json.dumps(candidate_ids, separators=(',', ':'))}."
    )


def _final_json_text(raw_completion: str) -> str:
    text = str(raw_completion)
    marker = "<|channel|>final<|message|>"
    if marker in text:
        text = text.rsplit(marker, maxsplit=1)[1]
    text = text.strip()
    for token in ("<|return|>", "<|end|>", "<|endoftext|>"):
        while text.endswith(token):
            text = text[: -len(token)].rstrip()
    if not text.startswith("{") or not text.endswith("}"):
        raise ValueError("Final channel is not a single JSON object")
    return text


def parse_strict_ranking(raw_completion: str, request: Mapping[str, Any]) -> dict[str, Any]:
    try:
        payload = json.loads(_final_json_text(raw_completion))
    except json.JSONDecodeError as error:
        raise ValueError(f"Invalid JSON: {error.msg} at position {error.pos}") from error
    expected = [row["candidate_id"] for row in request["candidates"]]
    validate_llm_response(payload, expected_candidate_ids=expected)
    if payload["pathway"] != request["pathway"]["name"] or payload["run_id"] != request["run_id"]:
        raise ValueError("Response pathway/run_id differs from request")
    if payload["backend"] != TransformersGPTOSSBackend.name:
        raise ValueError("Response backend differs from frozen Transformers backend")
    return payload


class TransformersGPTOSSBackend(LLMBackend):
    """Production adapter; construction and import never load model dependencies."""

    name = "transformers_gpt_oss_mxfp4_v1"
    version = "1"
    scientific_evaluation_allowed = True

    def __init__(
        self,
        protocol: Mapping[str, Any],
        *,
        model_path: str | Path | None = None,
        cache_path: str | Path | None = None,
        execution_authorized: bool | None = None,
        runtime_report: Mapping[str, Any] | None = None,
    ) -> None:
        self.protocol = protocol
        self.inference = protocol["inference"]
        self.model_path = Path(model_path).expanduser() if model_path else None
        configured_cache = self.inference["loading"].get("cache_path")
        self.cache_path = Path(cache_path or configured_cache).expanduser() if (cache_path or configured_cache) else None
        self.execution_authorized = (
            bool(protocol["execution_gate"]["real_llm_inference_allowed"])
            if execution_authorized is None else bool(execution_authorized)
        )
        self._runtime_report = dict(runtime_report) if runtime_report is not None else None
        self._model: Any = None
        self._tokenizer: Any = None
        self._torch: Any = None
        self._model_revision: str | None = self.inference.get("model_revision")
        self._tokenizer_revision: str | None = self.inference.get("tokenizer_revision")

    @property
    def is_loaded(self) -> bool:
        return self._model is not None or self._tokenizer is not None

    def model_load_kwargs(self) -> dict[str, Any]:
        loading = self.inference["loading"]
        if loading.get("local_files_only") is not True:
            raise ValueError("gpt-oss loading must remain local_files_only")
        if loading.get("allow_mxfp4_dequantization") is not False:
            raise ValueError("MXFP4 dequantization is forbidden")
        kwargs: dict[str, Any] = {
            "local_files_only": True,
            "device_map": {"": "cuda:0"},
            "dtype": "auto",
            "use_kernels": False,
        }
        if self.cache_path is not None:
            kwargs["cache_dir"] = str(self.cache_path)
        if self.inference.get("model_revision"):
            kwargs["revision"] = self.inference["model_revision"]
        return kwargs

    def config_load_kwargs(self) -> dict[str, Any]:
        kwargs: dict[str, Any] = {"local_files_only": True}
        if self.cache_path is not None:
            kwargs["cache_dir"] = str(self.cache_path)
        if self.inference.get("model_revision"):
            kwargs["revision"] = self.inference["model_revision"]
        return kwargs

    def tokenizer_load_kwargs(self) -> dict[str, Any]:
        kwargs: dict[str, Any] = {"local_files_only": True}
        if self.cache_path is not None:
            kwargs["cache_dir"] = str(self.cache_path)
        if self.inference.get("tokenizer_revision"):
            kwargs["revision"] = self.inference["tokenizer_revision"]
        return kwargs

    def _ensure_loaded(self) -> None:
        if self.is_loaded:
            return
        if not self.execution_authorized:
            raise RuntimeError("Real gpt-oss inference remains locked by the Phase 7 execution gate")
        report = self._runtime_report or build_runtime_report(self.protocol, self.cache_path)
        if report["classification"]["status"] != SUPPORTED_PRIMARY:
            reasons = report["classification"]["primary_failure_reasons"]
            raise RuntimeError(f"Primary gpt-oss runtime unsupported: {reasons}")
        transformers = importlib.import_module("transformers")
        torch = importlib.import_module("torch")
        model_source = str(self.model_path or self.inference["model"])
        config = transformers.AutoConfig.from_pretrained(
            model_source, **self.config_load_kwargs()
        )
        quantization = getattr(config, "quantization_config", None)
        quant_method = (
            quantization.get("quant_method") if isinstance(quantization, Mapping)
            else getattr(quantization, "quant_method", None)
        )
        dequantize = (
            quantization.get("dequantize", False) if isinstance(quantization, Mapping)
            else getattr(quantization, "dequantize", False)
        )
        if getattr(config, "model_type", None) != "gpt_oss" or quant_method != "mxfp4" or dequantize:
            raise RuntimeError("Local model is not the frozen prequantized MXFP4 gpt-oss runtime")
        tokenizer = transformers.AutoTokenizer.from_pretrained(
            model_source, **self.tokenizer_load_kwargs()
        )
        model = transformers.AutoModelForCausalLM.from_pretrained(
            model_source, config=config, **self.model_load_kwargs()
        )
        device_map = getattr(model, "hf_device_map", {"": str(getattr(model, "device", ""))})
        if any(str(device) not in {"0", "cuda:0"} for device in device_map.values()):
            raise RuntimeError(f"CPU/offload device detected in model map: {device_map}")
        loaded_quantization = getattr(model.config, "quantization_config", None)
        loaded_dequantize = (
            loaded_quantization.get("dequantize", False)
            if isinstance(loaded_quantization, Mapping)
            else getattr(loaded_quantization, "dequantize", False)
        )
        if loaded_dequantize:
            raise RuntimeError("Loaded model silently dequantized MXFP4")
        self._model = model
        self._tokenizer = tokenizer
        self._torch = torch
        self._runtime_report = report
        self._model_revision = getattr(config, "_commit_hash", None) or self._model_revision
        self._tokenizer_revision = (
            getattr(tokenizer, "_commit_hash", None)
            or getattr(tokenizer, "init_kwargs", {}).get("_commit_hash")
            or self._tokenizer_revision
        )

    def _messages(self, request: Mapping[str, Any], retry_number: int, failure_reason: str | None) -> list[dict[str, str]]:
        developer = exact_ranking_developer_instruction(request, self.name)
        if retry_number:
            developer += (
                f" This is fixed formatting retry {retry_number}; the previous response failed strict "
                f"validation for this non-scientific reason: {failure_reason}. Correct only formatting/schema."
            )
        scientific_payload = json.dumps(request, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
        return [
            {"role": "developer", "content": developer},
            {"role": "user", "content": scientific_payload},
        ]

    def _render_prompt(
        self, request: Mapping[str, Any], retry_number: int, failure_reason: str | None
    ) -> dict[str, Any]:
        messages = self._messages(request, retry_number, failure_reason)
        inputs = self._tokenizer.apply_chat_template(
            messages,
            tokenize=True,
            add_generation_prompt=True,
            return_tensors="pt",
            return_dict=True,
            reasoning_effort=self.inference["reasoning_effort"],
        )
        input_ids = inputs["input_ids"]
        prompt_hash = prompt_token_hash(input_ids.detach().cpu().tolist())
        model_device = getattr(self._model, "device", "cuda:0")
        inputs = {key: value.to(model_device) for key, value in inputs.items()}
        return {
            "messages": messages,
            "model_inputs": inputs,
            "prompt_hash": prompt_hash,
            "input_token_count": int(input_ids.shape[-1]),
        }

    def _generate_once(self, rendered: Mapping[str, Any]) -> dict[str, Any]:
        generation = deterministic_decoding_config(self.inference)
        started = time.perf_counter()
        with self._torch.inference_mode():
            generated = self._model.generate(**rendered["model_inputs"], **generation)
        latency = time.perf_counter() - started
        input_count = int(rendered["input_token_count"])
        output_ids = generated[0, input_count:]
        raw = self._tokenizer.decode(output_ids, skip_special_tokens=False)
        return {
            "raw_completion": raw,
            "wall_clock_latency_seconds": float(latency),
            "input_token_count": input_count,
            "output_token_count": int(output_ids.shape[-1]),
        }

    def _provenance(self) -> dict[str, Any]:
        packages = (self._runtime_report or {}).get("snapshot", {}).get("packages", {})
        snapshot = (self._runtime_report or {}).get("snapshot", {})
        return {
            "model_name": self.inference["model"],
            "model_revision": self._model_revision,
            "tokenizer_revision": self._tokenizer_revision,
            "transformers_version": packages.get("transformers") or importlib.metadata.version("transformers"),
            "torch_version": packages.get("torch") or importlib.metadata.version("torch"),
            "cuda_runtime": snapshot.get("cuda", {}).get("runtime"),
            "gpu_model": (snapshot.get("gpu") or {}).get("name"),
            "backend_name": self.name,
            "backend_version": self.version,
            "reasoning_effort": self.inference["reasoning_effort"],
            "decoding_configuration": deterministic_decoding_config(self.inference),
        }

    def rank_with_trace(
        self, request: Mapping[str, Any], *, invalid_raw_directory: str | Path | None = None
    ) -> dict[str, Any]:
        validate_llm_request(request)
        if request["backend"] != self.name:
            raise ValueError("Request backend does not match TransformersGPTOSSBackend")
        self._ensure_loaded()
        max_retries = int(
            self.inference["invalid_output"]["max_retries_after_initial_attempt"]
        )
        scientific_hash = _canonical_hash(request)
        order_hash = candidate_order_hash(request)
        attempts: list[dict[str, Any]] = []
        failure_reason: str | None = None
        parsed: dict[str, Any] | None = None
        for retry_number in range(max_retries + 1):
            rendered = self._render_prompt(request, retry_number, failure_reason)
            generated = self._generate_once(rendered)
            attempt: dict[str, Any] = {
                "retry_number": retry_number,
                "prompt_hash": rendered["prompt_hash"],
                "scientific_content_hash": scientific_hash,
                **generated,
                "parsed_completion": None,
                "valid": False,
                "failure_reason": None,
            }
            try:
                parsed = parse_strict_ranking(generated["raw_completion"], request)
                attempt["parsed_completion"] = parsed
                attempt["valid"] = True
            except (ValueError, TypeError) as error:
                failure_reason = f"{type(error).__name__}: {error}"
                attempt["failure_reason"] = failure_reason
                if invalid_raw_directory is not None:
                    safe_run_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", request["run_id"])
                    write_json_atomic(
                        {**self._provenance(), "run_id": request["run_id"], **attempt},
                        Path(invalid_raw_directory) / f"{safe_run_id}__retry{retry_number:02d}.json",
                    )
            attempts.append(attempt)
            if parsed is not None:
                break
        return {
            "status": "PASS" if parsed is not None else "INVALID_OUTPUT_EXHAUSTED",
            **self._provenance(),
            "run_id": request["run_id"],
            "candidate_order_seed": request["candidate_order_seed"],
            "candidate_order_hash": order_hash,
            "scientific_content_hash": scientific_hash,
            "attempt_count": len(attempts),
            "attempts": attempts,
            "response": parsed,
        }

    def rank(self, request: Mapping[str, Any]) -> dict[str, Any]:
        trace = self.rank_with_trace(request)
        if trace["response"] is None:
            raise InvalidLLMOutputError(
                f"gpt-oss returned invalid output after {trace['attempt_count']} attempts"
            )
        return trace["response"]
