import importlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np
import yaml

from gene_embedding_project.genept_scpa.phase7.gpt_oss_backend import (
    TransformersGPTOSSBackend,
    candidate_order_hash,
    deterministic_decoding_config,
    prompt_token_hash,
)
from gene_embedding_project.genept_scpa.phase7.llm_backend import MockLLMBackend
from gene_embedding_project.genept_scpa.phase7.llm_prompts import build_llm_request
from gene_embedding_project.genept_scpa.phase7.schemas import validate_llm_response


ROOT = Path(__file__).resolve().parents[1]


class FakeTransformersBackend(TransformersGPTOSSBackend):
    def __init__(self, protocol, raw_outputs):
        super().__init__(protocol, execution_authorized=True)
        self.raw_outputs = list(raw_outputs)
        self._runtime_report = {
            "snapshot": {
                "packages": {"transformers": "4.57.6", "torch": "2.8.0"},
                "cuda": {"runtime": "12.8"},
                "gpu": {"name": "fake GPU"},
            }
        }

    def _ensure_loaded(self):
        return None

    def _render_prompt(self, request, retry_number, failure_reason):
        scientific = json.dumps(request, sort_keys=True, separators=(",", ":"))
        token_ids = [[ord(character) for character in scientific] + [retry_number]]
        return {
            "messages": self._messages(request, retry_number, failure_reason),
            "model_inputs": {},
            "prompt_hash": prompt_token_hash(token_ids),
            "input_token_count": len(token_ids[0]),
        }

    def _generate_once(self, rendered):
        raw = self.raw_outputs.pop(0)
        return {
            "raw_completion": raw,
            "wall_clock_latency_seconds": 0.01,
            "input_token_count": rendered["input_token_count"],
            "output_token_count": len(raw),
        }


class Phase7GPTOSSBackendTests(unittest.TestCase):
    def setUp(self):
        self.protocol = yaml.safe_load(
            (ROOT / "config/phase7_gpt_oss_synthetic.yaml").read_text(encoding="utf-8")
        )
        genes = ["G1", "G2", "G3"]
        a = np.arange(18, dtype=float).reshape(6, 3) / 10
        b = a + np.array([0.1, 0.3, 0.2])
        descriptions = {gene: f"Gene Symbol {gene} description" for gene in genes}
        self.request = build_llm_request(
            experiment_id="E1", run_id="E1_stats_order1", pathway="P1",
            source_database="KEGG", genes=genes, condition_a=a, condition_b=b,
            descriptions=descriptions, prompt_condition="stats_only",
            candidate_order_seed=7, description_shuffle_seed=8,
            backend=TransformersGPTOSSBackend.name,
        ).request

    def valid_raw(self):
        response = {
            "schema_version": "phase7.llm-ranking.v1",
            "pathway": "P1", "run_id": "E1_stats_order1",
            "backend": TransformersGPTOSSBackend.name,
            "ranking": [
                {"candidate_id": row["candidate_id"], "rank": index + 1}
                for index, row in enumerate(self.request["candidates"])
            ],
        }
        return "<|channel|>final<|message|>" + json.dumps(response) + "<|return|>"

    def test_constructor_does_not_import_or_load_runtime_modules(self):
        with mock.patch("importlib.import_module") as importer:
            backend = TransformersGPTOSSBackend(self.protocol)
        importer.assert_not_called()
        self.assertFalse(backend.is_loaded)

    def test_local_files_only_and_no_cpu_or_precision_fallback(self):
        backend = TransformersGPTOSSBackend(self.protocol)
        model_kwargs = backend.model_load_kwargs()
        self.assertTrue(model_kwargs["local_files_only"])
        self.assertEqual(model_kwargs["device_map"], {"": "cuda:0"})
        self.assertFalse(model_kwargs["use_kernels"])
        self.assertNotIn("offload_folder", model_kwargs)
        self.assertTrue(backend.tokenizer_load_kwargs()["local_files_only"])
        self.assertFalse(self.protocol["inference"]["loading"]["allow_mxfp4_dequantization"])

    def test_deterministic_decoding_and_prompt_hash(self):
        decoding = deterministic_decoding_config(self.protocol["inference"])
        self.assertEqual(decoding, {
            "do_sample": False, "num_beams": 1,
            "max_new_tokens": 2048, "use_cache": True,
        })
        self.assertEqual(prompt_token_hash([[1, 2, 3]]), prompt_token_hash([[1, 2, 3]]))
        self.assertNotEqual(prompt_token_hash([[1, 2, 3]]), prompt_token_hash([[3, 2, 1]]))
        self.assertEqual(candidate_order_hash(self.request), candidate_order_hash(self.request))

    def test_response_schema_is_equivalent_to_mock_backend_contract(self):
        real = FakeTransformersBackend(self.protocol, [self.valid_raw()]).rank(self.request)
        mock_request = dict(self.request)
        mock_request["backend"] = MockLLMBackend.name
        mock_response = MockLLMBackend().rank(mock_request)
        self.assertEqual(set(real), set(mock_response))
        self.assertEqual(set(real["ranking"][0]), set(mock_response["ranking"][0]))
        validate_llm_response(real)

    def test_invalid_json_retries_are_bounded_and_raw_is_preserved(self):
        backend = FakeTransformersBackend(
            self.protocol, ["not json", "{bad json", self.valid_raw()]
        )
        with tempfile.TemporaryDirectory() as directory:
            trace = backend.rank_with_trace(
                self.request, invalid_raw_directory=Path(directory)
            )
            invalid_files = sorted(Path(directory).glob("*.json"))
            self.assertEqual(len(invalid_files), 2)
            self.assertEqual(trace["attempt_count"], 3)
            self.assertEqual(trace["status"], "PASS")
            self.assertEqual([attempt["raw_completion"] for attempt in trace["attempts"][:2]], ["not json", "{bad json"])
            self.assertEqual(
                len({attempt["scientific_content_hash"] for attempt in trace["attempts"]}), 1
            )
            saved = json.loads(invalid_files[0].read_text(encoding="utf-8"))
            self.assertIn("failure_reason", saved)
            self.assertIn("raw_completion", saved)

    def test_invalid_retry_limit_stops_after_three_total_attempts(self):
        backend = FakeTransformersBackend(self.protocol, ["bad1", "bad2", "bad3"])
        trace = backend.rank_with_trace(self.request)
        self.assertEqual(trace["status"], "INVALID_OUTPUT_EXHAUSTED")
        self.assertEqual(trace["attempt_count"], 3)
        self.assertIsNone(trace["response"])


if __name__ == "__main__":
    unittest.main()
