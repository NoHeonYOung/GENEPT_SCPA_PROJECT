# Phase 7 GPT-OSS synthetic benchmark protocol

Status: **FROZEN / IMPLEMENTATION SCAFFOLD / MOCK ONLY**

Decision: D-0031

Configuration source of truth: `config/phase7_gpt_oss_synthetic.yaml`

## Research questions

- RQ1: Can Vanilla SCPA, non-L2 GenePT-informed SCPA, and a context-aware LLM
  recover known synthetically perturbed pathway genes?
- RQ2: Does correct biological gene-description correspondence improve LLM
  recovery over stats-only and within-pathway shuffled descriptions?

## Frozen benchmark

- Source: GSE212270 naïve CD4 0h.
- Cohort: 1,000 unique source cells; disjoint pseudo A/B groups of 500; PCG64
  seed 20260814; original cell IDs retained.
- Preprocessing: full-transcriptome total-count normalization to 10,000, then
  `log1p`; perturb in normalized log1p space; no post-injection renormalization.
- Pathways: frozen Phase 4 paired universe only; no Phase 4/5 outcome score may
  enter selection; 15–60 analysis genes; usable sanitized NCBI description for
  every gene; frozen seed 20260815.
- Selection: 6 KEGG and 5 REACTOME pathways across size bins. REACTOME large has
  only one eligible pathway; no HALLMARK pathway is eligible. Criteria are not
  relaxed to fill these strata.
- Truth count: `max(3, min(10, ceil(0.15 * pathway_gene_count)))`.
- Scenarios: null, mean shift, seeded 30% cell subset, and mixed direction.
- Non-null strengths: 0.5 and 1.0 pooled-baseline SD, using
  `scale_g = max(pooled_baseline_sd_g, 0.1)`.
- The same pathway/scenario truth set is reused across strengths.

## Compared rankings

- Vanilla SCPA: zero the masked expression column.
- GenePT-informed SCPA: non-L2 `Z = X @ E`, then subtract the masked gene outer
  product.
- SCPA ranks retain the Phase 5 average-rank rule for exact delta ties; metric
  ordering resolves a remaining tie lexicographically by gene symbol.
- LLM: stats-only, true sanitized descriptions, and a within-pathway derangement
  of the exact same description multiset; three prompt-order repeats.
- Primary ranking metrics: Recall@truth-K, Average Precision, and NDCG@N.
- Secondary checks: NDCG@truth-K and prompt-order Spearman stability.

## Leakage boundary

Method runners may read expression, pathway metadata, descriptions, embeddings,
and an internal opaque-ID mapping. They must not read `phase7_ground_truth.csv`.
Only the evaluation entry point may read the ground-truth table. LLM-visible
requests contain opaque IDs and frozen summaries, not gene symbols, perturbation
metadata, precomputed tests, or truth labels.

## Current execution gate

Allowed now: unit tests and tiny toy smoke using `MockLLMBackend` and a stated
non-SCPA score surrogate. Mock outputs are never scientific results.

Locked until another explicit approval: production synthetic generation,
production SCPA/MCM, gpt-oss-20b download, and real LLM inference.

## Frozen gpt-oss inference

- Model/backend: `openai/gpt-oss-20b`, Hugging Face Transformers, pretrained
  MXFP4 without dequantization or backend/precision fallback.
- Formatting: official `tokenizer.apply_chat_template`; raw Harmony strings are
  never manually assembled; reasoning effort is `low`.
- Decoding: greedy, `do_sample=false`, `num_beams=1`, `max_new_tokens=2048`;
  temperature/top-p sampling is absent.
- Loading: local-files-only and lazy on a single CUDA device. Primary scientific
  production never falls back to CPU/offload.
- Invalid JSON: two fixed retries after the initial attempt. Scientific content
  and candidate order remain identical; every invalid raw output and reason is
  retained.
- Primary capability minimum: CC 7.5, 16GiB VRAM, 32GiB system RAM, 30GiB free
  disk, and all declared MXFP4 runtime packages.

## Split replication and runtime-only pilot

Prompt-order repeats measure LLM ordering stability. Independent pseudo-split
seeds measure biological/sampling robustness. Pilot seed 20260814 is frozen;
production seeds/count remain empty until runtime measurements are collected,
and must be frozen before any scientific metric is inspected.

The future pilot uses one split, the smallest and largest frozen eligible
pathways, stats-only, and one candidate order. It may report only loading,
resource use, latency, throughput, JSON validity and retry/parser behavior.
Recovery metrics, truth positions, method comparisons and ranking interpretation
remain hidden.
