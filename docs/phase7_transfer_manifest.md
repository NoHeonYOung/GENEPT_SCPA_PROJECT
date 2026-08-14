# Phase 7 transfer manifest

Generated: `2026-08-14T03:40:08.384765+00:00`

Schema: `phase7.transfer.v1`

This inventory does not copy or download any resource. On the destination machine,
restore each file at its exact repository-relative path and run the integrity checker.

## Summary

- Resources: 64
- Required: 61
- Required separate-transfer files: 2
- Required separate-transfer bytes: 832572752
- Optional source RDS is provenance only; the checksum-validated sparse export is the Phase 7 input.
- The historical export manifest contains lab absolute paths. Phase 7 uses the exact
  repo-relative `source.counts_export_files` map and verifies each original manifest hash.

## Resources

| Required | Separate | Git | Size (bytes) | Repository-relative path | SHA256 | Purpose |
|---|---:|---:|---:|---|---|---|
| required | false | true | 13772 | `config/genept_scpa.yaml` | `077d0e078b145bebc96049caae6d0c71a211b13d2e9887bd159a43891b564463` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 8489 | `config/phase7_gpt_oss_synthetic.yaml` | `21d7ab27ccd8550ba2ade6cd288c8af4b31b9171e4637a8d23d2a00aba2af044` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 36294 | `docs/genept_scpa_decision_log.md` | `00fd492ce709687f7018a392b65b0cd09f532ff517ac11f421a02ad9d78bd3b0` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 4361 | `docs/phase7_gpt_oss_synthetic_protocol.md` | `b3af550121d8d5f6f06cf2b3f643441368026f5f34dde5ff389a4d36a6d738c5` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 36837 | `genept_scpa_experiment_plan.md` | `6d26274465a611404704d78b1644e50336c8226ae5667717f44f61c73320e40a` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 420 | `pyproject.toml` | `a6a49a2d9e7af02f011359f2072faf9e130ff477fd0b19347d92cedc5ffd8c0e` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 26499 | `scripts/genept/build_genept_w.py` | `0add70bd055635f45d811534e23444149104e1e35b1baa0048a021a9df4f9a84` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 19791 | `scripts/phase3/run_cd4_cd8_benchmark.py` | `25f336ff4566b1e6ef7e8bf75976d28251a0e8e1f657be1c2eeb09d4b3641220` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 30797 | `scripts/phase4/run_pathway_comparison.py` | `09a3cfdb5ec3ec8873b1ffa42f0706be55628f95dbd12a3a8fabeddc00803557` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 48737 | `scripts/phase4/run_timecourse_validation.py` | `2a759c9f39e03e98a084a3a1dfb1e73567e5356901f7603689935d5af718d09c` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 7004 | `scripts/phase7/build_environment_snapshot.py` | `b5e430530f96ef4f08a5fbcf94dfcdd9c7803c88eb46eb808114831c72a01d9c` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 4365 | `scripts/phase7/build_llm_requests.py` | `74c03c6d7ffc46ab690abc4c99bbf1ed64b62c38603b24a2010f10dbea387179` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 9501 | `scripts/phase7/build_transfer_manifest.py` | `d903147c46938ecdd4e3905ceb13335e3961e2e9a54d8a62854b87b04bf4b26d` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 1221 | `scripts/phase7/check_gpt_oss_runtime.py` | `d8fc985aa2a2e147e6aecc07990d3e5c62165f1de1b3b6d8c46be7419ca024db` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 1222 | `scripts/phase7/check_transfer_integrity.py` | `148789db19ffc6be2c43ef3d1ab6caf2c37caf9fb8a5d08baa0b6d65f1ca2211` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 7386 | `scripts/phase7/evaluate_rankings.py` | `d1d792bfb026c8cc7e9ab72db168a8469bc83b59df9b0e463e9a06320541d7a5` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 15350 | `scripts/phase7/prepare_synthetic_benchmark.py` | `ad12911f9778a536d6c564290815b8f01c91f433c752dcc65c2b587c938f0ecd` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 3441 | `scripts/phase7/run_gpt_oss_inference.py` | `d4a94c6e4d5cc38de90c7a95ccba56dc8ec1693e189486a8d6b8c9efb6b82b09` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 1623 | `scripts/phase7/run_mock_llm.py` | `3b4bdd4443e0f8ccb8ddbea843c7815d8046167c2c1422bbd3fe59cbf6b2e2fc` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 1510 | `scripts/phase7/run_scpa_masking.py` | `2308b678fa1c00edceb5508a8180606c269e45c85a03abb2216c2c85f4dc6f1c` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 4174 | `scripts/phase7/run_toy_smoke.py` | `0d9c5475e3483a7befff4eeaf9d982173052d766377fd4ff8c4b8e1fb0bd04c5` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 1227 | `scripts/phase7/validate_llm_rankings.py` | `46ddb9d3b631fd46d84f6844c50df6923def031ca9c60db547963b09588e448c` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 5910 | `scripts/scpa/run_phase7_synthetic_masking_core.R` | `41f5f00786fc785cedc81cd2ab4bf28b6ffea04b5150945671a908fd4ca47ab9` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 5126 | `scripts/scpa/scpa_core_adapter.R` | `1052629d6885a393b8077213d3aafd0b462a3f18ee7ccc8167ce9bf55ee71162` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 12366 | `src/gene_embedding_project/genept_scpa/config.py` | `aa2e12bc904ac162f5d27ee6c099283dda29b36d2d478922c75f7405097a1762` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 4799 | `src/gene_embedding_project/genept_scpa/gene_mapping.py` | `b2969ff5b232be9f9de786ff7c667004fd5b33dc7e859f1cf7c6a6f5ac531c7d` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 8388 | `src/gene_embedding_project/genept_scpa/genept_projection.py` | `81c3c09e5b9997667c050c86717c43be1851ffd9893d04cad8706b43420e6ced` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 3939 | `src/gene_embedding_project/genept_scpa/io.py` | `2a4ab6d026d95ad437e56a4fd8833e5246abf6c5d99f67f17f3876d31cb1217d` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 8535 | `src/gene_embedding_project/genept_scpa/pathway_projection.py` | `8a24a41c1edbe7cecbe4c684aeab00e30891e1b48a8bc8748ba5c58335b750be` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 441 | `src/gene_embedding_project/genept_scpa/phase7/__init__.py` | `3abcd6fce4fb8219ff2456c47eeb07c4227ff499fac9e2ad732b48fd5089a30b` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 1775 | `src/gene_embedding_project/genept_scpa/phase7/cohort.py` | `02ec1b9193c9b4605050c8a0979a5a1092559c108ccbf74cde59f3a273e98bb1` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 3454 | `src/gene_embedding_project/genept_scpa/phase7/evaluation.py` | `876dbac85d1dc04bef23f40a2122df41a79af95ae0b2c70a06f30eff9e0bccfa` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 17065 | `src/gene_embedding_project/genept_scpa/phase7/gpt_oss_backend.py` | `7dfc03f056ea78dd73cab2add2d3d4f30a289f6155e388b008f850af2d1912b2` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 1820 | `src/gene_embedding_project/genept_scpa/phase7/llm_backend.py` | `847912972bc98eee00eed83704b0d897375cb1385666ce0a61052dbe826b50a4` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 5179 | `src/gene_embedding_project/genept_scpa/phase7/llm_prompts.py` | `4eba5b5ca31b929c7d6ef7657a8656d56abc1f891be0fe651c2a020f604c07b7` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 4299 | `src/gene_embedding_project/genept_scpa/phase7/pathway_selection.py` | `b064bc6225a75b2929a18532d0b3a28c4304761345681bdfdf5631352716e1b7` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 6390 | `src/gene_embedding_project/genept_scpa/phase7/ranking.py` | `fe0adaa6d3ccc2c056ef3675ed1804189c4b0c43f33da8d2949d6231057934ca` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 8386 | `src/gene_embedding_project/genept_scpa/phase7/runtime.py` | `5841b1535bffd689a4c9652df12034b366ec9fba31117f5db0dc931ae2739e91` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 6669 | `src/gene_embedding_project/genept_scpa/phase7/schemas.py` | `2f1a71b045fea0e8d6c0bf3bfe6e37f76d6448183cd98d1e52229cb280357eff` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 7087 | `src/gene_embedding_project/genept_scpa/phase7/synthetic_perturbation.py` | `af70f252a7ab4c945f78bb6ce4838abae9c973d0c1287b6141d735bcb35f97e3` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 2951 | `src/gene_embedding_project/genept_scpa/phase7/transfer.py` | `102dfd6b60b49f81a0fea5e9fa6e0a05f631a909bfc72673163b9a202d5e73e5` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 9259 | `tests/test_config.py` | `8c6d06bf1c3ce03b1d26c632e883a71fdd496c9c1d854f1c81a4c17626a188a4` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 1135 | `tests/test_phase7_cohort.py` | `7a0908d983fe265a0a3ac051acb38a7f7965bc001a663696b3f0c47fbb2c0d28` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 1585 | `tests/test_phase7_execution_gates.py` | `ff22003a73ac5edf59d62503dfdbc86e52250de987eff758c5e96f6efbadc9b9` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 6654 | `tests/test_phase7_gpt_oss_backend.py` | `9fb7fd870c27ff4f591f023dacd1139863f864d78ad054bafeaf26feb45ad1ee` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 2142 | `tests/test_phase7_gpt_oss_runtime.py` | `93298f93b2b042ff029d283fa6fbbf5a9aeba4c9775b315ec30bf0365d3825c0` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 4599 | `tests/test_phase7_llm_pipeline.py` | `f66e0a5af3b905999b9c86e2285552ee2df9351fc72ef84a7b186bdd9aa551af` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 5488 | `tests/test_phase7_ranking_metrics.py` | `be07179bf406061aba8167ef9dad7b2e942c8c035fca42ca7300b1ca4004e07c` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 800 | `tests/test_phase7_scpa_masking.R` | `e9fc6bbaff72f6d481a1a09ea864f191714895ce18ea7798f13540efcf5b1b30` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 3826 | `tests/test_phase7_synthetic.py` | `d18113db156eabce268f1d8e4a9312ead051c6b5b02a3baa86d54c909afc98fd` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 810 | `tests/test_phase7_toy_smoke.py` | `95540591521f1e3c9a7b00f63c831e54e1c27839a8be5a72c9ae963524c78ba1` | Phase 7 implementation, protocol dependency, or regression test |
| required | false | true | 4868 | `tests/test_phase7_transfer_integrity.py` | `3e16129d998f226b0eefda554cd5a2d1a93fe96a806a1222c0691bd6b6c342b3` | Phase 7 implementation, protocol dependency, or regression test |
| optional | false | true | 4670 | `data/interim/genept_scpa/phase1_dataset_qc.json` | `d8936266b3dbcbd732e9874c3fa952d53c4389a618b061097d7e5f7af9902b72` | Optional source-dataset validation provenance |
| optional | false | true | 414 | `data/interim/genept_scpa/phase1_download_metadata.json` | `51b402e940a29f166d903be585855955adb5ed21c16ca9809f5b327202934ced` | Optional GEO acquisition provenance |
| required | false | true | 312774 | `data/interim/genept_scpa/phase2_export/naive_cd4/naive_cd4_cell_ids.txt` | `670a1680d316e1f664daa88bfd550f719462d8bb9e7d6df00c1ce79ae54277ef` | Original source cell IDs |
| required | false | true | 3390 | `data/interim/genept_scpa/phase2_export/naive_cd4/naive_cd4_export_manifest.json` | `8b18bf6ae66f8a970ee913d8fff9ac73fcec0bd744fa0d6da4daeb890d9342af` | Frozen RNA/counts export dimensions and component hashes |
| required | false | true | 131560 | `data/interim/genept_scpa/phase2_export/naive_cd4/naive_cd4_gene_ids.txt` | `7a75f8cf23c63a87ca1ec5092a6aa7449c589d10f9cf28cf41ea3e5fa16ddf17` | RNA/counts gene axis |
| required | false | true | 559735 | `data/interim/genept_scpa/phase2_export/naive_cd4/naive_cd4_metadata.csv` | `9a36ce3d806727c559f4cb142e54b41ac2d4dc60ac983f17b5c5211a835b9ab7` | Cell metadata including Hour and Cell_Type |
| required | true | false | 371775504 | `data/interim/genept_scpa/phase2_export/naive_cd4/naive_cd4_rna_counts_genes_by_cells.mtx` | `6ea1626a0610d701fd23ae330ae384cfcbc90d013a458568f24d686d50ad9e88` | GSE212270 naïve CD4 RNA/counts sparse matrix |
| required | false | true | 850546 | `data/processed/genept_scpa/phase4/pathway_projection_manifest.json` | `d7d7ffcbdf71649b96600997251dc1db307dfded862babbb82f2514290015611` | Frozen Phase 4 paired pathway universe |
| optional | true | false | 2531273585 | `data/raw/genept_scpa/GSE212270_integrated_naive_cd4.rds` | `ee5e52b26ed39611007fd5522cc52e2aec23c83cf0c613032f814eef4140c358` | Optional source-object provenance; validated export is sufficient for Phase 7 |
| required | false | true | 245222 | `data/reference/genept_scpa/combined_metabolic_pathways.csv` | `6bc5977da3fa60f86d5ffb59fc938740bf418fa4d976182a314d65479eb8b744` | Source pathway/gene collection and provenance |
| required | true | false | 460797248 | `data/reference/genept_scpa/genept_ada002/GenePT_gene_embedding_ada_text.pickle` | `fd297510ddd3040744033fde0b0f2cf15a40ac8b2fd2fb02f10667295e55c862` | Official GenePT ada-002 gene embeddings |
| required | false | true | 10982188 | `data/reference/genept_scpa/genept_ada002/NCBI_summary_of_genes.json` | `3db721b7cfbc35795428c6a7ab4e8a9c59fb5f5ccb8c64af50dd3c97c2d4de2e` | Primary gene keys and sanitized LLM descriptions |

## Verification

```bash
PYTHONPATH=src python scripts/phase7/check_transfer_integrity.py \
  --manifest artifacts/phase7_transfer_manifest.json \
  --root "$PWD"
```

The checker never writes, downloads, or searches alternate paths.
