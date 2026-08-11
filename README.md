# GenePT × SCPA

GenePT의 문헌 기반 gene representation을 SCPA에 통합할 수 있는지 단계별로
검증하는 연구 프로젝트입니다. 현재 단계는 **Phase 4 — pathway-specific Vanilla vs
GenePT-informed SCPA comparison**입니다.

## 현재 상태

- Phase 0 PASS
- Phase 1A full dataset acquisition/validation PASS
- Phase 1B Hour-only 네 분석과 official reference reproduction PASS
- Phase 2 published GenePT-w reproduction PASS
- Phase 3 whole-cell GenePT-w CD4/CD8 SCPA-core feasibility PASS
- Phase 4A exploratory lineage comparison 완료, historical output 보존
- Phase 4B Naive CD4 activation 3-comparison production PASS
- Phase 4C optional CD8 generalization NOT SCHEDULED
- Phase 5 pathway-internal paired gene-masking sensitivity IN PROGRESS
- Phase 6 이후 NOT STARTED

설정 검증:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

전체 데이터 준비 및 검증 절차는
[`docs/phase1_data_download.md`](docs/phase1_data_download.md)를 따릅니다.

Phase 1B 고정 프로토콜은
[`docs/phase1b_protocol.md`](docs/phase1b_protocol.md)에 기록되어 있습니다. Hour-only
네 분석과 reference는 이미 실행 및 검토가 끝났으며, Phase 1을 다시 열지 않는 한
재실행하지 않습니다.
Phase 2 고정 protocol과 실제 PASS 근거는
[`docs/phase2_genept_w_protocol.md`](docs/phase2_genept_w_protocol.md)에 있습니다.
Phase 3 protocol과 qval audit은
[`docs/phase3_cd4_cd8_protocol.md`](docs/phase3_cd4_cd8_protocol.md) 및
[`docs/phase3_qval_implementation_audit.md`](docs/phase3_qval_implementation_audit.md)에
기록되어 있습니다.

Phase 4의 결과 확인 전 동결된 방법은
[`docs/phase4_pathway_method_decision.md`](docs/phase4_pathway_method_decision.md)에
있습니다. 아래 명령은 이미 완료된 historical Phase 4A lineage 분석의 runner입니다.

```bash
PYTHONPATH=src python scripts/phase4/run_pathway_comparison.py
```

이 명령은 Phase 3의 canonical CD4/CD8 0 h 500 cells/group을 재사용하며 full
transcriptome normalization 뒤 동일한 paired pathway genes로 Vanilla와
GenePT-informed branch를 실행합니다. 기본 primary는 non-L2 projection입니다. 사전
선언된 L2 sensitivity까지 함께 실행하려면 `--run-l2-sensitivity`를 추가합니다.
Phase 4A historical output은 보존됩니다. 현재 primary Phase 4B는 lineage를 CD4로
고정한 activation benchmark이며, full production 명령은 다음 하나입니다.

```bash
PYTHONPATH=src python scripts/phase4/run_timecourse_validation.py \
  --comparison-set cd4_activation
```

이 명령은 frozen CD4 0h/12h/24h 각 500 cells와 동일 123 paired pathways에서
Vanilla, GenePT non-L2 primary 및 L2 sensitivity를 정확히 3 comparisons에
실행합니다(1,107 MCM calls). 이전 9-comparison 기능은 explicit `all_9` option으로
보존하지만 현재 실행 대상이 아닙니다. CD8 generalization, Phase 5 gene contribution,
semantic controls와 classifier는 실행하지 않습니다.

Phase 5는 Phase 4B의 Vanilla-only/GenePT-only 30개 pathway-comparison pair에서
동일 gene을 Vanilla와 GenePT non-L2에 masking합니다. Full production 명령은 다음
하나이며 pathway-comparison 단위 checkpoint를 재사용합니다.

```bash
PYTHONPATH=src python scripts/phase5/run_gene_contribution.py
```

이 분석은 gene masking sensitivity를 측정하며 gene importance, biological
correctness 또는 GenePT superiority를 주장하지 않습니다. GenePT L2 gene-level
analysis, Phase 6 controls, CD8 generalization과 classifier는 실행하지 않습니다.
