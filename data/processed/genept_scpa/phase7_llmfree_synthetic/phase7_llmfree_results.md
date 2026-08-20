# Phase 7 LLM-free synthetic benchmark 결과

비-null synthetic benchmark에서 평균 AP는 genept_scpa가 더 높았고 (Vanilla=0.4820, GenePT=0.5541), paired Wilcoxon은 유의했다 (raw p=5.723e-14, rank-biserial=0.338).

- 상태: **PASS_WITH_NULL_WARNING**
- 완료 experiment: 1540 / 1540
- 방법: Vanilla SCPA zero-mask vs GenePT non-L2 exact subtraction mask
- 순위 신호: `-log10(raw_p_full) - -log10(raw_p_masked)` 내림차순
- 추론 단위: pathway × scenario × draw에서 두 strength AP를 평균한 paired unit
- Null 경고: **True** — 경고가 있으면 비-null 해석보다 방법론 점검이 우선
- Truth-pool fallback experiment: 560; FALLBACK/NO_FALLBACK 집계를 함께 확인

## 해석 제한

- 이 결과는 Naïve CD4 0h 한 집단과 동결된 synthetic perturbation에서의 injected-signal recovery만 뜻한다.
- GenePT의 biological superiority, 일반적 우월성, 인과성을 주장하지 않는다.
- Phase 6은 correct gene↔embedding correspondence의 특이성을 검정한 별도 질문이며 결론을 합치지 않는다.
- GenePT 열세가 나오더라도 pathway K genes와 1,536D projection의 geometry 차이가 confound일 수 있다.
