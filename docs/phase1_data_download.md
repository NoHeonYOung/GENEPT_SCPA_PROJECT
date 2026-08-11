# Phase 1A — full dataset acquisition and validation

현재 Phase 1A dataset gate는 실제 전체 객체 validation 후 `PASS`다. Phase 1
전체는 여전히 `in_progress`이며 Phase 1B SCPA reproduction은 별도 단계로
잠겨 있다. 아래 명령과 설명은 acquisition/validation 재현용으로 유지한다.

## 고정한 공식 입력

```text
GEO series: GSE212270
File: GSE212270_integrated_naive_cd4.rds.gz
GEO 표시 크기: 약 2.4 GB
범위: full processed Seurat object, naïve CD4, 0 h / 12 h / 24 h
```

공식 근거:

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212270>
- 논문 코드: <https://github.com/jackbibby1/scpa_paper>
- SCPA tutorial:
  <https://jackbibby1.github.io/SCPA/articles/comparing_two_populations.html>

`GSE212270_RAW.tar`와 SRA FASTQ는 Cell Ranger부터 다시 처리할 때 필요한
입력이다. Phase 1A/1B는 저자가 공개한 full processed Seurat object를 그대로
검사하고 사용하는 것이므로 raw matrix/FASTQ를 받지 않는다.

## 사용자 실행 명령

압축 archive와 해제된 RDS를 함께 저장할 수 있도록 최소 12 GiB 이상의 여유
공간을 권장한다.

```bash
cd /home/node00/nhy_python/GenePT_SCPA
bash scripts/data/download_phase1_data.sh
Rscript scripts/data/validate_phase1_data.R
```

## Download safeguards

다운로드 스크립트는 다음을 수행한다.

1. config의 Phase 1 gate와 `in_progress` 상태를 확인한다.
2. 공식 NCBI GEO HTTPS 주소만 사용한다.
3. destination directory를 자동 생성한다.
4. `.part` 파일과 HTTP range를 사용해 중단된 다운로드를 이어받는다.
5. 완료된 archive가 이미 있으면 다시 받지 않는다.
6. gzip stream을 끝까지 읽어 CRC/truncation 오류를 검사한다.
7. 실제 archive의 SHA-256과 byte size를 계산한다.
8. provenance JSON을 atomic write로 생성한다.
9. 실패 시 non-zero로 종료한다.

생성 파일:

```text
data/raw/genept_scpa/GSE212270_integrated_naive_cd4.rds.gz
data/raw/genept_scpa/GSE212270_integrated_naive_cd4.rds
data/interim/genept_scpa/phase1_download_metadata.json
```

`phase1_download_metadata.json`의 size/hash는 실행 시 관찰한 값만 기록한다.
사전에 가짜 hash나 정확한 cell/feature 수를 넣지 않는다.

## Dataset validation

validation은 알려진 column/assay 이름을 정답으로 고정하지 않는다.
`Hour` 같은 알려진 후보는 우선순위로만 사용하고, 모든 metadata column의 실제
값에서 0/12/24 h가 함께 검출되는 column을 선택한다.

검사 항목:

- archive 존재, gzip integrity, SHA-256과 download metadata 일치
- RDS read 성공 및 Seurat class
- serialized object/installed Seurat 버전, cell/feature 수, assay와 active assay
- 구형 Seurat 객체이면 공식 `UpdateSeuratObject()`를 메모리상 적용하고 성공 여부
  기록(원본 RDS는 수정하지 않음)
- 전체 metadata column과 실제 time column/label, 시간대별 cell 수
- 공식 filename 및 object/metadata에서 얻은 naïve CD4 identity evidence
- assay별 layer/slot, counts 및 normalized expression 존재 여부
- feature identifier 형식 추정, missing/duplicate feature names
- critical failure와 warning 분리

결과:

```text
data/interim/genept_scpa/phase1_dataset_qc.json
```

validation stdout의 마지막 줄은 다음 형태의 summary다.

```text
PHASE1A_VALIDATION_SUMMARY status=... cells=... features=... time_column=... failed_checks=... warnings=... qc_json=...
```

## Phase 1A PASS

다음 조건이 모두 실제 파일에서 확인되어야 `gate.status = PASS`다.

- 공식 GSE212270 filename/source와 download metadata 일치
- gzip integrity와 SHA-256 기록
- RDS read 및 Seurat class 확인
- 0/12/24 h 모두 존재하고 각 조건의 cell 수가 1 이상
- normalized expression layer/data 존재
- 사용 가능한 고유 gene identifiers 존재
- naïve CD4 identity evidence 존재
- critical failed check 없음

counts layer 부재나 복수 time-column 후보처럼 분석 전에 검토 가능한 사항은
warning으로 기록될 수 있다. warning과 failure는 자동으로 혼합하지 않는다.

## GPT에게 전달할 결과

- `data/interim/genept_scpa/phase1_download_metadata.json`
- `data/interim/genept_scpa/phase1_dataset_qc.json`
- validation stdout의 마지막 `PHASE1A_VALIDATION_SUMMARY` 줄
- 실패했다면 전체 error message와 `gate.failed_checks`
