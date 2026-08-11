# GSE212270 naïve CD8 acquisition and validation

Status: `READY_FOR_USER_RUN`. 이 pipeline은 향후 Phase 3 CD4-vs-CD8 benchmark를
준비할 뿐 Phase 1B, Phase 2 또는 Phase 3를 시작하거나 PASS 처리하지 않는다.

## Official input

NCBI GEO GSE212270 Series record는 다음 full processed Seurat object를 제공한다.

```text
File: GSE212270_integrated_naive_cd8.rds.gz
GEO displayed size: approximately 1.1 GB
Samples represented by the Series: naïve CD8 0 h, 12 h, 24 h
```

Official sources:

- GEO Series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212270>
- NCBI GEO FTP file:
  <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212270/suppl/GSE212270_integrated_naive_cd8.rds.gz>

제3자 mirror, raw tar 또는 SRA FASTQ를 사용하지 않는다. GEO 표시 크기 1.1 GB는
rounded display value이며 exact byte size와 SHA-256은 다운로드 후 로컬 파일에서
측정한다.

## User commands

프로젝트 루트에서 순서대로 실행한다.

```bash
bash scripts/data/download_naive_cd8.sh
Rscript scripts/data/validate_naive_cd8.R
```

압축 archive와 extracted RDS를 함께 저장하므로 표시 archive 크기보다 충분히 큰
여유 공간이 필요하다.

## Download safeguards

- official NCBI GEO HTTPS URL 고정
- `.part`와 HTTP range를 이용한 resume
- completed archive가 있으면 gzip integrity 확인 후 재사용
- invalid existing archive는 자동 overwrite하지 않고 non-zero 종료
- 다운로드 완료 전에 final filename으로 이동하지 않음
- gzip CRC/EOF 확인
- 실제 byte size, SHA-256, recorded time 기록
- metadata JSON atomic write
- extracted RDS도 `.part` 완료 후 final filename으로 이동
- 기존 naïve CD4 archive/RDS/metadata를 사용하거나 덮어쓰지 않음

생성 파일:

```text
data/raw/genept_scpa/GSE212270_integrated_naive_cd8.rds.gz
data/raw/genept_scpa/GSE212270_integrated_naive_cd8.rds
data/interim/genept_scpa/naive_cd8_download_metadata.json
```

## Validation and identity evidence

Validator는 사전 cell/feature 수나 metadata 구조를 가정하지 않고 실제 객체에서
동적으로 검사한다.

- gzip/RDS/Seurat/compatibility와 source non-modification
- cells/features, assays, active assay, RNA layers
- counts와 normalized expression
- metadata columns와 0/12/24 h detection/counts
- gene identifier type와 missing/duplicate names
- official GEO source/filename evidence와 internal naïve CD8 evidence 분리
- `Cell_Type` 존재 여부, unique values와 counts

Internal metadata에 명시적 naïve CD8 label이 없어도 official GEO source와 filename이
확인되면 identity evidence는 성립하되 warning을 남긴다.

## Read-only CD4 compatibility comparison

CD8 inspection 후 CD8 객체를 메모리에서 해제하고 기존 CD4 RDS를 읽어 다음만
기록한다.

- CD4/CD8 feature counts
- exact shared, CD4-only, CD8-only genes
- gene-order identity
- 양쪽 RNA assay/normalized expression
- 양쪽 0/12/24 h presence
- feature naming compatibility
- 양쪽 `Cell_Type` values/counts
- obvious representation/preprocessing mismatch와 이유

Intersection matrix, merged object 또는 GenePT input을 만들지 않으며 source objects를
수정하지 않는다. Gene matching policy는 Phase 2/3에서 별도로 결정한다.

## Gate interpretation

- `PASS`: CD8 critical checks가 통과하고 obvious CD4 compatibility mismatch 없음
- `NEEDS_REVIEW`: CD8 critical checks는 통과하지만 RNA/normalized representation,
  feature naming 또는 exact shared genes에 근본적 mismatch가 있음
- `FAIL`: official identity/download integrity/RDS/Seurat/RNA/time/features 같은
  critical check 실패

결과:

```text
data/interim/genept_scpa/naive_cd8_dataset_qc.json
```

마지막 stdout 줄은 `NAIVE_CD8_VALIDATION_SUMMARY`로 시작한다. 실행 후 GPT에게
download metadata JSON, dataset QC JSON, 이 summary 줄을 전달한다.
