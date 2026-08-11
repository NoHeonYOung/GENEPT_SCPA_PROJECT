#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/../.." && pwd)"
raw_dir="${project_root}/data/raw/genept_scpa"
interim_dir="${project_root}/data/interim/genept_scpa"
filename="GSE212270_integrated_naive_cd8.rds.gz"
archive="${raw_dir}/${filename}"
rds_file="${raw_dir}/${filename%.gz}"
metadata_file="${interim_dir}/naive_cd8_download_metadata.json"
dataset_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212270/suppl/${filename}"

# Official GEO GSE212270 lists this processed Seurat object at approximately 1.1 GB.
# The exact downloaded byte size and SHA-256 are observed locally, never predeclared.
PROJECT_ROOT="${project_root}" PYTHONPATH="${project_root}/src" python - <<'PY'
import os
from pathlib import Path
from gene_embedding_project.genept_scpa.config import load_config

root = Path(os.environ["PROJECT_ROOT"])
config = load_config(root / "config/genept_scpa.yaml")
if config.values["phase0"]["status"] != "passed":
    raise SystemExit("Refusing CD8 preparation: Phase 0 is not passed")
if config.values["phase1"]["dataset_gate_status"] != "passed":
    raise SystemExit("Refusing CD8 preparation: the naïve CD4 dataset gate is not passed")
PY

mkdir -p "${raw_dir}" "${interim_dir}"

if [[ -s "${archive}" ]]; then
    if gzip -t -- "${archive}"; then
        echo "Valid archive already present: ${archive}"
    else
        echo "Existing archive failed gzip integrity; refusing to overwrite: ${archive}" >&2
        exit 1
    fi
else
    partial="${archive}.part"
    if command -v curl >/dev/null 2>&1; then
        curl -fL --retry 5 --retry-delay 5 --continue-at - \
            --output "${partial}" "${dataset_url}"
    elif command -v wget >/dev/null 2>&1; then
        wget --continue --output-document="${partial}" "${dataset_url}"
    else
        echo "curl or wget is required" >&2
        exit 1
    fi
    gzip -t -- "${partial}"
    mv -- "${partial}" "${archive}"
fi

PROJECT_ROOT="${project_root}" PYTHONPATH="${project_root}/src" python -m \
    gene_embedding_project.genept_scpa.io download-metadata \
    --input "${archive}" \
    --output "${metadata_file}" \
    --geo-accession "GSE212270" \
    --download-source "${dataset_url}"

if [[ ! -s "${rds_file}" ]]; then
    extraction_partial="${rds_file}.part"
    gzip -dc -- "${archive}" > "${extraction_partial}"
    mv -- "${extraction_partial}" "${rds_file}"
else
    echo "Already extracted: ${rds_file}"
fi

echo "Naïve CD8 acquisition complete"
echo "Archive: ${archive}"
echo "Extracted RDS: ${rds_file}"
echo "Metadata: ${metadata_file}"
