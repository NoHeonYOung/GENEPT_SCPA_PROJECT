#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/../.." && pwd)"
raw_dir="${project_root}/data/raw/genept_scpa"
interim_dir="${project_root}/data/interim/genept_scpa"
archive="${raw_dir}/GSE212270_integrated_naive_cd4.rds.gz"
rds_file="${raw_dir}/GSE212270_integrated_naive_cd4.rds"
metadata_file="${interim_dir}/phase1_download_metadata.json"
dataset_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212270/suppl/GSE212270_integrated_naive_cd4.rds.gz"

PROJECT_ROOT="${project_root}" PYTHONPATH="${project_root}/src" python - <<'PY'
import os
from pathlib import Path

from gene_embedding_project.genept_scpa.config import load_config

root = Path(os.environ["PROJECT_ROOT"])
config = load_config(root / "config/genept_scpa.yaml")
if config.max_phase_allowed < 1:
    raise SystemExit(
        "Refusing download: Phase 1 is locked; set max_phase_allowed=1 only "
        "after recording Phase 0 PASS in the decision log"
    )
if config.active_phase != 1:
    raise SystemExit(
        f"Refusing download: active_phase={config.active_phase}; expected 1"
    )
if config.values["phase1"]["status"] != "in_progress":
    raise SystemExit(
        "Refusing download: phase1.status must be in_progress"
    )
PY

mkdir -p "${raw_dir}" "${interim_dir}"

download() {
    local url="$1"
    local destination="$2"
    local partial="${destination}.part"

    if [[ -s "${destination}" ]]; then
        echo "Already present: ${destination}"
        return
    fi

    if command -v curl >/dev/null 2>&1; then
        curl -fL --retry 5 --retry-delay 5 --continue-at - \
            --output "${partial}" "${url}"
    elif command -v wget >/dev/null 2>&1; then
        wget --continue --output-document="${partial}" "${url}"
    else
        echo "curl or wget is required" >&2
        exit 1
    fi
    mv "${partial}" "${destination}"
}

download "${dataset_url}" "${archive}"

# This reads the gzip stream through EOF, then records observed size and SHA-256.
# It exits non-zero and does not write valid-looking metadata if integrity fails.
PROJECT_ROOT="${project_root}" PYTHONPATH="${project_root}/src" python -m \
    gene_embedding_project.genept_scpa.io download-metadata \
    --input "${archive}" \
    --output "${metadata_file}" \
    --geo-accession "GSE212270" \
    --download-source "${dataset_url}"

if [[ ! -s "${rds_file}" ]]; then
    gzip -dk "${archive}"
else
    echo "Already extracted: ${rds_file}"
fi

echo "Phase 1A acquisition complete"
echo "Archive: ${archive}"
echo "Extracted RDS: ${rds_file}"
echo "Metadata: ${metadata_file}"
