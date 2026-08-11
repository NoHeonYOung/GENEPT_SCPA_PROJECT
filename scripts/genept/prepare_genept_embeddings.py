#!/usr/bin/env python3
"""Acquire and validate the pinned official GenePT ada-002 artifact."""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from gene_embedding_project.genept_scpa.gene_mapping import (  # noqa: E402
    load_official_genept_embeddings,
    load_primary_gene_keys,
)
from gene_embedding_project.genept_scpa.io import (  # noqa: E402
    sha256_file,
    write_json_atomic,
)


ZENODO_RECORD = "https://zenodo.org/records/10833191"
ARCHIVE_URL = (
    "https://zenodo.org/api/records/10833191/files/"
    "GenePT_emebdding_v2.zip/content"
)
ARCHIVE_NAME = "GenePT_emebdding_v2.zip"
ARCHIVE_SIZE = 574_395_233
ARCHIVE_MD5 = "3f6ce4317e3a0091978ae5cb8fbf05a3"
EMBEDDING_NAME = "GenePT_gene_embedding_ada_text.pickle"
SUMMARY_NAME = "NCBI_summary_of_genes.json"


def md5_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.md5()  # noqa: S324 - official Zenodo integrity checksum
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_with_curl(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "curl",
            "--fail",
            "--location",
            "--continue-at",
            "-",
            "--output",
            str(destination),
            url,
        ],
        check=True,
    )


def validate_archive(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Official GenePT archive not found: {path}")
    if path.stat().st_size != ARCHIVE_SIZE:
        raise ValueError(
            f"Archive size mismatch: {path.stat().st_size} != {ARCHIVE_SIZE}"
        )
    observed_md5 = md5_file(path)
    if observed_md5 != ARCHIVE_MD5:
        raise ValueError(f"Archive MD5 mismatch: {observed_md5} != {ARCHIVE_MD5}")
    if not zipfile.is_zipfile(path):
        raise ValueError("Official GenePT archive is not a valid ZIP file")


def extract_selected(archive: Path, output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    wanted = {EMBEDDING_NAME, SUMMARY_NAME}
    selected: dict[str, zipfile.ZipInfo] = {}
    with zipfile.ZipFile(archive) as handle:
        for member in handle.infolist():
            basename = Path(member.filename).name
            if basename in wanted:
                if basename in selected:
                    raise ValueError(f"Duplicate archive member basename: {basename}")
                selected[basename] = member
        missing = wanted - set(selected)
        if missing:
            raise ValueError(
                "Official archive is missing required files: "
                + ", ".join(sorted(missing))
            )
        extracted: dict[str, Path] = {}
        for basename, member in selected.items():
            destination = output_dir / basename
            temporary_name: str | None = None
            try:
                with handle.open(member) as source, tempfile.NamedTemporaryFile(
                    mode="wb",
                    dir=output_dir,
                    prefix=f".{basename}.",
                    suffix=".tmp",
                    delete=False,
                ) as target:
                    temporary_name = target.name
                    while True:
                        chunk = source.read(1024 * 1024)
                        if not chunk:
                            break
                        target.write(chunk)
                    target.flush()
                    os.fsync(target.fileno())
                os.replace(temporary_name, destination)
                extracted[basename] = destination
            finally:
                if temporary_name and os.path.exists(temporary_name):
                    os.unlink(temporary_name)
    return extracted


def inspect_artifact(embedding_path: Path, summary_path: Path) -> dict[str, object]:
    embeddings = load_official_genept_embeddings(
        embedding_path, expected_dimension=1536
    )
    primary_keys = load_primary_gene_keys(summary_path)
    first_vector = next(iter(embeddings.values()))
    alias_key_count = len(set(embeddings) - primary_keys)
    return {
        "embedding_gene_key_count": len(embeddings),
        "primary_ncbi_gene_key_count": len(primary_keys),
        "official_alias_key_count": alias_key_count,
        "dimension": int(first_vector.size),
        "dtype_after_load": str(first_vector.dtype),
        "finite_embeddings": True,
        "embedding_key_schema": "uppercase human gene symbols plus official HGNC alias keys",
        "identifier_fields": {
            "artifact": "dictionary key (gene symbol or official HGNC alias)",
            "primary_reference": "NCBI_summary_of_genes.json object key",
        },
        "embedding_side_duplicate_keys": 0,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--download",
        action="store_true",
        help="download/resume the pinned 574 MB Zenodo archive if needed",
    )
    parser.add_argument(
        "--archive",
        type=Path,
        default=PROJECT_ROOT / "data/reference/genept_scpa" / ARCHIVE_NAME,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "data/reference/genept_scpa/genept_ada002",
    )
    parser.add_argument(
        "--provenance-output",
        type=Path,
        default=PROJECT_ROOT
        / "data/interim/genept_scpa/phase2_genept_embedding_provenance.json",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    archive = args.archive.resolve()
    archive_incomplete = not archive.exists() or archive.stat().st_size < ARCHIVE_SIZE
    if archive_incomplete:
        if not args.download:
            raise FileNotFoundError(
                f"{archive} is missing or incomplete; rerun with --download "
                "to acquire/resume it"
            )
        download_with_curl(ARCHIVE_URL, archive)
    validate_archive(archive)
    extracted = extract_selected(archive, args.output_dir.resolve())
    embedding_path = extracted[EMBEDDING_NAME]
    summary_path = extracted[SUMMARY_NAME]
    inspection = inspect_artifact(embedding_path, summary_path)
    payload = {
        "source": {
            "repository": "https://github.com/yiqunchen/GenePT",
            "zenodo_record": ZENODO_RECORD,
            "doi": "10.5281/zenodo.10833191",
            "archive_url": ARCHIVE_URL,
            "archive_filename": ARCHIVE_NAME,
            "archive_size_bytes": archive.stat().st_size,
            "archive_md5": md5_file(archive),
            "archive_sha256": sha256_file(archive),
        },
        "primary_embedding": {
            "model": "text-embedding-ada-002",
            "filename": EMBEDDING_NAME,
            "path": str(embedding_path),
            "size_bytes": embedding_path.stat().st_size,
            "sha256": sha256_file(embedding_path),
            **inspection,
        },
        "primary_gene_summaries": {
            "filename": SUMMARY_NAME,
            "path": str(summary_path),
            "size_bytes": summary_path.stat().st_size,
            "sha256": sha256_file(summary_path),
        },
        "openai_api_called": False,
    }
    write_json_atomic(payload, args.provenance_output)
    print(
        "GENEPT_EMBEDDING_PREP status=PASS "
        f"genes={inspection['embedding_gene_key_count']} "
        f"dimension={inspection['dimension']} "
        f"provenance={args.provenance_output.resolve()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
