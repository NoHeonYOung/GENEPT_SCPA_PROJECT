"""Reproducible metadata helpers for versioned experiment artifacts."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Any


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 digest of a file without loading it into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gzip_is_valid(path: str | Path, chunk_size: int = 1024 * 1024) -> bool:
    """Read a gzip stream through EOF so CRC/truncation errors are detected."""
    try:
        with gzip.open(path, "rb") as handle:
            for _ in iter(lambda: handle.read(chunk_size), b""):
                pass
    except (OSError, EOFError):
        return False
    return True


def build_download_metadata(
    path: str | Path,
    *,
    geo_accession: str,
    download_source: str,
) -> dict[str, Any]:
    """Inspect a completed GEO download and build its provenance record."""
    archive = Path(path)
    if not archive.is_file():
        raise FileNotFoundError(f"Downloaded archive not found: {archive}")

    gzip_integrity = gzip_is_valid(archive)
    if not gzip_integrity:
        raise ValueError(f"Gzip integrity check failed: {archive}")

    recorded_at = datetime.now(timezone.utc).isoformat()
    return {
        "geo_accession": geo_accession,
        "filename": archive.name,
        "download_source": download_source,
        "file_size_bytes": archive.stat().st_size,
        "sha256": sha256_file(archive),
        "gzip_integrity": True,
        "recorded_at": recorded_at,
        "recorded_at_utc": recorded_at,
    }


def write_json_atomic(payload: dict[str, Any], output: str | Path) -> None:
    """Atomically write JSON so interrupted runs cannot leave a valid-looking file."""
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            json.dump(payload, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, destination)
    finally:
        if temporary_name is not None and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def _download_metadata_command(args: argparse.Namespace) -> int:
    metadata = build_download_metadata(
        args.input,
        geo_accession=args.geo_accession,
        download_source=args.download_source,
    )
    write_json_atomic(metadata, args.output)
    print(f"Download metadata written: {Path(args.output).resolve()}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    metadata_parser = subparsers.add_parser(
        "download-metadata",
        help="validate a gzip download and write SHA-256/size provenance JSON",
    )
    metadata_parser.add_argument("--input", required=True)
    metadata_parser.add_argument("--output", required=True)
    metadata_parser.add_argument("--geo-accession", required=True)
    metadata_parser.add_argument("--download-source", required=True)
    metadata_parser.set_defaults(handler=_download_metadata_command)
    args = parser.parse_args(argv)
    return int(args.handler(args))


if __name__ == "__main__":
    raise SystemExit(main())
