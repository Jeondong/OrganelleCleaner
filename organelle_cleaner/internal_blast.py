#!/usr/bin/env python3
"""Run internal BLAST+ searches and return intermediate TSV paths."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess


BLAST_OUTFMT_FIELDS = (
    "6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send"
)
BLAST_INTERMEDIATES_DIRNAME = "blast_intermediates"
BLAST_DB_PREFIX_NAME = "assembly_blast_db"
PLASTID_BLAST_TSV_NAME = "plastid_hits.tsv"
MITOCHONDRIAL_BLAST_TSV_NAME = "mitochondrial_hits.tsv"


@dataclass(frozen=True)
class InternalBlastOutputs:
    plastid_tsv: Path | None
    mit_tsv: Path | None
    database_prefix: Path
    intermediates_dir: Path


def ensure_blast_executables_available() -> None:
    missing: list[str] = []
    for executable in ("makeblastdb", "blastn"):
        if shutil.which(executable) is None:
            missing.append(executable)
    if missing:
        missing_text = ", ".join(missing)
        raise ValueError(
            f"BLAST+ executables are required but were not found on PATH: {missing_text}"
        )


def run_internal_blast(
    *,
    assembly_fasta: Path,
    plastid_fasta: Path | None,
    mit_fasta: Path | None,
    output_dir: Path,
    threads: int,
) -> InternalBlastOutputs:
    if threads < 1:
        raise ValueError("Internal BLAST requires --threads to be at least 1")

    ensure_blast_executables_available()
    intermediates_dir = output_dir / BLAST_INTERMEDIATES_DIRNAME
    intermediates_dir.mkdir(parents=True, exist_ok=True)
    database_prefix = intermediates_dir / BLAST_DB_PREFIX_NAME

    _run_makeblastdb(assembly_fasta=assembly_fasta, database_prefix=database_prefix)

    plastid_tsv = None
    mit_tsv = None
    if plastid_fasta is not None:
        plastid_tsv = intermediates_dir / PLASTID_BLAST_TSV_NAME
        _run_blastn(
            query_fasta=plastid_fasta,
            database_prefix=database_prefix,
            output_tsv=plastid_tsv,
            threads=threads,
        )
    if mit_fasta is not None:
        mit_tsv = intermediates_dir / MITOCHONDRIAL_BLAST_TSV_NAME
        _run_blastn(
            query_fasta=mit_fasta,
            database_prefix=database_prefix,
            output_tsv=mit_tsv,
            threads=threads,
        )

    return InternalBlastOutputs(
        plastid_tsv=plastid_tsv,
        mit_tsv=mit_tsv,
        database_prefix=database_prefix,
        intermediates_dir=intermediates_dir,
    )


def _run_makeblastdb(*, assembly_fasta: Path, database_prefix: Path) -> None:
    _run_blast_command(
        [
            "makeblastdb",
            "-in",
            str(assembly_fasta),
            "-dbtype",
            "nucl",
            "-out",
            str(database_prefix),
        ]
    )


def _run_blastn(
    *,
    query_fasta: Path,
    database_prefix: Path,
    output_tsv: Path,
    threads: int,
) -> None:
    _run_blast_command(
        [
            "blastn",
            "-query",
            str(query_fasta),
            "-db",
            str(database_prefix),
            "-out",
            str(output_tsv),
            "-num_threads",
            str(threads),
            "-outfmt",
            BLAST_OUTFMT_FIELDS,
        ]
    )


def _run_blast_command(command: list[str]) -> None:
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode == 0:
        return

    stderr = completed.stderr.strip()
    stdout = completed.stdout.strip()
    details = stderr or stdout or "no output captured"
    raise ValueError(f"BLAST command failed: {' '.join(command)}\n{details}")
