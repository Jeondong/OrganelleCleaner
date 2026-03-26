#!/usr/bin/env python3
"""Run internal BLAST+ searches and return intermediate TSV paths."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess
import tempfile


BLAST_OUTFMT_FIELDS = (
    "6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send"
)
BLAST_INTERMEDIATES_DIRNAME = "blast_intermediates"
BLAST_DB_PREFIX_NAME = "assembly_blast_db"
PLASTID_BLAST_TSV_NAME = "plastid_hits.tsv"
MITOCHONDRIAL_BLAST_TSV_NAME = "mitochondrial_hits.tsv"
GENERATED_SUBJECT_FASTA_NAME = "assembly_subject.fa"
_BLAST_TEMP_DIRS: list[tempfile.TemporaryDirectory[str]] = []


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
    assembly_fasta: Path | None,
    contigs: dict[str, dict[str, object]] | None = None,
    plastid_fasta: Path | None,
    mit_fasta: Path | None,
    output_dir: Path,
    threads: int,
) -> InternalBlastOutputs:
    if threads < 1:
        raise ValueError("Internal BLAST requires --threads to be at least 1")

    ensure_blast_executables_available()
    output_dir = output_dir.resolve()
    blast_work_dir = tempfile.TemporaryDirectory(prefix="organelle_cleaner_blast_")
    _BLAST_TEMP_DIRS.append(blast_work_dir)
    intermediates_dir = Path(blast_work_dir.name).resolve()
    intermediates_dir.mkdir(parents=True, exist_ok=True)
    database_prefix = intermediates_dir / BLAST_DB_PREFIX_NAME
    # BLAST subjects are the assembly contigs. Prefer sequences already present in
    # the parsed GFA; fall back to --assembly-fasta only when no usable GFA
    # sequences are available.
    subject_fasta = _resolve_blast_subject_fasta(
        assembly_fasta=assembly_fasta.resolve() if assembly_fasta is not None else None,
        contigs=contigs,
        intermediates_dir=intermediates_dir,
    )

    _run_makeblastdb(assembly_fasta=subject_fasta, database_prefix=database_prefix)

    plastid_tsv = None
    mit_tsv = None
    if plastid_fasta is not None:
        plastid_tsv = intermediates_dir / PLASTID_BLAST_TSV_NAME
        _run_blastn(
            query_fasta=plastid_fasta.resolve(),
            database_prefix=database_prefix,
            output_tsv=plastid_tsv,
            threads=threads,
        )
    if mit_fasta is not None:
        mit_tsv = intermediates_dir / MITOCHONDRIAL_BLAST_TSV_NAME
        _run_blastn(
            query_fasta=mit_fasta.resolve(),
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


def _resolve_blast_subject_fasta(
    *,
    assembly_fasta: Path | None,
    contigs: dict[str, dict[str, object]] | None,
    intermediates_dir: Path,
) -> Path:
    subject_fasta = (intermediates_dir / GENERATED_SUBJECT_FASTA_NAME).resolve()
    if contigs is not None and _write_subject_fasta_from_contigs(contigs, subject_fasta) > 0:
        return subject_fasta
    if assembly_fasta is not None:
        return assembly_fasta.resolve()
    raise ValueError(
        "Internal BLAST requires assembly contig sequences from GFA S lines or --assembly-fasta; "
        "no usable GFA sequences were available"
    )


def _write_subject_fasta_from_contigs(
    contigs: dict[str, dict[str, object]],
    subject_fasta: Path,
) -> int:
    written = 0
    with subject_fasta.open("w", encoding="utf-8", newline="\n") as handle:
        for contig_id in sorted(contigs):
            sequence_value = contigs[contig_id].get("sequence", "")
            sequence = sequence_value if isinstance(sequence_value, str) else ""
            if not sequence:
                continue
            handle.write(f">{contig_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(f"{sequence[start:start + 80]}\n")
            written += 1
    return written


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
