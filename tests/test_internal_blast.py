from __future__ import annotations

from pathlib import Path
import subprocess

import pytest

from organelle_cleaner import internal_blast


def test_run_internal_blast_uses_makeblastdb_then_blastn_with_cli_threads(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    assembly_fasta = tmp_path / "assembly.fa"
    plastid_fasta = tmp_path / "plastid.fa"
    mit_fasta = tmp_path / "mit.fa"
    output_dir = tmp_path / "out"
    assembly_fasta.write_text(">contigA\nACGT\n", encoding="utf-8")
    plastid_fasta.write_text(">plastid_ref\nACGT\n", encoding="utf-8")
    mit_fasta.write_text(">mit_ref\nACGT\n", encoding="utf-8")

    monkeypatch.setattr(internal_blast.shutil, "which", lambda name: f"/usr/bin/{name}")
    commands: list[list[str]] = []

    def fake_run(command: list[str], check: bool, capture_output: bool, text: bool):
        assert check is False
        assert capture_output is True
        assert text is True
        commands.append(command)
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(internal_blast.subprocess, "run", fake_run)

    outputs = internal_blast.run_internal_blast(
        assembly_fasta=assembly_fasta,
        plastid_fasta=plastid_fasta,
        mit_fasta=mit_fasta,
        output_dir=output_dir,
        threads=65,
    )

    db_prefix = output_dir / internal_blast.BLAST_INTERMEDIATES_DIRNAME / internal_blast.BLAST_DB_PREFIX_NAME
    assert commands == [
        [
            "makeblastdb",
            "-in",
            str(assembly_fasta),
            "-dbtype",
            "nucl",
            "-out",
            str(db_prefix),
        ],
        [
            "blastn",
            "-query",
            str(plastid_fasta),
            "-db",
            str(db_prefix),
            "-out",
            str(output_dir / internal_blast.BLAST_INTERMEDIATES_DIRNAME / internal_blast.PLASTID_BLAST_TSV_NAME),
            "-num_threads",
            "65",
            "-outfmt",
            internal_blast.BLAST_OUTFMT_FIELDS,
        ],
        [
            "blastn",
            "-query",
            str(mit_fasta),
            "-db",
            str(db_prefix),
            "-out",
            str(output_dir / internal_blast.BLAST_INTERMEDIATES_DIRNAME / internal_blast.MITOCHONDRIAL_BLAST_TSV_NAME),
            "-num_threads",
            "65",
            "-outfmt",
            internal_blast.BLAST_OUTFMT_FIELDS,
        ],
    ]
    assert outputs.plastid_tsv == output_dir / internal_blast.BLAST_INTERMEDIATES_DIRNAME / internal_blast.PLASTID_BLAST_TSV_NAME
    assert outputs.mit_tsv == output_dir / internal_blast.BLAST_INTERMEDIATES_DIRNAME / internal_blast.MITOCHONDRIAL_BLAST_TSV_NAME


def test_missing_blast_executable_fails_clearly(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        internal_blast.shutil,
        "which",
        lambda name: None if name == "blastn" else f"/usr/bin/{name}",
    )

    with pytest.raises(ValueError, match="BLAST\\+ executables are required but were not found on PATH: blastn"):
        internal_blast.ensure_blast_executables_available()
