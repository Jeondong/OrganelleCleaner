from __future__ import annotations

from pathlib import Path
import subprocess

import pytest

from organelle_cleaner import internal_blast


class _FakeTemporaryDirectory:
    def __init__(self, path: Path):
        self.name = str(path)


def test_run_internal_blast_uses_makeblastdb_then_blastn_with_cli_threads(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    assembly_fasta = tmp_path / "assembly.fa"
    plastid_fasta = tmp_path / "plastid.fa"
    mit_fasta = tmp_path / "mit.fa"
    output_dir = tmp_path / "out"
    work_dir = tmp_path / "blast-work"
    assembly_fasta.write_text(">contigA\nACGT\n", encoding="utf-8")
    plastid_fasta.write_text(">plastid_ref\nACGT\n", encoding="utf-8")
    mit_fasta.write_text(">mit_ref\nACGT\n", encoding="utf-8")

    monkeypatch.setattr(internal_blast.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(
        internal_blast.tempfile,
        "TemporaryDirectory",
        lambda prefix: _FakeTemporaryDirectory(work_dir),
    )
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
        contigs=None,
        plastid_fasta=plastid_fasta,
        mit_fasta=mit_fasta,
        output_dir=output_dir,
        threads=65,
    )

    resolved_work_dir = work_dir.resolve()
    db_prefix = resolved_work_dir / internal_blast.BLAST_DB_PREFIX_NAME
    assert commands == [
        [
            "makeblastdb",
            "-in",
            str(assembly_fasta.resolve()),
            "-dbtype",
            "nucl",
            "-out",
            str(db_prefix),
        ],
        [
            "blastn",
            "-query",
            str(plastid_fasta.resolve()),
            "-db",
            str(db_prefix),
            "-out",
            str(resolved_work_dir / internal_blast.PLASTID_BLAST_TSV_NAME),
            "-num_threads",
            "65",
            "-outfmt",
            internal_blast.BLAST_OUTFMT_FIELDS,
        ],
        [
            "blastn",
            "-query",
            str(mit_fasta.resolve()),
            "-db",
            str(db_prefix),
            "-out",
            str(resolved_work_dir / internal_blast.MITOCHONDRIAL_BLAST_TSV_NAME),
            "-num_threads",
            "65",
            "-outfmt",
            internal_blast.BLAST_OUTFMT_FIELDS,
        ],
    ]
    assert outputs.plastid_tsv == resolved_work_dir / internal_blast.PLASTID_BLAST_TSV_NAME
    assert outputs.mit_tsv == resolved_work_dir / internal_blast.MITOCHONDRIAL_BLAST_TSV_NAME
    assert outputs.intermediates_dir == resolved_work_dir
    assert not str(outputs.intermediates_dir).startswith(str(output_dir.resolve()))


def test_run_internal_blast_generates_subject_fasta_from_contigs_when_assembly_fasta_missing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    plastid_fasta = tmp_path / "plastid.fa"
    output_dir = tmp_path / "out"
    work_dir = tmp_path / "blast-work"
    plastid_fasta.write_text(">plastid_ref\nACGT\n", encoding="utf-8")

    monkeypatch.setattr(internal_blast.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(
        internal_blast.tempfile,
        "TemporaryDirectory",
        lambda prefix: _FakeTemporaryDirectory(work_dir),
    )
    commands: list[list[str]] = []

    def fake_run(command: list[str], check: bool, capture_output: bool, text: bool):
        commands.append(command)
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(internal_blast.subprocess, "run", fake_run)

    outputs = internal_blast.run_internal_blast(
        assembly_fasta=None,
        contigs={
            "contigA": {"sequence": "ACGTACGT", "length": 8, "coverage": None},
            "contigB": {"sequence": "", "length": 0, "coverage": None},
        },
        plastid_fasta=plastid_fasta,
        mit_fasta=None,
        output_dir=output_dir,
        threads=4,
    )

    resolved_work_dir = work_dir.resolve()
    subject_fasta = (
        resolved_work_dir
        / internal_blast.GENERATED_SUBJECT_FASTA_NAME
    )
    assert subject_fasta.exists()
    assert subject_fasta.read_text(encoding="utf-8") == ">contigA\nACGTACGT\n"
    assert commands[0][2] == str(subject_fasta)
    assert outputs.plastid_tsv == resolved_work_dir / internal_blast.PLASTID_BLAST_TSV_NAME


def test_run_internal_blast_prefers_gfa_contig_sequences_over_assembly_fasta(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    assembly_fasta = tmp_path / "assembly.fa"
    plastid_fasta = tmp_path / "plastid.fa"
    output_dir = tmp_path / "out"
    work_dir = tmp_path / "blast-work"
    assembly_fasta.write_text(">contigA\nTTTT\n", encoding="utf-8")
    plastid_fasta.write_text(">plastid_ref\nACGT\n", encoding="utf-8")

    monkeypatch.setattr(internal_blast.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(
        internal_blast.tempfile,
        "TemporaryDirectory",
        lambda prefix: _FakeTemporaryDirectory(work_dir),
    )
    commands: list[list[str]] = []

    def fake_run(command: list[str], check: bool, capture_output: bool, text: bool):
        commands.append(command)
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(internal_blast.subprocess, "run", fake_run)

    internal_blast.run_internal_blast(
        assembly_fasta=assembly_fasta,
        contigs={
            "contigA": {"sequence": "ACGTACGT", "length": 8, "coverage": None},
        },
        plastid_fasta=plastid_fasta,
        mit_fasta=None,
        output_dir=output_dir,
        threads=4,
    )

    subject_fasta = (
        work_dir.resolve()
        / internal_blast.GENERATED_SUBJECT_FASTA_NAME
    )
    assert subject_fasta.read_text(encoding="utf-8") == ">contigA\nACGTACGT\n"
    assert commands[0][2] == str(subject_fasta)


def test_run_internal_blast_falls_back_to_assembly_fasta_when_gfa_sequences_missing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    assembly_fasta = tmp_path / "assembly.fa"
    plastid_fasta = tmp_path / "plastid.fa"
    output_dir = tmp_path / "out"
    work_dir = tmp_path / "blast-work"
    assembly_fasta.write_text(">contigA\nTTTT\n", encoding="utf-8")
    plastid_fasta.write_text(">plastid_ref\nACGT\n", encoding="utf-8")

    monkeypatch.setattr(internal_blast.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(
        internal_blast.tempfile,
        "TemporaryDirectory",
        lambda prefix: _FakeTemporaryDirectory(work_dir),
    )
    commands: list[list[str]] = []

    def fake_run(command: list[str], check: bool, capture_output: bool, text: bool):
        commands.append(command)
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(internal_blast.subprocess, "run", fake_run)

    internal_blast.run_internal_blast(
        assembly_fasta=assembly_fasta,
        contigs={
            "contigA": {"sequence": "", "length": 4, "coverage": None},
            "contigB": {"sequence": "", "length": 5, "coverage": None},
        },
        plastid_fasta=plastid_fasta,
        mit_fasta=None,
        output_dir=output_dir,
        threads=4,
    )

    subject_fasta = (
        work_dir.resolve()
        / internal_blast.GENERATED_SUBJECT_FASTA_NAME
    )
    assert subject_fasta.exists()
    assert subject_fasta.read_text(encoding="utf-8") == ""
    assert commands[0][2] == str(assembly_fasta.resolve())


def test_run_internal_blast_requires_assembly_fasta_when_gfa_sequences_are_unusable(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    plastid_fasta = tmp_path / "plastid.fa"
    work_dir = tmp_path / "blast-work"
    plastid_fasta.write_text(">plastid_ref\nACGT\n", encoding="utf-8")

    monkeypatch.setattr(internal_blast.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(
        internal_blast.tempfile,
        "TemporaryDirectory",
        lambda prefix: _FakeTemporaryDirectory(work_dir),
    )

    with pytest.raises(
        ValueError,
        match="Internal BLAST requires assembly contig sequences from GFA S lines or --assembly-fasta; no usable GFA sequences were available",
    ):
        internal_blast.run_internal_blast(
            assembly_fasta=None,
            contigs={
                "contigA": {"sequence": "", "length": 4, "coverage": None},
                "contigB": {"sequence": "", "length": 5, "coverage": None},
            },
            plastid_fasta=plastid_fasta,
            mit_fasta=None,
            output_dir=tmp_path / "out",
            threads=4,
        )


def test_run_internal_blast_resolves_relative_output_dir_for_makeblastdb(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    assembly_fasta = tmp_path / "assembly.fa"
    plastid_fasta = tmp_path / "plastid.fa"
    work_dir = tmp_path / "blast-work"
    assembly_fasta.write_text(">contigA\nACGT\n", encoding="utf-8")
    plastid_fasta.write_text(">plastid_ref\nACGT\n", encoding="utf-8")

    monkeypatch.setattr(internal_blast.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(
        internal_blast.tempfile,
        "TemporaryDirectory",
        lambda prefix: _FakeTemporaryDirectory(work_dir),
    )
    monkeypatch.chdir(tmp_path)
    commands: list[list[str]] = []

    def fake_run(command: list[str], check: bool, capture_output: bool, text: bool):
        assert check is False
        assert capture_output is True
        assert text is True
        commands.append(command)
        if command[0] == "makeblastdb":
            db_prefix = Path(command[command.index("-out") + 1])
            assert db_prefix.is_absolute()
            assert db_prefix.parent.exists()
            for suffix in (".nhr", ".nin", ".nsq"):
                db_prefix.with_suffix(suffix).write_text("", encoding="utf-8")
        if command[0] == "blastn":
            output_tsv = Path(command[command.index("-out") + 1])
            assert output_tsv.is_absolute()
            output_tsv.write_text("", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(internal_blast.subprocess, "run", fake_run)

    outputs = internal_blast.run_internal_blast(
        assembly_fasta=Path("assembly.fa"),
        contigs=None,
        plastid_fasta=Path("plastid.fa"),
        mit_fasta=None,
        output_dir=Path("hybrid2"),
        threads=2,
    )

    resolved_work_dir = work_dir.resolve()
    db_prefix = resolved_work_dir / internal_blast.BLAST_DB_PREFIX_NAME
    assert commands[0] == [
        "makeblastdb",
        "-in",
        str(assembly_fasta.resolve()),
        "-dbtype",
        "nucl",
        "-out",
        str(db_prefix),
    ]
    assert db_prefix.parent.exists()
    assert outputs.database_prefix == db_prefix
    assert outputs.intermediates_dir == db_prefix.parent


def test_run_internal_blast_uses_temp_work_dir_separate_from_output_dir(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    assembly_fasta = tmp_path / "assembly.fa"
    plastid_fasta = tmp_path / "plastid.fa"
    output_dir = tmp_path / "final-out"
    work_dir = tmp_path / "blast-work"
    assembly_fasta.write_text(">contigA\nACGT\n", encoding="utf-8")
    plastid_fasta.write_text(">plastid_ref\nACGT\n", encoding="utf-8")

    monkeypatch.setattr(internal_blast.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(
        internal_blast.tempfile,
        "TemporaryDirectory",
        lambda prefix: _FakeTemporaryDirectory(work_dir),
    )

    def fake_run(command: list[str], check: bool, capture_output: bool, text: bool):
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(internal_blast.subprocess, "run", fake_run)

    outputs = internal_blast.run_internal_blast(
        assembly_fasta=assembly_fasta,
        contigs=None,
        plastid_fasta=plastid_fasta,
        mit_fasta=None,
        output_dir=output_dir,
        threads=2,
    )

    assert outputs.intermediates_dir == work_dir.resolve()
    assert outputs.plastid_tsv == work_dir.resolve() / internal_blast.PLASTID_BLAST_TSV_NAME
    assert not output_dir.exists()


def test_missing_blast_executable_fails_clearly(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        internal_blast.shutil,
        "which",
        lambda name: None if name == "blastn" else f"/usr/bin/{name}",
    )

    with pytest.raises(ValueError, match="BLAST\\+ executables are required but were not found on PATH: blastn"):
        internal_blast.ensure_blast_executables_available()
