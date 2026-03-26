from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

from organelle_cleaner import cli, internal_blast


REPO_ROOT = Path(__file__).resolve().parents[1]

SYNTHETIC_GFA = """\
S\tmissing_seq\t*\tLN:i:10
S\tnuclear_a\tACGTACGT\tLN:i:8\tRC:i:8
S\tnuclear_b\tATATATAT\tLN:i:8\tRC:i:8
"""


def _run_pipeline(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "organelle_cleaner.cli", *args],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )


def _write_synthetic_gfa(tmp_path: Path) -> Path:
    input_path = tmp_path / "synthetic.gfa"
    input_path.write_text(SYNTHETIC_GFA, encoding="utf-8")
    return input_path


def _write_synthetic_fasta(tmp_path: Path) -> Path:
    assembly_path = tmp_path / "synthetic.fa"
    assembly_path.write_text(
        ">missing_seq\nACGTACGTAC\n>nuclear_a\nACGTACGT\n>nuclear_b\nATATATAT\n",
        encoding="utf-8",
    )
    return assembly_path


def _read_report(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _read_fasta_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].strip().split()[0])
    return ids


def _read_id_list(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


def test_default_cli_invocation_succeeds_without_blast_and_keeps_outputs_consistent(tmp_path):
    input_path = _write_synthetic_gfa(tmp_path)
    output_dir = tmp_path / "out"

    result = _run_pipeline(str(input_path), "--output-dir", str(output_dir))

    assert result.returncode == 0, result.stderr
    assert (output_dir / "report.tsv").exists()
    assert (output_dir / "nuclear_contigs.txt").exists()
    assert (output_dir / "cleaned_assembly.fa").exists()

    report_rows = _read_report(output_dir / "report.tsv")
    report_by_id = {row["contig_id"]: row for row in report_rows}
    fasta_ids = _read_fasta_ids(output_dir / "cleaned_assembly.fa")
    nuclear_ids = _read_id_list(output_dir / "nuclear_contigs.txt")

    assert report_by_id["missing_seq"]["has_sequence"] == "False"
    assert report_by_id["nuclear_a"]["mode"] == "graph-only"
    assert fasta_ids == nuclear_ids
    assert "missing_seq" not in fasta_ids
    assert "missing_seq" not in nuclear_ids


def test_cli_writes_candidate_report_when_all_candidates_name_is_requested(tmp_path):
    input_path = _write_synthetic_gfa(tmp_path)
    output_dir = tmp_path / "out"

    result = _run_pipeline(
        str(input_path),
        "--output-dir",
        str(output_dir),
        "--all-candidates-name",
        "custom_candidates.tsv",
    )

    assert result.returncode == 0, result.stderr
    assert (output_dir / "custom_candidates.tsv").exists()


def test_help_hides_advanced_fallback_options_but_parser_still_accepts_them():
    parser = cli.build_parser()
    help_text = parser.format_help()

    assert "--assembly-fasta" not in help_text
    assert "--all-candidates-name" not in help_text

    args = parser.parse_args(
        [
            "assembly.gfa",
            "--assembly-fasta",
            "assembly.fa",
            "--all-candidates-name",
            "custom_candidates.tsv",
        ]
    )
    assert args.assembly_fasta == Path("assembly.fa")
    assert args.all_candidates_name == "custom_candidates.tsv"


def test_hybrid_and_blast_only_require_blast_inputs(tmp_path):
    input_path = _write_synthetic_gfa(tmp_path)

    hybrid_result = _run_pipeline(str(input_path), "--mode", "hybrid")
    blast_only_result = _run_pipeline(str(input_path), "--mode", "blast-only")

    assert hybrid_result.returncode != 0
    assert "requires at least one organelle FASTA input" in hybrid_result.stderr
    assert blast_only_result.returncode != 0
    assert "requires at least one organelle FASTA input" in blast_only_result.stderr


def test_hybrid_accepts_organelle_fasta_without_assembly_fasta(tmp_path, monkeypatch):
    input_path = _write_synthetic_gfa(tmp_path)
    plastid_fasta = tmp_path / "plastid.fa"
    plastid_fasta.write_text(">plastid_ref\nACGTACGT\n", encoding="utf-8")
    output_dir = tmp_path / "hybrid_out"

    captured: dict[str, Path | None] = {}

    def fake_run_internal_blast(**kwargs):
        captured["assembly_fasta"] = kwargs["assembly_fasta"]
        captured["plastid_fasta"] = kwargs["plastid_fasta"]
        captured["output_dir"] = kwargs["output_dir"]
        return internal_blast.InternalBlastOutputs(
            plastid_tsv=None,
            mit_tsv=None,
            database_prefix=kwargs["output_dir"] / "db",
            intermediates_dir=kwargs["output_dir"] / "blast_intermediates",
        )

    monkeypatch.setattr(cli, "run_internal_blast", fake_run_internal_blast)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "organelle-cleaner",
            str(input_path),
            "--mode",
            "hybrid",
            "--plastid-fasta",
            str(plastid_fasta),
            "--output-dir",
            str(output_dir),
        ],
    )

    assert cli.main() == 0
    assert captured["assembly_fasta"] is None
    assert captured["plastid_fasta"] == plastid_fasta
    assert (output_dir / "report.tsv").exists()


def test_blast_only_accepts_assembly_fasta_without_gfa(tmp_path, monkeypatch):
    assembly_fasta = _write_synthetic_fasta(tmp_path)
    plastid_fasta = tmp_path / "plastid.fa"
    plastid_fasta.write_text(">plastid_ref\nACGTACGT\n", encoding="utf-8")
    output_dir = tmp_path / "blast_only_out"

    captured: dict[str, Path | None] = {}

    def fake_run_internal_blast(**kwargs):
        captured["assembly_fasta"] = kwargs["assembly_fasta"]
        captured["plastid_fasta"] = kwargs["plastid_fasta"]
        return internal_blast.InternalBlastOutputs(
            plastid_tsv=None,
            mit_tsv=None,
            database_prefix=kwargs["output_dir"] / "db",
            intermediates_dir=kwargs["output_dir"] / "blast_intermediates",
        )

    monkeypatch.setattr(cli, "run_internal_blast", fake_run_internal_blast)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "organelle-cleaner",
            "--assembly-fasta",
            str(assembly_fasta),
            "--mode",
            "blast-only",
            "--plastid-fasta",
            str(plastid_fasta),
            "--output-dir",
            str(output_dir),
        ],
    )

    assert cli.main() == 0
    assert captured["assembly_fasta"] == assembly_fasta
    assert captured["plastid_fasta"] == plastid_fasta
    assert (output_dir / "report.tsv").exists()
    report_rows = _read_report(output_dir / "report.tsv")
    assert {row["contig_id"] for row in report_rows} == {
        "missing_seq",
        "nuclear_a",
        "nuclear_b",
    }


def test_blast_only_requires_a_contig_source(tmp_path):
    plastid_fasta = tmp_path / "plastid.fa"
    plastid_fasta.write_text(">plastid_ref\nACGTACGT\n", encoding="utf-8")

    result = _run_pipeline(
        "--mode",
        "blast-only",
        "--plastid-fasta",
        str(plastid_fasta),
    )

    assert result.returncode != 0
    assert "blast-only mode requires either input_gfa or --assembly-fasta" in result.stderr


def test_graph_only_ignores_organelle_fasta_without_requiring_internal_blast(tmp_path):
    input_path = _write_synthetic_gfa(tmp_path)
    plastid_fasta = tmp_path / "plastid.fa"
    plastid_fasta.write_text(">plastid_ref\nACGTACGT\n", encoding="utf-8")
    output_dir = tmp_path / "graph_only_out"

    result = _run_pipeline(
        str(input_path),
        "--mode",
        "graph-only",
        "--plastid-fasta",
        str(plastid_fasta),
        "--output-dir",
        str(output_dir),
    )

    assert result.returncode == 0, result.stderr
    assert (output_dir / "report.tsv").exists()
    assert not (output_dir / "blast_intermediates").exists()
