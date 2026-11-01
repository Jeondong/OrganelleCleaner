from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

from organelle_cleaner import unitig


REPO_ROOT = Path(__file__).resolve().parents[1]

# org_ctg reads sit in two isolated unitigs; nuc_ctg reads sit in linked unitigs;
# big_ctg is isolated but longer than the 2 Mb cap used without references.
CONTIG_GFA = """\
S\torg_ctg\tACGTACGTAC\tLN:i:10\trd:i:1
A\torg_ctg\t0\t+\tread_o1\t0\t5\tid:i:1\tHG:A:a
A\torg_ctg\t2\t+\tread_o2\t0\t5\tid:i:2\tHG:A:a
A\torg_ctg\t4\t+\tread_o3\t0\t5\tid:i:3\tHG:A:a
A\torg_ctg\t6\t+\tread_n9\t0\t4\tid:i:4\tHG:A:a
S\tnuc_ctg\tATATATATAT\tLN:i:10\trd:i:20
A\tnuc_ctg\t0\t+\tread_n1\t0\t5\tid:i:5\tHG:A:a
A\tnuc_ctg\t5\t+\tread_n2\t0\t5\tid:i:6\tHG:A:a
S\tbig_ctg\tGGGGCCCCGG\tLN:i:3000000\trd:i:1
A\tbig_ctg\t0\t+\tread_b1\t0\t5\tid:i:7\tHG:A:a
"""

UTG_GFA = """\
S\tutg_o1\tACGT\tLN:i:4
A\tutg_o1\t0\t+\tread_o1\t0\t4\tid:i:1\tHG:A:a
A\tutg_o1\t1\t+\tread_o2\t0\t4\tid:i:2\tHG:A:a
S\tutg_o2\tACGT\tLN:i:4
A\tutg_o2\t0\t+\tread_o3\t0\t4\tid:i:3\tHG:A:a
S\tutg_n1\tATAT\tLN:i:4
A\tutg_n1\t0\t+\tread_n1\t0\t4\tid:i:5\tHG:A:a
A\tutg_n1\t1\t+\tread_n9\t0\t4\tid:i:4\tHG:A:a
S\tutg_n2\tATAT\tLN:i:4
A\tutg_n2\t0\t+\tread_n2\t0\t4\tid:i:6\tHG:A:a
S\tutg_b1\tGGGG\tLN:i:4
A\tutg_b1\t0\t+\tread_b1\t0\t4\tid:i:7\tHG:A:a
L\tutg_n1\t+\tutg_n2\t+\t0M
"""


def _write(tmp_path: Path) -> tuple[Path, Path]:
    contig_gfa = tmp_path / "asm.p_ctg.gfa"
    utg_gfa = tmp_path / "asm.p_utg.gfa"
    contig_gfa.write_text(CONTIG_GFA, encoding="utf-8")
    utg_gfa.write_text(UTG_GFA, encoding="utf-8")
    return contig_gfa, utg_gfa


def test_iso_frac_counts_reads_in_degree_zero_unitigs(tmp_path):
    contig_gfa, utg_gfa = _write(tmp_path)
    isolated_reads, linked_reads = unitig.read_isolated_unitig_reads(utg_gfa)
    iso = unitig.calculate_iso_frac(contig_gfa, isolated_reads, linked_reads)
    assert iso["org_ctg"] == 0.75
    assert iso["nuc_ctg"] == 0.0
    assert iso["big_ctg"] == 1.0


def test_max_len_is_twice_longest_reference_capped_at_two_megabases(tmp_path):
    plastid = tmp_path / "pt.fa"
    mit = tmp_path / "mt.fa"
    plastid.write_text(">pt\n" + "A" * 150 + "\n", encoding="utf-8")
    mit.write_text(">mt_a\n" + "C" * 300 + "\n>mt_b\n" + "G" * 100 + "\n", encoding="utf-8")
    assert unitig.max_contig_length([plastid, mit]) == 800
    assert unitig.max_contig_length([None, None]) == 2_000_000
    big_mit = tmp_path / "big_mt.fa"
    big_mit.write_text(">mt\n" + "A" * 1_500_000 + "\n", encoding="utf-8")
    assert unitig.max_contig_length([plastid, big_mit]) == 2_000_000


def test_decide_applies_length_cap_and_blast_evidence():
    weak = unitig.BlastEvidence(weak=True)
    strong = unitig.BlastEvidence(strong=True)
    none = unitig.BlastEvidence()
    assert unitig.decide(mode="hybrid", length=10, iso_frac=0.6, blast=weak, max_len=100) == unitig.RULE_WEAK_BLAST_ISOLATED
    assert unitig.decide(mode="hybrid", length=1000, iso_frac=0.6, blast=weak, max_len=100) == ""
    assert unitig.decide(mode="hybrid", length=1000, iso_frac=0.0, blast=strong, max_len=100) == unitig.RULE_STRONG_BLAST
    assert unitig.decide(mode="hybrid", length=10, iso_frac=0.9, blast=none, max_len=100) == ""
    assert unitig.decide(mode="graph-only", length=10, iso_frac=0.5, blast=none, max_len=100) == unitig.RULE_ISOLATED
    assert unitig.decide(mode="graph-only", length=10, iso_frac=None, blast=none, max_len=100) == ""


def test_read_blast_evidence_uses_identity_coverage_length_and_evalue(tmp_path):
    tsv = tmp_path / "hits.tsv"
    tsv.write_text(
        # strong: 85% identity over 30% of a 10 kb contig
        "pt\tc_strong\t85.0\t3000\t10000\t0\t0\t1\t3000\t1\t3000\t0.0\n"
        # weak only: 600 bp at 82% on a 100 kb contig
        "pt\tc_weak\t82.0\t600\t100000\t0\t0\t1\t600\t5000\t5599\t1e-50\n"
        # rejected: e-value above 1e-20
        "pt\tc_evalue\t99.0\t600\t1000\t0\t0\t1\t600\t1\t600\t1e-5\n"
        # rejected: identity below 80
        "pt\tc_lowid\t75.0\t5000\t5000\t0\t0\t1\t5000\t1\t5000\t0.0\n",
        encoding="utf-8",
    )
    evidence = unitig.read_blast_evidence([tsv, None])
    assert evidence["c_strong"] == unitig.BlastEvidence(strong=True, weak=True)
    assert evidence["c_weak"] == unitig.BlastEvidence(strong=False, weak=True)
    assert "c_evalue" not in evidence
    assert "c_lowid" not in evidence


def _run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "organelle_cleaner.cli", *args],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )


def _report(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {row["contig_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def test_graph_only_with_utg_gfa_classifies_isolated_contigs_under_the_cap(tmp_path):
    contig_gfa, utg_gfa = _write(tmp_path)
    output_dir = tmp_path / "out"
    result = _run_cli(str(contig_gfa), "--utg-gfa", str(utg_gfa), "--output-dir", str(output_dir))
    assert result.returncode == 0, result.stderr
    report = _report(output_dir / "report.tsv")
    assert report["org_ctg"]["classification"] == "organelle"
    assert report["org_ctg"]["utg_rule"] == unitig.RULE_ISOLATED
    assert report["org_ctg"]["iso_frac"] == "0.75"
    assert report["nuc_ctg"]["classification"] == "nuclear"
    assert report["big_ctg"]["classification"] == "nuclear"
    assert report["big_ctg"]["utg_max_len"] == "2000000"
    organelle_ids = (output_dir / "organelle_contigs.txt").read_text(encoding="utf-8").split()
    assert organelle_ids == ["org_ctg"]


def test_without_utg_gfa_classification_is_unchanged(tmp_path):
    contig_gfa, utg_gfa = _write(tmp_path)
    plain = _run_cli(str(contig_gfa), "--output-dir", str(tmp_path / "plain"))
    assert plain.returncode == 0, plain.stderr
    report = _report(tmp_path / "plain" / "report.tsv")
    assert all(row["utg_used"] == "False" for row in report.values())
    assert all(row["iso_frac"] == "" and row["utg_rule"] == "" for row in report.values())


def test_utg_gfa_is_rejected_in_blast_only_mode(tmp_path):
    contig_gfa, utg_gfa = _write(tmp_path)
    plastid = tmp_path / "pt.fa"
    plastid.write_text(">pt\nACGT\n", encoding="utf-8")
    result = _run_cli(
        str(contig_gfa), "--mode", "blast-only", "--plastid-fasta", str(plastid),
        "--utg-gfa", str(utg_gfa), "--output-dir", str(tmp_path / "out"),
    )
    assert result.returncode != 0
    assert "--utg-gfa applies to graph-only and hybrid" in result.stderr
