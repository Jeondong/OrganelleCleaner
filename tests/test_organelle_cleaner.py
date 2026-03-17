#!/usr/bin/env python3
"""Regression tests for organelle_cleaner output consistency."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_PATH = REPO_ROOT / "organelle_cleaner.py"

SYNTHETIC_GFA = """\
S\tmissing_seq\t*\tLN:i:10
S\tnuclear_a\tACGTACGT\tLN:i:8\tRC:i:8
S\tnuclear_b\tATATATAT\tLN:i:8\tRC:i:8
"""


def fasta_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].strip().split()[0])
    return ids


def contig_id_list(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


def read_report(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class OrganelleCleanerContractTests(unittest.TestCase):
    def run_pipeline(self) -> tuple[Path, subprocess.CompletedProcess[str]]:
        temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(temp_dir.cleanup)
        temp_path = Path(temp_dir.name)
        input_path = temp_path / "synthetic.gfa"
        output_dir = temp_path / "out"
        input_path.write_text(SYNTHETIC_GFA, encoding="utf-8")

        result = subprocess.run(
            [sys.executable, str(PIPELINE_PATH), str(input_path), "--output-dir", str(output_dir)],
            check=True,
            capture_output=True,
            text=True,
        )
        return output_dir, result

    def test_fasta_matches_nuclear_contig_ids(self) -> None:
        output_dir, _ = self.run_pipeline()

        self.assertEqual(
            fasta_ids(output_dir / "cleaned_assembly.fa"),
            contig_id_list(output_dir / "nuclear_contigs.txt"),
        )

    def test_missing_sequence_contig_is_reported_but_excluded_from_nuclear_outputs(self) -> None:
        output_dir, _ = self.run_pipeline()
        report_rows = {row["contig_id"]: row for row in read_report(output_dir / "report.tsv")}

        self.assertIn("missing_seq", report_rows)
        self.assertEqual(report_rows["missing_seq"]["has_sequence"], "False")
        self.assertNotIn("missing_seq", fasta_ids(output_dir / "cleaned_assembly.fa"))
        self.assertNotIn("missing_seq", contig_id_list(output_dir / "nuclear_contigs.txt"))

    def test_all_nuclear_contigs_with_sequence_appear_in_fasta(self) -> None:
        output_dir, _ = self.run_pipeline()
        report_rows = read_report(output_dir / "report.tsv")
        fasta_contigs = fasta_ids(output_dir / "cleaned_assembly.fa")

        expected_ids = {
            row["contig_id"]
            for row in report_rows
            if row["classification"] == "nuclear" and row["has_sequence"] == "True"
        }

        self.assertTrue(expected_ids)
        self.assertTrue(expected_ids.issubset(fasta_contigs))


if __name__ == "__main__":
    unittest.main()
