#!/usr/bin/env python3
"""Unitig-isolation evidence from a hifiasm unitig graph (p_utg.gfa or r_utg.gfa).

Organelle reads are sequenced far deeper than nuclear reads, and hifiasm tends to
build them into many small unitigs that have no links to the rest of the unitig
graph. For each contig, iso_frac is the fraction of its reads (A-lines of the
contig GFA) that sit in unitigs with no L-line at all (degree 0).

Decision rules used when --utg-gfa is given (thresholds fixed in advance):
    hybrid:     strong BLAST hit
                OR (weak BLAST hit AND iso_frac >= 0.5 AND length <= max_len)
    graph-only: iso_frac >= 0.5 AND length <= max_len
where a strong hit is identity >= 80 and per-HSP contig coverage >= 0.20, a weak
hit is identity >= 80 and alignment length >= 500 (both with e-value <= 1e-20),
and max_len = min(2 x the longest organelle reference genome, 2 Mb); 2 Mb without references.
The 2 Mb ceiling keeps a reference with a very large mitogenome (e.g. Amborella,
3.9 Mb) from admitting multi-megabase nuclear contigs through the weak rule.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

ISO_FRAC_THRESHOLD = 0.5
STRONG_MIN_IDENTITY = 80.0
STRONG_MIN_CONTIG_COVERAGE = 0.20
WEAK_MIN_IDENTITY = 80.0
WEAK_MIN_ALIGNMENT_LENGTH = 500
MAX_EVALUE = 1e-20
MAX_LEN_REFERENCE_MULTIPLIER = 2
DEFAULT_REFERENCE_GENOME_LENGTH = 1_000_000
MAX_LEN_CEILING = 2_000_000

RULE_STRONG_BLAST = "strong_blast"
RULE_WEAK_BLAST_ISOLATED = "weak_blast+isolated_unitigs"
RULE_ISOLATED = "isolated_unitigs"


@dataclass(frozen=True)
class BlastEvidence:
    strong: bool = False
    weak: bool = False


def read_isolated_unitig_reads(utg_gfa: Path) -> tuple[set[str], set[str]]:
    """Return (reads in degree-0 unitigs, reads in linked unitigs) from a unitig GFA."""
    linked_unitigs: set[str] = set()
    read_to_unitig: dict[str, str] = {}
    with utg_gfa.open("r", encoding="utf-8") as handle:
        for line in handle:
            record_type = line[:2]
            if record_type == "L\t":
                fields = line.split("\t", 4)
                linked_unitigs.add(fields[1])
                linked_unitigs.add(fields[3])
            elif record_type == "A\t":
                fields = line.split("\t", 5)
                read_to_unitig[fields[4]] = fields[1]
    isolated_reads: set[str] = set()
    linked_reads: set[str] = set()
    for read_name, unitig in read_to_unitig.items():
        if unitig in linked_unitigs:
            linked_reads.add(read_name)
        else:
            isolated_reads.add(read_name)
    return isolated_reads, linked_reads


def calculate_iso_frac(
    contig_gfa: Path,
    isolated_reads: set[str],
    linked_reads: set[str],
) -> dict[str, float]:
    """iso_frac per contig over its reads that are present in the unitig graph."""
    mapped: dict[str, int] = {}
    isolated: dict[str, int] = {}
    with contig_gfa.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("A\t"):
                continue
            fields = line.split("\t", 5)
            contig_id, read_name = fields[1], fields[4]
            if read_name in isolated_reads:
                isolated[contig_id] = isolated.get(contig_id, 0) + 1
            elif read_name not in linked_reads:
                continue
            mapped[contig_id] = mapped.get(contig_id, 0) + 1
    return {contig_id: isolated.get(contig_id, 0) / count for contig_id, count in mapped.items()}


def reference_genome_length(fasta: Path | None) -> int:
    """Total sequence length of one organelle reference FASTA (all records = one genome)."""
    if fasta is None:
        return 0
    total = 0
    with fasta.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def max_contig_length(reference_fastas: list[Path | None]) -> int:
    longest = max((reference_genome_length(path) for path in reference_fastas), default=0)
    if longest <= 0:
        longest = DEFAULT_REFERENCE_GENOME_LENGTH
    return min(MAX_LEN_REFERENCE_MULTIPLIER * longest, MAX_LEN_CEILING)


def read_blast_evidence(blast_tsvs: list[Path | None]) -> dict[str, BlastEvidence]:
    """Strong/weak hit flags per contig from internal BLAST TSVs (organelle query, contig subject).

    Columns: qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send evalue.
    """
    evidence: dict[str, BlastEvidence] = {}
    for path in blast_tsvs:
        if path is None or not path.exists():
            continue
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12 or fields[0].startswith("#"):
                    continue
                try:
                    identity = float(fields[2])
                    alignment_length = int(fields[3])
                    subject_length = int(fields[4])
                    subject_span = abs(int(fields[10]) - int(fields[9])) + 1
                    evalue = float(fields[11])
                except ValueError:
                    continue
                if evalue > MAX_EVALUE or subject_length <= 0:
                    continue
                coverage = min(1.0, max(subject_span, alignment_length) / subject_length)
                strong = identity >= STRONG_MIN_IDENTITY and coverage >= STRONG_MIN_CONTIG_COVERAGE
                weak = identity >= WEAK_MIN_IDENTITY and alignment_length >= WEAK_MIN_ALIGNMENT_LENGTH
                if not (strong or weak):
                    continue
                previous = evidence.get(fields[1], BlastEvidence())
                evidence[fields[1]] = BlastEvidence(
                    strong=previous.strong or strong,
                    weak=previous.weak or weak,
                )
    return evidence


def decide(
    *,
    mode: str,
    length: int,
    iso_frac: float | None,
    blast: BlastEvidence,
    max_len: int,
) -> str:
    """Return the rule that classifies the contig as organelle, or "" for nuclear."""
    isolated = iso_frac is not None and iso_frac >= ISO_FRAC_THRESHOLD and length <= max_len
    if mode == "hybrid":
        if blast.strong:
            return RULE_STRONG_BLAST
        if blast.weak and isolated:
            return RULE_WEAK_BLAST_ISOLATED
        return ""
    if mode == "graph-only":
        return RULE_ISOLATED if isolated else ""
    raise ValueError(f"Unitig rules support graph-only and hybrid modes, not {mode!r}")


def apply_unitig_rules(
    report_rows: list[dict[str, object]],
    *,
    mode: str,
    contig_gfa: Path,
    utg_gfa: Path,
    blast_tsvs: list[Path | None],
    reference_fastas: list[Path | None],
) -> None:
    """Replace classification with the unitig rules and record the evidence in each row."""
    isolated_reads, linked_reads = read_isolated_unitig_reads(utg_gfa)
    iso_by_contig = calculate_iso_frac(contig_gfa, isolated_reads, linked_reads)
    blast_by_contig = read_blast_evidence(blast_tsvs) if mode == "hybrid" else {}
    max_len = max_contig_length(reference_fastas)
    for row in report_rows:
        contig_id = str(row["contig_id"])
        length_value = row.get("contig_length")
        length = length_value if isinstance(length_value, int) else 0
        iso_frac = iso_by_contig.get(contig_id)
        rule = decide(
            mode=mode,
            length=length,
            iso_frac=iso_frac,
            blast=blast_by_contig.get(contig_id, BlastEvidence()),
            max_len=max_len,
        )
        row["utg_used"] = True
        row["iso_frac"] = None if iso_frac is None else round(iso_frac, 4)
        row["utg_max_len"] = max_len
        row["utg_rule"] = rule
        row["classification"] = "organelle" if rule else "nuclear"
