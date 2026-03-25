from organelle_cleaner.blast_features import (
    BLAST_SUPPORT_MODERATE,
    BLAST_SUPPORT_STRONG,
    load_blast_support_by_contig,
)


def test_blast_support_attaches_to_subject_contig_id(tmp_path):
    blast_tsv = tmp_path / "plastid.tsv"
    blast_tsv.write_text(
        "qseqid\tsseqid\tpident\tlength\tslen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\n"
        "plastid_ref\tcontigA\t99.5\t100\t100\t0\t0\t1\t100\t1\t100\n",
        encoding="utf-8",
    )

    rows = load_blast_support_by_contig(plastid_blast_tsv=blast_tsv)

    assert "contigA" in rows
    assert "plastid_ref" not in rows
    assert rows["contigA"]["blast_support_level"] == BLAST_SUPPORT_STRONG
    assert rows["contigA"]["blast_selected_source"] == "plastid"


def test_moderate_support_uses_90_percent_identity_and_80_percent_coverage_defaults(tmp_path):
    blast_tsv = tmp_path / "plastid.tsv"
    blast_tsv.write_text(
        "qseqid\tsseqid\tpident\tlength\tslen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\n"
        "plastid_ref\tcontigA\t90.0\t80\t100\t0\t0\t1\t80\t1\t80\n",
        encoding="utf-8",
    )

    rows = load_blast_support_by_contig(plastid_blast_tsv=blast_tsv)

    assert rows["contigA"]["blast_support_level"] == BLAST_SUPPORT_MODERATE
