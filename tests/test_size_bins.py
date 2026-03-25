from organelle_cleaner.organelle_scoring import (
    SIZE_BIN_INTERMEDIATE,
    SIZE_BIN_LARGE,
    ScoringConfig,
    determine_size_bin,
    score_contigs,
)


def test_large_contig_boundary_respects_configured_lower_bound():
    config = ScoringConfig(
        intermediate_max_length=200_000,
        large_contig_length=300_000,
        very_large_contig_length=500_000,
    )

    assert determine_size_bin(200_000, config) == SIZE_BIN_INTERMEDIATE
    assert determine_size_bin(200_001, config) == SIZE_BIN_INTERMEDIATE
    assert determine_size_bin(300_000, config) == SIZE_BIN_LARGE
    assert determine_size_bin(300_001, config) == SIZE_BIN_LARGE


def test_score_contigs_defaults_to_graph_only_mode():
    rows = score_contigs(
        {
            "contigA": {
                "sequence": "ACGTACGT",
                "length": 8,
                "coverage": 8,
            }
        },
        graph=None,
    )

    assert len(rows) == 1
    assert rows[0]["mode"] == "graph-only"
