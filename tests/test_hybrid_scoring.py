from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from organelle_scoring import (
    BLAST_SUPPORT_NONE,
    HIGH_CONFIDENCE,
    MEDIUM_CONFIDENCE,
    MODE_HYBRID,
    ScoringConfig,
    SIZE_BIN_PREFERRED_RANGE,
    determine_graph_only_tier,
    determine_hybrid_tier,
    evaluate_contig,
)


def _preferred_range_singleton_graph_features():
    return {
        "graph_score": 5,
        "support_signal_count": 4,
        "strong_signal_count": 0,
        "topology_signal_count": 2,
        "topology_strong_signal_count": 0,
        "size_bin": SIZE_BIN_PREFERRED_RANGE,
        "is_compact_component": True,
        "has_gc_signal": True,
        "has_coverage_signal": False,
        "has_preferred_length_signal": True,
        "has_preferred_range_promotion_pattern": True,
        "is_large": False,
    }


def _no_blast_features():
    return {
        "blast_support_level": BLAST_SUPPORT_NONE,
        "blast_best_identity": None,
        "blast_best_aligned_fraction": None,
        "blast_merged_aligned_bp": 0,
        "blast_merged_coverage_fraction": 0.0,
        "blast_support_sources": "",
        "blast_selected_source": None,
        "plastid_support_level": BLAST_SUPPORT_NONE,
        "plastid_best_identity": None,
        "plastid_merged_aligned_bp": 0,
        "plastid_merged_coverage_fraction": 0.0,
        "mit_support_level": BLAST_SUPPORT_NONE,
        "mit_best_identity": None,
        "mit_merged_aligned_bp": 0,
        "mit_merged_coverage_fraction": 0.0,
        "blast_score": 0,
        "has_strong_blast_support": False,
        "has_moderate_blast_support": False,
        "has_exceptional_blast_support": False,
        "blast_support_reasons": [],
    }


def test_graph_only_keeps_preferred_range_singleton_promotion():
    config = ScoringConfig()

    tier, _ = determine_graph_only_tier(_preferred_range_singleton_graph_features(), config)

    assert tier == HIGH_CONFIDENCE


def test_hybrid_blocks_blast_negative_singleton_promotion_without_strong_topology():
    config = ScoringConfig()
    graph_features = _preferred_range_singleton_graph_features()
    graph_only_tier, _ = determine_graph_only_tier(graph_features, config)

    tier, notes = determine_hybrid_tier(
        size_bin=SIZE_BIN_PREFERRED_RANGE,
        graph_features=graph_features,
        blast_features=_no_blast_features(),
        graph_only_tier=graph_only_tier,
        config=config,
    )

    assert graph_only_tier == HIGH_CONFIDENCE
    assert tier == MEDIUM_CONFIDENCE
    assert any("lacked strong topology" in note for note in notes)


def test_hybrid_keeps_blast_supported_preferred_range_recovery():
    config = ScoringConfig()
    row = evaluate_contig(
        contig_id="ptg-test",
        mode=MODE_HYBRID,
        feature_row={"length": 120000, "coverage": 30.0, "gc_content": 40.0},
        component_row={
            "component_size": 1,
            "component_edge_count": 0,
            "component_branching_nodes": 0,
            "component_total_length": 120000,
            "is_circular": False,
            "is_isolated": True,
            "has_self_loop": False,
            "is_branching": False,
        },
        blast_row={
            "blast_support_level": "strong_blast_support",
            "blast_best_identity": 99.5,
            "blast_best_aligned_fraction": 0.95,
            "blast_merged_aligned_bp": 114000,
            "blast_merged_coverage_fraction": 0.95,
            "blast_support_sources": "plastid",
            "blast_selected_source": "plastid",
            "plastid_support_level": "strong_blast_support",
            "plastid_best_identity": 99.5,
            "plastid_merged_aligned_bp": 114000,
            "plastid_merged_coverage_fraction": 0.95,
            "mit_support_level": BLAST_SUPPORT_NONE,
            "mit_best_identity": None,
            "mit_merged_aligned_bp": 0,
            "mit_merged_coverage_fraction": 0.0,
            "blast_score": 4,
            "has_strong_blast_support": True,
            "has_moderate_blast_support": True,
            "has_exceptional_blast_support": False,
            "blast_supporting_evidence": "plastid BLAST support",
        },
        nuclear_coverage=10.0,
        gc_baseline=35.0,
        config=config,
    )

    assert row["confidence_tier"] == HIGH_CONFIDENCE
    assert row["classification"] == "organelle"
