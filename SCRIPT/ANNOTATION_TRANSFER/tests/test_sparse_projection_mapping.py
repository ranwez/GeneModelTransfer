from __future__ import annotations

import pytest

from locus_alignment_transfer import (
    CdsFeature,
    ProjectionError,
    map_requested_model_positions,
    projection_candidate_positions,
)


def test_projection_candidates_are_only_endpoints_and_in_frame_fallbacks():
    features = [CdsFeature(1, 20, "+", "0"), CdsFeature(30, 35, "+", "0")]

    assert projection_candidate_positions(features) == {
        1,
        4,
        7,
        10,
        11,
        14,
        17,
        20,
        30,
        32,
        33,
        35,
    }


def test_projection_candidates_never_cross_a_short_cds_boundary():
    assert projection_candidate_positions([CdsFeature(10, 14, "+", "0")]) == {
        10,
        11,
        13,
        14,
    }


def test_mapper_does_not_retain_unrequested_model_positions():
    requested = {1, 500, 1000}
    model_alignment = "A" * 1000
    target_alignment = "A" * 1000

    mapping = map_requested_model_positions(
        model_alignment,
        target_alignment,
        requested,
    )

    assert mapping == {1: 1, 500: 500, 1000: 1000}
    assert mapping.keys() == requested


def test_sparse_mapper_records_target_gap_only_at_requested_position():
    mapping = map_requested_model_positions("ACGT", "A-GT", {2, 4})

    assert mapping == {2: None, 4: 3}


@pytest.mark.parametrize("requested", [{0}, {5}, {-1, 1}])
def test_requested_positions_must_be_inside_ungapped_model(requested):
    with pytest.raises(ProjectionError, match="requested"):
        map_requested_model_positions("ACGT", "ACGT", requested)
