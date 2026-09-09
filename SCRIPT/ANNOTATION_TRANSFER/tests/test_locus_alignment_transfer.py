from __future__ import annotations

from pathlib import Path

import pytest

from locus_alignment_transfer import (
    AlignmentError,
    BlastHsp,
    CdsFeature,
    InputValidationError,
    ProjectionError,
    assess_eligibility,
    map_requested_model_positions,
    normalize_model_gff,
    parse_blast_tabular,
    parse_model_gff,
    project_cds_features,
    read_single_fasta,
    render_draft_gff,
    reverse_coordinate,
)


def hsp(**overrides) -> BlastHsp:
    values = {
        "query_id": "model",
        "subject_id": "target",
        "query_length": 100,
        "alignment_length": 85,
        "query_start": 1,
        "query_end": 85,
        "subject_start": 11,
        "subject_end": 95,
        "identical": 81,
        "percent_identity": 95.0,
        "bit_score": 200.0,
        "evalue": 1e-50,
    }
    values.update(overrides)
    return BlastHsp(**values)


def write(path: Path, content: str) -> Path:
    path.write_text(content)
    return path


class TestBlastEligibility:
    def test_inclusive_thresholds_pass_without_rounding(self):
        decision = assess_eligibility([hsp()])
        assert decision.eligible is True
        assert decision.coverage == pytest.approx(0.85)
        assert decision.identity == pytest.approx(0.95)
        assert decision.orientation == "forward"

    @pytest.mark.parametrize(
        "candidate,reason",
        [
            (hsp(alignment_length=84), "coverage_below_threshold"),
            (hsp(percent_identity=94.999999), "identity_below_threshold"),
        ],
    )
    def test_value_immediately_below_either_threshold_fails(self, candidate, reason):
        decision = assess_eligibility([candidate])
        assert decision.eligible is False
        assert reason in decision.reasons

    def test_no_significant_hit_is_ineligible(self):
        decision = assess_eligibility([])
        assert decision.eligible is False
        assert decision.selected_hsp is None
        assert decision.reasons == ("no_significant_hit",)

    def test_best_hsp_ranking_is_deterministic(self):
        candidates = [
            hsp(bit_score=199, evalue=0.0, alignment_length=100),
            hsp(bit_score=200, evalue=1e-10, alignment_length=84),
            hsp(bit_score=200, evalue=1e-20, alignment_length=83),
            hsp(bit_score=200, evalue=1e-20, alignment_length=85, percent_identity=96),
        ]
        decision = assess_eligibility(candidates)
        assert decision.selected_hsp == candidates[3]

    def test_final_tie_break_uses_coordinates_not_input_order(self):
        later = hsp(query_start=2, query_end=86, subject_start=20, subject_end=104)
        earlier = hsp(query_start=1, query_end=85, subject_start=30, subject_end=114)
        assert assess_eligibility([later, earlier]).selected_hsp == earlier
        assert assess_eligibility([earlier, later]).selected_hsp == earlier

    def test_reverse_subject_coordinates_set_reverse_orientation(self):
        decision = assess_eligibility([hsp(subject_start=95, subject_end=11)])
        assert decision.orientation == "reverse"

    def test_blast_parser_uses_v3_column_order_and_rejects_bad_rows(self):
        text = "model\ttarget\t100\t85\t1\t85\t11\t95\t81\t95.0\t200\t1e-50\n"
        assert parse_blast_tabular(text) == [hsp()]
        with pytest.raises(AlignmentError, match="BLAST"):
            parse_blast_tabular("model\ttarget\ttoo-few-columns\n")


class TestCoordinateProjection:
    def test_no_gap_alignment_maps_only_requested_positions(self):
        assert map_requested_model_positions("ACGT", "ACGT", {1, 4}) == {1: 1, 4: 4}

    def test_target_gap_is_explicitly_unmappable(self):
        assert map_requested_model_positions("ACGT", "A-GT", {1, 2, 4}) == {1: 1, 2: None, 4: 3}

    def test_model_gap_advances_only_target_coordinate(self):
        assert map_requested_model_positions("A-CGT", "ATCGT", {1, 2, 4}) == {1: 1, 2: 3, 4: 5}

    def test_alignment_lengths_and_ungapped_inputs_are_validated(self):
        with pytest.raises(AlignmentError, match="same length"):
            map_requested_model_positions("ACGT", "ACG", {1, 4})
        with pytest.raises(AlignmentError, match="model"):
            map_requested_model_positions("ACGT", "ACGT", {1, 4}, expected_model="ACGA")
        with pytest.raises(AlignmentError, match="target"):
            map_requested_model_positions("ACGT", "ACGT", {1, 4}, expected_target="ACGA")

    def test_first_last_and_one_base_cds_have_no_off_by_one_error(self):
        mapping = map_requested_model_positions("ACGT", "ACGT", {1, 4})
        projected = project_cds_features(
            [CdsFeature(1, 1, "+", "0"), CdsFeature(4, 4, "+", "0")],
            mapping,
            target_length=4,
            orientation="forward",
        )
        assert [(item.start, item.end) for item in projected] == [(1, 1), (4, 4)]

    def test_internal_target_gap_does_not_change_mapped_boundaries(self):
        mapping = map_requested_model_positions("ACGTAC", "AC-TAC", {2, 5})
        projected = project_cds_features(
            [CdsFeature(2, 5, "+", "0")], mapping, target_length=5, orientation="forward"
        )
        assert (projected[0].start, projected[0].end) == (2, 4)

    def test_start_boundary_gap_snaps_inward_by_whole_codons_only(self):
        mapping = {1: None, 4: 3, 7: 6, 10: 9}
        projected = project_cds_features(
            [CdsFeature(1, 10, "+", "0")], mapping, target_length=9, orientation="forward"
        )
        assert (projected[0].start, projected[0].end) == (3, 9)

    def test_end_boundary_gap_snaps_inward_by_whole_codons_only(self):
        mapping = {1: 1, 4: 4, 7: 7, 10: None}
        projected = project_cds_features(
            [CdsFeature(1, 10, "+", "0")], mapping, target_length=9, orientation="forward"
        )
        assert (projected[0].start, projected[0].end) == (1, 7)

    @pytest.mark.parametrize("missing", [[1, 4, 7, 10], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]])
    def test_boundary_gap_rejects_when_no_same_frame_base_within_three_codons(self, missing):
        mapping = {position: (None if position in missing else position) for position in (1, 4, 7, 10)}
        with pytest.raises(ProjectionError, match="boundary"):
            project_cds_features(
                [CdsFeature(1, 10, "+", "0")], mapping, target_length=10, orientation="forward"
            )

    def test_reverse_coordinate_formula_and_projection_apply_exactly_once(self):
        assert reverse_coordinate(1, 10) == 10
        assert reverse_coordinate(10, 10) == 1
        projected = project_cds_features(
            [CdsFeature(2, 5, "+", "0")],
            {2: 2, 5: 5},
            target_length=10,
            orientation="reverse",
        )
        assert (projected[0].start, projected[0].end, projected[0].strand) == (6, 9, "-")

    def test_multi_cds_output_is_genomically_ordered_on_both_orientations(self):
        features = [CdsFeature(7, 9, "+", "0"), CdsFeature(1, 3, "+", "0")]
        mapping = {1: 1, 3: 3, 7: 7, 9: 9}
        forward = project_cds_features(features, mapping, target_length=9, orientation="forward")
        reverse = project_cds_features(features, mapping, target_length=9, orientation="reverse")
        assert [(x.start, x.end, x.strand) for x in forward] == [(1, 3, "+"), (7, 9, "+")]
        assert [(x.start, x.end, x.strand) for x in reverse] == [(1, 3, "-"), (7, 9, "-")]

    def test_overlapping_or_reordered_projection_is_rejected(self):
        with pytest.raises(ProjectionError):
            project_cds_features(
                [CdsFeature(1, 3, "+", "0"), CdsFeature(4, 6, "+", "0")],
                {1: 1, 3: 4, 4: 4, 6: 6},
                target_length=6,
                orientation="forward",
            )


class TestInputAndGffContracts:
    def test_single_fasta_accepts_lowercase_and_iupac(self, tmp_path):
        record = read_single_fasta(write(tmp_path / "one.fa", ">model\nacgtryswkmbdhvn\n"))
        assert record.identifier == "model"
        assert record.sequence == "ACGTRYSWKMBDHVN"

    @pytest.mark.parametrize(
        "content,message",
        [
            ("", "one record"),
            (">a\nACG\n>b\nACG\n", "one record"),
            (">a\n", "non-empty"),
            (">a\nACGZ\n", "IUPAC"),
        ],
    )
    def test_invalid_fasta_is_rejected(self, tmp_path, content, message):
        with pytest.raises(InputValidationError, match=message):
            read_single_fasta(write(tmp_path / "bad.fa", content))

    def test_negative_genome_annotation_is_normalized_to_transcription_orientation(self, tmp_path):
        source = write(
            tmp_path / "model.gff",
            "Chr4\tsrc\tgene\t100\t199\t.\t-\t.\tID=model\n"
            "Chr4\tsrc\tmRNA\t100\t199\t.\t-\t.\tID=model_rna;Parent=model\n"
            "Chr4\tsrc\tCDS\t100\t129\t.\t-\t0\tID=c2;Parent=model_rna\n"
            "Chr4\tsrc\tCDS\t170\t199\t.\t-\t0\tID=c1;Parent=model_rna\n",
        )
        normalized = normalize_model_gff(source, model_id="model", locus_length=100)
        assert [(cds.start, cds.end, cds.strand) for cds in normalized.cds] == [
            (1, 30, "+"),
            (71, 100, "+"),
        ]

    def test_positive_genome_annotation_uses_gene_start_as_locus_origin(self, tmp_path):
        source = write(
            tmp_path / "model.gff",
            "Chr1\tsrc\tgene\t100\t199\t.\t+\t.\tID=model\n"
            "Chr1\tsrc\tmRNA\t100\t199\t.\t+\t.\tID=model_rna;Parent=model\n"
            "Chr1\tsrc\tCDS\t110\t139\t.\t+\t0\tID=c1;Parent=model_rna\n",
        )
        normalized = normalize_model_gff(source, model_id="model", locus_length=100)
        assert [(cds.start, cds.end, cds.strand) for cds in normalized.cds] == [(11, 40, "+")]

    @pytest.mark.parametrize(
        "body",
        [
            "seq\tsrc\tgene\t1\t12\t.\t+\t.\tID=g\nseq\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=t;Parent=g\n",
            "seq\tsrc\tgene\t1\t12\t.\t+\t.\tID=g\nseq\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g\nseq\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=t2;Parent=g\nseq\tsrc\tCDS\t1\t3\t.\t+\t0\tParent=t1\n",
            "seq\tsrc\tgene\t1\t12\t.\t+\t.\tID=g\nseq\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=t;Parent=g\nseq\tsrc\tCDS\t1\t6\t.\t+\t0\tParent=t\nseq\tsrc\tCDS\t6\t9\t.\t+\t0\tParent=t\n",
            "seq\tsrc\tgene\t1\t12\t.\t+\t.\tID=g\nseq\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=t;Parent=g\nseq\tsrc\tCDS\t1\t3\t.\t-\t0\tParent=t\n",
            "seq\tsrc\tgene\t1\t12\t.\t+\t.\tID=g\nseq\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=t;Parent=g\nseq\tsrc\tCDS\t10\t13\t.\t+\t0\tParent=t\n",
            "seq\tsrc\tgene\t1\t12\t.\t+\t.\tID=g\nseq\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=t;Parent=g\nseq\tsrc\tCDS\t1\t3\t.\t+\t0\tParent=other\n",
            "seq\tsrc\tgene\t1\t12\t.\t+\t.\tID=g\nseq\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=t;Parent=g\nseq\tsrc\tCDS\t1\t3\t.\t+\t0\tParent=t\nseq\tsrc\tCDS\t7\t9\t.\t+\t0\tParent=other\n",
        ],
    )
    def test_invalid_model_gff_is_rejected(self, tmp_path, body):
        with pytest.raises(InputValidationError):
            parse_model_gff(write(tmp_path / "bad.gff", body), sequence_length=12)

    def test_draft_gff_is_minimal_deterministic_and_genomically_sorted(self):
        features = [CdsFeature(10, 12, "+", "0"), CdsFeature(1, 3, "+", "0")]
        rendered = render_draft_gff("target", features, source="locusAlignment")
        assert rendered == (
            "target\tlocusAlignment\tgene\t1\t12\t.\t+\t.\tID=target\n"
            "target\tlocusAlignment\tCDS\t1\t3\t.\t+\t0\tID=target_cds1;Parent=target\n"
            "target\tlocusAlignment\tCDS\t10\t12\t.\t+\t0\tID=target_cds2;Parent=target\n"
        )
