import io
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from CANDIDATE_LOCI.blastn_utils import (
    BlastnHit,
    ParametersBlastn,
    blastn_desired_expansion,
    blastn_new_loci_discovery,
    prepare_blastn,
    select_same_model_candidate,
)
from CANDIDATE_LOCI.bounds import Bounds
from CANDIDATE_LOCI.candidate_loci import (
    CandidateLocus,
    ExpansionPair,
    ParametersCandidateLoci,
    ParametersExpansion,
    ParametersLociScoring,
    find_candidate_loci,
    find_candidate_loci_from_file,
    reconcile_expansions,
)
from CANDIDATE_LOCI.gff_utils import gff_to_geneInfo


FIXTURE = (
    Path(__file__).parent
    / "data"
    / "OsjNip_HIPP_Nip2Nip_Chr04_19814005"
)


def _blastn_row(
    query="model1",
    chromosome="chr1",
    qlen=100,
    length=100,
    qstart=1,
    qend=100,
    sstart=100,
    send=199,
    nident=100,
    pident=100,
    width=12,
):
    fields = [
        query,
        chromosome,
        str(qlen),
        str(length),
        str(qstart),
        str(qend),
        str(sstart),
        str(send),
        str(nident),
        str(pident),
        "0",
        "0",
    ]
    fields.extend(str(index) for index in range(12, width))
    return "\t".join(fields)


def _write_gff(tmp_path, cds=((100, 399),), mrnas=("model1_rna",), cds_parent=None):
    path = tmp_path / "models.gff"
    lines = ["chrR\ttest\tgene\t100\t999\t.\t+\t.\tID=model1"]
    for mrna in mrnas:
        lines.append(
            f"chrR\ttest\tmRNA\t100\t999\t.\t+\t.\tID={mrna};Parent=model1"
        )
    parent = cds_parent if cds_parent is not None else (mrnas[0] if mrnas else "model1")
    for index, (start, end) in enumerate(cds):
        lines.append(
            f"chrR\ttest\tCDS\t{start}\t{end}\t.\t+\t0\t"
            f"ID=cds{index};Parent={parent}"
        )
    path.write_text("\n".join(lines) + "\n")
    return path


def _hit(query="model1", start=100, end=199, quality=1.0, strand=1):
    return BlastnHit(
        query_id=query,
        chr_id="chr1",
        query_length=100,
        alignment_length=100,
        nident=100,
        query_start=1,
        query_end=100,
        start=start,
        end=end,
        strand=strand,
        coverage=quality,
        identity=1.0,
        quality=quality,
    )


def _candidate(model, start, end, score=100.0):
    return SimpleNamespace(
        prot_id=model,
        score=score,
        chr_bounds=Bounds(start, end),
    )


@pytest.mark.parametrize("width", [12, 13, 16])
def test_parser_accepts_consistent_supported_widths(tmp_path, width):
    table = tmp_path / f"blastn_{width}.tsv"
    table.write_text(_blastn_row(width=width) + "\n")
    prepared = prepare_blastn(table, {"model1": object()})
    assert prepared.stats.raw_rows == 1
    assert prepared.stats.qualifying_rows == 1
    assert len(prepared.hits_by_chr["chr1"]) == 1


def test_parser_threshold_boundaries_and_just_below(tmp_path):
    table = tmp_path / "thresholds.tsv"
    table.write_text(
        "\n".join(
            [
                _blastn_row(length=85, qend=85, nident=81, pident=95),
                _blastn_row(length=84, qend=84, nident=84, pident=100, sstart=300, send=383),
                _blastn_row(length=100, pident=94.999, sstart=500, send=599),
            ]
        )
        + "\n"
    )
    prepared = prepare_blastn(table, {"model1": object()})
    assert prepared.stats.qualifying_rows == 1
    assert prepared.hits_by_chr["chr1"][0].coverage == pytest.approx(0.85)
    assert prepared.hits_by_chr["chr1"][0].identity == pytest.approx(0.95)


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"qlen": "bad"}, "qlen"),
        ({"length": 0}, "positive"),
        ({"qstart": 0}, "positive"),
        ({"qend": 101}, "outside"),
        ({"pident": "nan"}, "finite"),
        ({"pident": 101}, "pident"),
    ],
)
def test_parser_rejects_malformed_required_fields(tmp_path, kwargs, message):
    table = tmp_path / "bad.tsv"
    table.write_text(_blastn_row(**kwargs) + "\n")
    with pytest.raises(ValueError, match=message):
        prepare_blastn(table, {"model1": object()})


def test_parser_rejects_short_inconsistent_and_unknown_rows(tmp_path):
    short = tmp_path / "short.tsv"
    short.write_text("\t".join(["x"] * 11) + "\n")
    with pytest.raises(ValueError, match=r"short.tsv:1.*at least 12"):
        prepare_blastn(short, {"model1": object()})

    inconsistent = tmp_path / "inconsistent.tsv"
    inconsistent.write_text(_blastn_row(width=12) + "\n" + _blastn_row(width=13) + "\n")
    with pytest.raises(ValueError, match=r"inconsistent.tsv:2.*expected 12"):
        prepare_blastn(inconsistent, {"model1": object()})

    unknown = tmp_path / "unknown.tsv"
    unknown.write_text(_blastn_row(query="missing") + "\n")
    with pytest.raises(ValueError, match=r"unknown.tsv:1.*missing.*reference GFF"):
        prepare_blastn(unknown, {"model1": object()})


@pytest.mark.parametrize("value", [-0.01, 1.01, float("nan"), float("inf")])
def test_invalid_thresholds_fail_before_rows(value):
    with pytest.raises(ValueError, match="fraction"):
        ParametersBlastn(min_coverage=value)
    with pytest.raises(ValueError, match="fraction"):
        ParametersBlastn(min_identity=value)


def test_parser_normalizes_reverse_coordinates_deduplicates_and_orders(tmp_path):
    table = tmp_path / "ranked.tsv"
    reverse = _blastn_row(qstart=100, qend=1, sstart=300, send=201)
    lower = _blastn_row(sstart=500, send=599, pident=96, nident=96)
    table.write_text(lower + "\n" + reverse + "\n" + reverse + "\n")
    prepared = prepare_blastn(table, {"model1": object()})
    hits = prepared.hits_by_chr["chr1"]
    assert prepared.stats.qualifying_rows == 3
    assert prepared.stats.retained_hits == 2
    assert [(hit.start, hit.end) for hit in hits] == [(201, 300), (500, 599)]
    assert hits[0].query_start == 1 and hits[0].query_end == 100
    assert hits[0].strand == -1


def test_empty_and_no_qualifying_input(tmp_path):
    empty = tmp_path / "empty.tsv"
    empty.write_text("\n")
    assert prepare_blastn(empty, {"model1": object()}).hits_by_chr == {}
    low = tmp_path / "low.tsv"
    low.write_text(_blastn_row(length=10, qend=10, pident=10) + "\n")
    prepared = prepare_blastn(low, {"model1": object()})
    assert prepared.stats.raw_rows == 1
    assert prepared.stats.qualifying_rows == 0
    assert prepared.hits_by_chr == {}


def test_same_model_selection_and_expansion_are_deterministic():
    hit = _hit(start=150, end=250)
    lower_score = _candidate("model1", 100, 200, score=50)
    higher_score = _candidate("model1", 200, 300, score=75)
    other = _candidate("other", 140, 260, score=1000)
    assert select_same_model_candidate(hit, [lower_score, higher_score, other]) is higher_score
    candidate, left, right = blastn_desired_expansion(
        _hit(start=180, end=340), [higher_score]
    )
    assert candidate is higher_score
    assert (left, right) == (20, 40)


def test_discovery_overlap_contract_covers_models_priority_and_adjacency():
    protein = _candidate("other", 100, 200)
    accepted = _candidate("model1", 300, 400)
    assert not blastn_new_loci_discovery(_hit(start=200, end=250), [protein], [])
    assert blastn_new_loci_discovery(_hit(start=201, end=299), [protein], [accepted])
    assert not blastn_new_loci_discovery(_hit(start=250, end=300), [], [accepted])
    assert not blastn_new_loci_discovery(_hit(start=320, end=350), [], [accepted])


def test_reference_protein_lengths_single_and_multi_cds(tmp_path):
    single = _write_gff(tmp_path, cds=((100, 399),))
    info, _ = gff_to_geneInfo(single, 0.5, require_protein_length=True)
    assert info["model1"].protein_length == 100

    multi_dir = tmp_path / "multi"
    multi_dir.mkdir()
    multi = _write_gff(multi_dir, cds=((100, 248), (300, 450)))
    info, _ = gff_to_geneInfo(multi, 0.5, require_protein_length=True)
    assert info["model1"].protein_length == 100


@pytest.mark.parametrize(
    "mrnas, cds, parent, message",
    [
        ((), ((100, 399),), "model1", "exactly one mRNA"),
        (("rna1", "rna2"), ((100, 399),), "rna1", "exactly one mRNA"),
        (("rna1",), (), None, "no CDS"),
        (("rna1",), ((100, 399),), "model1", "inconsistent parentage"),
        (("rna1",), ((100, 398),), "rna1", "divisible by three"),
    ],
)
def test_reference_model_validation_errors(tmp_path, mrnas, cds, parent, message):
    gff = _write_gff(tmp_path, cds=cds, mrnas=mrnas, cds_parent=parent)
    with pytest.raises(ValueError, match=message):
        gff_to_geneInfo(gff, 0.5, require_protein_length=True)


def test_blastn_only_scoring_source_serialization_and_query_target():
    loci = find_candidate_loci(
        FIXTURE / "ref.gff",
        FIXTURE / "blastp_default.tsv",
        blastn_file=FIXTURE / "blastn.tsv",
    )
    candidate = next(
        locus for locus in loci["Chr4"] if locus.prot_id == "OsjNip_Chr04_19836437"
    )
    assert candidate.origin == "blastn"
    assert candidate.prot_len == 267
    assert candidate.pc_similarity == 1.0
    assert candidate.score == 267.0
    gff = candidate.as_gff()
    assert gff.count("\n") == 1
    assert "\tcandidateLoci_blastn\tgene\t" in gff
    assert "\tCDS\t" not in gff
    assert candidate.prot_id in candidate.as_query_target()
    protein = next(locus for locus in loci["Chr4"] if locus.origin == "blastp")
    assert "\tcandidateLoci_blastp\tgene\t" in protein.as_gff()


def test_HIPP_default_blastp_recovers_one_true_blastn_only_locus(capsys):
    """Lock discovery, diagnostics, origin, and feature shape for the V2 fixture."""
    missing_model = "OsjNip_Chr04_19836437"
    expansion = ParametersExpansion(
        nb_aa_for_missing_part=10,
        nb_nt_default=0,
        nb_nt_when_missing_part=0,
    )
    params = ParametersCandidateLoci(expansion=expansion)

    protein_only = find_candidate_loci(
        FIXTURE / "ref.gff",
        FIXTURE / "blastp_default.tsv",
        params=params,
    )
    assert all(locus.prot_id != missing_model for locus in protein_only["Chr4"])
    capsys.readouterr()

    combined = find_candidate_loci(
        FIXTURE / "ref.gff",
        FIXTURE / "blastp_default.tsv",
        params=params,
        blastn_file=FIXTURE / "blastn.tsv",
    )
    diagnostic = capsys.readouterr().out
    assert "expanded=1, new=1" in diagnostic

    recovered = [
        locus for locus in combined["Chr4"] if locus.prot_id == missing_model
    ]
    assert len(recovered) == 1
    candidate = recovered[0]
    assert candidate.origin == "blastn"
    assert candidate.chr_bounds == Bounds(19835501, 19836437)
    assert sum(locus.origin == "blastn" for locus in combined["Chr4"]) == 1
    assert sum(locus.origin == "blastp" for locus in combined["Chr4"]) == 6

    features = [line.split("\t") for line in candidate.as_gff().splitlines()]
    assert [fields[1] for fields in features] == [
        "candidateLoci_blastn",
        "candidateLoci_blastn",
    ]
    assert [fields[2] for fields in features] == ["gene", "mRNA"]


def test_evidence_expansion_runs_with_expansion_none_and_unions_hits(tmp_path):
    gff = _write_gff(tmp_path)
    blastp = tmp_path / "blastp.tsv"
    blastp.write_text("model1\tchr1\t100\t34\t1\t34\t100\t200\t34\t100\t0\t0\t100\n")
    blastn = tmp_path / "blastn.tsv"
    blastn.write_text(
        _blastn_row(qlen=101, length=101, qend=101, sstart=80, send=180)
        + "\n"
        + _blastn_row(qlen=121, length=121, qend=121, sstart=120, send=240)
        + "\n"
    )
    params = ParametersCandidateLoci(
        expansion=None,
        loci_scoring=ParametersLociScoring(min_score=0, min_similarity=0),
    )
    loci = find_candidate_loci(gff, blastp, params=params, blastn_file=blastn)
    assert loci["chr1"][0].chr_bounds == Bounds(80, 240)


def test_new_chromosome_filtering_and_final_score_filters(tmp_path):
    gff = _write_gff(tmp_path)
    empty_protein = io.StringIO("")
    blastn = tmp_path / "blastn.tsv"
    blastn.write_text(_blastn_row(chromosome="chrNew") + "\n")
    params = ParametersCandidateLoci(
        expansion=None,
        loci_scoring=ParametersLociScoring(min_score=101, min_similarity=0),
    )
    assert find_candidate_loci_from_file(
        gff, empty_protein, params=params, blastn_file=blastn
    ) == {}

    params.loci_scoring.min_score = 0
    params.loci_scoring.min_similarity = 1.01
    assert find_candidate_loci_from_file(
        gff, io.StringIO(""), params=params, blastn_file=blastn
    ) == {}
    params.loci_scoring.min_similarity = 0
    loci = find_candidate_loci_from_file(
        gff, io.StringIO(""), params=params, blastn_file=blastn
    )
    assert list(loci) == ["chrNew"]
    assert find_candidate_loci_from_file(
        gff,
        io.StringIO(""),
        params=params,
        chr="other",
        blastn_file=blastn,
    ) == {}


def test_higher_ranked_blastn_hit_blocks_lower_ranked_overlap(tmp_path):
    gff = _write_gff(tmp_path)
    blastn = tmp_path / "blastn.tsv"
    blastn.write_text(
        _blastn_row(sstart=100, send=199, pident=96, nident=96)
        + "\n"
        + _blastn_row(sstart=150, send=249, pident=100, nident=100)
        + "\n"
    )
    params = ParametersCandidateLoci(
        expansion=None,
        loci_scoring=ParametersLociScoring(min_score=0, min_similarity=0),
    )
    loci = find_candidate_loci_from_file(
        gff, io.StringIO(""), params=params, blastn_file=blastn
    )
    assert [locus.chr_bounds for locus in loci["chr1"]] == [Bounds(150, 249)]


def test_neighbor_clipping_and_genomic_start_remain_non_overlapping():
    left = SimpleNamespace(
        chr_id="chr1",
        chr_bounds=Bounds(2, 100),
        expansion=ExpansionPair(20, 100),
        shrink_info=(False, False),
    )
    right = SimpleNamespace(
        chr_id="chr1",
        chr_bounds=Bounds(151, 200),
        expansion=ExpansionPair(100, 0),
        shrink_info=(False, False),
    )
    reconcile_expansions([right, left])
    assert left.chr_bounds.start == 1
    assert left.chr_bounds.end < right.chr_bounds.start
    assert right.chr_bounds.end == 200


def test_blastn_parameter_without_file_and_missing_file_errors(tmp_path):
    gff = _write_gff(tmp_path)
    with pytest.raises(ValueError, match="requires blastn_file"):
        find_candidate_loci_from_file(
            gff, io.StringIO(""), blastn_params=ParametersBlastn()
        )
    with pytest.raises(FileNotFoundError, match="missing or unreadable"):
        find_candidate_loci_from_file(
            gff, io.StringIO(""), blastn_file=tmp_path / "missing.tsv"
        )


def test_no_blastn_preserves_legacy_source_and_structure():
    loci = find_candidate_loci(FIXTURE / "ref.gff", FIXTURE / "blastp_default.tsv")
    gff = loci["Chr4"][0].as_gff()
    assert "\tcandidateLoci\tgene\t" in gff
    assert "candidateLoci_blast" not in gff
    assert "\tCDS\t" in gff


@pytest.mark.parametrize("protein_table", ["blastp_softmax.tsv", "blastp_default.tsv"])
def test_cli_blastn_fixture_outputs_are_consistent(tmp_path, protein_table):
    output_gff = tmp_path / f"{protein_table}.gff"
    output_list = tmp_path / f"{protein_table}.list"
    command = [
        sys.executable,
        "SCRIPT/candidate_loci_VR.py",
        "-g",
        str(FIXTURE / "ref.gff"),
        "-t",
        str(FIXTURE / protein_table),
        "--blastn_table",
        str(FIXTURE / "blastn.tsv"),
        "-o",
        str(output_gff),
        "-l",
        str(output_list),
        "--ignore_gff_for_expansion",
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    gene_lines = [line for line in output_gff.read_text().splitlines() if "\tgene\t" in line]
    associations = output_list.read_text().splitlines()
    assert len(gene_lines) == len(associations) == 7
    if protein_table == "blastp_default.tsv":
        assert any("candidateLoci_blastn" in line for line in gene_lines)
    else:
        assert all("candidateLoci_blastp" in line for line in gene_lines)


def test_cli_missing_blastn_does_not_replace_existing_outputs(tmp_path):
    output_gff = tmp_path / "existing.gff"
    output_list = tmp_path / "existing.list"
    output_gff.write_text("old gff\n")
    output_list.write_text("old list\n")
    command = [
        sys.executable,
        "SCRIPT/candidate_loci_VR.py",
        "-g",
        str(FIXTURE / "ref.gff"),
        "-t",
        str(FIXTURE / "blastp_default.tsv"),
        "--blastn_table",
        str(tmp_path / "missing.tsv"),
        "-o",
        str(output_gff),
        "-l",
        str(output_list),
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode != 0
    assert output_gff.read_text() == "old gff\n"
    assert output_list.read_text() == "old list\n"



def test_cli_staging_failure_does_not_replace_first_output(tmp_path):
    output_gff = tmp_path / "existing.gff"
    output_gff.write_text("old gff\n")
    output_list = tmp_path / "missing_directory" / "result.list"
    command = [
        sys.executable,
        "SCRIPT/candidate_loci_VR.py",
        "-g",
        str(FIXTURE / "ref.gff"),
        "-t",
        str(FIXTURE / "blastp_default.tsv"),
        "-o",
        str(output_gff),
        "-l",
        str(output_list),
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode != 0
    assert output_gff.read_text() == "old gff\n"
    assert not list(tmp_path.glob(".candidate_loci_*"))
