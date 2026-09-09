from __future__ import annotations

import json

import pytest

from locus_alignment_transfer import AlignmentError, BlastHsp, run_transfer


def _write(path, content):
    path.write_text(content)
    return path


def _valid_inputs(tmp_path):
    model_fasta = _write(tmp_path / "model.fa", ">model\nACGTACGTACGT\n")
    target_fasta = _write(tmp_path / "target.fa", ">target\nACGTACGTACGT\n")
    model_gff = _write(
        tmp_path / "model.gff",
        "model\tsrc\tgene\t1\t12\t.\t+\t.\tID=model\n"
        "model\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=model_rna;Parent=model\n"
        "model\tsrc\tCDS\t1\t12\t.\t+\t0\tID=model_cds;Parent=model_rna\n",
    )
    return model_fasta, model_gff, target_fasta


def _passing_hsp():
    return BlastHsp(
        query_id="model",
        subject_id="target",
        query_length=12,
        alignment_length=12,
        query_start=1,
        query_end=12,
        subject_start=1,
        subject_end=12,
        identical=12,
        percent_identity=100.0,
        bit_score=24.0,
        evalue=0.0,
    )


def test_ineligible_pair_creates_empty_output_and_never_calls_mafft(tmp_path):
    model_fasta, model_gff, target_fasta = _valid_inputs(tmp_path)
    output = tmp_path / "draft.gff"
    diagnostics = tmp_path / "status.json"
    mafft_calls = []

    result = run_transfer(
        model_fasta=model_fasta,
        model_gff=model_gff,
        target_fasta=target_fasta,
        output_gff=output,
        diagnostics_json=diagnostics,
        blast_runner=lambda *_: [],
        mafft_runner=lambda *_: mafft_calls.append(True),
    )

    assert result.status == "ineligible"
    assert result.exit_code == 0
    assert output.exists() and output.stat().st_size == 0
    assert mafft_calls == []
    assert json.loads(diagnostics.read_text())["reason"] == "no_significant_hit"


def test_success_writes_complete_draft_and_structured_diagnostics(tmp_path):
    model_fasta, model_gff, target_fasta = _valid_inputs(tmp_path)
    output = tmp_path / "draft.gff"
    diagnostics = tmp_path / "status.json"

    result = run_transfer(
        model_fasta=model_fasta,
        model_gff=model_gff,
        target_fasta=target_fasta,
        output_gff=output,
        diagnostics_json=diagnostics,
        blast_runner=lambda *_: [_passing_hsp()],
        mafft_runner=lambda model, target: (model, target),
    )

    assert result.status == "success"
    assert result.exit_code == 0
    assert output.read_text() == (
        "target\tlocusAlignment\tgene\t1\t12\t.\t+\t.\tID=target\n"
        "target\tlocusAlignment\tCDS\t1\t12\t.\t+\t0\tID=target_cds1;Parent=target\n"
    )
    status = json.loads(diagnostics.read_text())
    assert status["status"] == "success"
    assert status["coverage"] == 1.0
    assert status["identity"] == 1.0


@pytest.mark.parametrize("failing_boundary", ["blast", "mafft"])
def test_tool_failure_is_error_and_never_leaves_partial_draft(tmp_path, failing_boundary):
    model_fasta, model_gff, target_fasta = _valid_inputs(tmp_path)
    output = _write(tmp_path / "draft.gff", "stale partial data\n")
    diagnostics = tmp_path / "status.json"

    def fail(*_):
        raise AlignmentError(f"{failing_boundary} failed")

    with pytest.raises(AlignmentError):
        run_transfer(
            model_fasta=model_fasta,
            model_gff=model_gff,
            target_fasta=target_fasta,
            output_gff=output,
            diagnostics_json=diagnostics,
            blast_runner=fail if failing_boundary == "blast" else lambda *_: [_passing_hsp()],
            mafft_runner=fail if failing_boundary == "mafft" else lambda model, target: (model, target),
        )

    assert not output.exists() or output.stat().st_size == 0
    assert json.loads(diagnostics.read_text())["status"] == "error"
