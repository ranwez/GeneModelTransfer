from __future__ import annotations

from locus_alignment_transfer import BlastHsp, run_transfer


def test_reverse_hsp_prepares_reverse_complement_and_projects_negative_strand(tmp_path):
    model_fasta = tmp_path / "model.fa"
    model_fasta.write_text(">model\nAACCGGTTACGA\n")
    target_fasta = tmp_path / "target.fa"
    target_fasta.write_text(">target\nTCGTAACCGGTT\n")
    model_gff = tmp_path / "model.gff"
    model_gff.write_text(
        "model\tsrc\tgene\t1\t12\t.\t+\t.\tID=model\n"
        "model\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=model_rna;Parent=model\n"
        "model\tsrc\tCDS\t1\t12\t.\t+\t0\tParent=model_rna\n"
    )
    hsp = BlastHsp(
        "model", "target", 12, 12, 1, 12, 12, 1, 12, 100.0, 24.0, 0.0
    )

    def aligned(model_sequence, prepared_target):
        assert prepared_target == model_sequence == "AACCGGTTACGA"
        return model_sequence, prepared_target

    output = tmp_path / "draft.gff"
    result = run_transfer(
        model_fasta,
        model_gff,
        target_fasta,
        output,
        tmp_path / "diagnostics.json",
        blast_runner=lambda *_: [hsp],
        mafft_runner=aligned,
    )

    assert result.status == "success"
    assert "\tCDS\t1\t12\t.\t-\t0\t" in output.read_text()
