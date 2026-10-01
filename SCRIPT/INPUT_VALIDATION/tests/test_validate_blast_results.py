"""Tests for external BLAST result validation."""

from pathlib import Path
from typing import Optional

import pytest

from validate_blast_results import BlastResultsValidator, ValidationError


def write_blast(
    tmp_path: Path,
    content: str,
    filename: str = "blast.tsv",
) -> Path:
    """Write a temporary BLAST result table."""
    path = tmp_path / filename
    path.write_text(content, encoding="utf-8")
    return path


def validate(
    gff_path: Path,
    target_fai_path: Path,
    tblastn_path: Optional[Path] = None,
    blastn_path: Optional[Path] = None,
) -> BlastResultsValidator:
    """Run BLAST result validation."""
    validator = BlastResultsValidator(
        reference_gff_path=gff_path,
        target_fai_path=target_fai_path,
        tblastn_path=tblastn_path,
        blastn_path=blastn_path,
    )
    validator.validate()
    return validator


def test_valid_tblastn_passes(
    valid_gff: Path,
    valid_target_fai: Path,
    valid_tblastn: Path,
) -> None:
    """Accept valid external TBLASTN results."""
    validate(
        valid_gff,
        valid_target_fai,
        tblastn_path=valid_tblastn,
    )


def test_valid_blastn_passes(
    valid_gff: Path,
    valid_target_fai: Path,
    valid_blastn: Path,
) -> None:
    """Accept valid external BLASTN results."""
    validate(
        valid_gff,
        valid_target_fai,
        blastn_path=valid_blastn,
    )


def test_both_blast_files_pass(
    valid_gff: Path,
    valid_target_fai: Path,
    valid_tblastn: Path,
    valid_blastn: Path,
) -> None:
    """Accept compatible TBLASTN and BLASTN results together."""
    validate(
        valid_gff,
        valid_target_fai,
        tblastn_path=valid_tblastn,
        blastn_path=valid_blastn,
    )


def test_no_blast_file_fails(
    valid_gff: Path,
    valid_target_fai: Path,
) -> None:
    """Require at least one external BLAST result file."""
    with pytest.raises(
        ValidationError,
        match="At least one BLAST result file",
    ):
        validate(
            valid_gff,
            valid_target_fai,
        )


@pytest.mark.parametrize(
    "blast_type",
    [
        "tblastn",
        "blastn",
    ],
)
def test_unknown_query_id_fails(
    valid_gff: Path,
    valid_target_fai: Path,
    tmp_path: Path,
    blast_type: str,
) -> None:
    """Reject BLAST query IDs absent from the reference GFF."""
    blast_path = write_blast(
        tmp_path,
        "unknown_gene\tchrTarget\t100\t90\n",
    )

    kwargs = {f"{blast_type}_path": blast_path}

    with pytest.raises(
        ValidationError,
        match="absent from the reference GFF",
    ):
        validate(
            valid_gff,
            valid_target_fai,
            **kwargs,
        )


@pytest.mark.parametrize(
    "blast_type",
    [
        "tblastn",
        "blastn",
    ],
)
def test_unknown_subject_seqid_fails(
    valid_gff: Path,
    valid_target_fai: Path,
    tmp_path: Path,
    blast_type: str,
) -> None:
    """Reject BLAST subject seqids absent from the target genome."""
    blast_path = write_blast(
        tmp_path,
        "gene1\tunknown_chr\t100\t90\n",
    )

    kwargs = {f"{blast_type}_path": blast_path}

    with pytest.raises(
        ValidationError,
        match="absent from the target genome",
    ):
        validate(
            valid_gff,
            valid_target_fai,
            **kwargs,
        )


def test_multiple_valid_query_and_subject_ids_pass(
    valid_gff: Path,
    tmp_path: Path,
) -> None:
    """Accept BLAST tables containing multiple valid identifiers."""
    with valid_gff.open("a", encoding="utf-8") as handle:
        handle.write(
            "chr1\ttest\tgene\t1000\t1500\t.\t+\t.\tID=gene2\n"
        )

    target_fai = tmp_path / "target.fasta.fai"
    target_fai.write_text(
        "chrTarget1\t5000\t11\t80\t81\n"
        "chrTarget2\t6000\t5080\t80\t81\n",
        encoding="utf-8",
    )

    blast_path = write_blast(
        tmp_path,
        "gene1\tchrTarget1\t100\t90\n"
        "gene2\tchrTarget2\t100\t90\n",
    )

    validate(
        valid_gff,
        target_fai,
        tblastn_path=blast_path,
    )


def test_extra_reference_gene_ids_are_allowed(
    valid_gff: Path,
    valid_target_fai: Path,
    valid_tblastn: Path,
) -> None:
    """Allow reference genes that are absent from BLAST results."""
    with valid_gff.open("a", encoding="utf-8") as handle:
        handle.write(
            "chr1\ttest\tgene\t1000\t1500\t.\t+\t.\tID=gene2\n"
        )

    validate(
        valid_gff,
        valid_target_fai,
        tblastn_path=valid_tblastn,
    )


def test_extra_target_seqids_are_allowed(
    valid_gff: Path,
    tmp_path: Path,
    valid_tblastn: Path,
) -> None:
    """Allow target sequences that are absent from BLAST results."""
    target_fai = tmp_path / "target.fasta.fai"
    target_fai.write_text(
        "chrTarget\t5000\t11\t80\t81\n"
        "chrUnused\t6000\t5080\t80\t81\n",
        encoding="utf-8",
    )

    validate(
        valid_gff,
        target_fai,
        tblastn_path=valid_tblastn,
    )


def test_blank_lines_are_ignored(
    valid_gff: Path,
    valid_target_fai: Path,
    tmp_path: Path,
) -> None:
    """Ignore blank lines in BLAST result tables."""
    blast_path = write_blast(
        tmp_path,
        "\n"
        "gene1\tchrTarget\t100\t90\n"
        "\n",
    )

    validate(
        valid_gff,
        valid_target_fai,
        tblastn_path=blast_path,
    )


def test_malformed_blast_row_fails(
    valid_gff: Path,
    valid_target_fai: Path,
    tmp_path: Path,
) -> None:
    """Reject BLAST rows without both qseqid and sseqid."""
    blast_path = write_blast(
        tmp_path,
        "gene1\n",
    )

    with pytest.raises(
        ValidationError,
        match="at least qseqid and sseqid",
    ):
        validate(
            valid_gff,
            valid_target_fai,
            tblastn_path=blast_path,
        )


@pytest.mark.parametrize(
    "row",
    [
        "\tchrTarget\t100\t90\n",
        "gene1\t\t100\t90\n",
    ],
)
def test_empty_blast_identifiers_fail(
    valid_gff: Path,
    valid_target_fai: Path,
    tmp_path: Path,
    row: str,
) -> None:
    """Reject empty qseqid or sseqid values."""
    blast_path = write_blast(
        tmp_path,
        row,
    )

    with pytest.raises(
        ValidationError,
        match="must not be empty",
    ):
        validate(
            valid_gff,
            valid_target_fai,
            tblastn_path=blast_path,
        )


def test_empty_blast_file_fails(
    valid_gff: Path,
    valid_target_fai: Path,
    tmp_path: Path,
) -> None:
    """Reject an empty BLAST result file."""
    blast_path = write_blast(
        tmp_path,
        "",
    )

    with pytest.raises(
        ValidationError,
        match="no BLAST result rows found",
    ):
        validate(
            valid_gff,
            valid_target_fai,
            tblastn_path=blast_path,
        )


def test_duplicate_blast_rows_are_allowed(
    valid_gff: Path,
    valid_target_fai: Path,
    tmp_path: Path,
) -> None:
    """Allow repeated BLAST hits with valid identifiers."""
    blast_path = write_blast(
        tmp_path,
        "gene1\tchrTarget\t100\t90\n"
        "gene1\tchrTarget\t100\t90\n",
    )

    validator = validate(
        valid_gff,
        valid_target_fai,
        tblastn_path=blast_path,
    )

    assert validator.summaries["TBLASTN"] == (2, 1, 1)


def test_malformed_target_fai_fails(
    valid_gff: Path,
    valid_tblastn: Path,
    tmp_path: Path,
) -> None:
    """Reject malformed target FASTA index rows."""
    fai_path = tmp_path / "target.fasta.fai"
    fai_path.write_text(
        "chrTarget\n",
        encoding="utf-8",
    )

    with pytest.raises(
        ValidationError,
        match="malformed FAI row",
    ):
        validate(
            valid_gff,
            fai_path,
            tblastn_path=valid_tblastn,
        )


def test_empty_target_fai_fails(
    valid_gff: Path,
    valid_tblastn: Path,
    tmp_path: Path,
) -> None:
    """Reject a target FASTA index without sequence identifiers."""
    fai_path = tmp_path / "target.fasta.fai"
    fai_path.write_text("", encoding="utf-8")

    with pytest.raises(
        ValidationError,
        match="no sequence identifiers found",
    ):
        validate(
            valid_gff,
            fai_path,
            tblastn_path=valid_tblastn,
        )


def test_validation_report_with_both_files(
    valid_gff: Path,
    valid_target_fai: Path,
    valid_tblastn: Path,
    valid_blastn: Path,
    tmp_path: Path,
) -> None:
    """Write validation statistics for both BLAST inputs."""
    report_path = tmp_path / "blast_validation.tsv"

    validator = validate(
        valid_gff,
        valid_target_fai,
        tblastn_path=valid_tblastn,
        blastn_path=valid_blastn,
    )

    validator.write_report(report_path)

    report = report_path.read_text(encoding="utf-8")

    assert "status\tOK" in report
    assert "tblastn\tOK" in report
    assert "tblastn_rows\t1" in report
    assert "tblastn_query_ids\t1" in report
    assert "tblastn_subject_seqids\t1" in report
    assert "blastn\tOK" in report
    assert "blastn_rows\t1" in report


def test_validation_report_marks_missing_input_as_skipped(
    valid_gff: Path,
    valid_target_fai: Path,
    valid_tblastn: Path,
    tmp_path: Path,
) -> None:
    """Mark a BLAST input as skipped when it was not supplied."""
    report_path = tmp_path / "blast_validation.tsv"

    validator = validate(
        valid_gff,
        valid_target_fai,
        tblastn_path=valid_tblastn,
    )

    validator.write_report(report_path)

    report = report_path.read_text(encoding="utf-8")

    assert "tblastn\tOK" in report
    assert "blastn\tSKIPPED" in report
