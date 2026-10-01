"""Tests for reference GFF validation."""

from pathlib import Path

import pytest

from validate_ref_gff import ReferenceGffValidator, ValidationError


GENE = "chr1\ttest\tgene\t100\t900\t.\t+\t.\tID=gene1"
MRNA = "chr1\ttest\tmRNA\t150\t850\t.\t+\t.\tID=mrna1;Parent=gene1"
CDS1 = "chr1\ttest\tCDS\t200\t300\t.\t+\t0\tParent=mrna1"
CDS2 = "chr1\ttest\tCDS\t400\t500\t.\t+\t0\tID=cds2;Parent=mrna1"


def write_gff(tmp_path: Path, lines: list[str]) -> Path:
    """Write GFF lines to a temporary test file."""
    path = tmp_path / "test.gff"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def validate(
    gff_path: Path,
    ref_fai_path: Path | None = None,
    ref_locus_info_path: Path | None = None,
) -> None:
    """Run reference GFF validation."""
    validator = ReferenceGffValidator(
        gff_path=gff_path,
        ref_fai_path=ref_fai_path,
        ref_locus_info_path=ref_locus_info_path,
    )
    validator.validate()


def test_valid_gff_passes(
    valid_gff: Path,
    valid_fai: Path,
    valid_info_locus: Path,
) -> None:
    """Accept a valid GFF and matching reference inputs."""
    validate(valid_gff, valid_fai, valid_info_locus)


def test_optional_reference_inputs_can_be_omitted(
    valid_gff: Path,
) -> None:
    """Accept a valid GFF without optional reference inputs."""
    validate(valid_gff)


def test_unsorted_gff_passes(tmp_path: Path) -> None:
    """Accept valid features regardless of their input order."""
    gff_path = write_gff(
        tmp_path,
        [CDS2, MRNA, CDS1, GENE],
    )

    validate(gff_path)


def test_cds_without_id_passes(tmp_path: Path) -> None:
    """Accept CDS features without an ID attribute."""
    cds2_without_id = (
        "chr1\ttest\tCDS\t400\t500\t.\t+\t0\tParent=mrna1"
    )

    gff_path = write_gff(
        tmp_path,
        [GENE, MRNA, CDS1, cds2_without_id],
    )

    validate(gff_path)


def test_extra_info_locus_ids_are_allowed(
    valid_gff: Path,
    tmp_path: Path,
) -> None:
    """Allow locus information to contain additional genes."""
    info_path = tmp_path / "info_locus.tsv"
    info_path.write_text(
        "gene1\tNLR\tCanonical\n"
        "gene2\tNLR\tCanonical\n",
        encoding="utf-8",
    )

    validate(valid_gff, ref_locus_info_path=info_path)


def test_extra_fasta_seqids_are_allowed(
    valid_gff: Path,
    tmp_path: Path,
) -> None:
    """Allow the reference FASTA to contain additional sequences."""
    fai_path = tmp_path / "reference.fasta.fai"
    fai_path.write_text(
        "chr1\t2000\t6\t80\t81\n"
        "chr2\t3000\t2032\t80\t81\n",
        encoding="utf-8",
    )

    validate(valid_gff, ref_fai_path=fai_path)


@pytest.mark.parametrize(
    "line",
    [
        "chr1\ttest\tgene\t100\t900\t.\t+\t.",
        "chr1\ttest\tgene\tfoo\t900\t.\t+\t.\tID=gene1",
        "chr1\ttest\tgene\t0\t900\t.\t+\t.\tID=gene1",
        "chr1\ttest\tgene\t900\t100\t.\t+\t.\tID=gene1",
    ],
)
def test_invalid_gff_rows_fail(
    tmp_path: Path,
    line: str,
) -> None:
    """Reject malformed GFF rows and invalid coordinates."""
    gff_path = write_gff(tmp_path, [line])

    with pytest.raises(ValidationError):
        validate(gff_path)


@pytest.mark.parametrize(
    "gene",
    [
        "chr1\ttest\tgene\t100\t900\t.\t+\t.\tName=gene1",
        "chr1\ttest\tgene\t100\t900\t.\t.\t.\tID=gene1",
    ],
)
def test_invalid_gene_features_fail(
    tmp_path: Path,
    gene: str,
) -> None:
    """Reject genes without a valid ID or strand."""
    gff_path = write_gff(
        tmp_path,
        [gene, MRNA, CDS1, CDS2],
    )

    with pytest.raises(ValidationError):
        validate(gff_path)


def test_duplicate_gene_ids_fail(tmp_path: Path) -> None:
    """Reject duplicate gene IDs."""
    duplicate_gene = (
        "chr1\ttest\tgene\t1000\t1500\t.\t+\t.\tID=gene1"
    )

    gff_path = write_gff(
        tmp_path,
        [GENE, duplicate_gene, MRNA, CDS1, CDS2],
    )

    with pytest.raises(ValidationError):
        validate(gff_path)


def test_gene_without_mrna_fails(tmp_path: Path) -> None:
    """Reject genes without an mRNA."""
    gff_path = write_gff(tmp_path, [GENE])

    with pytest.raises(ValidationError):
        validate(gff_path)


def test_gene_with_multiple_mrnas_fails(tmp_path: Path) -> None:
    """Reject genes containing more than one mRNA."""
    second_mrna = (
        "chr1\ttest\tmRNA\t150\t850\t.\t+\t.\t"
        "ID=mrna2;Parent=gene1"
    )

    gff_path = write_gff(
        tmp_path,
        [GENE, MRNA, second_mrna, CDS1, CDS2],
    )

    with pytest.raises(ValidationError):
        validate(gff_path)


@pytest.mark.parametrize(
    "mrna",
    [
        "chr1\ttest\tmRNA\t150\t850\t.\t+\t.\tParent=gene1",
        "chr1\ttest\tmRNA\t150\t850\t.\t+\t.\tID=mrna1",
        (
            "chr1\ttest\tmRNA\t150\t850\t.\t+\t.\t"
            "ID=mrna1;Parent=gene1,gene2"
        ),
        (
            "chr1\ttest\tmRNA\t150\t850\t.\t+\t.\t"
            "ID=mrna1;Parent=missing_gene"
        ),
        (
            "chr2\ttest\tmRNA\t150\t850\t.\t+\t.\t"
            "ID=mrna1;Parent=gene1"
        ),
        (
            "chr1\ttest\tmRNA\t150\t850\t.\t-\t.\t"
            "ID=mrna1;Parent=gene1"
        ),
        (
            "chr1\ttest\tmRNA\t50\t850\t.\t+\t.\t"
            "ID=mrna1;Parent=gene1"
        ),
    ],
)
def test_invalid_mrna_features_fail(
    tmp_path: Path,
    mrna: str,
) -> None:
    """Reject invalid mRNA identifiers, parentage and coordinates."""
    gff_path = write_gff(
        tmp_path,
        [GENE, mrna, CDS1, CDS2],
    )

    with pytest.raises(ValidationError):
        validate(gff_path)


def test_duplicate_mrna_ids_fail(tmp_path: Path) -> None:
    """Reject duplicate mRNA IDs."""
    second_gene = (
        "chr1\ttest\tgene\t1000\t1800\t.\t+\t.\tID=gene2"
    )
    duplicate_mrna = (
        "chr1\ttest\tmRNA\t1100\t1700\t.\t+\t.\t"
        "ID=mrna1;Parent=gene2"
    )

    gff_path = write_gff(
        tmp_path,
        [GENE, second_gene, MRNA, duplicate_mrna, CDS1, CDS2],
    )

    with pytest.raises(ValidationError):
        validate(gff_path)


def test_mrna_without_cds_fails(tmp_path: Path) -> None:
    """Reject mRNAs without CDS features."""
    gff_path = write_gff(
        tmp_path,
        [GENE, MRNA],
    )

    with pytest.raises(ValidationError):
        validate(gff_path)


@pytest.mark.parametrize(
    "cds",
    [
        "chr1\ttest\tCDS\t200\t300\t.\t+\t0\tID=cds1",
        (
            "chr1\ttest\tCDS\t200\t300\t.\t+\t0\t"
            "Parent=mrna1,mrna2"
        ),
        (
            "chr1\ttest\tCDS\t200\t300\t.\t+\t0\t"
            "Parent=missing_mrna"
        ),
        "chr2\ttest\tCDS\t200\t300\t.\t+\t0\tParent=mrna1",
        "chr1\ttest\tCDS\t200\t300\t.\t-\t0\tParent=mrna1",
        "chr1\ttest\tCDS\t100\t300\t.\t+\t0\tParent=mrna1",
    ],
)
def test_invalid_cds_features_fail(
    tmp_path: Path,
    cds: str,
) -> None:
    """Reject invalid CDS parentage, sequence, strand and coordinates."""
    gff_path = write_gff(
        tmp_path,
        [GENE, MRNA, cds, CDS2],
    )

    with pytest.raises(ValidationError):
        validate(gff_path)


def test_overlapping_cds_fail(tmp_path: Path) -> None:
    """Reject overlapping CDS features from the same mRNA."""
    overlapping_cds = (
        "chr1\ttest\tCDS\t250\t400\t.\t+\t0\tParent=mrna1"
    )

    gff_path = write_gff(
        tmp_path,
        [GENE, MRNA, CDS1, overlapping_cds],
    )

    with pytest.raises(ValidationError):
        validate(gff_path)


def test_missing_gene_from_info_locus_fails(
    valid_gff: Path,
    tmp_path: Path,
) -> None:
    """Reject reference genes absent from locus information."""
    info_path = tmp_path / "info_locus.tsv"
    info_path.write_text(
        "other_gene\tNLR\tCanonical\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError):
        validate(
            valid_gff,
            ref_locus_info_path=info_path,
        )


def test_missing_seqid_from_reference_fai_fails(
    valid_gff: Path,
    tmp_path: Path,
) -> None:
    """Reject GFF seqids absent from the reference FASTA."""
    fai_path = tmp_path / "reference.fasta.fai"
    fai_path.write_text(
        "chr2\t2000\t6\t80\t81\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError):
        validate(
            valid_gff,
            ref_fai_path=fai_path,
        )


def test_validation_report_is_written(
    valid_gff: Path,
    valid_fai: Path,
    valid_info_locus: Path,
    tmp_path: Path,
) -> None:
    """Write a validation report after successful validation."""
    report_path = tmp_path / "validation.tsv"

    validator = ReferenceGffValidator(
        gff_path=valid_gff,
        ref_fai_path=valid_fai,
        ref_locus_info_path=valid_info_locus,
    )

    validator.validate()
    validator.write_report(report_path)

    report = report_path.read_text(encoding="utf-8")

    assert "status\tOK" in report
    assert "genes\t1" in report
    assert "mrnas\t1" in report
    assert "cds\t2" in report
    assert "ref_genome_seqids\tOK" in report
    assert "ref_locus_info_gene_ids\tOK" in report
