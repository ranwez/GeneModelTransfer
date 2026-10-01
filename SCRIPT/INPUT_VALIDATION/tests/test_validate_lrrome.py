"""Tests for external LRRome validation."""

from pathlib import Path

import pytest

from validate_lrrome import LrromeValidator, ValidationError


def validate(
    lrrome_path: Path,
    gff_path: Path,
) -> None:
    """Run LRRome validation."""
    validator = LrromeValidator(
        lrrome_path=lrrome_path,
        reference_gff_path=gff_path,
    )
    validator.validate()


def test_valid_lrrome_passes(
    valid_lrrome: Path,
    valid_gff: Path,
) -> None:
    """Accept a complete and consistent LRRome."""
    validate(valid_lrrome, valid_gff)


@pytest.mark.parametrize(
    "name",
    [
        "REF_proteins.fasta",
        "REF_loci.fasta",
    ],
)
def test_missing_required_file_fails(
    valid_lrrome: Path,
    valid_gff: Path,
    name: str,
) -> None:
    """Reject an LRRome missing a required file."""
    (valid_lrrome / name).unlink()

    with pytest.raises(ValidationError):
        validate(valid_lrrome, valid_gff)


@pytest.mark.parametrize(
    "name",
    [
        "REF_PEP",
        "REF_EXONS",
        "REF_cDNA",
        "REF_LOCI",
        "REF_LOCI_GFF",
    ],
)
def test_missing_required_directory_fails(
    valid_lrrome: Path,
    valid_gff: Path,
    name: str,
) -> None:
    """Reject an LRRome missing a required directory."""
    directory = valid_lrrome / name

    for path in directory.iterdir():
        path.unlink()

    directory.rmdir()

    with pytest.raises(ValidationError):
        validate(valid_lrrome, valid_gff)


@pytest.mark.parametrize(
    "name",
    [
        "REF_proteins.fasta",
        "REF_loci.fasta",
    ],
)
def test_empty_required_file_fails(
    valid_lrrome: Path,
    valid_gff: Path,
    name: str,
) -> None:
    """Reject an empty required LRRome file."""
    (valid_lrrome / name).write_text("", encoding="utf-8")

    with pytest.raises(ValidationError):
        validate(valid_lrrome, valid_gff)


def test_ref_loci_ids_must_match_protein_ids(
    valid_lrrome: Path,
    valid_gff: Path,
) -> None:
    """Require REF_loci identifiers to match reference protein models."""
    (valid_lrrome / "REF_loci.fasta").write_text(
        ">other_gene\nATGCGTACGTAG\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError):
        validate(valid_lrrome, valid_gff)


@pytest.mark.parametrize(
    "dirname",
    [
        "REF_PEP",
        "REF_cDNA",
        "REF_LOCI",
    ],
)
def test_missing_model_file_fails(
    valid_lrrome: Path,
    valid_gff: Path,
    dirname: str,
) -> None:
    """Reject missing per-model FASTA files."""
    (valid_lrrome / dirname / "gene1").unlink()

    with pytest.raises(ValidationError):
        validate(valid_lrrome, valid_gff)


@pytest.mark.parametrize(
    "dirname",
    [
        "REF_PEP",
        "REF_cDNA",
        "REF_LOCI",
    ],
)
def test_extra_model_file_fails(
    valid_lrrome: Path,
    valid_gff: Path,
    dirname: str,
) -> None:
    """Reject unexpected per-model FASTA files."""
    (valid_lrrome / dirname / "gene2").write_text(
        ">gene2\nATGC\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError):
        validate(valid_lrrome, valid_gff)


@pytest.mark.parametrize(
    "dirname",
    [
        "REF_PEP",
        "REF_cDNA",
        "REF_LOCI",
    ],
)
def test_per_model_fasta_id_must_match_filename(
    valid_lrrome: Path,
    valid_gff: Path,
    dirname: str,
) -> None:
    """Require a per-model FASTA ID to match its filename."""
    (valid_lrrome / dirname / "gene1").write_text(
        ">wrong_id\nATGC\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError):
        validate(valid_lrrome, valid_gff)


def test_missing_model_gff_fails(
    valid_lrrome: Path,
    valid_gff: Path,
) -> None:
    """Reject an LRRome missing a model GFF."""
    (valid_lrrome / "REF_LOCI_GFF" / "gene1.gff").unlink()

    with pytest.raises(ValidationError):
        validate(valid_lrrome, valid_gff)


def test_extra_model_gff_fails(
    valid_lrrome: Path,
    valid_gff: Path,
) -> None:
    """Reject unexpected model GFF files."""
    (valid_lrrome / "REF_LOCI_GFF" / "gene2.gff").write_text(
        "gene2\ttest\tgene\t1\t10\t.\t+\t.\tID=gene2\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError):
        validate(valid_lrrome, valid_gff)


def test_model_without_exon_fails(
    valid_lrrome: Path,
    valid_gff: Path,
) -> None:
    """Require at least one REF_EXONS file per model."""
    (valid_lrrome / "REF_EXONS" / "gene1_CDS_1").unlink()

    with pytest.raises(ValidationError):
        validate(valid_lrrome, valid_gff)


def test_orphan_exon_file_fails(
    valid_lrrome: Path,
    valid_gff: Path,
) -> None:
    """Reject exon files that cannot be assigned to a model."""
    (valid_lrrome / "REF_EXONS" / "other_gene_CDS_1").write_text(
        ">other_gene_CDS_1\nATGC\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError):
        validate(valid_lrrome, valid_gff)


@pytest.mark.parametrize(
    "exon_name",
    [
        "gene1_CDS_1",
        "gene1:CDS_1",
        "gene1-CDS_1",
    ],
)
def test_supported_exon_separators_pass(
    valid_lrrome: Path,
    valid_gff: Path,
    exon_name: str,
) -> None:
    """Accept exon filenames using separators supported by the workflow."""
    exon_dir = valid_lrrome / "REF_EXONS"

    for path in exon_dir.iterdir():
        path.unlink()

    (exon_dir / exon_name).write_text(
        f">{exon_name}\nATGC\n",
        encoding="utf-8",
    )

    validate(valid_lrrome, valid_gff)


def test_exon_fasta_id_must_match_filename(
    valid_lrrome: Path,
    valid_gff: Path,
) -> None:
    """Require REF_EXONS FASTA IDs to match their filenames."""
    exon_path = valid_lrrome / "REF_EXONS" / "gene1_CDS_1"
    exon_path.write_text(
        ">wrong_id\nATGC\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError):
        validate(valid_lrrome, valid_gff)


def test_reference_gff_missing_lrrome_model_fails(
    valid_lrrome: Path,
    tmp_path: Path,
) -> None:
    """Require every LRRome model to exist in the reference GFF."""
    gff_path = tmp_path / "reference.gff"
    gff_path.write_text(
        "chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=gene2\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError):
        validate(valid_lrrome, gff_path)


def test_reference_gff_extra_model_fails(
    valid_lrrome: Path,
    valid_gff: Path,
) -> None:
    """Reject reference GFF genes absent from the LRRome."""
    with valid_gff.open("a", encoding="utf-8") as handle:
        handle.write(
            "chr1\ttest\tgene\t1000\t1200\t.\t+\t.\tID=gene2\n"
        )

    with pytest.raises(ValidationError):
        validate(valid_lrrome, valid_gff)


def test_hidden_files_are_ignored(
    valid_lrrome: Path,
    valid_gff: Path,
) -> None:
    """Ignore hidden metadata files in LRRome directories."""
    for dirname in [
        "REF_PEP",
        "REF_EXONS",
        "REF_cDNA",
        "REF_LOCI",
        "REF_LOCI_GFF",
    ]:
        (valid_lrrome / dirname / ".snakemake_timestamp").touch()

    validate(valid_lrrome, valid_gff)


def test_duplicate_reference_protein_ids_fail(
    valid_lrrome: Path,
    valid_gff: Path,
) -> None:
    """Reject duplicate IDs in REF_proteins.fasta."""
    (valid_lrrome / "REF_proteins.fasta").write_text(
        ">gene1\nMPEPTIDE\n"
        ">gene1\nMPEPTIDE\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError):
        validate(valid_lrrome, valid_gff)


def test_validation_report_is_written(
    valid_lrrome: Path,
    valid_gff: Path,
    tmp_path: Path,
) -> None:
    """Write a validation report after successful validation."""
    report_path = tmp_path / "lrrome_validation.tsv"

    validator = LrromeValidator(
        lrrome_path=valid_lrrome,
        reference_gff_path=valid_gff,
    )

    validator.validate()
    validator.write_report(report_path)

    report = report_path.read_text(encoding="utf-8")

    assert "status\tOK" in report
    assert "models\t1" in report
    assert "exon_files\t1" in report
    assert "ref_gff_gene_ids\tOK" in report
