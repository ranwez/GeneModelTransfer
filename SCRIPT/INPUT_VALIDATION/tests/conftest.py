"""Shared pytest fixtures for input validation tests."""

import sys
from pathlib import Path

import pytest


VALID_GFF = """\
chr1	test	gene	100	900	.	+	.	ID=gene1
chr1	test	mRNA	150	850	.	+	.	ID=mrna1;Parent=gene1
chr1	test	CDS	200	300	.	+	0	Parent=mrna1
chr1	test	CDS	400	500	.	+	0	ID=cds2;Parent=mrna1
"""


TEST_DIR = Path(__file__).resolve().parent
VALIDATION_DIR = TEST_DIR.parent

sys.path.insert(0, str(VALIDATION_DIR))


@pytest.fixture
def valid_gff(tmp_path: Path) -> Path:
    """Create a valid reference GFF."""
    path = tmp_path / "reference.gff"
    path.write_text(VALID_GFF, encoding="utf-8")
    return path


@pytest.fixture
def valid_fai(tmp_path: Path) -> Path:
    """Create a FASTA index containing the reference seqid."""
    path = tmp_path / "reference.fasta.fai"
    path.write_text(
        "chr1\t2000\t6\t80\t81\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture
def valid_info_locus(tmp_path: Path) -> Path:
    """Create locus information containing the reference gene."""
    path = tmp_path / "info_locus.tsv"
    path.write_text(
        "gene1\tNLR\tCanonical\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture
def valid_lrrome(tmp_path: Path) -> Path:
    """Create a minimal valid LRRome."""
    lrrome = tmp_path / "LRRome"
    lrrome.mkdir()

    for dirname in [
        "REF_PEP",
        "REF_EXONS",
        "REF_cDNA",
        "REF_LOCI",
        "REF_LOCI_GFF",
    ]:
        (lrrome / dirname).mkdir()

    (lrrome / "REF_proteins.fasta").write_text(
        ">gene1\nMPEPTIDE\n",
        encoding="utf-8",
    )
    (lrrome / "REF_loci.fasta").write_text(
        ">gene1\nATGCGTACGTAG\n",
        encoding="utf-8",
    )

    (lrrome / "REF_PEP" / "gene1").write_text(
        ">gene1\nMPEPTIDE\n",
        encoding="utf-8",
    )
    (lrrome / "REF_cDNA" / "gene1").write_text(
        ">gene1\nATGCGTACGTAG\n",
        encoding="utf-8",
    )
    (lrrome / "REF_LOCI" / "gene1").write_text(
        ">gene1\nATGCGTACGTAG\n",
        encoding="utf-8",
    )

    (lrrome / "REF_EXONS" / "gene1_CDS_1").write_text(
        ">gene1_CDS_1\nATGCGT\n",
        encoding="utf-8",
    )

    (lrrome / "REF_LOCI_GFF" / "gene1.gff").write_text(
        "gene1\ttest\tgene\t1\t12\t.\t+\t.\tID=gene1\n",
        encoding="utf-8",
    )

    return lrrome


@pytest.fixture
def valid_target_fai(tmp_path: Path) -> Path:
    """Create a target FASTA index containing one sequence."""
    path = tmp_path / "target.fasta.fai"
    path.write_text(
        "chrTarget\t5000\t11\t80\t81\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture
def valid_tblastn(tmp_path: Path) -> Path:
    """Create a valid TBLASTN result table."""
    path = tmp_path / "tblastn.tsv"
    path.write_text(
        "gene1\tchrTarget\t100\t90\t1\t90\t100\t369\t80\t88.9\t0\t1e-20\t100\t85\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture
def valid_blastn(tmp_path: Path) -> Path:
    """Create a valid BLASTN result table."""
    path = tmp_path / "blastn.tsv"
    path.write_text(
        "gene1\tchrTarget\t700\t680\t1\t680\t500\t1179\t650\t95.6\t0\t1e-50\t500\t670\n",
        encoding="utf-8",
    )
    return path
