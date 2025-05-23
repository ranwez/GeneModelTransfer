import pytest
from CANDIDATE_LOCI.blast_utils import  blast_to_HSPs, blast_to_sortedHSPs
from pathlib import Path

def test_blast_tsv2tuples():
    test_data_path = Path(__file__).parent / "data" / "ENSG00000169598_tblastn_2chr.tsv"
    tuples = blast_to_HSPs(test_data_path)
    assert len(tuples) == 2

def test_blast_tsv2file():
    test_data_path = Path(__file__).parent / "data" / "ENSG00000169598_tblastn_2chr.tsv"
    blast_to_sortedHSPs("/Users/ranwez/Desktop/TEST_REGION/blast_refProt.tsv","__test_output_blast.tsv")
