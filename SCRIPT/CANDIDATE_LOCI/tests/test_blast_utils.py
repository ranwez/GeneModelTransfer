import pytest
import polars as pl
from CANDIDATE_LOCI.blast_utils import parse_blast_results, blast_to_HSPs, HSP, HSP_chr
from pathlib import Path

def test_blast_tsv2tuples():
    test_data_path = Path(__file__).parent / "data" / "ENSG00000169598_tblastn_2chr.tsv"
    tuples = blast_to_HSPs(test_data_path)
    assert len(tuples) == 2