import pytest
import polars as pl
from CANDIDATE_LOCI.candidate_loci import  CandidateLocus , find_candidate_loci
from pathlib import Path

def test_find_candidate_loci():
    ENSG00000169598_tsv = Path(__file__).parent / "data" / "ENSG00000169598_tblastn_2chr.tsv"
    ENSG00000169598_gff = Path(__file__).parent / "data" / "ENSG00000169598.gff"
    CandidateLocus = find_candidate_loci(ENSG00000169598_gff,ENSG00000169598_tsv)
    chr1_bounds = CandidateLocus["1"][0].chr_bounds
    assert (chr1_bounds.start == 3857604)
    assert (chr1_bounds.end == 3883741)
    chr2_bounds = CandidateLocus["2"][0].chr_bounds
    assert (chr2_bounds.start == 3857604)
    assert (chr2_bounds.end == 3883741)
def test_find_candidate_loci_numeroushsps():
    ENSG00000081870_tsv = Path(__file__).parent / "data" / "ENSG00000081870_tblastn.tsv"
    ENSG00000081870_gff = Path(__file__).parent / "data" / "ENSG00000081870.gff"
    CandidateLocus = find_candidate_loci(ENSG00000081870_gff,ENSG00000081870_tsv)
    best_locus = CandidateLocus["1"][0]
    assert best_locus.chr_bounds.start <= 53916790
    assert best_locus.chr_bounds.end >= 53940082