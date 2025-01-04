import pytest
import polars as pl
from CANDIDATE_LOCI.gff_utils import get_coding_regions, get_longest_intron, parse_gff, GeneInfo, gff_to_geneInfo
from pathlib import Path

def test_gff_to_geneInfo():
    test_data_path = Path(__file__).parent / "data" / "ENSG00000160679.gff"
    (prot_dict, def_intron_lg) = gff_to_geneInfo(test_data_path,0.7)
    exp_info= GeneInfo (
        gene_id="gene:ENSG00000160679",
        chr_id="1",
        coding_start=153634301,
        coding_end=153645269,
        strand=1,
        longest_intron=3797
    )
    assert prot_dict["gene:ENSG00000160679"] == exp_info
    assert len(prot_dict) == 1
    assert def_intron_lg == 3797

