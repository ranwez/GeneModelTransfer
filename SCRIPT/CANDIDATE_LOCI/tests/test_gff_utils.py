import pytest
import polars as pl
from CANDIDATE_LOCI.gff_utils import get_coding_regions, get_longest_intron, parse_gff, GeneInfo, gff_to_geneInfo, sort_gff
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

def test_gff_sort():
    test_data_path = Path(__file__).parent / "data" / "ENSG000_unsorted.gff"
    df = parse_gff(test_data_path)
    df = sort_gff(df)
    assert df["seqid"].is_sorted()
    last_gene = None
    last_mNRA = None
    last_other = None
    for feature in df.to_dicts():
        if feature["type"]=="gene":
            if last_gene is not None:
                assert last_gene["start"] <= last_gene["end"]
            last_gene = feature
        elif feature["type"]=="mRNA":
            assert feature["gene"] == last_gene["gene"]
            if last_mNRA is not None:
                assert last_mNRA["start"] <= last_mNRA["end"]
            last_mNRA = feature
        else:
            if last_other is not None:
                assert last_other["start"] <= last_other["end"]
            assert feature["gene"] == last_gene["gene"]
            assert feature["mRNA"] == last_mNRA["mRNA"]
            last_other = feature