import pytest
import polars as pl
from CANDIDATE_LOCI.candidate_loci import  CandidateLocus, ParametersExpansion , find_candidate_loci, Bounds
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

def test_find_candidate_loci_exp_def():
    ENSG00000169598_tsv = Path(__file__).parent / "data" / "ENSG00000169598_tblastn_2chr.tsv"
    ENSG00000169598_gff = Path(__file__).parent / "data" / "ENSG00000169598.gff"
    CandidateLocus = find_candidate_loci(ENSG00000169598_gff,ENSG00000169598_tsv,  ParametersExpansion(nb_aa_for_missing_part=10, nb_nt_default=300, nb_nt_when_missing_part=3000))
    chr1_bounds = CandidateLocus["1"][0].chr_bounds
    assert (chr1_bounds.start == (3857604-300))
    assert (chr1_bounds.end == (3883741+300))
    chr2_bounds = CandidateLocus["2"][0].chr_bounds
    assert (chr2_bounds.start == (3857604-300))
    assert (chr2_bounds.end == (3883741+300))
   
def test_find_candidate_loci_numeroushsps():
    ENSG00000081870_tsv = Path(__file__).parent / "data" / "ENSG00000081870_tblastn.tsv"
    ENSG00000081870_gff = Path(__file__).parent / "data" / "ENSG00000081870.gff"
    CandidateLocus = find_candidate_loci(ENSG00000081870_gff,ENSG00000081870_tsv)
    best_locus = CandidateLocus["1"][0]
    assert best_locus.chr_bounds.start <= 53916790
    assert best_locus.chr_bounds.end >= 53940082

def test_cluster_Chr2B_0013899XXX():
    param_ext= ParametersExpansion(nb_aa_for_missing_part=10, nb_nt_default=300, nb_nt_when_missing_part=3000)
    Chr2B_0013899XXX_tsv = Path(__file__).parent / "data" / "clusters"/"Chr2B_0013899XXX_tblastn.tsv"
    Chr2B_0013899XXX_gff = Path(__file__).parent / "data" / "IRGSP_SVEVO_JULY_LRR.gff"
    CandidateLocus = find_candidate_loci(Chr2B_0013899XXX_gff,Chr2B_0013899XXX_tsv, param_ext)
    expBounds=( Bounds(13899741,13903702), Bounds(13923253,13927271), Bounds(13946742,13950765), Bounds(13970328,13974345), Bounds(13993920,13997937), Bounds(14017525,14021500))
    assert len(CandidateLocus["Chr2B"]) == 6
    predicted_bounds = [locus.chr_bounds for locus in CandidateLocus["Chr2B"]]
    predicted_bounds.sort(key=lambda bound: bound.start)
    for i, predBound in enumerate(predicted_bounds):
        assert predBound.start <= expBounds[i].start
        assert predBound.end >= expBounds[i].end

def test_cluster_Chr7A_0752599XXX():
    param_ext= ParametersExpansion(nb_aa_for_missing_part=10, nb_nt_default=300, nb_nt_when_missing_part=3000)
    Chr7A_0752599XXX_tsv = Path(__file__).parent / "data" / "clusters"/"Chr7A_0752599XXX_tblastn.tsv"
    Chr7A_0752599XXX_gff = Path(__file__).parent / "data"  / "IRGSP_SVEVO_JULY_LRR.gff"
    CandidateLocus = find_candidate_loci(Chr7A_0752599XXX_gff, Chr7A_0752599XXX_tsv, param_ext)
    expBounds=( Bounds(752599479,752606036), Bounds(752615666,752622223), Bounds(752631527,752638408),Bounds(752647965,752654608), Bounds(753266855,753274882))
    predicted_bounds = [locus.chr_bounds for locus in CandidateLocus["Chr7A"]]
    chrPath=CandidateLocus["Chr7A"][2].chr_path
    predicted_score = [locus.score for locus in CandidateLocus["Chr7A"]]
    predicted_bounds.sort(key=lambda bound: bound.start)
    assert len(CandidateLocus["Chr7A"]) == len(expBounds)
    for i, predBound in enumerate(predicted_bounds):
        assert predBound.start <= expBounds[i].start
        assert predBound.end >= expBounds[i].end

# HERE
def test_cluster_Chr2B_0004452XXX():
    param_ext= ParametersExpansion(nb_aa_for_missing_part=10, nb_nt_default=300, nb_nt_when_missing_part=3000)
    Chr2B_0004452XXX_tsv = Path(__file__).parent / "data" / "clusters"/"Chr2B_0004452XXX_tblastn.tsv"
    Chr2B_0004452XXX_gff = Path(__file__).parent / "data"  / "IRGSP_SVEVO_JULY_LRR.gff"
    CandidateLocus = find_candidate_loci(Chr2B_0004452XXX_gff, Chr2B_0004452XXX_tsv, param_ext)
    expBounds=( Bounds(752599479,752606036), Bounds(752615666,752622223), Bounds(752631527,752638408),Bounds(752647965,752654608), Bounds(753266855,753274882))
    predicted_bounds = [locus.chr_bounds for locus in CandidateLocus["Chr2B"]]
    predicted_score = [locus.score for locus in CandidateLocus["Chr2B"]]
    predicted_bounds.sort(key=lambda bound: bound.start)
    assert len(CandidateLocus["Chr2B"]) == len(expBounds)
    for i, predBound in enumerate(predicted_bounds):
        assert predBound.start <= expBounds[i].start
        assert predBound.end >= expBounds[i].end

def Xtest_bug_Chr6A_0017833XXX():
    param_ext= ParametersExpansion(nb_aa_for_missing_part=10, nb_nt_default=300, nb_nt_when_missing_part=3000)
    Chr6A_0017833XXX_tsv = Path(__file__).parent / "data" / "bug"/"Chr6A_0017833XXX_tblastn.tsv"
    Chr6A_0017833XXX_gff = Path(__file__).parent / "data" / "bug"/"IRGSP_DWSvevo3January_LRR.gff"
    CandidateLocus = find_candidate_loci(Chr6A_0017833XXX_gff, Chr6A_0017833XXX_tsv, param_ext)
    expBounds=( Bounds(13899741,13903702), Bounds(13923253,13927271), Bounds(13946742,13950765), Bounds(13970328,13974345), Bounds(13993920,13997937), Bounds(14017525,14021500))
    assert len(CandidateLocus["Chr6A"]) == len( expBounds)
    predicted_bounds = [locus.chr_bounds for locus in CandidateLocus["Chr6A"]]
    predicted_bounds.sort(key=lambda bound: bound.start)
    for i, predBound in enumerate(predicted_bounds):
        assert predBound.start <= expBounds[i].start
        assert predBound.end >= expBounds[i].en
