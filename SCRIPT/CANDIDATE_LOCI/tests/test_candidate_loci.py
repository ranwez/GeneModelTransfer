import pytest
import polars as pl
from CANDIDATE_LOCI.candidate_loci import  CandidateLocus, ParametersCandidateLoci , find_candidate_loci, Bounds, find_candidate_loci_from_file, ParametersExpansion
from CANDIDATE_LOCI.blast_utils import blast_to_sortedHSPs
from pathlib import Path

def test_find_candidate_loci():
    ENSG00000169598_tsv = Path(__file__).parent / "data" / "ENSG00000169598_tblastn_2chr.tsv"
    ENSG00000169598_gff = Path(__file__).parent / "data" / "ENSG00000169598.gff"
    params = ParametersCandidateLoci(expansion=None)
    CandidateLocus = find_candidate_loci(ENSG00000169598_gff,ENSG00000169598_tsv, params)
    chr1_bounds = CandidateLocus["1"][0].chr_bounds
    assert (chr1_bounds.start == 3857604)
    assert (chr1_bounds.end == 3883741)
    chr2_bounds = CandidateLocus["2"][0].chr_bounds
    assert (chr2_bounds.start == 3857604)
    assert (chr2_bounds.end == 3883741)

def test_find_candidate_loci_exp_def():
    ENSG00000169598_tsv = Path(__file__).parent / "data" / "ENSG00000169598_tblastn_2chr.tsv"
    ENSG00000169598_gff = Path(__file__).parent / "data" / "ENSG00000169598.gff"
    CandidateLocus = find_candidate_loci(ENSG00000169598_gff,ENSG00000169598_tsv)
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
    # The best locus is the one with the highest score need to sort the loci by score to assert
    CandidateLocus["1"].sort(key=lambda locus: locus.score, reverse=True)
    best_locus = CandidateLocus["1"][0]
    assert best_locus.chr_bounds.start <= 53916790
    assert best_locus.chr_bounds.end >= 53940082

def test_cluster_Chr2B_0013899XXX():
    Chr2B_0013899XXX_tsv = Path(__file__).parent / "data" / "clusters"/"Chr2B_0013899XXX_tblastn.tsv"
    Chr2B_0013899XXX_gff = Path(__file__).parent / "data" / "IRGSP_SVEVO_JULY_LRR.gff"
    CandidateLocus = find_candidate_loci(Chr2B_0013899XXX_gff,Chr2B_0013899XXX_tsv)
    expBounds=( Bounds(13899741,13903702), Bounds(13923253,13927271), Bounds(13946742,13950765), Bounds(13970328,13974345), Bounds(13993920,13997937), Bounds(14017525,14021500))
    assert len(CandidateLocus["Chr2B"]) == 6
    predicted_bounds = [locus.chr_bounds for locus in CandidateLocus["Chr2B"]]
    predicted_bounds.sort(key=lambda bound: bound.start)
    for i, predBound in enumerate(predicted_bounds):
        assert predBound.start <= expBounds[i].start
        assert predBound.end >= expBounds[i].end

def test_cluster_Chr7A_0752599XXX():
    Chr7A_0752599XXX_tsv = Path(__file__).parent / "data" / "clusters"/"Chr7A_0752599XXX_tblastn.tsv"
    Chr7A_0752599XXX_gff = Path(__file__).parent / "data"  / "IRGSP_SVEVO_JULY_LRR.gff"
    CandidateLocus = find_candidate_loci(Chr7A_0752599XXX_gff, Chr7A_0752599XXX_tsv)
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
    Chr2B_0004452XXX_tsv = Path(__file__).parent / "data" / "clusters"/"Chr2B_0004452XXX_tblastn_debug.tsv"
    Chr2B_0004452XXX_gff = Path(__file__).parent / "data"  / "IRGSP_SVEVO_JULY_LRR.gff"
    CandidateLocus = find_candidate_loci(Chr2B_0004452XXX_gff, Chr2B_0004452XXX_tsv,ParametersCandidateLoci(expansion=None))
    expBounds=( Bounds(4452160,4456183), Bounds(4456186,4457979), Bounds(4457982,4459778),Bounds(4459781,4461991))
    predicted_bounds = [locus.chr_bounds for locus in CandidateLocus["Chr2B"]]
    predicted_bounds.sort(key=lambda bound: bound.start)
    assert len(CandidateLocus["Chr2B"]) == len(expBounds)
    for i, predBound in enumerate(predicted_bounds):
        assert abs(predBound.start - expBounds[i].start)<100
        assert abs(predBound.end -expBounds[i].end)<100

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

### test the version using intermediate file
def test_find_candidate_loci_inter(tmp_path):
    ENSG00000169598_tsv = Path(__file__).parent / "data" / "ENSG00000169598_tblastn_2chr.tsv"
    ENSG00000169598_gff = Path(__file__).parent / "data" / "ENSG00000169598.gff"
    params = ParametersCandidateLoci(expansion=None)
    # Create a temporary output file in the pytest temporary directory
    temp_sorted_blast = tmp_path / "test_output_blast.tsv"
    blast_to_sortedHSPs(ENSG00000169598_tsv,temp_sorted_blast)
    CandidateLocus = find_candidate_loci_from_file(ENSG00000169598_gff,temp_sorted_blast, params)
    chr1_bounds = CandidateLocus["1"][0].chr_bounds
    assert (chr1_bounds.start == 3857604)
    assert (chr1_bounds.end == 3883741)
    chr2_bounds = CandidateLocus["2"][0].chr_bounds
    assert (chr2_bounds.start == 3857604)
    assert (chr2_bounds.end == 3883741)

def test_find_candidate_loci_exp_def_inter(tmp_path):
    ENSG00000169598_tsv = Path(__file__).parent / "data" / "ENSG00000169598_tblastn_2chr.tsv"
    ENSG00000169598_gff = Path(__file__).parent / "data" / "ENSG00000169598.gff"
    CandidateLocus = find_candidate_loci(ENSG00000169598_gff,ENSG00000169598_tsv)
    temp_sorted_blast = tmp_path / "test_output_blast.tsv"
    blast_to_sortedHSPs(ENSG00000169598_tsv,temp_sorted_blast)
    chr1_bounds = CandidateLocus["1"][0].chr_bounds
    assert (chr1_bounds.start == (3857604-300))
    assert (chr1_bounds.end == (3883741+300))
    chr2_bounds = CandidateLocus["2"][0].chr_bounds
    assert (chr2_bounds.start == (3857604-300))
    assert (chr2_bounds.end == (3883741+300))
   
def test_find_candidate_loci_numeroushsps_inter(tmp_path):
    ENSG00000081870_tsv = Path(__file__).parent / "data" / "ENSG00000081870_tblastn.tsv"
    ENSG00000081870_gff = Path(__file__).parent / "data" / "ENSG00000081870.gff"
    temp_sorted_blast = tmp_path / "test_output_blast.tsv"
    blast_to_sortedHSPs(ENSG00000081870_tsv,temp_sorted_blast)
    CandidateLocus = find_candidate_loci(ENSG00000081870_gff,temp_sorted_blast)
    # The best locus is the one with the highest score need to sort the loci by score to assert
    CandidateLocus["1"].sort(key=lambda locus: locus.score, reverse=True)
    best_locus = CandidateLocus["1"][0]
    assert best_locus.chr_bounds.start <= 53916790
    assert best_locus.chr_bounds.end >= 53940082

def test_cluster_Chr2B_0013899XXX_inter(tmp_path):
    Chr2B_0013899XXX_tsv = Path(__file__).parent / "data" / "clusters"/"Chr2B_0013899XXX_tblastn.tsv"
    Chr2B_0013899XXX_gff = Path(__file__).parent / "data" / "IRGSP_SVEVO_JULY_LRR.gff"
    temp_sorted_blast = tmp_path / "test_output_blast.tsv"
    blast_to_sortedHSPs(Chr2B_0013899XXX_tsv,temp_sorted_blast)
    CandidateLocus = find_candidate_loci(Chr2B_0013899XXX_gff,temp_sorted_blast)
    expBounds=( Bounds(13899741,13903702), Bounds(13923253,13927271), Bounds(13946742,13950765), Bounds(13970328,13974345), Bounds(13993920,13997937), Bounds(14017525,14021500))
    assert len(CandidateLocus["Chr2B"]) == 6
    predicted_bounds = [locus.chr_bounds for locus in CandidateLocus["Chr2B"]]
    predicted_bounds.sort(key=lambda bound: bound.start)
    for i, predBound in enumerate(predicted_bounds):
        assert predBound.start <= expBounds[i].start
        assert predBound.end >= expBounds[i].end

def test_cluster_Chr7A_0752599XXX_inter(tmp_path):
    Chr7A_0752599XXX_tsv = Path(__file__).parent / "data" / "clusters"/"Chr7A_0752599XXX_tblastn.tsv"
    Chr7A_0752599XXX_gff = Path(__file__).parent / "data"  / "IRGSP_SVEVO_JULY_LRR.gff"
    temp_sorted_blast = tmp_path / "test_output_blast.tsv"
    blast_to_sortedHSPs(Chr7A_0752599XXX_tsv,temp_sorted_blast)
    CandidateLocus = find_candidate_loci(Chr7A_0752599XXX_gff, temp_sorted_blast)
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
    Chr2B_0004452XXX_tsv = Path(__file__).parent / "data" / "clusters"/"Chr2B_0004452XXX_tblastn_debug.tsv"
    Chr2B_0004452XXX_gff = Path(__file__).parent / "data"  / "IRGSP_SVEVO_JULY_LRR.gff"
    CandidateLocus = find_candidate_loci(Chr2B_0004452XXX_gff, Chr2B_0004452XXX_tsv,ParametersCandidateLoci(expansion=None))
    expBounds=( Bounds(4452160,4456183), Bounds(4456186,4457979), Bounds(4457982,4459778),Bounds(4459781,4461991))
    predicted_bounds = [locus.chr_bounds for locus in CandidateLocus["Chr2B"]]
    predicted_bounds.sort(key=lambda bound: bound.start)
    assert len(CandidateLocus["Chr2B"]) == len(expBounds)
    for i, predBound in enumerate(predicted_bounds):
        assert abs(predBound.start - expBounds[i].start)<100
        assert abs(predBound.end -expBounds[i].end)<100

def test_cluster_Chr2B_0004452XXX_gff():
    Chr2B_0004452XXX_tsv = Path(__file__).parent / "data" / "clusters"/"Chr2B_0004452XXX_tblastn_debug.tsv"
    Chr2B_0004452XXX_gff = Path(__file__).parent / "data"  / "IRGSP_SVEVO_JULY_LRR.gff"
    paramExp=ParametersExpansion(nb_aa_for_missing_part=10, nb_nt_default=300, nb_nt_when_missing_part=3000, template_gff=Chr2B_0004452XXX_gff)
    CandidateLocus = find_candidate_loci(Chr2B_0004452XXX_gff, Chr2B_0004452XXX_tsv,ParametersCandidateLoci(expansion=paramExp))
    expBounds=( Bounds(4452160,4456183), Bounds(4456186,4457979), Bounds(4457982,4459778),Bounds(4459781,4461991))
    predicted_bounds = [locus.chr_bounds for locus in CandidateLocus["Chr2B"]]
    predicted_bounds.sort(key=lambda bound: bound.start)
    assert len(CandidateLocus["Chr2B"]) == len(expBounds)
    for i, predBound in enumerate(predicted_bounds):
        assert abs(predBound.start - expBounds[i].start)<=300
        assert abs(predBound.end -expBounds[i].end)<=300

def test_missing_hsp_Chr1_ENSMMUG00000004466 ():
    ENSMMUG00000004466_tsv = Path(__file__).parent / "data" / "missing_hsp"/ "ENSMMUG00000004466.tsv"
    ENSMMUG00000004466_gff = Path(__file__).parent / "data"  /  "missing_hsp"/ "ENSMMUG00000004466.gff"
    paramExp=ParametersExpansion(nb_aa_for_missing_part=10, nb_nt_default=300, nb_nt_when_missing_part=3000, template_gff=ENSMMUG00000004466_gff)
    CandidateLocus = find_candidate_loci(ENSMMUG00000004466_gff, ENSMMUG00000004466_tsv,ParametersCandidateLoci(expansion=paramExp))
    expBounds=Bounds(24419540, 24470086)
    
    CandidateLocus["1"].sort(key=lambda locus: locus.score, reverse=True)
    best_locus = CandidateLocus["1"][0]
    print (best_locus.chr_bounds)
    assert best_locus.chr_bounds.start <= expBounds.start
    assert best_locus.chr_bounds.end >= expBounds.end

    assert abs (best_locus.chr_bounds.start - expBounds.start) < 500
    assert abs (best_locus.chr_bounds.end - expBounds.end) < 500
###
def test_blast_tsv2file():
    #test_data_path = Path(__file__).parent / "data" / "ENSG00000169598_tblastn_2chr.tsv"
    #blast_to_sortedHSPs("/Users/ranwez/Desktop/TEST_REGION/blast_refProt.tsv","__test_output_blast.tsv")
    IRGSP_SVEVO_gff = Path(__file__).parent / "data"  / "IRGSP_SVEVO_JULY_LRR.gff"
    param_ext= ParametersExpansion(nb_aa_for_missing_part=10, nb_nt_default=300, nb_nt_when_missing_part=3000)
    candidateLoci=find_candidate_loci_from_file(IRGSP_SVEVO_gff,"__test_output_blast.tsv",ParametersCandidateLoci(expansion=param_ext))
    # Write results to the output GFF file
    with open("__test_output_blast_july.gff", "w") as out_file:
        for chr in candidateLoci:
            for locus in candidateLoci[chr]:
                out_file.write(locus.as_gff() + "\n")

    # Write results to the list query/target file
    with open("list_query_target_july.output_list", "w") as out_file:
        for chr in candidateLoci:
            for locus in candidateLoci[chr]:
                out_file.write(locus.as_query_target() + "\n")

def test_genomic_coordinate_issue ():
    OSJnip_Chr01_02834335_tsv = Path(__file__).parent / "data" / "OSJnip_Chr01_02834335_tblastn.tsv"
    OSJnip_Chr01_02834335_gff = Path(__file__).parent / "data"  /  "OSJnip_Chr01_02834335_tblastn.gff"
    paramExp=ParametersExpansion(nb_aa_for_missing_part=10, nb_nt_default=300, nb_nt_when_missing_part=3000, template_gff=OSJnip_Chr01_02834335_gff)
    CandidateLocus = find_candidate_loci(OSJnip_Chr01_02834335_gff, OSJnip_Chr01_02834335_tsv,ParametersCandidateLoci(expansion=paramExp))

    
    CandidateLocus["Chr1"].sort(key=lambda locus: locus.score, reverse=True)
    best_locus = CandidateLocus["Chr1"][0]
    print (best_locus.chr_bounds)