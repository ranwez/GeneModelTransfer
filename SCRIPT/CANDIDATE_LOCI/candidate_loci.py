import argparse
import csv
import math
from typing import Optional, Union, TextIO
from io import StringIO
from pathlib import Path
from attrs import define, field
from CANDIDATE_LOCI.gff_utils import GeneInfo, gff_to_geneInfo, gff_to_cdsInfo, CdsInfo
from CANDIDATE_LOCI.interlap import InterLap, Interval
from CANDIDATE_LOCI.blast_utils import HSP, HSP_chr, blast_to_sortedHSPs
from CANDIDATE_LOCI.bounds import Bounds, RangeCoverage


##############################################################################
# Candidate loci finding
# The goal of this module is to find candidate loci in the genome that are homologous to a protein.
# The input is a GFF file with the gene annotations and a BLAST file with the protein homology information.
# The output is a list of candidate loci : a pair  (genomic_region, model protein).
# The goal is to identify the best protein to annotate each candidate loci.
# ##############################################################################

# avoid importing numpy just for this
def argmax(x):
    return max(range(len(x)), key=lambda i: x[i])

# by default the region is expand over HSP coordinates by 300 nt on each side, 
# but if more than 10 AA are missing on one side the extenstion is 3000 on this side
# if the template gff is provided the extension is based on the genomic size of the corresponding missing part in the template protein
@define
class ParametersExpansion:
    nb_aa_for_missing_part: int = field(default=10)
    nb_nt_default: int = field(default=300)
    nb_nt_when_missing_part: int = field(default=3000)
    template_gff: Optional[str] = field(default=None)

# HSPs of the sam protein are allowd to be merged in a single candidat loci 
# if they are separated by less than the maximum allowed distance (max intron length)
# this distance could be a fixed number (4000) or a quantile (0.5 for median) of the intron length distribution of the input GFF 
@define
class ParametersHspClustering:
    maxIntronLength : int = field(default=4000)
    quantileForMaxIntronLength : float = field(default=0.5)
    useQuantile : bool = field(default=False)

# Candidate loci are penalized if their length is too different from the genomic length of the protein
# the penalty is expressed as a number of extra intron length using a len_penalty_percentage parameter
# the higher the len_penalty_percentage the more the penalty is important default is 0.1
# length_penalty_in_percent = max(-1, len_penalty_percentage * (1 - max(1, math.exp(length_deviation))))
# similarity is a weighted average of pc_identity and pc_coverage obtained using
# similarity = (identityWeigthvsCoverage * self.nident + max_homology)/(identityWeigthvsCoverage+1); pc_similarity = similarity/ali_lg
# Candidate loci with low score or similarity (percentage of identity with template protein) are discared   
# finally loci kept cannot overlap; a shrink of nt_shrink nucleotide is tested to try to keep loci with small overlap  
@define
class ParametersLociScoring:
    length_penalty_percentage: float = field(default=0.1)
    identityWeigthvsCoverage: int = field(default=3)
    min_similarity: float = field(default=0.25)
    min_score: float = field(default=50)
    nt_shrink: int = field(default=60)
    identity_vs_missed_importance: int = field(default=4) # not used hard coded for the moment
    
@define
class ParametersCandidateLoci:
    expansion: ParametersExpansion = field(default=ParametersExpansion())
    hsp_clustering: ParametersHspClustering = field(default=ParametersHspClustering())
    loci_scoring: ParametersLociScoring = field(default=ParametersLociScoring())
    skip_neighborhood_dist: Optional[int] = field(default=None) # for evaluation purpose not reannotating a locus with itself or its neighbors

@define(slots=True,frozen=True)
class HspCompatibleOverlap:
    compatible: bool
    max_id_overlap: int

# 
@ define(slots=True)
class HspOverlapCacher:
    hsp_overlaps: list[HspCompatibleOverlap]
    nb_hsp: int
    @classmethod
    def from_hsp(cls, hsp:HSP) -> "HspOverlapCacher":
        return cls(
            hsp_overlaps = cls.pre_compute_overlap(hsp),
            nb_hsp = len(hsp)
        ) 
    @classmethod
    def pre_compute_overlap(cls,hsps:list[HSP])-> list[HspCompatibleOverlap]:
        n = len(hsps)
        # store precomputed overlap between each pair of HSPs in 1D, lower triangular matrix
        hsp_overlaps = [HspCompatibleOverlap(False, 0)] * (n*(n-1)//2)
        for i in range(1,n):
            row_idx = (i*(i-1))//2
            hspi=hsps[i]
            for j in range(i):
                hspj=hsps[j]
                if ((hspj.prot_bounds.start> hspi.prot_bounds.start) and hspj.locS_bounds.start > hspi.locS_bounds.start) \
                or ((hspj.prot_bounds.start< hspi.prot_bounds.start) and hspj.locS_bounds.start < hspi.locS_bounds.start):
                    loc_overlap = math.ceil(hspi.locS_bounds.overlap(hspj.locS_bounds)/3)
                    prot_overlap = hspi.prot_bounds.overlap(hspj.prot_bounds)
                    # if the two overlap size are too different compare to the hsp size, the hsp are not really compatible 
                    # the same region is used twice (duplication)
                    diff_overlap = abs(loc_overlap - prot_overlap)
                    hspij_min_length = min( hspi.prot_bounds.length(),hspj.prot_bounds.length(), hspi.locS_bounds.length()//3, hspj.locS_bounds.length()//3)
                    if diff_overlap > 0.5 * hspij_min_length:
                        continue
                    max_id_overlap = min(
                        max(loc_overlap, prot_overlap),  # max hsp overlap with locus and prot
                        hspj.nident,                   # Identities in hsp_j
                        hspi.nident                    # Identities in hsp_i
                    )
                    hsp_overlaps[ row_idx+ j]= HspCompatibleOverlap(True, max_id_overlap)
        return hsp_overlaps
    def max_overlap(self, i:int, j:int) -> HspCompatibleOverlap:
        if i < j:
            i, j = j, i  # Ensure j <= i for lower triangular part
        return self.hsp_overlaps[(i * (i - 1)) // 2 + j]


# use left and right instead of start and end to avoid confusion with gene orientation
# left is the start of the genomic region, right is the end of the genomic region
@define(slots=True)
class ExpansionPair:
    left: int
    right: int
    def __str__(self):
        return f"({self.left}, {self.right})"
    def update_left(self, new_left: int):
        # raise an error if new_left is negative
        if new_left < 0:
            raise ValueError(f"new_start cannot be negative: {new_left}")
        if new_left > self.left:
            self.left = new_left
    def update_right(self, new_right: int):
        # raise an error if new_right is negative
        if new_right < 0:
            raise ValueError(f"new_right cannot be negative: {new_right}")
        if new_right > self.right:
            self.right = new_right  
    def swap(self):
        # swap start and end
        self.left, self.right = self.right, self.left

@define(slots=True)
class CandidateLocus:
    chr_id: str
    strand: int
    prot_id: str
    prot_len: int
    chr_path: list[Bounds]
    prot_path: list[Bounds]
    chr_bounds: Bounds
    prot_bounds: Bounds
    nident: int
    nhomol_loc: int
    nhomol_prot: int
    pc_similarity: Optional[float]= field(default=None)
    score: Optional[float]= field(default=None)
    expansion: Optional[ExpansionPair]= field(default=None)
    shrink_info: Optional[tuple[bool,bool]]= field(default=None)

    @classmethod
    def from_hsp(cls, hsp:HSP) -> "CandidateLocus":
        return cls(
            chr_id = hsp.chr_id,
            strand = hsp.strand,
            prot_id = hsp.prot_id,
            prot_len = hsp.prot_len,
            chr_path = [Bounds.clone(hsp.loc_bounds)],
            chr_bounds = Bounds.clone(hsp.loc_bounds),
            prot_path = [Bounds.clone(hsp.prot_bounds)],
            prot_bounds = Bounds.clone(hsp.prot_bounds),
            nident = hsp.nident,
            nhomol_loc = hsp.loc_bounds.length(),
            nhomol_prot = hsp.prot_bounds.length()
            
        )
    @classmethod
    def from_hsp_path(cls, hsp_path:list[HSP], nident:int) -> "CandidateLocus":
        chr_path = [hsp.loc_bounds for hsp in hsp_path]
        prot_path = [hsp.prot_bounds for hsp in hsp_path]
        homolog_fraction_prot = Interval()# end is exclude so add 1
        for hsp_prot in prot_path:
            homolog_fraction_prot.add([(hsp_prot.start, (hsp_prot.end+1))])
        
        homolog_fraction_loc = Interval()# end is exclude so add 1
        for hsp_loc in chr_path:
            homolog_fraction_loc.add([(hsp_loc.start, (hsp_loc.end+1))])

        return cls( 
            chr_id = hsp_path[0].chr_id,
            strand = hsp_path[0].strand,
            prot_id= hsp_path[0].prot_id,
            prot_len = hsp_path[0].prot_len,
            chr_path = chr_path,
            prot_path = prot_path,
            chr_bounds = path_bound(chr_path),
            prot_bounds = path_bound(prot_path),
            nident = nident,
            nhomol_prot=homolog_fraction_prot.coverage_VR(),
            nhomol_loc = homolog_fraction_loc.coverage_VR()
            
        ) 
    def compute_score(self, protInfo:GeneInfo, max_intron_len:int, params : ParametersLociScoring) -> float:
        self.nhomol_loc= int(self.nhomol_loc/3)
        max_homology= min(self.nhomol_prot, self.nhomol_loc)
        ali_lg= max(self.prot_len, self.nhomol_loc)
        #pc_similarity is a weighted average of pc_identity and pc_coverage
        similarity = (params.identityWeigthvsCoverage * self.nident + max_homology)/(params.identityWeigthvsCoverage+1)
        self.pc_similarity =similarity/ali_lg
        if( protInfo != None): 
            prot_genomic_len = protInfo.coding_region.length()
            # how many extra intron length do we have
            len_penalty_percentage= params.length_penalty_percentage
            length_deviation=((self.chr_bounds).length() - prot_genomic_len) / max_intron_len
            length_penalty_in_percent = max(-1, len_penalty_percentage * (1 - max(1, math.exp(length_deviation))))
        else:
            length_penalty_in_percent=-0.1
        # score favor long stretch of high percentage similarity; and having a somehow similar genomic length
        self.score = similarity * self.pc_similarity * (1-length_penalty_in_percent)
        return
    def as_query_target(self)-> str:
         # Determine strand representation
        strand = '+' if self.strand == 1 else '-' if self.strand == -1 else '.'
        gene_id = self.build_id()
        return f"{gene_id}\t{self.prot_id}\t{strand}"

    def build_id(self)-> str:
        if (self.strand == 1):
            gene_id = f"{self.chr_id}_{self.chr_bounds.start:010}"
        else:
            gene_id = f"{self.chr_id}_{self.chr_bounds.end:010}"
        return gene_id
     
    def as_gff(self) -> str:
        # Determine strand representation
        strand = '+' if self.strand == 1 else '-' if self.strand == -1 else '.'
        gff_lines = []

        # Generate detailed prot_path info
        prot_path_info = ",".join(f"{bound.start}-{bound.end}" for bound in self.prot_path) if self.prot_path else "None"

        # Gene line with prot_path comment
        gene_id = self.build_id()
        gene_note = f"protId={self.prot_id};protLg={self.prot_len};prot_path={prot_path_info};score={self.score:.2f};nident={self.nident};pc_sym={self.pc_similarity:.2f}"
        gff_lines.append(
            f"{self.chr_id}\tcandidateLoci\tgene\t{self.chr_bounds.start}\t{self.chr_bounds.end}\t.\t{strand}\t.\tID={gene_id};{gene_note}"
        )

        # mRNA line
        mrna_id = f"{gene_id}_mRNA"
        gff_lines.append(
            f"{self.chr_id}\tcandidateLoci\tmRNA\t{self.chr_bounds.start}\t{self.chr_bounds.end}\t.\t{strand}\t.\tID={mrna_id};Parent={gene_id}"
        )

        # Loop over chr_path to add CDS lines (sorted by start)
        sorted_chr_path = sorted(self.chr_path, key=lambda b: b.start)
        for idx, bounds in enumerate(sorted_chr_path):
            cds_id = f"{mrna_id}_CDS_{idx + 1}"
            gff_lines.append(
                f"{self.chr_id}\tcandidateLoci\tCDS\t{bounds.start}\t{bounds.end}\t.\t{strand}\t0\tID={cds_id};Parent={mrna_id}"
            )

        # Combine all lines into a single GFF-formatted string
        return "\n".join(gff_lines)

@define(slots=True)
class PathScoring:
    nb_ident: int
    nb_missed: int
    def score(self) -> int:
        return self.nb_ident - int(self.nb_missed/12)
    @classmethod
    def eval(cls, nb_ident:int, nb_missed:int) -> int:
        return nb_ident - int(nb_missed/12)
     
@define(slots=True)
class MergeableHSP:
    prot_id: str 
    hsps: list[HSP]
    max_coord: int
    max_intron_len: int
    hsp_overlap_cache: Optional[HspOverlapCacher]=None
    covered_loci_cache: Optional[RangeCoverage]=None

    def dist_to_protein(self, protInfos:dict[GeneInfo]) -> Optional[int]:
        protInfo:GeneInfo = None
        if(protInfos.keys().__contains__(self.prot_id)):
            protInfo =protInfos[self.prot_id]
        if protInfo == None or (self.hsps[0].chr_id != protInfo.chr_id) or (protInfo == None):
            return None
        hsp_region = Bounds(self.hsps[0].loc_bounds.start, self.hsps[-1].loc_bounds.end)
        return hsp_region.distance(protInfo.coding_region)
        
    def mergeHSP(self, hsp:HSP) -> bool:
        if self.prot_id == hsp.prot_id and hsp.locS_bounds.start - self.max_coord <= self.max_intron_len:
            self.hsps.append(hsp)
            self.max_coord = max(self.max_coord, hsp.locS_bounds.end)
            return True
        return False
    
    def max_possible_coverage(self) -> float:
        if len(self.hsps) == 1:
            coverage = self.hsps[0].nident
        else:
            homolog_fraction = Interval()# end is exclude so add 1
            for hsp in self.hsps:
                homolog_fraction.add([(hsp.prot_bounds.start, (hsp.prot_bounds.end+1))])
            coverage = homolog_fraction.coverage_VR()
        return coverage/self.hsps[0].prot_len

    def compute_candidate_loci_rec(self, param: ParametersLociScoring)-> list[CandidateLocus]:
        candidateLocus=self.compute_candidate_loci()
        res=[candidateLocus]
        chr_bounds = candidateLocus.chr_bounds
        HPSr=[]
        HSPl=[]
        for hsp in self.hsps :
            if (hsp.loc_bounds.end < chr_bounds.start ) or \
              ((hsp.loc_bounds.end - chr_bounds.start ) < param.nt_shrink and hsp.loc_bounds.length() > 10 * param.nt_shrink):
                HPSr.append(hsp)
            if (chr_bounds.end < hsp.loc_bounds.start )  or \
              ((chr_bounds.end - hsp.loc_bounds.start )  < param.nt_shrink and hsp.loc_bounds.length() > 10 * param.nt_shrink):
                HSPl.append(hsp)
        for hsp_list in (HSPl, HPSr) :
            if (len(hsp_list) >0):
                sub_region= MergeableHSP(self.prot_id,hsp_list,self.max_coord,self.max_intron_len)
                res.extend(sub_region.compute_candidate_loci_rec(param))
        return res

    def compute_candidate_loci(self)->CandidateLocus:
        if len(self.hsps) == 1:
            return CandidateLocus.from_hsp(self.hsps[0]) 
        # Initialize scores and bounds, to store best resuts ending at each HSP
        self.covered_loci_cache = RangeCoverage([hsp.locS_bounds for hsp in self.hsps])
        self.hsps.sort(key=lambda hsp: hsp.locS_bounds.end)
        self.hsp_overlap_cache = HspOverlapCacher.from_hsp(self.hsps)
        
        best_prev_nident = [self.hsps[i].nident for i in range(len(self.hsps))]
        nident = [self.hsps[i].nident for i in range(len(self.hsps))]
        best_scoring_path_to =[PathScoring(nident[i], self.covered_loci_cache.get_coverage_up_to(self.hsps[i].locS_bounds.start)) for i in range(len(self.hsps))]
        paths = [[i] for i in range(len(self.hsps))]
        #return CandidateLocus.from_hsp(self.hsps[0]) 
        # Compute scores and update bounds
        for j, hspj in enumerate(self.hsps):
            # Only consider previous HSPs.
            # Reversed order for better shortcuts, the further away the hsp the most unlikely it is to get a good score
            # (since this likely lead to missing large part of the protein)
            best_prev = -1
            for i in reversed(range(j)):
                # to update with coverage
                #if (nident[i]+nident[j] < best_ident): # overlap only decrease score so no need to compute them
                #    if (best_prev_nident[i]+nident[j] < best_ident):
                #        break # nothing better earlier so stop the loop
                #    else:
                #        continue # something worth checking earlier so continue the loop 
                path_compatible_overlap = self.max_path_overlap(paths[i], j)
                nident_ij = best_scoring_path_to[i].nb_ident+nident[j]-path_compatible_overlap.max_id_overlap
                # add penalty for missing genomic region with HSP
                missed_ij = best_scoring_path_to[i].nb_missed + self.covered_loci_cache.get_coverage(self.hsps[i].locS_bounds.end, hspj.locS_bounds.start )
                score_ij = PathScoring.eval(nident_ij, missed_ij)
                # prefer solution with shorter path limit hsp included in hsp with the same scoring
                if path_compatible_overlap.compatible and score_ij > best_scoring_path_to[j].score() : #nident_ij > best_ident:
                    best_scoring_path_to[j].nb_ident = nident_ij
                    best_scoring_path_to[j].nb_missed = missed_ij
                    best_scoring_path_to[j].nb_missed = missed_ij
                    best_prev = i
            if best_prev != -1:
                paths[j] = paths[best_prev] + [j]
            best_prev_nident[j] = max(best_prev_nident[j-1], nident[j]) if j >= 1 else nident[j]
        # Find the index of the best scoring HSP
        best_idx = argmax([scoring.score() for scoring in best_scoring_path_to])
        candidate_locus= CandidateLocus.from_hsp_path([self.hsps[i] for i in paths[best_idx]], best_scoring_path_to[best_idx].nb_ident)
        return candidate_locus

    def max_path_overlap(self, prev_path:list[int], j:int)-> HspCompatibleOverlap:
        sum_overlap = 0
        for i in prev_path:
            resi = self.hsp_overlap_cache.max_overlap(i, j)
            if (not resi.compatible):
                return HspCompatibleOverlap(False, 0)
            else:
                sum_overlap += resi.max_id_overlap
        return HspCompatibleOverlap(True, sum_overlap)

def path_bound(path:list[Bounds]) -> Bounds:
    return Bounds(min(path[0].start, path[-1].start),
                 max( path[0].end, path[-1].end))

def small_shrinking(locus:CandidateLocus, covered:InterLap, ntShrink:int)-> tuple[bool,bool]:
    # if the locus does not intesect previous ones no need to shrink
    if not covered.__contains__((locus.chr_bounds.start, locus.chr_bounds.end)):
        return (False, False)
    # else can it be kept by shrinking it slightly ?
    if(locus.chr_bounds.length() < 10*ntShrink):
        return (False, False)
    if not covered.__contains__((locus.chr_bounds.start +ntShrink, locus.chr_bounds.end-ntShrink)):
        if not covered.__contains__((locus.chr_bounds.start +ntShrink, locus.chr_bounds.end)):
            locus.chr_bounds.start += ntShrink
            return (True, False)
        if not covered.__contains__((locus.chr_bounds.start, locus.chr_bounds.end-ntShrink)):
            locus.chr_bounds.end -= ntShrink
            return (False, True)
        if(locus.chr_bounds.length() >= 10*(2*ntShrink)):
            locus.chr_bounds.start += ntShrink
            locus.chr_bounds.end -= ntShrink
            return (True, True)
    return (False, False)


def keep_best_non_overlaping_loci(candidate_loci:list[CandidateLocus], params:ParametersLociScoring) -> list[CandidateLocus]:
    # Sort candidate loci by score
    candidate_loci.sort(key=lambda x: x.score, reverse=True)
    # Filter candidate loci
    filtered_candidate_loci = []
    covered = InterLap()
    for locus in candidate_loci:
        #print(locus.as_gff)
        skrink_info=small_shrinking(locus, covered, params.nt_shrink)
        if covered.__contains__((locus.chr_bounds.start, locus.chr_bounds.end)):
            continue
        else:
            covered.add((locus.chr_bounds.start, locus.chr_bounds.end))
            locus.shrink_info=skrink_info
            filtered_candidate_loci.append(locus)
    return filtered_candidate_loci

# define the expension based on the genomic size of the corresponding missing part in the template protein
def set_desired_expansion(candidate_loci:list[CandidateLocus], expand_params:ParametersExpansion, protInfo:dict) -> None:
    gff_file = expand_params.template_gff
    if gff_file is not None:
        prots = set([locus.prot_id for locus in candidate_loci])
        cds_infos = gff_to_cdsInfo(gff_file, prots)
    # set the expansion based on the genomic size of the corresponding missing part in the template protein 
    for locus in candidate_loci:
        # set default expansion
        expansion: ExpansionPair = ExpansionPair(expand_params.nb_nt_default, expand_params.nb_nt_default)
        # adjust for missing protein part using default or specific prot intron 
        expansion_for_missing_part = expand_params.nb_nt_when_missing_part
        if gff_file is not None :
            template_prot_info = protInfo.get(locus.prot_id) if gff_file is not None else None
            expansion_for_missing_part= max(expand_params.nb_nt_when_missing_part, int(template_prot_info.longest_intron * 1.1)) 
        if locus.prot_bounds.start > expand_params.nb_aa_for_missing_part:
            expansion.update_left(expansion_for_missing_part)
        if locus.prot_len - locus.prot_bounds.end > expand_params.nb_aa_for_missing_part:
            expansion.update_right(expansion_for_missing_part)
        if locus.strand == -1:
            expansion.swap()  
        # adjust for missing protein using specific prot genomic information (useful when several introns are missed, could replace previous logic when gff_file is provided)
        if gff_file is not None:
            # now adjust the expansion based on the genomic size of the corresponding missing part in the template protein
            cds_info = cds_infos[locus.prot_id]
            if locus.prot_bounds.start > 1:
                # blast can overextends HPS over non similar region we thus move forward 6 AA to handle this
                adjusted_start = min(locus.prot_len, locus.prot_bounds.start + 6)
                (gene_start, locus_start, gene_end)= cds_info.get_genomic_coord(adjusted_start, locus.strand)
                if( locus.strand==1):    
                    missing_part = int((locus_start - gene_start +1) * 1.1)
                    expansion.update_left( missing_part)
                else:
                    missing_part = int((gene_start - locus_start +1) * 1.1)
                    expansion.update_right(missing_part)
            if locus.prot_len - locus.prot_bounds.end > 0:
                # blast can overextends HPS over non similar region we thus move back 6 AA to handle this
                adjusted_end = max(1, locus.prot_bounds.end - 6)
                (_gene_start, locus_end, gene_end)= cds_info.get_genomic_coord(adjusted_end, locus.strand)
                if( locus.strand==1):    
                    missing_part = int((gene_end - locus_end +1) * 1.1)
                    expansion.right = max(expansion.right, missing_part)
                else:
                    missing_part = int((locus_end - gene_end +1) * 1.1)
                    expansion.left = max(expansion.left, missing_part)    
        locus.expansion = expansion
        
#start=4452163, end=4461979
def expands(candidate_loci:list[CandidateLocus], expand_params:ParametersExpansion, protInfo:dict):
    if (len(candidate_loci))==0:
        return 

    # first start by sorting by start of the locus (rather than score)
    candidate_loci.sort(key=lambda locus: locus.chr_bounds.start)
    set_desired_expansion(candidate_loci, expand_params, protInfo)
    # handle start of first locus
    if (candidate_loci[0].expansion.left > 0):
        candidate_loci[0].chr_bounds.start = max(0, candidate_loci[0].chr_bounds.start - candidate_loci[0].expansion.left)
    # handle end of last locus (coord out of chromosome should be taken care during region extraction)
    if (candidate_loci[-1].expansion.right > 0):
        candidate_loci[-1].chr_bounds.end += candidate_loci[-1].expansion.right
    # now consider the space between loci to expand them
    for loc_id in range(1,len(candidate_loci)):
        (prev_locus, locus) = (candidate_loci[loc_id-1], candidate_loci[loc_id])
        available_nt_for_expansion = locus.chr_bounds.start - prev_locus.chr_bounds.end -1
        if available_nt_for_expansion >= prev_locus.expansion.right + locus.expansion.left:
            prev_locus.chr_bounds.end += prev_locus.expansion.right
            locus.chr_bounds.start -= locus.expansion.left
        elif prev_locus.expansion.right ==0:
            locus.chr_bounds.start -= available_nt_for_expansion
        elif locus.expansion.left == 0:
            prev_locus.chr_bounds.end += available_nt_for_expansion
        else: # split the available nt between the two loci; adjusting shrinking or proportionnally to their demands
            if( prev_locus.shrink_info[1]):
                prev_locus.chr_bounds.end += available_nt_for_expansion
            elif(locus.shrink_info[0]): # only the locus with the lowest score can be shrinked on this interval
                locus.chr_bounds.start -= available_nt_for_expansion
            else:           
                nt_for_prev = math.floor(available_nt_for_expansion *  prev_locus.expansion.right / (locus.expansion.left+prev_locus.expansion.right))
                prev_locus.chr_bounds.end += nt_for_prev
                locus.chr_bounds.start -= (available_nt_for_expansion - nt_for_prev)




def add_loci_from_mergeableHSPs(mergeableHSP: MergeableHSP, protInfos : dict, params: ParametersCandidateLoci, candidateLoci:list[CandidateLocus])-> None:
    if(mergeableHSP is None):
        return
    protInfo = None
    if(protInfos.keys().__contains__(mergeableHSP.prot_id)):
        protInfo =protInfos[mergeableHSP.prot_id]
    if(protInfo != None and params.skip_neighborhood_dist != None):
        dist=mergeableHSP.dist_to_protein(protInfos) 
        if(dist != None and dist < params.skip_neighborhood_dist):
            return
    
    candidate_loci = mergeableHSP.compute_candidate_loci_rec(params.loci_scoring)
    for candidate_locus in candidate_loci:
        candidate_locus.compute_score(protInfo, mergeableHSP.max_intron_len, params.loci_scoring)
        if(candidate_locus.score > params.loci_scoring.min_score and (candidate_locus.pc_similarity >= params.loci_scoring.min_similarity )) :#or candidate_locus.nhomol_prot>500)):
            candidateLoci.append(candidate_locus)

def init_mergeagleHSPs(hsp:HSP, protInfo:dict, default_intron_lg:int, params:ParametersCandidateLoci) -> MergeableHSP:
    if(protInfo.keys().__contains__(hsp.prot_id)):
        longest_allowed_intron = max (default_intron_lg, int(protInfo[hsp.prot_id].longest_intron * 1.1))
    else:
        longest_allowed_intron = default_intron_lg
    return MergeableHSP(hsp.prot_id, [hsp], hsp.locS_bounds.end, longest_allowed_intron )

def find_candidate_loci_from_file(gff_file:str, sorted_blast_file:Union[str, TextIO, StringIO], params=None, chr=None) -> dict:
    if(params is None):
        candideLociParams=ParametersCandidateLoci()
    else:
        candideLociParams=params    
    (protInfo, def_intron_lg) = gff_to_geneInfo(gff_file, candideLociParams.hsp_clustering.quantileForMaxIntronLength)
    if(not candideLociParams.hsp_clustering.useQuantile):
        def_intron_lg=candideLociParams.hsp_clustering.maxIntronLength
    print("finding candidate loci")
    
    if isinstance(sorted_blast_file, str) or isinstance(sorted_blast_file, Path):
        file_handle = open(sorted_blast_file, "r")
    else:
        # If it's a StringIO or TextIO object, use it directl
        file_handle = sorted_blast_file
        
    with file_handle:
        reader = csv.reader(file_handle, delimiter="\t") 
        candidate_loci_per_chr={}
        prev_chr : Optional[str] = None
        prev_strand : Optional[int] = None
        mergeableHSP: Optional[MergeableHSP]=None
        chr_candidate_loci=[]
        for row in reader:
            hsp = HSP(*row)
            # if its the first HSP initialiaze and continue
            if(prev_chr is None): 
                mergeableHSP=init_mergeagleHSPs(hsp, protInfo, def_intron_lg, candideLociParams)
                (prev_chr, prev_strand) = (hsp.chr_id, hsp.strand)
                continue
            # if the HSP is compatible with previous ones, accumulate it and continue
            if (prev_chr, prev_strand)== (hsp.chr_id, hsp.strand) and mergeableHSP.mergeHSP(hsp):
                continue
            # else infer candidate loci for accumulated HSP and initiate the next set of HSP
            add_loci_from_mergeableHSPs(mergeableHSP, protInfo, candideLociParams, chr_candidate_loci)
            mergeableHSP=init_mergeagleHSPs(hsp, protInfo, def_intron_lg, candideLociParams)
            # then check, if we are moving to a new chromosome, if so handle accumulated candidate loci for the previous chromosome
            if  prev_chr != hsp.chr_id : 
                candidate_loci_per_chr[prev_chr]=keep_best_non_overlaping_loci(chr_candidate_loci,candideLociParams.loci_scoring)
                if(candideLociParams.expansion is not None):
                    expands(candidate_loci_per_chr[prev_chr], candideLociParams.expansion, protInfo)
                chr_candidate_loci=[]
            (prev_chr, prev_strand)= (hsp.chr_id, hsp.strand)
            
        # handle last mergeableHSP and last chromosome 
        add_loci_from_mergeableHSPs(mergeableHSP, protInfo, candideLociParams, chr_candidate_loci)
        if(prev_chr is not None): 
            candidate_loci_per_chr[prev_chr]=keep_best_non_overlaping_loci(chr_candidate_loci,candideLociParams.loci_scoring)
            if(candideLociParams.expansion is not None):
                expands(candidate_loci_per_chr[prev_chr], candideLociParams.expansion, protInfo)
            
    return candidate_loci_per_chr

def find_candidate_loci(gff_file:str, blast_file:str, params=None, chr=None) -> dict:
    """ 
    Find candidate loci from a GFF file and a blast file. 
    Avoid intermediate sorted blast file by storing it in memory (mainly for testing).
    """
    memory_blast_sorted_file = StringIO()
    blast_to_sortedHSPs(blast_file, memory_blast_sorted_file, chr)
    # Set the cursor to the beginning of the StringIO object and use it as input for find_candidate_loci_from_file
    memory_blast_sorted_file.seek(0)
    return find_candidate_loci_from_file(gff_file, memory_blast_sorted_file, params, chr)
