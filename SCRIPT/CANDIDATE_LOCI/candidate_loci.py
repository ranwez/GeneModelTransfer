import argparse
import csv
import math
from typing import Optional
from attrs import define, field
from CANDIDATE_LOCI.gff_utils import GeneInfo, gff_to_geneInfo, gff_to_cdsInfo, CdsInfo
from CANDIDATE_LOCI.interlap import InterLap, Interval
from CANDIDATE_LOCI.blast_utils import HSP, HSP_chr, blast_to_HSPs, Bounds


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
@ define
class ParametersLociScoring:
    length_penalty_percentage: float = field(default=0.1)
    identityWeigthvsCoverage: int = field(default=3)
    min_similarity: float = field(default=0.4)
    min_score: float = field(default=0.0)
    nt_shrink: int = field(default=60)
    
@define
class ParametersCandidateLoci:
    expansion: ParametersExpansion = field(default=ParametersExpansion())
    hsp_clustering: ParametersHspClustering = field(default=ParametersHspClustering())
    loci_scoring: ParametersLociScoring = field(default=ParametersLociScoring())

@define(slots=True,frozen=True)
class HspCompatibleOverlap:
    compatible: bool
    max_id_overlap: int

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
    expansion: Optional[Bounds]= field(default=None)
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
            prot_genomic_len = protInfo.coding_region_length()
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
        gene_note = f"protId={self.prot_id};protLg={self.prot_len};prot_path={prot_path_info};score={self.score};nident={self.nident}"
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
class MergeableHSP:
    prot_id: str 
    hsps: list[HSP]
    max_coord: int
    max_intron_len: int
    hsp_overlap_cache: Optional[HspOverlapCacher]=None

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

    def compute_candidate_loci_rec(self)-> list[CandidateLocus]:
        candidateLocus=self.compute_candidate_loci()
        res=[candidateLocus]
        chr_bounds = candidateLocus.chr_bounds
        HPSr=[]
        HSPl=[]
        for hsp in self.hsps :
            if (hsp.loc_bounds.end < chr_bounds.start):
                HPSr.append(hsp)
            if(hsp.loc_bounds.start > chr_bounds.end):
                HSPl.append(hsp)
        for hsp_list in (HSPl, HPSr) :
            if (len(hsp_list) >0):
                sub_region= MergeableHSP(self.prot_id,hsp_list,self.max_coord,self.max_intron_len)
                res.extend(sub_region.compute_candidate_loci_rec())
        return res

    def compute_candidate_loci(self)->CandidateLocus:
        if len(self.hsps) == 1:
            return CandidateLocus.from_hsp(self.hsps[0]) 
        # Initialize scores and bounds, to store best resuts ending at each HSP
        self.hsps.sort(key=lambda hsp: hsp.locS_bounds.end)
        self.hsp_overlap_cache = HspOverlapCacher.from_hsp(self.hsps)

        nident = [self.hsps[i].nident for i in range(len(self.hsps))]
        best_prev_nident = [self.hsps[i].nident for i in range(len(self.hsps))]
        paths = [[i] for i in range(len(self.hsps))]
        #return CandidateLocus.from_hsp(self.hsps[0]) 
        # Compute scores and update bounds
        for j, hspj in enumerate(self.hsps):
            best_prev = -1
            best_ident = nident[j]
            # Only consider previous HSPs.
            # Reversed order for better shortcuts, the further away the hsp the most unlikely it is to get a good score
            # (since this likely lead to missing large part of the protein)
            for i in reversed(range(j)):
               # print("i", i, "j", j, nident[i], nident[j], " best_ident", best_ident, best_prev_nident[i])
                if (nident[i]+nident[j] < best_ident): # overlap only decrease score so no need to compute them
                    if (best_prev_nident[i]+nident[j] < best_ident):
                        #print("break")
                        break # nothing better earlier so stop the loop
                    else:
                        #print("continue")
                        continue # something worth checking earlier so continue the loop 
                path_compatible_overlap = self.max_path_overlap(paths[i], j)
                nident_ij = nident[i]+nident[j]-path_compatible_overlap.max_id_overlap
                # prefer solution with shorter path limit hsp included in hsp with the same scoring
                if path_compatible_overlap.compatible and nident_ij > best_ident:
                    best_ident = nident_ij
                    best_prev = i
            if best_prev != -1:
                paths[j] = paths[best_prev] + [j]
            nident[j] = best_ident
            best_prev_nident[j] = max(best_prev_nident[j-1], nident[j]) if j >= 1 else nident[j]
        # Find the index of the best scoring HSP
        best_idx = argmax(nident)
        candidate_locus= CandidateLocus.from_hsp_path([self.hsps[i] for i in paths[best_idx]], nident[best_idx])
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
        # get the list of relevant proteins present in candidate loci
        prots = set([locus.prot_id for locus in candidate_loci])
        # get the CDS coordinates of these proteins
        cds_infos = gff_to_cdsInfo(gff_file, prots)
    # set the expansion based on the genomic size of the corresponding missing part in the template protein
    for locus in candidate_loci:
        # set default expansion
        locus.expansion = Bounds(expand_params.nb_nt_default, expand_params.nb_nt_default)
        if locus.prot_bounds.start > expand_params.nb_aa_for_missing_part:
            locus.expansion.start = expand_params.nb_nt_when_missing_part
            if gff_file is not None:
                template_prot_info=protInfo[locus.prot_id]
                locus.expansion.start= max(locus.expansion.start, int(template_prot_info.longest_intron * 1.1))
        if locus.prot_len - locus.prot_bounds.end > expand_params.nb_aa_for_missing_part:
            locus.expansion.end = expand_params.nb_nt_when_missing_part
            if gff_file is not None:
                template_prot_info=protInfo[locus.prot_id]
                locus.expansion.start= max(locus.expansion.start, int(template_prot_info.longest_intron * 1.1) )
        if gff_file is not None:
            # now adjust the expansion based on the genomic size of the corresponding missing part in the template protein
            cds_info = cds_infos[locus.prot_id]
            if locus.prot_bounds.start > 0:
                (gene_start, locus_start, _gene_end)= cds_info.get_genomic_coord(locus.prot_bounds.start)
                missing_part = int((locus_start - gene_start +1) * 1.1)
                locus.expansion.start = max(locus.expansion.start, missing_part)
            if locus.prot_len - locus.prot_bounds.end > 0:
                coords = cds_info.get_genomic_coord(locus.prot_bounds.end)
                if coords is None:
                    raise RuntimeError(
                        f"No genomic coordinate found for protein end: {locus.prot_bounds.end} "
                        f"in locus:\n {locus}\n"
                        f"with cds_info:\n {cds_info}\n"
                    )
                _gene_start, locus_end, gene_end = coords
                #(_gene_start, locus_end, gene_end)= cds_info.get_genomic_coord(locus.prot_bounds.end)
                missing_part = int((gene_end - locus_end + 1) * 1.1) 
                locus.expansion.end = max(locus.expansion.end, missing_part)

#start=4452163, end=4461979
def expands(candidate_loci:list[CandidateLocus], expand_params:ParametersExpansion, protInfo:dict):
    if (len(candidate_loci))==0:
        return 

    # first start by sorting by start of the locus (rather than score)
    candidate_loci.sort(key=lambda locus: locus.chr_bounds.start)
    set_desired_expansion(candidate_loci, expand_params, protInfo)
    # handle start of first locus
    if (candidate_loci[0].expansion.start > 0):
        candidate_loci[0].chr_bounds.start = max(0, candidate_loci[0].chr_bounds.start - candidate_loci[0].expansion.start)
    # handle end of last locus (coord out of chromosome should be taken care during region extraction)
    if (candidate_loci[-1].expansion.end > 0):
        candidate_loci[-1].chr_bounds.end += candidate_loci[-1].expansion.end
    # now consider the space between loci to expand them
    for loc_id in range(1,len(candidate_loci)):
        (prev_locus, locus) = (candidate_loci[loc_id-1], candidate_loci[loc_id])
        available_nt_for_expansion = locus.chr_bounds.start - prev_locus.chr_bounds.end -1
        if available_nt_for_expansion >= prev_locus.expansion.end + locus.expansion.start:
            prev_locus.chr_bounds.end += prev_locus.expansion.end
            locus.chr_bounds.start -= locus.expansion.start
        elif prev_locus.expansion.end ==0:
            locus.chr_bounds.start -= available_nt_for_expansion
        elif locus.expansion.start == 0:
            prev_locus.chr_bounds.end += available_nt_for_expansion
        else: # split the available nt between the two loci; adjusting shrinking or proportionnally to their demands
            if( prev_locus.shrink_info[1]):
                prev_locus.chr_bounds.end += available_nt_for_expansion
            elif(locus.shrink_info[0]): # only the locus with the lowest score can be shrinked on this interval
                locus.chr_bounds.start -= available_nt_for_expansion
            else:           
                nt_for_prev = math.floor(available_nt_for_expansion *  prev_locus.expansion.end / (locus.expansion.start+prev_locus.expansion.end))
                prev_locus.chr_bounds.end += nt_for_prev
                locus.chr_bounds.start -= (available_nt_for_expansion - nt_for_prev)

def find_candidate_loci_from_hsps(list_hspchr:list[HSP_chr], protInfo:dict, default_intron_lg:int, params:ParametersCandidateLoci) -> dict:
    
    list_hspchr.sort(key=lambda hspchr: (hspchr.chr_id, hspchr.strand))
    list_hspchr.append(HSP_chr("dummmyyChr", "", []))  # add a dummy HSP_chr to make sure the last chromosome is processed
    prev_chr : Optional[str] = None
    candidate_loci_per_chr={}
    for hsp_chr in list_hspchr:
        print( "treating chromosome", hsp_chr.chr_id, "strand", hsp_chr.strand)
        # if the chromosome changes, process the candidate loci for the previous chromosome
        if prev_chr is None or prev_chr != hsp_chr.chr_id:
            if prev_chr is not None:
                print("start filtering")
                candidate_loci_per_chr[prev_chr]=keep_best_non_overlaping_loci(chr_candidate_loci,params.loci_scoring)
                if(params.expansion is not None):
                    print("start expanding")
                    expands(candidate_loci_per_chr[prev_chr], params.expansion, protInfo)
                print("end treating chromosome", prev_chr)
            chr_candidate_loci=[]
            prev_chr = hsp_chr.chr_id
        # now process the HSPs of the current chromosome / strand pair
        hsps=hsp_chr.HSP
        print("sorting HSPs")
        hsps.sort(key=lambda hsp: (hsp.prot_id, hsp.locS_bounds.start))
        hsps.append(HSP.build_dummy("dummmyProt"))  # add a dummy HSP to make sure the last HSP is processed
        # mergeable HSPs are HSPs of the same prot in the same region that can be merged
        print("start merging")
        mergeableHSP: Optional[MergeableHSP]=None
        for hsp in hsps:
            #print("treating HSP", hsp.prot_id, hsp.loc_startS, hsp.loc_endS)
            if mergeableHSP is None or not mergeableHSP.mergeHSP(hsp): 
                if mergeableHSP is not None:
                    candidate_loci = mergeableHSP.compute_candidate_loci_rec()
                    for candidate_locus in candidate_loci:
                        candidate_locus.compute_score(protInfo[mergeableHSP.prot_id], mergeableHSP.max_intron_len, params.loci_scoring)
                        if(candidate_locus.score > params.loci_scoring.min_score and (candidate_locus.pc_similarity >= params.loci_scoring.min_similarity )) :#or candidate_locus.nhomol_prot>500)):
                            chr_candidate_loci.append(candidate_locus)
                if(hsp.prot_id != "dummmyProt"):
                    if(protInfo.keys().__contains__(hsp.prot_id)):
                        longest_allowed_intron = max (default_intron_lg, int(protInfo[hsp.prot_id].longest_intron * 1.1))
                    else:
                        continue # skip HSPs with unknown protein for tests
                    mergeableHSP = MergeableHSP(hsp.prot_id, [hsp], hsp.locS_bounds.end, longest_allowed_intron )
        hsps.pop()  # remove the dummy HSP
    return candidate_loci_per_chr

def find_candidate_loci(gff_file:str, blast_file:str, params=None, chr=None) -> dict:
    if(params is None):
        candideLociParams=ParametersCandidateLoci()
    else:
        candideLociParams=params    
    print("parsing blast")
    list_hspchr = blast_to_HSPs(blast_file)
    print("parsing gff")
    (protInfo, def_intron_lg) = gff_to_geneInfo(gff_file, candideLociParams.hsp_clustering.quantileForMaxIntronLength)
    if(not candideLociParams.hsp_clustering.useQuantile):
        def_intron_lg=candideLociParams.hsp_clustering.maxIntronLength
    print("finding candidate loci")
    candidateLoci = find_candidate_loci_from_hsps(list_hspchr, protInfo, def_intron_lg, candideLociParams)
    return candidateLoci




#############################
def add_loci_from_mergeableHSPs(mergeableHSP: MergeableHSP, protInfos : dict, params: ParametersCandidateLoci, candidateLoci:[CandidateLocus])-> None:
    if(mergeableHSP is None):
        return
    candidate_loci = mergeableHSP.compute_candidate_loci_rec()
    protInfo = None
    if(protInfos.keys().__contains__(mergeableHSP.prot_id)):
        protInfo =protInfos[mergeableHSP.prot_id]
        
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

def find_candidate_loci_from_file(gff_file:str, sorted_blast_file:str, params=None, chr=None) -> dict:
    if(params is None):
        candideLociParams=ParametersCandidateLoci()
    else:
        candideLociParams=params    
    (protInfo, def_intron_lg) = gff_to_geneInfo(gff_file, candideLociParams.hsp_clustering.quantileForMaxIntronLength)
    if(not candideLociParams.hsp_clustering.useQuantile):
        def_intron_lg=candideLociParams.hsp_clustering.maxIntronLength
    print("finding candidate loci")
    
    with open(sorted_blast_file, "r") as file:
        reader = csv.reader(file, delimiter="\t")  # Use "," if it's a comma-separated file
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
                print("handle candidate loci :", prev_chr)
                candidate_loci_per_chr[prev_chr]=keep_best_non_overlaping_loci(chr_candidate_loci,candideLociParams.loci_scoring)
                if(candideLociParams.expansion is not None):
                    expands(candidate_loci_per_chr[prev_chr], candideLociParams.expansion, protInfo)
                chr_candidate_loci=[]
            (prev_chr, prev_strand)= (hsp.chr_id, hsp.strand)
            
        # handle last mergeableHSP and last chromosome 
        add_loci_from_mergeableHSPs(mergeableHSP, protInfo, candideLociParams, chr_candidate_loci)
        if(prev_chr is not None): 
            print("handle candidate loci :", prev_chr)
            candidate_loci_per_chr[prev_chr]=keep_best_non_overlaping_loci(chr_candidate_loci,candideLociParams.loci_scoring)
            if(candideLociParams.expansion is not None):
                expands(candidate_loci_per_chr[prev_chr], candideLociParams.expansion, protInfo)
    return candidate_loci_per_chr

