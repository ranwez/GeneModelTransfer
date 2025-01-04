import argparse
import math
from typing import Optional
from attrs import define
from CANDIDATE_LOCI.gff_utils import GeneInfo, gff_to_geneInfo
from CANDIDATE_LOCI.interlap import InterLap, Interval
from CANDIDATE_LOCI.blast_utils import HSP, HSP_chr, blast_to_HSPs, Bounds

# avoid importing numpy just for this
def argmax(x):
    return max(range(len(x)), key=lambda i: x[i])

@define(slots=True,frozen=True)
class HspCompatibleOverlap:
    compatible: bool
    max_overlap: int

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
    nhomol: int
    pc_similarity: Optional[float]
    score: Optional[float]

    @classmethod
    def from_hsp(cls, hsp:HSP) -> "CandidateLocus":
        return cls(
            chr_id = hsp.chr_id,
            strand = hsp.strand,
            prot_id = hsp.prot_id,
            prot_len = hsp.prot_len,
            chr_path = [hsp.locS_bounds],
            chr_bounds = [hsp.locS_bounds],
            prot_path = [hsp.prot_bounds],
            prot_bounds = [hsp.prot_bounds],
            nident = hsp.nident,
            nhomol = hsp.prot_bounds.length(),
            pc_similarity= None,
            score = None
        )
    @classmethod
    def from_hsp_path(cls, hsp_path:list[HSP], nident:int) -> "CandidateLocus":
        chr_path = [hsp.loc_bounds for hsp in hsp_path]
        prot_path = [hsp.prot_bounds for hsp in hsp_path]
        homolog_fraction = Interval()# end is exclude so add 1
        for hsp_prot in prot_path:
            homolog_fraction.add([(hsp_prot.start, (hsp_prot.end+1))])
        
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
            nhomol = homolog_fraction.coverage_VR(),
            pc_similarity= None,
            score =None
        ) 
    def compute_score(self, protInfo:GeneInfo, max_intron_len:int, len_penalty_percentage:float) -> float:
        # compute homology as the covered hsp length; score using similarity : 1 for id and 0.5 for non id
        similarity = (self.nhomol+self.nident)/2
        self.pc_similarity = similarity / self.prot_len
        # penalty for too large genomic region
        prot_genomic_len = protInfo.coding_region_length()
        # how many extra intron length do we have
        length_deviation=(self.chr_bounds.length() - prot_genomic_len) / max_intron_len
        length_penalty_in_percent = max(-1, len_penalty_percentage * (1 - max(1, math.exp(length_deviation))))
        # score favor high nident; covering most of the template protein; and having a somehow similar length
        self.score = similarity * self.pc_similarity * (1-length_penalty_in_percent)
        return 
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
                nident_ij = nident[i]+nident[j]-path_compatible_overlap.max_overlap
                # prefer solution with shorter path limit hsp included in hsp with the same scoring
                if path_compatible_overlap.compatible and nident_ij >= best_ident:
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
                sum_overlap += resi.max_overlap
        return HspCompatibleOverlap(True, sum_overlap)

def path_bound(path:list[Bounds]) -> Bounds:
    return Bounds(min(path[0].start, path[-1].start),
                 max( path[0].end, path[-1].end))


def keep_best_non_overlaping_loci(candidate_loci:list[CandidateLocus]) -> list[CandidateLocus]:
    # Sort candidate loci by score
    candidate_loci.sort(key=lambda x: x.score, reverse=True)
    # Filter candidate loci
    filtered_candidate_loci = []
    covered = InterLap()
    for locus in candidate_loci:
        if covered.__contains__((locus.chr_bounds.start, locus.chr_bounds.end)):
            continue
        else:
            covered.add((locus.chr_bounds.start, locus.chr_bounds.end))
            filtered_candidate_loci.append(locus)
    return filtered_candidate_loci

def find_candidate_loci_from_hsps(list_hspchr:list[HSP_chr], protInfo:dict, default_intron_lg:int) -> dict:
    list_hspchr.sort(key=lambda hspchr: (hspchr.chr_id, hspchr.strand))
    list_hspchr.append(HSP_chr("dummmyyChr", "", []))  # add a dummy HSP_chr to make sure the last chromosome is processed
    prev_chr : Optional[str] = None
    candidate_loci={}
    for hsp_chr in list_hspchr:
        print( "treating chromosome", hsp_chr.chr_id, "strand", hsp_chr.strand)
        # if the chromosome changes, process the candidate loci for the previous chromosome
        if prev_chr is None or prev_chr != hsp_chr.chr_id:
            if prev_chr is not None:
                print("start filtering")
                candidate_loci[prev_chr]=keep_best_non_overlaping_loci(chr_candidate_loci)
                print("end filtering")
                #return candidate_loci
                #print ("candidate loci", candidate_loci)
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
                    candidate_locus = mergeableHSP.compute_candidate_loci()
                    max_coverage = mergeableHSP.max_possible_coverage()
                    #print("computing candidate locus", mergeableHSP.prot_id, len(mergeableHSP.hsps), mergeableHSP.hsps, max_coverage)
                    if(max_coverage >= 0.5):
                        candidate_locus.compute_score(protInfo[mergeableHSP.prot_id], mergeableHSP.max_intron_len, 0.1)
                     #   print("candidate locus", candidate_locus)
                      #  print()
                        if(candidate_locus.score > 0 and candidate_locus.pc_similarity >= 0.5):
                            chr_candidate_loci.append(candidate_locus)
                if(hsp.prot_id != "dummmyProt"):
                    longest_allowed_intron = int(max (default_intron_lg, protInfo[hsp.prot_id].longest_intron) * 1.1)
                    mergeableHSP = MergeableHSP(hsp.prot_id, [hsp], hsp.locS_bounds.end, longest_allowed_intron )
        hsps.pop()  # remove the dummy HSP
    return candidate_loci

def find_candidate_loci(gff_file:str, blast_file:str) -> dict:
    print("parsing blast")
    list_hspchr = blast_to_HSPs(blast_file)
    print("parsing gff")
    (protInfo, def_intron_lg) = gff_to_geneInfo(gff_file, 0.5)
    print("finding candidate loci")
    candidateLoci = find_candidate_loci_from_hsps(list_hspchr, protInfo, 4000)
    return candidateLoci
