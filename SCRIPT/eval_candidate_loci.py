from typing import Dict, List, Tuple
from pathlib import Path
import argparse
import sys
from CANDIDATE_LOCI.gff_utils import gff_to_cdsInfo, GeneInfo, CdsInfo
from CANDIDATE_LOCI.bounds import Bounds
from CANDIDATE_LOCI.intervals_utils import OrderedIntervals
from attrs import define, field

@define 
class OverlapScore:
    loci_overlap : float
    CDS_overlap : float

@define 
class OverlapMatch:
    candidate_id: str
    ref_id:str
    score :OverlapScore 



def evaluate_candidate_loci(candidate_gff: str, reference_gff: str) -> Dict[str, List[OverlapMatch]]:
    """
    Evaluate candidate loci by comparing them to reference annotations.
    
    Args:
        candidate_gff: Path to GFF file containing candidate loci
        reference_gff: Path to GFF file containing reference annotations
        
    Returns:
        Dictionary mapping chromosome to list of (candidate_gene, ref_gene, overlap_score) tuples
    """
    # Parse both GFF files to get CDS information
    candidate_genes = gff_to_cdsInfo(candidate_gff,None)
    reference_genes = gff_to_cdsInfo(reference_gff, None)
    
    # Group genes by chromosome
    candidate_by_chr = _group_genes_by_chromosome(candidate_genes)
    reference_by_chr = _group_genes_by_chromosome(reference_genes)
    
    results = {}
    
    # Process each chromosome
    for chr_id in candidate_by_chr.keys():
        if chr_id not in reference_by_chr:
            continue
            
        chr_results = _evaluate_chromosome(
            candidate_by_chr[chr_id], 
            reference_by_chr[chr_id]
        )
        
        if chr_results:
            results[chr_id] = chr_results
    
    return results

def _group_genes_by_chromosome(genes: Dict[str, CdsInfo]) -> Dict[str, List[Tuple[str, CdsInfo]]]:
    """Group genes by chromosome ID."""
    by_chr = {}
    for gene_id, gene_info in genes.items():
        chr_id = gene_info.chr_id
        if chr_id not in by_chr:
            by_chr[chr_id] = []
        by_chr[chr_id].append((gene_id, gene_info))
    
    # Sort genes by start position for efficient merging
    for chr_id in by_chr:
        by_chr[chr_id].sort(key=lambda x: x[1].coding_region.start)
    
    return by_chr

def _evaluate_chromosome(candidate_loci: List[Tuple[str, CdsInfo]], 
                        reference_genes: List[Tuple[str, CdsInfo]]) -> List[OverlapMatch]:
    """
    Evaluate candidate genes against reference genes on a single chromosome.
    Uses merge-sort inspired algorithm for efficient overlap detection.
    """
    results = []
    ref_idx = 0
    
    for cand_id, cand_locus in candidate_loci:
        # Find all reference genes that could overlap with current candidate
        overlapping_refs = []
        
        # Start from current reference index and move forward
        temp_idx = ref_idx
        while temp_idx < len(reference_genes):
            ref_id, ref_gene = reference_genes[temp_idx]
            
            # If reference gene starts after candidate ends, no more overlaps possible
            if ref_gene.coding_region.start > cand_locus.gene_bounds.end:
                break
                
            # If reference gene ends before candidate starts, advance reference index
            if ref_gene.coding_region.end < cand_locus.gene_bounds.start:
                if temp_idx == ref_idx:
                    ref_idx += 1
                temp_idx += 1
                continue
                
            # Genes overlap, add to candidates for detailed evaluation
            overlapping_refs.append((ref_id, ref_gene))
            temp_idx += 1
        
        # Calculate overlap scores for all overlapping reference genes
        best_match = None
        best_score = OverlapScore(0,0)
        
        for ref_id, ref_gene in overlapping_refs:
            overlap_score = _calculate_cds_overlap_score(cand_locus, ref_gene)
            if overlap_score.loci_overlap > best_score.loci_overlap:
                best_score = overlap_score
                best_match = OverlapMatch(cand_id, ref_id, overlap_score)
        
        if best_match :
            results.append(best_match)
    
    return results

def _geneInfo_to_CDS_interval(gene_info: CdsInfo) -> OrderedIntervals:
    """
    Create an OrderedIntervals object from gene CDS bounds.
    
    Args:
        gene_info: GeneInfo object containing CDS bounds
        
    Returns:
        OrderedIntervals object containing all CDS intervals
    """
    intervals = []
    cds_bounds:list[Bounds]=gene_info.cds_bounds
    cds_bounds.sort(key=lambda x: x.start)
    if(len(cds_bounds)>0):
        (current_range_start, current_range_end) = (cds_bounds[0].start, cds_bounds[0].end)
        for cds_bound in gene_info.cds_bounds:
            if cds_bound.start > current_range_end:
                intervals.extend([current_range_start, current_range_end])
                (current_range_start, current_range_end) =(cds_bound.start, cds_bound.end)
            current_range_end = max(current_range_end, cds_bound.end)
        intervals.extend([current_range_start, current_range_end])
    return OrderedIntervals(intervals, True)

def _calculate_cds_overlap_score(candidate_locus: CdsInfo, reference_gene: CdsInfo) -> OverlapScore:
    """
    Calculate overlap score between candidate locus and CDS regions of the reference gene .
    Score = overlap_size / sum_of_reference_CDS_length
    """
    # Create OrderedIntervals objects for candidaate region and CDS regions of the reference gene
    candidate_loci_interval = OrderedIntervals([candidate_locus.gene_bounds.start, candidate_locus.gene_bounds.end], True)
    candidate_CDS_interval = _geneInfo_to_CDS_interval(candidate_locus)
    reference_interval = _geneInfo_to_CDS_interval(reference_gene)
    
    ref_total_length = reference_interval.total_length()
    if ref_total_length == 0:
        return OverlapScore(0,0)
    
    # Calculate intersection
    locus_intersection = candidate_loci_interval.intersection(reference_interval)
    locus_overlap_length = locus_intersection.total_length()

    CDS_intersection = candidate_CDS_interval.intersection(reference_interval)
    CDS_overlap_length = CDS_intersection.total_length()
    return OverlapScore(locus_overlap_length / ref_total_length, CDS_overlap_length/ref_total_length) 

def print_evaluation_results(results: Dict[str, List[OverlapMatch]], 
                           reference_genes: Dict[str, any],
                           min_score: float = 0.1):
    """
    Print evaluation results in TSV format.
    
    Args:
        results: Dictionary of evaluation results
        reference_genes: Dictionary of reference genes to get CDS bounds
        min_score: Minimum overlap score to display
    """
    #print("chr_id\tref_gene_id\tref_gene_CDS_start\tref_gene_CDS_end\tcand_id\toverlap_score")
    
    # Organize reference_genes by chromosome
    reference_genes_by_chr = {}
    for ref_id, ref_gene in reference_genes.items():
        chr_id = ref_gene.chr_id
        if chr_id not in reference_genes_by_chr:
            reference_genes_by_chr[chr_id] = {}
        reference_genes_by_chr[chr_id][ref_id] = ref_gene
    
    total_matches = 0
    total_missed = 0
    
    # Process each chromosome
    for chr_id in reference_genes_by_chr.keys():
        matched_genes = set()
        
        # Print matches for this chromosome
        if chr_id in results:
            significant_matches = [r for r in results[chr_id] if r.score.loci_overlap >= min_score]
            for match in significant_matches:
                ref_gene = reference_genes[match.ref_id]    
                print(f"{chr_id}\t{match.ref_id}\t{ref_gene.coding_region.start}\t{ref_gene.coding_region.end}\t{match.candidate_id}\t{match.score.loci_overlap:.3f}\t{match.score.CDS_overlap:.3f}")
                matched_genes.add(match.ref_id)
                total_matches += 1
        
        # Print missed reference genes for this chromosome
        nb_missed_genes = 0
        for ref_id, ref_gene in reference_genes_by_chr[chr_id].items():
            if ref_id not in matched_genes:
                nb_missed_genes += 1
                print(f"{chr_id}\t{ref_id}\t{ref_gene.coding_region.start}\t{ref_gene.coding_region.end}\tNA\t0.0\t0.0")
        
        total_missed += nb_missed_genes
        
        # Print summary to stderr
        print(f"# {chr_id}: matched: {len(matched_genes)}, missed: {nb_missed_genes}", file=sys.stderr)
    
    # Print global summary to stderr
    print(f"# Total: matched: {total_matches}, missed: {total_missed}", file=sys.stderr)

# Example usage
def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Evaluate candidate loci by comparing them to reference annotations"
    )
    
    parser.add_argument(
        "candidate_gff",
        help="Path to GFF file containing candidate loci"
    )
    
    parser.add_argument(
        "reference_gff", 
        help="Path to GFF file containing reference annotations"
    )
    
    parser.add_argument(
        "-m", "--min-score",
        type=float,
        default=0.1,
        help="Minimum overlap score to display (default: 0.1)"
    )
    
    parser.add_argument(
        "-o", "--output",
        type=str,
        help="Output file to write results (if not specified, prints to stdout)"
    )
    
    args = parser.parse_args()
    
    # Check if input files exist
    if not Path(args.candidate_gff).exists():
        print(f"Error: Candidate GFF file '{args.candidate_gff}' not found")
        return 1
        
    if not Path(args.reference_gff).exists():
        print(f"Error: Reference GFF file '{args.reference_gff}' not found")
        return 1
    
    # Evaluate candidate loci
    print("Evaluating candidate loci...", file=sys.stderr)
    results = evaluate_candidate_loci(args.candidate_gff, args.reference_gff)
    
    # Get reference genes for CDS bounds information
    from CANDIDATE_LOCI.gff_utils import parse_gff
    reference_df = parse_gff(args.reference_gff)
    reference_gene_ids = reference_df.select("gene").unique().to_series().to_list()
    reference_genes = gff_to_cdsInfo(args.reference_gff, reference_gene_ids)
    
    # Output results
    if args.output:
        # Redirect output to file
        with open(args.output, 'w') as f:
            old_stdout = sys.stdout
            sys.stdout = f
            print_evaluation_results(results, reference_genes, args.min_score)
            sys.stdout = old_stdout
        print(f"Results written to {args.output}", file=sys.stderr)
    else:
        # Print to stdout
        print_evaluation_results(results, reference_genes, args.min_score)
    
    return 0

if __name__ == "__main__":
    main()