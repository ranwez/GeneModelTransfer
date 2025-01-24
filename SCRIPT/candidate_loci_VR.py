
import argparse
import os

from CANDIDATE_LOCI.blast_utils import blast_to_HSPs
from CANDIDATE_LOCI.candidate_loci import ParametersExpansion, find_candidate_loci
from CANDIDATE_LOCI.gff_utils import gff_to_geneInfo


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-g", "--gff_file", type=str,
                        help="Exact path of gff file")
    parser.add_argument("-t", "--table", type=str,
                        help="Exact path of blast tabular output file")
    parser.add_argument("-o", "--output_gff", type=str, required=True,
                        help="Path to the output GFF file to store the results")
    parser.add_argument("-l", "--output_list", type=str, required=True,
                        help="Path to the output file of candidate loci /model prot association")
    parser.add_argument("--chr", type=str, default=None,
                        help="Chromosome identifier to filter the data (optional)")
    
    args = parser.parse_args()
    

    if args.gff_file is None:
        print("Please provide the path to the gff file use, -h for help")
        return
    if args.table is None:
        print("Please provide the path to the blast tabular output, use -h for help")
        return
    if args.output_gff is None:
        print("Please provide the path to the output gff file, use -h for help")
        return
    if args.output_list is None:
        print("Please provide the path to the output list of query/target, use -h for help")
        return
    
    # Check if file exists
    if not os.path.isfile(args.gff_file):
        raise FileNotFoundError(f"GFF file not found: {args.gff_file}")
    if not os.path.isfile(args.table):
        raise FileNotFoundError(f"BLAST table file not found: {args.table}")
    
    param_ext= ParametersExpansion(nb_aa_for_missing_part=10, nb_nt_default=300, nb_nt_when_missing_part=3000)
    candidateLoci = find_candidate_loci(args.gff_file, args.table, param_ext, args.chr)

    # Write results to the output GFF file
    with open(args.output_gff, "w") as out_file:
        for chr in candidateLoci:
            for locus in candidateLoci[chr]:
                out_file.write(locus.as_gff() + "\n")

    # Write results to the list query/target file
    with open(args.output_list, "w") as out_file:
        for chr in candidateLoci:
            for locus in candidateLoci[chr]:
                out_file.write(locus.as_query_target() + "\n")

if __name__ == "__main__":
    main()
