
import argparse

from CANDIDATE_LOCI.blast_utils import blast_to_HSPs
from CANDIDATE_LOCI.candidate_loci import find_candidate_loci
from CANDIDATE_LOCI.gff_utils import gff_to_geneInfo


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-g", "--gff_file", type=str,
                        help="Exact path of gff file")
    parser.add_argument("-t", "--table", type=str,
                        help="Exact path of blast tabular output file")

    args = parser.parse_args()
    if args.gff_file is None:
        print("Please provide the path to the gff file using -g or --gff_file")
        print("Use -h or --help for more information")
        return
    if args.table is None:
        print("Please provide the path to the blast tabulart output using -t or --table")
        print("Use -h or --help for more information")
        return
    candidateLoci = find_candidate_loci(args.gff_file, args.table)
    for chr in candidateLoci:
        print(f"Chromosome {chr}")
        for locus in candidateLoci[chr]:
            print(locus)
    #print(candidateLoci)
if __name__ == "__main__":
    main()
