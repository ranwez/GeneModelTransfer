import sys
import os
import argparse
import csv
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from CANDIDATE_LOCI.gff_utils import parse_gff, sort_gff


def sort_and_write_gff(gff_input, output):
    gff = parse_gff(gff_input)
    gff = sort_gff(gff)
    gff = gff.drop(["mRNA", "gene", "ID", "ParentID"])
    gff.write_csv(output, separator="\t", include_header=False)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-g", "--gff", type=str, help="Exact path of gff file")
    parser.add_argument("-o", "--output", type=str, help="Output file name")
    args = parser.parse_args()

    sort_and_write_gff(args.gff, args.output)


if __name__ == "__main__":
    main()
