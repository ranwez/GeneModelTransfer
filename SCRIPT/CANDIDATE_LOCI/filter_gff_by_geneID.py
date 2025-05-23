import sys
import argparse
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from CANDIDATE_LOCI.gff_utils import filter_feature_by_geneID


def main():
    parser = argparse.ArgumentParser(
        description="Filter a GFF based on a list of gene IDs. Child features of the genes are kept."
    )
    parser.add_argument(
        "-g",
        "--gff",
        metavar="input_ut8.gff",
        type=str,
        required=True,
        help="Path to the input GFF file.",
    )
    parser.add_argument(
        "-l",
        "--id_list",
        metavar="input_ut8.gff",
        type=str,
        required=True,
        help="File listing gene IDs to extract from the gff (one per line).",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="output_filtered.gff",
        type=str,
        required=True,
        help="Name of the output filtered GFF file.",
    )

    args = parser.parse_args()

    filter_feature_by_geneID(args.gff, args.id_list, args.output)


if __name__ == "__main__":
    main()
