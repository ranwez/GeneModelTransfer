import sys
import argparse
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from gff_utils import gff_to_geneInfo


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-g",
        "--gff",
        metavar="input_ut8.gff",
        type=str,
        required=True,
        help="Path to the input GFF file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="CDS.bed",
        type=str,
        required=True,
        help="Name of the output bed.",
    )

    args = parser.parse_args()

    (gene_info, intron_lg) = gff_to_geneInfo(args.gff, 0.5)
    # print(gene_info)
    with open(args.output, "w") as out_file:
        for gene in gene_info.values():
            out_file.write(
                f"{gene.chr_id}\t{gene.coding_start}\t{gene.coding_end}\t{gene.gene_id}\n"
            )


if __name__ == "__main__":
    main()
