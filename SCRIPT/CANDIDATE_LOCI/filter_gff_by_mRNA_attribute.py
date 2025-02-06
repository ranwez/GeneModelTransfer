import argparse
from gff_utils import filter_mRNA_by_attribute


def main():
    parser = argparse.ArgumentParser(description="Filter a GFF based on the value of a mRNA attribute. Child features of the selected mRNA are kept. Gene features are kept if at least one of their child mRNA is kept.")
    parser.add_argument("-g", "--gff", metavar='input_ut8.gff', type=str, required=True, help="Path to the input GFF file.")
    parser.add_argument("-a", "--attribute", metavar='representative', type=str, required=True, help="Name of the mRNA attribute to use to filter the GFF file.")
    parser.add_argument("-v", "--value", metavar='True', type=str, required=True, help="Value of the mRNA attribute to keep the feature and its child features.")
    parser.add_argument("-o", "--output", metavar='output_filtered.gff', type=str, required=True, help="Name of the output filtered GFF file.")
    parser.add_argument("-u", "--output_unresolved", metavar='svevo1_HC_without_repr.gff', type=str, required=False, help="Name of the output GFF file with all genes that had no mRNA passing the attribute filter.")


    args = parser.parse_args()

    filter_mRNA_by_attribute(args.gff, args.attribute, args.value, args.output, args.output_unresolved)

if __name__ == "__main__":
    main()
