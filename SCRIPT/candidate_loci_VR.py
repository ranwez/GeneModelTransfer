import polars as pl
import math
import argparse

# Helper function to compute max intron length


def max_intron_length(mRNA_rows):
    starts = mRNA_rows["start"].cast(int)
    ends = mRNA_rows["end"].cast(int)
    intron_lengths = (starts.shift(-1) - ends).drop_nulls()
    return intron_lengths.max() if len(intron_lengths) > 0 else 0


def max_intron_from_gff_v3(gff_file):
    """
    Process a GFF file to compute maximum intron length and genomic length for genes and mRNAs.
    """
    # Load the GFF file without headers and assign column names
    gff = pl.read_csv(
        gff_file, separator='\t', has_header=False, comment_prefix='#', new_columns=[
            "seqid", "source", "type", "start", "end",
            "score", "strand", "phase", "attributes"]
    )

    gff = gff.filter(gff["type"].is_in(["gene", "mRNA", "CDS"]))

    # Extract `type` and `attributes` into Python lists
    types = gff["type"].to_list()
    attributes = gff["attributes"].to_list()

    # Parse gene_id and mRNA_id from attributes based on type
    gene_ids, mRNA_ids = [], []

    for t, attr in zip(types, attributes):
        if t == "gene":
            gene_ids.append(attr.split("ID=")[1].split(";")[
                            0] if "ID=" in attr else None)
            mRNA_ids.append(None)
        elif t == "mRNA":
            gene_ids.append(attr.split("Parent=")[1].split(";")[
                            0] if "Parent=" in attr else None)
            mRNA_ids.append(attr.split("ID=")[1].split(";")[
                            0] if "ID=" in attr else None)
        elif t == "CDS":
            gene_ids.append(None)
            mRNA_ids.append(attr.split("Parent=")[1].split(";")[
                            0] if "Parent=" in attr else None)

    # Add gene_id and mRNA_id as new columns
    gff = gff.with_columns([
        pl.Series("prot_id", gene_ids),
        pl.Series("mRNA_id", mRNA_ids)
    ])

    # Compute intron length and genomic length for each mRNA
    cds = gff.filter(gff["type"] == "CDS")
    mRNA_info_df = cds.group_by("mRNA_id").agg([
        (pl.col("start").shift(-1) - pl.col("end"))
        .filter(pl.col("mRNA_id").is_not_null()).max().alias("maxIntronLength"),
        (pl.col("end").max() - pl.col("start").min() + 1).alias("genomicLength")
    ])

    # Add gene-level information
    mRNA_info_df = mRNA_info_df.join(
        gff.filter(gff["type"] == "mRNA").select(["prot_id", "mRNA_id"]),
        on="mRNA_id",
        how="inner"
    )

    # Compute gene-level metrics
    gene_info_df = mRNA_info_df.group_by("prot_id").agg([
        pl.col("maxIntronLength").max().fill_null(0).alias("maxIntronSize"),
        pl.col("genomicLength").max().alias("maxGenomicLength")
    ])

    return gene_info_df


def combined_score(hsp_i, hsp_j, score_i):
    # Define indices for relevant columns
    prot_start = 0
    prot_end = 1
    loc_startS = 2
    loc_endS = 3
    nident = 4
    # Check consistency: hsp_j must be after hsp_i in loc_startS
    if hsp_i[loc_startS] >= hsp_j[loc_startS]:
        return 0

    # Compute overlaps
    loc_overlap = max(0, hsp_i[loc_endS] - hsp_j[loc_startS])
    prot_overlap = max(0, hsp_i[prot_end] - hsp_j[prot_start])

    # Calculate the maximum number of identities in the overlap
    max_id_overlap = min(
        max(loc_overlap, prot_overlap),  # The overlap size
        hsp_j[nident],                   # Identities in hsp_j
        hsp_i[nident]                    # Identities in hsp_i
    )
    # Compute the combined pessimistic score
    combined_pessimistic_score = score_i + hsp_j[nident] - max_id_overlap
    # print(loc_overlap, prot_overlap, max_id_overlap, combined_pessimistic_score)
    return combined_pessimistic_score


def compute_class_score(class_rows):
    # Sort rows by prot_start
    class_rows = class_rows.sort("prot_start")
    # Extract relevant columns as a list of rows
    hsp = [[row[0], row[1], row[2], row[3], row[4], ]
           for row in zip(
        class_rows["prot_start"].to_list(),
        class_rows["prot_end"].to_list(),
        class_rows["loc_startS"].to_list(),
        class_rows["loc_endS"].to_list(),
        class_rows["nident"].to_list(),
    )
    ]
    nident = 4  # indices in hsp
    # Initialize score array
    score = [hsp[i][nident] for i in range(len(hsp))]
    # Compute scores
    for j in range(len(hsp)):
        for i in range(j):  # Only consider previous HSPs
            score_ij = combined_score(hsp[i], hsp[j], score[i])
            if score_ij > score[j]:
                score[j] = score_ij
            # print(hsp[i], hsp[j], score[j])
    return max(score)


def add_classification_with_lists(blast_df, threshold):
    """
    Add class_id, class_id_startS, and class_id_endS columns to the BLAST DataFrame.
    """
    # Sort and extract necessary columns
    chrS_ids, prot_ids, loc_startS, loc_endS, max_intron_sizes, prot_len, strand, chr_id = (
        blast_df.sort(["chrS_id", "prot_id", "loc_startS"])
        .select(["chrS_id", "prot_id", "loc_startS", "loc_endS", "maxIntronSize", "prot_len", "strand", "chr_id"]).to_numpy().T
    )

    # Initialize variables for classification
    hsp_class = []
    class_chrS, class_protId, class_startS, class_endS = [], [], [], []
    class_prot_lg, class_strand, class_chr = [], [], []
    # Iterate through rows and assign class_ids
    for i, (chrS_id, prot_id, startS, endS, max_intron, prot_lg, strand, chr) in enumerate(
        zip(chrS_ids, prot_ids, loc_startS, loc_endS,
            max_intron_sizes, prot_len, strand, chr_id)
    ):
        max_intron = math.ceil(1.1 * max(max_intron, threshold))

        # Start a new class if necessary
        if not class_startS or not (
            prot_id == class_protId[-1]
            and chrS_id == class_chrS[-1]
            and startS - class_endS[-1] <= max_intron
        ):
            # Start a new class
            class_chrS.append(chrS_id)
            class_protId.append(prot_id)
            class_startS.append(startS)
            class_endS.append(endS)
            class_prot_lg.append(prot_lg)
            class_strand.append(strand)
            class_chr.append(chr)
        else:
            # Update the current class's end position
            class_endS[-1] = max(class_endS[-1], endS)

        # Assign class_id to the row
        hsp_class.append(len(class_startS))

    # Create class_info DataFrame
    class_info_df = pl.DataFrame({
        "class_id": list(range(1, len(class_startS) + 1)),
        "chrS_id": class_chrS,
        "prot_id": class_protId,
        "startS": class_startS,
        "endS": class_endS,
        "prot_len": class_prot_lg,
        "strand": class_strand,
        "chr_id": class_chr,
    })

    # Add class_id to the main DataFrame
    return blast_df.with_columns(pl.Series("class_id", hsp_class)), class_info_df


def add_simplified_coord(blast_full):
    """
    Add extra columns to handle hits on the reverse complement as simply as possible.

    This function performs the following steps:
    1. To avoid mixing +/- hits that should be merged, they are assigned different
       chromosome identifiers in the `chrS_id` column:
       - Hits on the '+' strand retain their original `chr_id`.
       - Hits on the '-' strand are given a modified `chr_id` (e.g., "chr_id_revcomp").

    2. The loc_start and loc_end genomic positions of '-' strand hits are redefined:
       - `loc_startS` and `endS` are calculated relative to the chromosome's end,
         including a margin, for '-' strand hits.
       - For '+' strand hits, `loc_startS` and `loc_endS` remain unchanged.

    By organizing hits this way:
    - Two HSPs (hspi and hspj) can be merged if they are nearby on the same `chrS_id`
      and their `loc_startS` values are ordered according to their `prot_start`.

    Args:
        blast_full (polars.DataFrame): The input DataFrame containing blast results.

    Returns:
        polars.DataFrame: The modified DataFrame with additional columns for simplified
        handling of hits on the reverse complement.
    """
    # Step 1: Compute chr_max_coord
    chr_max_coord = (
        blast_full
        .with_columns(
            pl.max_horizontal("loc_start", "loc_end").alias("loc_max")
        )
        .group_by("chr_id")
        .agg(
            pl.max("loc_max").alias("max_loc_coord")
        )
    )

    # Step 2: Add strand column and chrS column
    blast_full = blast_full.with_columns(
        pl.when(pl.col("loc_start") < pl.col("loc_end"))
        .then(1)
        .otherwise(-1)
        .alias("strand")
    )
    blast_full = blast_full.with_columns(
        pl.when(pl.col("strand") == 1)
        .then(pl.col("chr_id").cast(str))
        .otherwise(pl.col("chr_id").cast(str) + "_revcomp")
        .alias("chrS_id")
    )

    # Step 3: Join chr_max_coord back to the main DataFrame
    blast_full = blast_full.join(chr_max_coord, on="chr_id")

    # Step 4: Add loc_startS and loc_endS based on strand
    blast_full = blast_full.with_columns(
        # Compute loc_startS
        pl.when(pl.col("strand") == 1)
        .then(pl.col("loc_start"))
        .otherwise(
            pl.col("max_loc_coord") + 1000000 - pl.col("loc_end")
        )
        .alias("loc_startS"),

        # Compute loc_endS
        pl.when(pl.col("strand") == 1)
        .then(pl.col("loc_end"))
        .otherwise(
            pl.col("max_loc_coord") + 1000000 - pl.col("loc_start")
        )
        .alias("loc_endS")
    )
    # Drop the helper column
    return blast_full.drop(["max_loc_coord"])


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-g", "--gff_file", type=str,
                        help="Exact path of gff file")
    parser.add_argument("-t", "--table", type=str,
                        help="Exact path of alignment res table")

    args = parser.parse_args()
    # debug
    # args.table = "/Users/ranwez/Desktop/TEST_REGION/blast_refProt_ENSG00000169598_full_pos.tsv"
    max_intron_lg = max_intron_from_gff_v3(args.gff_file)

    # qseqid sseqid qlen length qstart qend sstart send nident pident gapopen evalue bitscore ### blast model prot to genome using tblastn
    blast_full = pl.read_csv(args.table, separator='\t', has_header=False, new_columns=["prot_id", "chr_id", "prot_len", "length", "prot_start",
                                                                                        "prot_end", "loc_start", "loc_end", "nident", "pident",
                                                                                        "gapopen", "evalue", "bitscore"])
    blast_full = add_simplified_coord(blast_full)
   # perform the join separately else some rows are lost ???
    blast_full2 = blast_full.join(
        max_intron_lg, on="prot_id", how="left")

    rows_with_none = blast_full2.filter(pl.col("maxIntronSize").is_null())
    if (rows_with_none.shape[0] > 0):
        print("Some proteins are not found in the GFF :")
        print(rows_with_none["prot_id"].to_list()[
              1:min(10, rows_with_none.shape[0])])
        exit(1)
    blast_class, class_info = add_classification_with_lists(blast_full2, 4000)
    print(blast_class.head())
    # Group by class_id and compute scores
    scores = []
    for group in blast_class.group_by("class_id"):
        class_rows = group[1]  # Extract the DataFrame of the current group
        class_score = compute_class_score(class_rows)
        scores.append({"class_id": group[0][0], "class_score": class_score})
    # print(scores)
    scores_df = pl.DataFrame(scores).join(
        class_info, on="class_id", how="left")

    print(scores_df.head())


if __name__ == "__main__":
    main()
