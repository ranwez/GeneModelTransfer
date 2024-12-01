import polars as pl
import math
import argparse

# Helper function to compute max intron length


def max_intron_length(mRNA_rows):
    starts = mRNA_rows["start"].cast(int)
    ends = mRNA_rows["end"].cast(int)
    intron_lengths = (starts.shift(-1) - ends).drop_nulls()
    return intron_lengths.max() if len(intron_lengths) > 0 else 0


def max_intron_from_gff_v2(gff_file):
    # Load the GFF file
    gff = pl.read_csv(
        gff_file,
        separator='\t',
        has_header=False,
        comment_prefix='#'
    )

    # Add column names
    gff.columns = [
        "seqid", "source", "type", "start", "end",
        "score", "strand", "phase", "attributes"
    ]

    # Filter rows to keep only gene, mRNA, and CDS
    gff = gff.filter(gff["type"].is_in(["gene", "mRNA", "CDS"]))
    # Parse column 9 to add feature_id and parent_id columns
    gff = gff.with_columns([
        gff["attributes"].str.extract(r"ID=([^;]+)", 1).alias("feature_id"),
        gff["attributes"].str.extract(r"Parent=([^;]+)", 1).alias("parent_id"),
        (gff["end"] - gff["start"] + 1).alias("genomic_size")
    ])

# Populate temporary columns based on conditions
    gff = gff.with_columns([
        pl.when(gff["type"] == "gene").then(
            gff["feature_id"]).otherwise(None).alias("gene_id1"),
        pl.when(gff["type"] == "mRNA").then(
            gff["parent_id"]).otherwise(None).alias("gene_id2"),
        pl.when(gff["type"] == "mRNA").then(
            gff["feature_id"]).otherwise(None).alias("mRNA_id1"),
        pl.when(gff["type"] == "CDS").then(
            gff["parent_id"]).otherwise(None).alias("mRNA_id2")
    ])

# Create final gene_id and mRNA_id columns by combining the temporary columns
    gff = gff.with_columns([
        pl.coalesce([gff["gene_id1"], gff["gene_id2"]]).alias("gene_id"),
        pl.coalesce([gff["mRNA_id1"], gff["mRNA_id2"]]).alias("mRNA_id")
    ])

    # Drop the temporary columns
    gff = gff.drop(["gene_id1", "gene_id2", "mRNA_id1", "mRNA_id2"])

    # Group by mRNA_id and compute max intron length for each group
    mRNA_introns = []
    mRNA_genomic_lengths = []

    for group in gff.filter(gff["type"] == "CDS").group_by("mRNA_id"):
        mRNA_rows = group[1]
        intron_length = max_intron_length(mRNA_rows)
        mRNA_introns.append(
            {"mRNA_id": group[0][0], "maxIntronLength": intron_length})

    # Compute genomic length for each mRNA
    for group in gff.filter(gff["type"] == "CDS").group_by("mRNA_id"):
        mRNA_rows = group[1]
        genomic_length = (mRNA_rows["end"].max() -
                          mRNA_rows["start"].min() + 1)
        mRNA_genomic_lengths.append(
            {"mRNA_id": group[0][0], "genomicLength": genomic_length})

    # Convert introns and genomic lengths to DataFrames
    mRNA_introns_df = pl.DataFrame(mRNA_introns)
    mRNA_genomic_lengths_df = pl.DataFrame(mRNA_genomic_lengths)

    # Merge introns and genomic lengths into a single mRNA-level DataFrame
    mRNA_info_df = mRNA_introns_df.join(
        mRNA_genomic_lengths_df, on="mRNA_id", how="inner")

    # Merge max intron lengths and genomic lengths back with gene information
    gff_genes = gff.filter(gff["type"] == "mRNA").select(
        ["gene_id", "mRNA_id"])

    mRNA_info_df = mRNA_info_df.join(
        gff_genes, on="mRNA_id", how="inner")

    # Group by gene_id and compute max intron length and genomic length for each gene
    gene_info = []
    for group in mRNA_info_df.group_by("gene_id"):
        gene_rows = group[1]
        max_intron = gene_rows["maxIntronLength"].max()
        max_intron = max_intron if max_intron is not None else 0
        max_genomic_length = gene_rows["genomicLength"].max()
        gene_info.append(
            {"prot_id": group[0][0], "maxIntronSize": max_intron, "maxGenomicLength": max_genomic_length})

    # Build and return the result DataFrame
    return pl.DataFrame(gene_info)


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
    Add class_id, class_id_startS, and class_id_endS columns to the BLAST DataFrame,
    optimized by extracting relevant columns into Python lists.
    """
    # Step 1: Sort by chrS_id, prot_id and loc_startS
    sorted_df = blast_df.sort(["chrS_id", "prot_id", "loc_startS"])

    # Step 2: Extract relevant columns into Python lists
    chrS_ids = sorted_df["chrS_id"].to_list()
    prot_ids = sorted_df["prot_id"].to_list()
    loc_startS = sorted_df["loc_startS"].to_list()
    loc_endS = sorted_df["loc_endS"].to_list()
    max_intron_sizes = sorted_df["maxIntronSize"].to_list()

    # Initialize variables for dynamic class tracking
    class_ids = []
    class_chrS_list = []
    class_prot_list = []
    class_startS_list = []
    class_endS_list = []
    current_class_id = 1
    class_startS = None
    class_endS = -float("inf")
    class_prot_id = None
    class_chrS_id = None

    # Step 3: Process rows iteratively using lists
    for i in range(len(prot_ids)):
        row_prot_id = prot_ids[i]
        row_startS = loc_startS[i]
        row_endS = loc_endS[i]
        row_chrS_id = chrS_ids[i]
        max_intron_size = math.ceil(
            1.1*float(max(max_intron_sizes[i], threshold)))
        if class_startS is None:
            # Initialize the first class
            class_startS = row_startS
            class_endS = row_endS
            class_chrS_id = row_chrS_id
            class_prot_id = row_prot_id

        if (
            row_prot_id == class_prot_id and
            row_chrS_id == class_chrS_id and
            row_startS - class_endS <= max_intron_size
        ):
            class_endS = max(class_endS, row_endS)
        else:
            # Finalize the current class
            class_startS_list.append(class_startS)
            class_endS_list.append(class_endS)
            class_chrS_list.append(class_chrS_id)
            class_prot_list.append(class_prot_id)
            # Start a new class
            current_class_id += 1
            class_startS = row_startS
            class_endS = row_endS
            class_chrS_id = row_chrS_id
            class_prot_id = row_prot_id

        # Assign class_id to the row
        class_ids.append(current_class_id)

    # Step 4: Finalize the last class
    class_startS_list.append(class_startS)
    class_endS_list.append(class_endS)
    class_chrS_list.append(class_chrS_id)
    class_prot_list.append(class_prot_id)

    # Step 5: Add class_id, class_id_startS, and class_id_endS to the DataFrame
    class_info_df = pl.DataFrame({
        "class_id": list(range(1, len(class_startS_list) + 1)),
        "chrS_id": class_chrS_list,
        "class_id_startS": class_startS_list,
        "class_id_endS": class_endS_list,
        "prot_id": class_prot_list,
    })

    # Add class_id to the main DataFrame
    sorted_df = sorted_df.with_columns(pl.Series("class_id", class_ids))

    return sorted_df, class_info_df


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
    max_intron_lg = max_intron_from_gff_v2(args.gff_file)

    blast_full = pl.read_csv(args.table, separator='\t', has_header=False)
    # qseqid sseqid qlen length qstart qend sstart send nident pident gapopen evalue bitscore ### blast model prot to genome using tblastn
    blast_full.columns = ["prot_id", "chr_id", "prot_len", "length", "prot_start",
                          "prot_end", "loc_start", "loc_end", "nident", "pident", "gapopen", "evalue", "bitscore"]
    blast_full = add_simplified_coord(blast_full)

    blast_full2 = blast_full.join(
        max_intron_lg, on="prot_id", how="left")  # Perform the join

    rows_with_none = blast_full2.filter(pl.col("maxIntronSize").is_null())
    if (rows_with_none.shape[0] > 0):
        print("Some proteins are not found in the GFF :")
        # print first prot_ids that are not found in the GFF
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