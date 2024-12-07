import polars as pl
import math
import argparse
from VR.interlap import InterLap

# avoid importing numpy just for this


def argmax(x):
    return max(range(len(x)), key=lambda i: x[i])

# Helper function to compute max intron length


def max_intron_length(mRNA_rows):
    starts = mRNA_rows["start"].cast(int)
    ends = mRNA_rows["end"].cast(int)
    intron_lengths = (starts.shift(-1) - ends).drop_nulls()
    return intron_lengths.max() if len(intron_lengths) > 0 else 0


def get_prot_info(gff_file):
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
    prot_info_df = mRNA_info_df.group_by("prot_id").agg([
        pl.col("maxIntronLength").max().fill_null(0).alias("maxIntronSize"),
        pl.col("genomicLength").max().alias("maxGenomicLength")
    ])

    return prot_info_df


def overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def max_id_overlap(hsp_i, hsp_j):
    # Define indices for relevant columns
    prot_start = 0
    prot_end = 1
    loc_startS = 2
    loc_endS = 3
    interval_prot_i = [hsp_i[prot_start], hsp_i[prot_end]]
    interval_prot_j = [hsp_j[prot_start], hsp_j[prot_end]]
    interval_loc_i = [hsp_i[loc_startS], hsp_i[loc_endS]]
    interval_loc_j = [hsp_j[loc_startS], hsp_j[loc_endS]]
    nident = 4

    # Compute overlaps
    loc_overlap = math.ceil(overlap(interval_loc_i, interval_loc_j)/3)
    prot_overlap = overlap(interval_prot_i, interval_prot_j)

    # Calculate the maximum number of identities in the overlap
    max_id_overlap = min(
        max(loc_overlap, prot_overlap),  # max hsp overlap with locus and prot
        hsp_j[nident],                   # Identities in hsp_j
        hsp_i[nident]                    # Identities in hsp_i
    )
    return max_id_overlap


def combined_nident(hsp_i, hsp_j, score_i):
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
    loc_overlap = max(0, math.ceil((hsp_i[loc_endS] - hsp_j[loc_startS]) / 3))
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


def max_path_overlap(prev_path, hsp, hspj):
    # Define indices for relevant columns
    prot_start = 0
    loc_startS = 2
    nident = 4

    max_total_overlap = 0
    for hsp_prev_id in prev_path:
        hspi = hsp[hsp_prev_id]
        # hsp start should be in the same order to be compatible
        if (hspj[prot_start]-hspi[prot_start]) * (hspj[loc_startS]-hspi[loc_startS]) < 0:
            return (False, 0)
        max_total_overlap += max_id_overlap(hspi, hspj)
    return(True, max_total_overlap)


def compute_class_score(class_rows, genomic_length_penalty, pc_lenth_penalty):
    class_rows = class_rows.sort("prot_end")

    # Extract relevant columns as a list of rows
    hsp = [[row[0], row[1], row[2], row[3], row[4]]
           for row in zip(
               class_rows["prot_start"].to_list(),
               class_rows["prot_end"].to_list(),
               class_rows["loc_startS"].to_list(),
               class_rows["loc_endS"].to_list(),
               class_rows["nident"].to_list(),
    )]

    nident_idx = 4  # Index of the nident column in hsp
    start_idx, end_idx = 2, 3  # Indices for loc_startS and loc_endS

    # Initialize scores and bounds
    nident = [hsp[i][nident_idx] for i in range(len(hsp))]
    best_prev_nident = [hsp[i][nident_idx] for i in range(len(hsp))]
    bounds = [[hsp[i][start_idx], hsp[i][end_idx]] for i in range(len(hsp))]
    paths = [[i] for i in range(len(hsp))]

    # Compute scores and update bounds
    for j in range(len(hsp)):
        best_prev = -1
        best_ident = nident[j]
        # Only consider previous HSPs. 
        # Reversed for better shortcuts, the further away the hsp the most unlikely it is to get a good score
        # (since this likely lead to missing large part of the protein)
        for i in reversed(range(j)):
            if (nident[i]+nident[j] < best_ident):
                if (best_prev_nident[j]+nident[j] < best_ident): # nothing better earlier 
                    break
                else:
                    continue; # overlap only decrease score so no need to compute them
            (compatible, max_overlap) = max_path_overlap(paths[i], hsp,  hsp[j])
            nident_ij=nident[i]+nident[j]-max_overlap
            # prefer solution with shorter path limit hsp included in hsp with the same scoring
            if compatible and nident_ij >= best_ident:
                best_ident = nident_ij
                best_prev = i
                # Update bounds
                bounds[j][0] = min(bounds[i][0], bounds[j][0])
                bounds[j][1] = max(bounds[i][1], bounds[j][1])
        if best_prev != -1:
            paths[j] = paths[best_prev] + [j]
        nident[j] = best_ident
        best_prev_nident[j] = max(best_prev_nident[j-1], nident[j]) if j > 1 else nident[j]
    # Find the index of the best scoring HSP
    best_idx = argmax(nident)
    best_nident = nident[best_idx]
    best_bounds = bounds[best_idx]
    best_path_ids = paths[best_idx]
    best_path = [[hsp[i][start_idx], hsp[i][end_idx]] for i in best_path_ids]
    print(best_path)
    print(best_path_ids)
    
    # Compute best_hsp_total_len
    best_hsp_total_len = best_bounds[1] - best_bounds[0] + 1

    # Extract additional info
    prot_genomic_len = class_rows["prot_len"].to_list()[0]

    # Compute penalty
    penalty = max(-1, pc_lenth_penalty * (1 - max(1,
                  math.exp((best_hsp_total_len - prot_genomic_len) / genomic_length_penalty))))

    # Calculate adjusted score
    adjusted_score = best_nident * (1 + penalty)

    # Return results
    return best_nident, adjusted_score, best_bounds, best_path 


def compute_class_score_old(class_rows, genomic_length_penalty, pc_lenth_penalty):
    # Sort rows by prot_start
    class_rows=class_rows.sort("prot_start")

    # Extract relevant columns as a list of rows
    hsp=[[row[0], row[1], row[2], row[3], row[4]]
           for row in zip(
               class_rows["prot_start"].to_list(),
               class_rows["prot_end"].to_list(),
               class_rows["loc_startS"].to_list(),
               class_rows["loc_endS"].to_list(),
               class_rows["nident"].to_list(),
    )]

    nident_idx=4  # Index of the nident column in hsp
    start_idx, end_idx=2, 3  # Indices for loc_startS and loc_endS

    # Initialize scores and bounds
    nident=[hsp[i][nident_idx] for i in range(len(hsp))]
    bounds=[[hsp[i][start_idx], hsp[i][end_idx]] for i in range(len(hsp))]

    # Compute scores and update bounds
    for j in range(len(hsp)):
        for i in range(j):  # Only consider previous HSPs
            nident_ij=combined_nident(hsp[i], hsp[j], nident[i])
            if nident_ij > nident[j]:
                nident[j]=nident_ij
                # Update bounds
                bounds[j][0]=min(bounds[i][0], bounds[j][0])
                bounds[j][1]=max(bounds[i][1], bounds[j][1])

    # Find the index of the best scoring HSP
    best_idx=argmax(nident)
    best_nident=nident[best_idx]
    best_bounds=bounds[best_idx]

    # Compute best_hsp_total_len
    best_hsp_total_len=best_bounds[1] - best_bounds[0] + 1

    # Extract additional info
    prot_genomic_len=class_rows["prot_len"].to_list()[0]

    # Compute penalty
    penalty=max(-1, pc_lenth_penalty * (1 - max(1,
                  math.exp((best_hsp_total_len - prot_genomic_len) / genomic_length_penalty))))

    # Calculate adjusted score
    adjusted_score=best_nident * (1 + penalty)

    # Return results
    return best_nident, adjusted_score, best_bounds


def add_classification_with_lists(blast_df, threshold):
    """
    Add class_id, class_id_startS, and class_id_endS columns to the BLAST DataFrame.
    """
    # Sort and extract necessary columns
    chrS_ids, prot_ids, loc_startS, loc_endS, max_intron_sizes, prot_len, strand, chr_id=(
        blast_df.sort(["chrS_id", "prot_id", "loc_startS"])
        .select(["chrS_id", "prot_id", "loc_startS", "loc_endS", "maxIntronSize", "prot_len", "strand", "chr_id"]).to_numpy().T
    )

    # Initialize variables for classification
    hsp_class=[]
    class_chrS, class_protId, class_startS, class_endS=[], [], [], []
    class_prot_lg, class_strand, class_chr=[], [], []
    # Iterate through rows and assign class_ids
    for i, (chrS_id, prot_id, startS, endS, max_intron, prot_lg, strand, chr) in enumerate(
        zip(chrS_ids, prot_ids, loc_startS, loc_endS,
            max_intron_sizes, prot_len, strand, chr_id)
    ):
        max_intron=math.ceil(1.1 * max(max_intron, threshold))

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
            class_endS[-1]=max(class_endS[-1], endS)

        # Assign class_id to the row
        hsp_class.append(len(class_startS))

    # Create class_info DataFrame
    class_info_df=pl.DataFrame({
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
    chr_max_coord=(
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
    blast_full=blast_full.with_columns(
        pl.when(pl.col("loc_start") < pl.col("loc_end"))
        .then(1)
        .otherwise(-1)
        .alias("strand")
    )
    blast_full=blast_full.with_columns(
        pl.when(pl.col("strand") == 1)
        .then(pl.col("chr_id").cast(str))
        .otherwise(pl.col("chr_id").cast(str) + "_revcomp")
        .alias("chrS_id")
    )

    # Step 3: Join chr_max_coord back to the main DataFrame
    blast_full=blast_full.join(chr_max_coord, on="chr_id")

    # Step 4: Add loc_startS and loc_endS based on strand
    blast_full=blast_full.with_columns(
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


def filter_class(scored_class, min_pident):
    max_class_id = int(scored_class.select(pl.col("class_id").max()).to_numpy()[0, 0])
    keep_class = [False] * (max_class_id+1)
    for group in scored_class.group_by("chr_id"):
        class_rows = group[1]  # Extract the DataFrame of the current group
        sorted_class = class_rows.filter((pl.col("score") / pl.col("prot_len")) > min_pident).sort("score", reverse=True)
        candidate_loci = InterLap()
        # Extract relevant columns as a list of rows
        loci = [[row[0], row[1], row[2]]
           for row in zip(
               sorted_class["loc_startS"].to_list(),
               sorted_class["loc_endS"].to_list(),
            sorted_class["class_i"].to_list(),
    )]
        for i in range(len(loci)):
            if not candidate_loci.__contains__((loci[i][0], loci[i][1])):
                candidate_loci.add((loci[i][0], loci[i][1]))
                keep_class[loci[i][2]][1]= True
    return pl.DataFrame({ "class_id": len(keep_class), "keep": keep_class })

def main():
    parser=argparse.ArgumentParser()

    parser.add_argument("-g", "--gff_file", type=str,
                        help="Exact path of gff file")
    parser.add_argument("-t", "--table", type=str,
                        help="Exact path of alignment res table")

    args=parser.parse_args()
    # debug
    # args.table = "/Users/ranwez/Desktop/TEST_REGION/blast_refProt_ENSG00000169598_full_pos.tsv"
    protein_info=get_prot_info(args.gff_file)

    # qseqid sseqid qlen length qstart qend sstart send nident pident gapopen evalue bitscore ### blast model prot to genome using tblastn
    blast_full=pl.read_csv(args.table, separator='\t', has_header=False, new_columns=["prot_id", "chr_id", "prot_len", "length", "prot_start",
                                                                                        "prot_end", "loc_start", "loc_end", "nident", "pident",
                                                                                        "gapopen", "evalue", "bitscore"])
    blast_full=add_simplified_coord(blast_full)
   # perform the join separately else some rows are lost ???
    blast_full2=blast_full.join(
        protein_info, on="prot_id", how="left")

    rows_with_none=blast_full2.filter(pl.col("maxIntronSize").is_null())
    if (rows_with_none.shape[0] > 0):
        print("Some proteins are not found in the GFF :")
        print(rows_with_none["prot_id"].to_list()[
              1:min(10, rows_with_none.shape[0])])
        exit(1)

    default_authorized_intron_length=int(protein_info.select(
        pl.col("maxIntronSize").quantile(0.7)).to_numpy()[0, 0])
    print(default_authorized_intron_length)
    blast_class, class_info=add_classification_with_lists(
        blast_full2, default_authorized_intron_length)
    print(blast_class.head())
    # Group by class_id and compute scores
    scores=[]
    for group in blast_class.group_by("class_id"):
        class_rows=group[1]  # Extract the DataFrame of the current group
        genomic_size_penalty=max(default_authorized_intron_length,
                                   class_rows["maxIntronSize"].to_list()[0])
        nident_penalty_per_size=0.01
        nident, score, bound,path=compute_class_score(
            class_rows, genomic_size_penalty, nident_penalty_per_size)
        scores.append({"class_id": group[0][0], "nident": nident, "score":
                       score, "startS": bound[0], "endS": bound[1], "path":path})
    # print(scores)
    scores_df=pl.DataFrame(scores).join(
        class_info, on="class_id", how="left")

    print(scores_df.head())
    scores_df.write_json("scores.json")

if __name__ == "__main__":
    main()
