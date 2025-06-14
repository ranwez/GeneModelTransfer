import sys
from typing import Optional
import polars as pl
from attrs import define
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from CANDIDATE_LOCI.bounds import Bounds


# from bounds import Bounds
@define(slots=True, eq=True)
class CdsInfo:
    gene_id: str
    chr_id: str
    cds_bounds: list[Bounds]
    coding_region: Bounds
    gene_bounds: Bounds

    def get_genomic_coord(self, prot_coordinate: int, strand :int) -> tuple[int, int, int]:
        """
        Convert the protein coordinate to genomic coordinate.
        Parameters:
            prot_coord (int): The protein coordinate.
        Returns:
            (int, int, int) : The genomic coordinate of (the_start_of_the_protein, the_protein_coordinate, the_end_of_the_protein)
        """
        cds_coordinate = prot_coordinate * 3
        gene_start = self.coding_region.start if strand == 1 else self.coding_region.end
        gene_end = self.coding_region.end if strand == 1 else self.coding_region.start
        cds_Sbounds = self.cds_bounds if strand == 1 else self.cds_bounds[::-1]
        prev_cds_lg = 0
        for cds_bound in cds_Sbounds:
            if prev_cds_lg <= cds_coordinate <= prev_cds_lg + cds_bound.length():
                genomic_coord = cds_bound.start + (cds_coordinate - prev_cds_lg)
                if strand == -1:
                    genomic_coord = cds_bound.end - (cds_coordinate - prev_cds_lg)
                return gene_start, genomic_coord, gene_end
            prev_cds_lg += cds_bound.length()
             
        raise RuntimeError(
            f"No genomic coordinate found for protein end: {self.gene_id} "
            f"for coordinate:\n {prot_coordinate}\n"
        )

    def __init__(self, gene_id: str, chr_id: str, cds_bounds: list[Bounds], gene_bounds: Bounds):
        self.gene_id = gene_id
        self.chr_id = chr_id
        self.cds_bounds = sorted(cds_bounds, key=lambda b: b.start)
        self.coding_region = Bounds(
            start=self.cds_bounds[0].start, end=self.cds_bounds[-1].end
        )
        self.gene_bounds = gene_bounds


@define(slots=True, eq=True)
class GeneInfo:
    gene_id: str
    chr_id: str
    coding_region:Bounds
    strand: int
    longest_intron: int


def filter_mRNA_by_attribute(
    file_path,
    mRNA_attribute,
    mRNA_value,
    output_filtered_gff,
    output_unresolved_gff=None,
):
    df = parse_gff(file_path)

    df = df.with_columns(
        [
            df["attributes"]
            .str.extract(rf"{mRNA_attribute}=([^;]+)", 1)
            .map_elements(lambda x: 1 if x == mRNA_value else 0, return_dtype=int)
            .alias("mRNA_to_keep")
        ]
    )

    genes = df.filter(df["type"] == "gene").select(["gene"]).unique()

    mRNA_rep = (
        df.filter(df["type"] == "mRNA")
        .select(["ID", "mRNA_to_keep"])
        .rename({"ID": "mRNA", "mRNA_to_keep": "to_keep"})
    )
    save_df = df
    df = df.join(mRNA_rep, on="mRNA", how="left")

    gene_with_mRNA = (
        df.filter((df["type"] == "mRNA") & (df["to_keep"] == 1))
        .select(["gene"])
        .unique()
        .with_columns(pl.lit(1).alias("has_mRNA"))
    )

    if output_unresolved_gff != None:
        genes_without_mrna = (
            genes.join(gene_with_mRNA, on="gene", how="anti").select(["gene"]).unique()
        )
        save_df = save_df.join(genes_without_mrna, on="gene", how="inner")
        save_df = sort_gff(save_df)
        save_df = save_df.drop(["gene", "mRNA", "mRNA_to_keep", "ID", "ParentID"])
        save_df = save_df.with_columns(
            [
                pl.col(col).fill_null(".")
                for col in ["source", "score", "strand", "phase"]
            ]
        )
        save_df.write_csv(output_unresolved_gff, separator="\t", include_header=False)
        genes_without_mrna = None
        save_df = None

    df = df.join(gene_with_mRNA, on="gene", how="left")

    df = df.with_columns(
        to_keep=pl.when(pl.col("type") == "gene")
        .then(pl.col("has_mRNA"))
        .otherwise(pl.col("to_keep"))
    )

    df = df.filter(df["to_keep"] == 1)
    df = sort_gff(df)
    df = df.drop(
        ["mRNA_to_keep", "mRNA", "gene", "ID", "ParentID", "to_keep", "has_mRNA"]
    )
    df = df.with_columns(
        [pl.col(col).fill_null(".") for col in ["source", "score", "strand", "phase"]]
    )
    df.write_csv(output_filtered_gff, separator="\t", include_header=False)


def filter_feature_by_geneID(file_path, file_ID_list, output_filtered_gff):
    df = parse_gff(file_path)
    ID_list = pl.read_csv(
        file_ID_list,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        new_columns=["gene"],
    )
    df = df.filter(df["gene"].is_in(ID_list["gene"]))
    df = df.drop(["gene", "ID", "mRNA", "ParentID"])
    df = df.with_columns(
        [pl.col(col).fill_null(".") for col in ["source", "score", "strand", "phase"]]
    )
    df.write_csv(output_filtered_gff, separator="\t", include_header=False)


def parse_gff(file_path):
    """
    Parse a GFF file and return a Polars DataFrame.

    Parameters:
        file_path (str): Path to the GFF file.

    Returns:
        pl.DataFrame: A Polars DataFrame containing the parsed GFF data with added ID, ParentID, mRNA and gene columns.
    """
    column_names = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]

    column_types = {
        "seqid": pl.Utf8,
        "source": pl.Utf8,
        "type": pl.Utf8,
        "start": pl.Int64,
        "end": pl.Int64,
        "score": pl.Utf8,
        "strand": pl.Utf8,
        "phase": pl.Int64,
        "attributes": pl.Utf8,
    }

    try:
        df = pl.read_csv(
            file_path,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            new_columns=column_names,
            dtypes=column_types,
            null_values=".",
        )

        # Extract IDs from attributes column
        df = df.with_columns(
            [
                df["attributes"].str.extract(r"ID=([^;]+)", 1).alias("ID"),
                df["attributes"].str.extract(r"Parent=([^;]+)", 1).alias("ParentID"),
            ]
        )

        # Populate mRNA column
        df = df.with_columns(
            [
                pl.when(df["type"] == "mRNA")
                .then(df["ID"])
                .otherwise(df["ParentID"])
                .alias("mRNA")
            ]
        )

        # Create mapping for mRNA to gene relationships
        mRNA_mapping = (
            df.filter(df["type"] == "mRNA")
            .select(["ID", "ParentID"])
            .rename({"ID": "mRNA", "ParentID": "gene"})
        )

        # Join to populate gene column based on mRNA
        df = df.join(mRNA_mapping, on="mRNA", how="left")

        # Assign gene IDs to gene features
        df = df.with_columns(
            [
                pl.when(df["type"] == "gene")
                .then(df["ID"])
                .otherwise(df["gene"])
                .alias("gene")
            ]
        )

        return df
    except Exception as e:
        print(f"Error parsing GFF file: {e}")
        exit(1)


def sort_gff(df):
    # sort genes by chr_id and start, add row number,
    genes_order = (
        df.filter(df["type"] == "gene")
        .sort(["seqid", "start"])
        .with_row_index("gene_order")
        .select(["gene", "gene_order"])
    )
    # sort genes by chr_id and start, add row number,
    mrna_order = (
        df.filter(df["type"] == "mRNA")
        .sort(["seqid", "start"])
        .with_row_index("mRNA_order", offset=2)
        .select(["mRNA", "mRNA_order"])
    )
    # add a feature ordering in "type_order" column : 1 for gene 2 for mRNA 3 for others
    df = df.with_columns(
        pl.when(pl.col("type") == "gene")
        .then(1)
        .when(pl.col("type") == "mRNA")
        .then(2)
        .otherwise(3)
        .alias("type_order")
    )
    df = df.join(genes_order, on="gene", how="left")
    df = df.join(mrna_order, on="mRNA", how="left")
    # set mRNA order for gene features to 0 so that they appear first
    df = df.with_columns(
        pl.when(pl.col("type") == "gene")
        .then(0)
        .otherwise(pl.col("mRNA_order"))
        .alias("mRNA_order")
    )
    df = df.sort(["gene_order", "mRNA_order", "type_order", "start"])
    df = df.drop(["gene_order", "mRNA_order", "type_order"])
    return df


def get_coding_regions(df):
    """
    Compute the min and max CDS coordinates for each gene.

    Parameters:
        df (pl.DataFrame): A Polars DataFrame containing parsed GFF data.

    Returns:
        pl.DataFrame: A DataFrame with the min and max CDS coordinates for each gene.
    """
    cds_coordinates = (
        df.filter(df["type"] == "CDS")
        .group_by("gene")
        .agg(
            [
                pl.col("start").min().alias("coding_start"),
                pl.col("end").max().alias("coding_end"),
                pl.col("strand").first(),
                pl.col("seqid").first().alias("chr_id"),
            ]
        )
    )

    return cds_coordinates


def get_longest_intron(df) -> pl.DataFrame:
    """
    Compute the length of the longest intron for each gene.

    Parameters:
        df (pl.DataFrame): A Polars DataFrame containing parsed GFF data.

    Returns:
        pl.DataFrame: A DataFrame with the longest intron length for each gene.
    """
    # Step 1: Filter CDS features, sort by mRNA and start, and calculate previous mRNA and stop
    cds_with_prev = (
        df.filter(df["type"] == "CDS")
        .sort(["mRNA", "start"])
        .with_columns(
            [
                pl.col("mRNA").shift(1).alias("prev_mRNA"),
                pl.col("end").shift(1).alias("prev_end"),
            ]
        )
    )

    # Step 2: Calculate intron lengths
    intron_lengths = cds_with_prev.with_columns(
        [
            pl.when(pl.col("mRNA") != pl.col("prev_mRNA"))
            .then(0)
            .otherwise(pl.col("start") - pl.col("prev_end") - 1)
            .map_elements(lambda x: max(0, x), return_dtype=pl.Int64)
            .alias("intron_length")
        ]
    )

    # Step 3: Ensure genes with single exons are included with a value of 0
    all_genes = df.select("gene").unique()

    intron_aggregated = (
        intron_lengths.filter(pl.col("intron_length") > 0)  # Keep only valid introns
        .join(df.select(["mRNA", "gene"]).unique(), on="mRNA", how="left")
        .group_by("gene")
        .agg([pl.col("intron_length").max().alias("longest_intron")])
    )

    longest_introns_by_gene = all_genes.join(
        intron_aggregated, on="gene", how="left"
    ).with_columns([pl.col("longest_intron").fill_null(0)])

    return longest_introns_by_gene


def gff_to_cdsInfo(gff_file: str, relevant_gene_ids:Optional[set[str]]) -> dict:
    df = parse_gff(gff_file)
    # filter to keep only feature of relevant genes
    if relevant_gene_ids != None:
        df = df.filter(df["gene"].is_in(list(relevant_gene_ids)))
    
    # Get CDS coordinates
    cds_coordinates = (
        df.filter(df["type"] == "CDS")
        .group_by("gene")
        .agg([
            pl.col("start").alias("start"), 
            pl.col("end").alias("end"),
            pl.col("seqid").first().alias("chr_id")
        ])
    )
    
    # Get gene bounds
    gene_coordinates = (
        df.filter(df["type"] == "gene")
        .group_by("gene")
        .agg([
            pl.col("start").first().alias("gene_start"),
            pl.col("end").first().alias("gene_end")
        ])
    )
    
    # Join CDS and gene information
    merged_data = cds_coordinates.join(gene_coordinates, on="gene", how="left")

    cds_dict = {
        row["gene"]: CdsInfo(
            gene_id=row["gene"],
            chr_id=row["chr_id"],
            cds_bounds=[
                Bounds(start=start, end=end)
                for start, end in zip(row["start"], row["end"])
            ],
            gene_bounds=Bounds(start=row["gene_start"], end=row["gene_end"])
        )
        for row in merged_data.to_dicts()
    }
    return cds_dict


def gff_to_geneInfo(gff_file: str, intron_quantile: float) -> tuple[dict, int]:
    """
    Parse a GFF file and return a dictionary of GeneInfo objects keyed by prot_id.

    Parameters:
        gff_file (str): Path to the GFF file.
        intron_lg_quantile (float): The quantile to use for the default intron length.
    Returns:
        tuple: A tuple containing the dictionary of GeneInfo objects and the default intron length.
        dict: A dictionary with gene_id as keys and GeneInfo objects as values.
    """
    df = parse_gff(gff_file)
    coding_regions_df = get_coding_regions(df)
    longest_intron_df = get_longest_intron(df)

    # Merge info from the three DataFrames
    merged_df = coding_regions_df.join(
        longest_intron_df, on="gene", how="left"
    ).fill_null(0)

    # Construct the dictionary of ProtInfo objects
    prot_dict = {
        row["gene"]: GeneInfo(
            gene_id=str(row["gene"]),
            chr_id=str(row["chr_id"]),
            coding_region=Bounds(row["coding_start"],row["coding_end"]),
            strand=1 if row["strand"] == "+" else -1,
            longest_intron=row["longest_intron"],
        )
        for row in merged_df.to_dicts()
    }
    default_intron_length = longest_intron_df.select(
        pl.col("longest_intron").quantile(intron_quantile)
    ).to_numpy()[0, 0]
    return (prot_dict, default_intron_length)


# Only expose ProtInfo and parse_gff_to_protinfo for external use
__all__ = ["ProtInfo", "gff_to_geneInfo"]


def main():
    # Example usage
    test_data_path = Path(__file__).parent / "tests" / "data" / "ENSG00000160679.gff"
    df = parse_gff(test_data_path)
    coding_regions = get_coding_regions(df)
    longest_introns = get_longest_intron(df)

    print(coding_regions)
    print(longest_introns)

    (prot_dict, def_intron_lg) = gff_to_geneInfo(test_data_path, 0.7)
    print(prot_dict)
    print(def_intron_lg)


if __name__ == "__main__":
    main()
