import sys
import polars as pl
from attrs import define
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from CANDIDATE_LOCI.gff_utils import GeneInfo
from CANDIDATE_LOCI.bounds import Bounds


@define(slots=True)
class HSP:
    prot_id: str
    chr_id: str
    prot_len: int
    length: int
    nident: int
    pident: float
    strand: int
    prot_bounds: Bounds
    loc_bounds: Bounds
    locS_bounds: Bounds

    def __init__(
        self,
        prot_id: str,
        chr_id: str,
        prot_len: int,
        length: int,
        prot_start: int,
        prot_end: int,
        loc_start: int,
        loc_end: int,
        nident: int,
        pident: float,
        strand: int,
        loc_startS: int,
        loc_endS: int,
    ):
        self.prot_id = prot_id
        self.chr_id = chr_id
        self.prot_len = int(prot_len)
        self.length = int(length)
        self.nident = int(nident)
        self.pident = float(pident)
        self.strand = int(strand)

        self.prot_bounds = Bounds(int(prot_start), int(prot_end))
        self.loc_bounds = Bounds(int(loc_start), int(loc_end))
        self.locS_bounds = Bounds(int(loc_startS), int(loc_endS))

    @classmethod
    def build_dummy(cls, prot_id="") -> "HSP":
        return cls(
            prot_id=prot_id,
            chr_id="",
            prot_len=0,
            length=0,
            nident=0,
            pident=0.0,
            strand=0,
            prot_start=0,
            prot_end=0,
            loc_start=0,
            loc_end=0,
            loc_startS=0,
            loc_endS=0,
        )


@define(slots=True)
class HSP_chr:
    chr_id: str
    strand: str
    HSP: list[HSP]


def parse_blast_results(blast_tab_file: str, columns: list[str] = None) -> pl.LazyFrame:
    """
    Parses a BLAST tabular output file into a Polars LazyFrame.

    Parameters:
        blast_tab_file (str): Path to the BLAST tabular output file.
        columns (list[str], optional): List of column names. Defaults to standard BLAST tabular columns.

    Returns:
        pl.LazyFrame: A Polars LazyFrame containing the parsed BLAST results.
    """
    # Default column names based on standard BLAST tabular format (default 12 columns)
    if columns is None:
        columns = [
            "qseqid",  # Query sequence ID
            "sseqid",  # Subject sequence ID
            "pident",  # Percentage identity
            "length",  # Alignment length
            "mismatch",  # Number of mismatches
            "gapopen",  # Number of gap openings
            "qstart",  # Start of alignment in query
            "qend",  # End of alignment in query
            "sstart",  # Start of alignment in subject
            "send",  # End of alignment in subject
            "evalue",  # Expect value
            "bitscore",  # Bit score
        ]
    # Define column types
    column_types = {
        "chr_id": pl.Utf8,  # Force chr_id to be a string
        "prot_len": pl.Int64,  # Ensure prot_len is an integer
        "length": pl.Int64,
        "prot_start": pl.Int64,
        "prot_end": pl.Int64,
        "loc_start": pl.Int64,
        "loc_end": pl.Int64,
        "nident": pl.Int64,
        "pident": pl.Float64,  # Float for percentage identity
        "gapopen": pl.Int64,
        "evalue": pl.Float64,  # Scientific notation for evalue
        "bitscore": pl.Float64,
    }

    # Read the BLAST tabular file using Polars LazyFrame
    # bitscore may, in rare case, have one decimal point, so we need to set the schema length to a large value
    blast_lf = pl.scan_csv(
        blast_tab_file,
        separator="\t",
        has_header=False,
        new_columns=columns,
        schema_overrides=column_types,
        infer_schema_length=1000,
    )

    return blast_lf


def blast_to_HSPs(blast_tab_file: str, chr: str = None) -> list[HSP_chr]:
    """
    Processes a BLAST tabular output file and returns a list of HSP_chr objects.

    Parameters:
        blast_tab_file (str): Path to the BLAST tabular output file.

    Returns:
        list[HSP_chr]: A list of HSP_chr objects with grouped HSPs.
    """
    # Define the desired column names
    columns = [
        "prot_id",
        "chr_id",
        "prot_len",
        "length",
        "prot_start",
        "prot_end",
        "loc_start",
        "loc_end",
        "nident",
        "pident",
        "gapopen",
        "evalue",
        "bitscore",
    ]

    # Parse the BLAST results into a LazyFrame
    blast_lf = parse_blast_results(blast_tab_file, columns)
    if chr != None:
        blast_lf = blast_lf.filter(pl.col("chr_id") == chr)
    blast_lf = blast_lf.drop(["gapopen", "evalue", "bitscore"])

    # Calculate max_coord using Python
    blast_df = blast_lf.collect()
    max_loc_start = blast_df["loc_start"].max()
    max_loc_end = blast_df["loc_end"].max()
    max_coord = max(max_loc_start, max_loc_end) + 1000

    # Add strand column
    blast_df = blast_df.with_columns(
        pl.when(pl.col("loc_start") < pl.col("loc_end"))
        .then(+1)
        .otherwise(-1)
        .alias("strand")
    )

    # Add loc_startS and loc_endS columns
    blast_df = blast_df.with_columns(
        [
            pl.when(pl.col("strand") == +1)
            .then(pl.col("loc_start"))
            .otherwise(max_coord - pl.col("loc_start"))
            .alias("loc_startS"),
            pl.when(pl.col("strand") == +1)
            .then(pl.col("loc_end"))
            .otherwise(max_coord - pl.col("loc_end"))
            .alias("loc_endS"),
        ]
    )

    # Group by chromosome and strand
    grouped = blast_df.group_by(["chr_id", "strand"]).agg(
        [
            pl.struct(
                [
                    "prot_id",
                    "chr_id",
                    "prot_len",
                    "length",
                    "prot_start",
                    "prot_end",
                    "loc_start",
                    "loc_end",
                    "nident",
                    "pident",
                    "strand",
                    "loc_startS",
                    "loc_endS",
                ]
            ).alias("HSP")
        ]
    )

    # Convert to a list of HSP_chr objects
    hsp_chr_list = [
        HSP_chr(
            str(row["chr_id"]), str(row["strand"]), [HSP(**hsp) for hsp in row["HSP"]]
        )
        for row in grouped.to_dicts()
    ]

    return hsp_chr_list


# def compute_loci(hsp_chr_list: list[HSP_chr], prot_infos : dict[]) -> list[HSP_chr]:


def blast_to_sortedHSPs(blast_tab_file: str, output_blast_file: str, chr: str = None):
    """
    Processes a BLAST tabular output file and returns a list of HSP_chr objects.

    Parameters:
        blast_tab_file (str): Path to the BLAST tabular output file.

    Returns:
        list[HSP_chr]: A list of HSP_chr objects with grouped HSPs.
    """
    # Define the desired column names
    columns = [
        "prot_id",
        "chr_id",
        "prot_len",
        "length",
        "prot_start",
        "prot_end",
        "loc_start",
        "loc_end",
        "nident",
        "pident",
        "gapopen",
        "evalue",
        "bitscore",
    ]

    # Parse the BLAST results into a LazyFrame
    blast_lf = parse_blast_results(blast_tab_file, columns)
    if chr != None:
        blast_lf = blast_lf.filter(pl.col("chr_id") == chr)

    # drop unecessary columns

    blast_lf = blast_lf.select(columns)
    blast_lf = blast_lf.drop(["gapopen", "evalue", "bitscore"])

    # Calculate max_coord using Python
    blast_df = blast_lf.collect()
    max_loc_start = blast_df["loc_start"].max()
    max_loc_end = blast_df["loc_end"].max()
    max_coord = max(max_loc_start, max_loc_end) + 1000

    # Add strand column
    blast_df = blast_df.with_columns(
        pl.when(pl.col("loc_start") < pl.col("loc_end"))
        .then(+1)
        .otherwise(-1)
        .alias("strand")
    )

    # Add loc_startS and loc_endS columns
    blast_df = blast_df.with_columns(
        [
            pl.when(pl.col("strand") == +1)
            .then(pl.col("loc_start"))
            .otherwise(max_coord - pl.col("loc_start"))
            .alias("loc_startS"),
            pl.when(pl.col("strand") == +1)
            .then(pl.col("loc_end"))
            .otherwise(max_coord - pl.col("loc_end"))
            .alias("loc_endS"),
        ]
    )

    blast_df = blast_df.sort("chr_id", "strand", "prot_id", "loc_startS")
    blast_df.write_csv(output_blast_file, separator="\t", include_header=False)


# Example usage
def main():
    test_data_path = (
        Path(__file__).parent / "tests" / "data" / "tblastn_ENSG00000169598_2chr.tsv"
    )
    print(test_data_path)
    blast_to_HSPs(test_data_path, "__test_output_blast.tsv")


if __name__ == "__main__":
    main()
