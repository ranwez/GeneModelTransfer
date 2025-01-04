import polars as pl
from attrs import define
from pathlib import Path
from CANDIDATE_LOCI.gff_utils import GeneInfo

@define(slots=True)
class Bounds:
    start:int
    end:int # inclusive
    def __init__(self, start:int, end:int):
        self.start = min(start,end)
        self.end = max(start,end)
    def length(self) -> int:
        return self.end - self.start + 1
    def overlap(self, other: "Bounds") -> int:
        return max(0, min(self.end, other.end) - max(self.start, other.start))

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
        loc_endS: int
    ):
        self.prot_id = prot_id
        self.chr_id = chr_id
        self.prot_len = prot_len
        self.length = length
        self.nident = nident
        self.pident = pident
        self.strand = strand

        self.prot_bounds = Bounds(prot_start, prot_end)
        self.loc_bounds = Bounds(loc_start, loc_end)
        self.locS_bounds = Bounds(loc_startS, loc_endS)

    @classmethod
    def build_dummy(cls,prot_id="") -> "HSP":
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
            loc_endS=0
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
            "qseqid",    # Query sequence ID
            "sseqid",    # Subject sequence ID
            "pident",    # Percentage identity
            "length",    # Alignment length
            "mismatch",  # Number of mismatches
            "gapopen",   # Number of gap openings
            "qstart",    # Start of alignment in query
            "qend",      # End of alignment in query
            "sstart",    # Start of alignment in subject
            "send",      # End of alignment in subject
            "evalue",    # Expect value
            "bitscore"   # Bit score
        ]

    # Read the BLAST tabular file using Polars LazyFrame
    # bitscore may, in rare case, have one decimal point, so we need to set the schema length to a large value
    blast_lf = pl.scan_csv(blast_tab_file, separator="\t", has_header=False, new_columns=columns, infer_schema_length=10000)

    return blast_lf

def blast_to_HSPs(blast_tab_file: str) -> list[HSP_chr]:
    """
    Processes a BLAST tabular output file and returns a list of HSP_chr objects.

    Parameters:
        blast_tab_file (str): Path to the BLAST tabular output file.

    Returns:
        list[HSP_chr]: A list of HSP_chr objects with grouped HSPs.
    """
    # Define the desired column names
    columns = [
        "prot_id", "chr_id", "prot_len", "length", "prot_start",
        "prot_end", "loc_start", "loc_end", "nident", "pident",
        "gapopen", "evalue", "bitscore"
    ]

    # Parse the BLAST results into a LazyFrame
    blast_lf = parse_blast_results(blast_tab_file, columns)

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
    blast_df = blast_df.with_columns([
        pl.when(pl.col("strand") == +1)
          .then(pl.col("loc_start"))
          .otherwise(max_coord - pl.col("loc_start"))
          .alias("loc_startS"),

        pl.when(pl.col("strand") == +1)
          .then(pl.col("loc_end"))
          .otherwise(max_coord - pl.col("loc_end"))
          .alias("loc_endS")
    ])

    # Group by chromosome and strand
    grouped = blast_df.group_by(["chr_id", "strand"]).agg(
        [pl.struct([
            "prot_id", "chr_id", "prot_len", "length", "prot_start", "prot_end", 
            "loc_start", "loc_end", "nident", "pident", "strand", "loc_startS", "loc_endS"
        ]).alias("HSP")]
    )

    # Convert to a list of HSP_chr objects
    hsp_chr_list = [
        HSP_chr(
            str(row["chr_id"]),
            str(row["strand"]),
            [HSP(**hsp) for hsp in row["HSP"]]
        )
        for row in grouped.to_dicts()
    ]

    return hsp_chr_list
#def compute_loci(hsp_chr_list: list[HSP_chr], prot_infos : dict[]) -> list[HSP_chr]:
    
# Example usage
def main():
    test_data_path = Path(__file__).parent / "tests" / "data" / "tblastn_ENSG00000169598_2chr.tsv"
    print (test_data_path)
    hsp_chr_list = blast_to_HSPs(test_data_path)
    # Sort and display the first few HSP_chr tuples
    for hsp_chr in hsp_chr_list:
        sorted_hsp = sorted(hsp_chr.HSP, key=lambda h: (h.prot_id, h.loc_startS))
        print(f"Chromosome: {hsp_chr.chr_id}, Strand: {hsp_chr.strand}")
        for hsp in sorted_hsp:
            print(hsp)

if __name__ == "__main__":
    main()