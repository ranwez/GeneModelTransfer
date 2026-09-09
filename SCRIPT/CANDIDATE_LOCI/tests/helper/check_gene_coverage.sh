#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage:
  $(basename "$0") \
    --ref_gff FILE \
    --candidate_gff FILE \
    --out FILE \
    [--chr CHROMOSOME]

Required arguments:
  --ref_gff FILE        Reference annotation GFF
  --candidate_gff FILE  Candidate loci GFF
  --out FILE            Output TSV file

Optional arguments:
  --chr CHROMOSOME      Chromosome name to retain [default: all]
  -h, --help            Show this help message

Examples:
  $(basename "$0") \
    --ref_gff DATA_IN/OsjNip_HIPPHPP_202606v1.gff \
    --candidate_gff _test_chr4.gff \
    --out chr4_coverage.tsv \
    --chr Chr4

  $(basename "$0") \
    --ref_gff reference.gff \
    --candidate_gff candidate_loci.gff \
    --out all_gene_coverage.tsv
EOF
}

REF_GFF=""
CANDIDATE_GFF=""
OUTPUT=""
CHR="all"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref_gff)
            REF_GFF="${2:?Missing value after --ref_gff}"
            shift 2
            ;;
        --candidate_gff)
            CANDIDATE_GFF="${2:?Missing value after --candidate_gff}"
            shift 2
            ;;
        --out)
            OUTPUT="${2:?Missing value after --out}"
            shift 2
            ;;
        --chr)
            CHR="${2:?Missing value after --chr}"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$REF_GFF" || -z "$CANDIDATE_GFF" || -z "$OUTPUT" ]]; then
    echo "Error: --ref_gff, --candidate_gff and --out are required." >&2
    usage >&2
    exit 1
fi

if [[ ! -f "$REF_GFF" ]]; then
    echo "Error: reference GFF not found: $REF_GFF" >&2
    exit 1
fi

if [[ ! -f "$CANDIDATE_GFF" ]]; then
    echo "Error: candidate GFF not found: $CANDIDATE_GFF" >&2
    exit 1
fi

if ! command -v bedtools >/dev/null 2>&1; then
    echo "Error: bedtools is not available in PATH." >&2
    exit 1
fi

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

filter_genes() {
    local input_gff="$1"
    local output_gff="$2"

    awk -F'\t' -v chr="$CHR" '
        BEGIN { OFS="\t" }

        $0 !~ /^#/ &&
        NF >= 9 &&
        $3 == "gene" &&
        (chr == "all" || $1 == chr) {
            print
        }
    ' "$input_gff" > "$output_gff"
}

filter_genes "$REF_GFF"       "$tmpdir/ref_genes.gff"
filter_genes "$CANDIDATE_GFF" "$tmpdir/candidate_genes.gff"

if [[ ! -s "$tmpdir/ref_genes.gff" ]]; then
    echo "Error: no reference genes found for chromosome filter '$CHR'." >&2
    exit 1
fi

bedtools coverage \
    -a "$tmpdir/ref_genes.gff" \
    -b "$tmpdir/candidate_genes.gff" |
awk -F'\t' '
    BEGIN {
        OFS="\t"
        print "chr", "start", "end", "gene_id",
              "gene_length", "covered_bases",
              "coverage_percent", "fully_covered"
    }

    {
        gene_id = "."
        n = split($9, attrs, ";")

        for (i = 1; i <= n; i++) {
            if (attrs[i] ~ /^ID=/) {
                gene_id = substr(attrs[i], 4)
                break
            }
        }

        covered_bases = $(NF - 2)
        gene_length   = $(NF - 1)
        coverage      = $NF * 100

        print $1, $4, $5, gene_id,
              gene_length, covered_bases,
              sprintf("%.2f", coverage),
              (covered_bases == gene_length ? "YES" : "NO")
    }
' > "$OUTPUT"

echo "Chromosome filter: $CHR"
echo "Result written to: $OUTPUT"
