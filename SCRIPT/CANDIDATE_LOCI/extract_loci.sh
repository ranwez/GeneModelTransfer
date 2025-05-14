#!/bin/bash
# extract_loci.sh -- Extract gene regions from a GFF file into individual FASTA files
# using the gene ID as the filename, and reverse‐complement sequences on the negative strand.
#
# Usage: extract_loci.sh <gff_file> <reference_fasta> <output_directory>
#
# Dependencies: awk, bedtools

#/storage/replicated/cirad/projects/GE2POP/2023_LRR/SVEVO3/Data_Package_01_12_22/DW_Svevo4.fasta

set -euo pipefail

extract_loci() {
  if [ "$#" -ne 3 ]; then
    echo "Usage: extract_loci.sh <gff_file> <reference_fasta> <output_directory>" >&2
    return 1
  fi

  local gff_file="$1"
  local ref_genome="$2"
  local OUT_DIR="$3"

  # If the reference ends with .fasta, get basename without extension.
  local ref_basename
  ref_basename=$(basename "$ref_genome" .fasta)

  # Create a temporary directory for intermediate files.
  local tmp_dir=$(mktemp -d -t EXTRACT_LOCI_$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXXXX)
  #echo $tmp_dir
  #trap 'rm -rf "$tmp_dir"' EXIT

  # 1. Convert GFF to BED.
  #    For each line with feature "gene", extract the gene ID from the attributes (assumes "ID=..."),
  #    adjust start coordinate to 0-based, and output:
  #      CHR, start, end, geneID, placeholder, strand.
  local tmp_bed="${tmp_dir}/genes_tmp.bed"
  awk 'BEGIN { OFS="\t" }
       $3=="gene" {
         split($9, infos, ";");
         id = substr(infos[1], 4);  # Assumes first attribute is like "ID=gene123"
         start = $4 - 1;
         print $1, start, $5, id, ".", $7
       }' "$gff_file" >"$tmp_bed"

  # 2. Ensure an index (.fai) for the reference exists.
  #    We check if "$ref_genome.fai" or "$ref_genome.fasta.fai" exists.
  local fai_file="${ref_genome}.fai"
  if [ ! -f "${ref_genome}.fai" ]; then
    # If not, use bedtools getfasta on a trivial region (1 10 on first chromosome) to force creation of the index.
    local chr=$(head -n 1 "$tmp_bed" | cut -f1)
    local trivial_bed="${tmp_dir}/__tmp.bed"
    echo -e "$chr\t1\t10\t+\tdumb" >"$trivial_bed"
    # Run bedtools getfasta; this should force creation of an index.
    bedtools getfasta -name -s -fi "$ref_genome" -bed "$trivial_bed" >/dev/null
    # Now try to set fai_file (it might be created as ref_genome.fai).
    if [ ! -f "${ref_genome}.fai" ]; then
      echo "Error: unable to create FAI index for $ref_genome" >&2
      return 1
    fi
  fi

  # 3. Use awk (with the NR==FNR idiom) to adjust BED coordinates based on the FAI.
  local corrected_bed="${tmp_dir}/genes_corrected.bed"
  awk 'BEGIN { OFS = "\t" }
      NR==FNR {
          # FAI: column1 = chrom, column2 = length.
          chr_size[$1] = $2;
          next
       }
       {
          if ($2 < 0) $2 = 0;
          if ($3 > chr_size[$1]) $3 = chr_size[$1];
          print
       }' "$fai_file" "$tmp_bed" >"$corrected_bed"

  # 4. Extract FASTA sequences for these regions.
  #    The -s flag makes bedtools reverse-complement negative strand sequences.
  local multi_fasta="${tmp_dir}/multi_fasta.fasta"
  bedtools getfasta -name -s -fi "$ref_genome" -bed "$corrected_bed" >"$multi_fasta"

  # 5. Split the multi-FASTA file into individual FASTA files,
  #    using the header (without the '>' and without parentheses) as the filename,
  #    and saving them in ${OUT_DIR}.
  mkdir -p "${OUT_DIR}"
  awk -v outdir="${OUT_DIR}" '
    /^>/ {
      if (out) close(out);
      header = substr($0, 2);
      sub(/\([+-]\)$/, "", header);
      sub(/::.*/, "", header);
      out = outdir "/" header;
      print ">" header > out;
      next;
    }
    { print >> out }
    END { if (out) close(out) }
	' "$multi_fasta"

  echo "Individual FASTA files have been created in ${OUT_DIR}"

  rm -rf "$tmp_dir"
}

# If the script is executed directly, run the function with the arguments.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  extract_loci "$@"
fi
