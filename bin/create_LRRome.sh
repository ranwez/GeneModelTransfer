#!/bin/bash
#========================================================
# PROJET : LRRtransfer
# SCRIPT : create_LRRome.sh
# AUTHOR : Celine Gottin & Thibaud Vicat & Vincent Ranwez
# CREATION : 2021.05.07
#========================================================
# DESCRIPTION
# Build an LRRome from a reference genome FASTA and its
# corresponding reference GFF annotation.
#
# The script extracts reference loci and generates the
# sequence resources required by the annotation transfer
# workflow:
#   - REF_LOCI/
#   - REF_LOCI_GFF/
#   - REF_PEP/
#   - REF_cDNA/
#   - REF_EXONS/
#   - REF_loci.fasta
#   - REF_proteins.fasta
#   - REF_cDNA.fasta
#   - REF_exons.fasta
#   - REF_LOCI_PROVENANCE.tsv
#
# ARGUMENTS
#   $1 : Reference genome FASTA
#   $2 : Reference GFF
#   $3 : Output LRRome directory
#   $4 : Path to the LRRtransfer SCRIPT directory
#
# DEPENDENCIES
#   - bash
#   - gawk
#   - python3
#   - bedtools
#========================================================

set -euo pipefail

#========================================================
#                Environment & variables
#========================================================
REF_GENOME=$1
REF_GFF=$2
LRROME_DIR=$3
LRR_SCRIPT=$4

#========================================================
#                        Functions
#========================================================

function extractSeq {
	##usage :: extractSeq multifasta.file
	##Extracting each sequence from a fasta in separate files
	gawk -F"[;]" '{
    if ($1~/>/) {
      line=$1
      gsub(">","")
      filename=$1
      print(line) > filename
    } else {
      print > filename
    }
  }' $1
}

export -f extractSeq


#========================================================
#                Script
#========================================================

mkdir -p "$LRROME_DIR"
cd "$LRROME_DIR"


mkdir -p REF_PEP
mkdir -p REF_EXONS
mkdir -p REF_cDNA
mkdir -p REF_LOCI
mkdir -p REF_LOCI_GFF

bash "${LRR_SCRIPT}/CANDIDATE_LOCI/extract_loci.sh" "${REF_GFF}" "${REF_GENOME}" REF_LOCI

cat REF_LOCI/* > REF_loci.fasta

python3 "${LRR_SCRIPT}/ANNOTATION_TRANSFER/prepare_reference_loci.py" "${REF_GFF}" REF_LOCI REF_LOCI_GFF

python3 ${LRR_SCRIPT}/Extract_sequences_from_genome.py -g ${REF_GFF} -f ${REF_GENOME} -o REF_proteins.fasta -t FSprot --no_FS_codon
python3 ${LRR_SCRIPT}/Extract_sequences_from_genome.py -g ${REF_GFF} -f ${REF_GENOME} -o REF_cDNA.fasta -t FScdna --no_FS_codon
python3 ${LRR_SCRIPT}/Extract_sequences_from_genome.py -g ${REF_GFF} -f ${REF_GENOME} -o REF_exons.fasta -t cds

cd REF_PEP
extractSeq ../REF_proteins.fasta
cd ../REF_cDNA
extractSeq ../REF_cDNA.fasta
cd ../REF_EXONS
extractSeq ../REF_exons.fasta
cd ..
