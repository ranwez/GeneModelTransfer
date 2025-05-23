#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 <input_fasta_file> <output_directory>"
    echo
    echo "This script splits a fasta file into individual files for each sequence."
    echo
    echo "Arguments:"
    echo "  <input_fasta_file>    Path to the input fasta file."
    echo "  <output_directory>    Directory to save the individual fasta files."
}

# Check the number of arguments
if [ "$#" -ne 2 ]; then
    usage
    exit 1
fi

# Assign arguments to variables
input_fasta_file=$1
output_directory=$2

# Ensure output directory exists
mkdir -p "$output_directory"

# AWK script to split fasta file
awk -v output_dir="$output_directory" '
BEGIN {
    seqname = ""
    seq = ""
}

/^>/ {
    if (seqname != "") {
        filename = output_dir "/" seqname ".fasta"
        print ">" seqname > filename
        print seq >> filename
        close(filename)
    }
    seqname = substr($1, 2)
    seq = ""
    next
}

{
    seq = seq $0
}

END {
    if (seqname != "") {
        filename = output_dir "/" seqname ".fasta"
        print ">" seqname > filename
        print seq >> filename
    }
}
' "$input_fasta_file"