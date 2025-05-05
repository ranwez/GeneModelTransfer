#!/bin/bash

set -euo pipefail

usage() {
    echo "Usage: $0 <input_GFF_folder> <input_jobs_folder> <output_stats_folder>"
    echo
    echo "This script computes a few statistics on LRRtransfer gff output and SLURM jobs log files (.err and .out) for the genePrediction step."
    echo
    echo "Arguments:"
    echo "  <input_GFF_folder>    Path to the folder containing the annot_best.gff produced by LRRtransfer."
    echo "  <input_jobs_folder>    Path to the folder containing the SLURM jobs log files (.err and .out)."
    echo "  <output_stats_folder>    Path to the output stat file folder."
}

printGFFstats(){
  local gff=$1
  local out_file=$2
  echo -n "Total number of genes: " > $out_file
  cut -f3 $gff | grep "gene" | uniq -c >> $out_file
  echo "Number of genes per prediction method:" >> $out_file
  grep gene $gff | sed -e 's/.*pred://' | cut -f1 -d ' ' | sort | uniq -c >> $out_file
  echo "Number of genes per chromsome:" >> $out_file
  grep gene $gff | cut -f1 | sort | uniq -c >> $out_file

}

_printJobStats1(){
  input_dir=$1
  ext=$2
  out_file=$3
  wc -l ${input_dir}/LRRtransfer.genePrediction.*${ext} | awk '{print $1}' | sort | uniq -c | sort -nr > ${out_file}
  nb_expected=$(awk '{if (NR ==1){print $2}}' ${out_file})
  echo -e "\nDetail: file with a number of lines different from the mod (${nb_expected})" >> ${out_file}
  wc -l ${input_dir}/LRRtransfer.genePrediction.*${ext} | grep -v ${nb_expected} >> ${out_file}
}

printJobStats(){
  local input_dir=$1
  local out_file_prefix=$2
  _printJobStats1 ${input_dir} "out" "${out_file_prefix}_out.txt"
  _printJobStats1 ${input_dir} "err" "${out_file_prefix}_err.txt"
}


## Main ##

# Check the number of arguments
if [ "$#" -ne 3 ]; then
    usage
    exit 1
fi


input_GFF_folder=$1
input_jobs_folder=$2
output_LRRtransfer_stats=$3
mkdir -p ${output_LRRtransfer_stats}
printGFFstats ${input_GFF_folder}/annot_best_chr.gff ${output_LRRtransfer_stats}/GFFstats.txt
printJobStats ${input_jobs_folder} ${output_LRRtransfer_stats}/jobsStats
