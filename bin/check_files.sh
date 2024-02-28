#!/bin/bash

set -euo pipefail

target_genome=$1
ref_genome=$2
ref_gff=$3
ref_locus_info=$4
lrrome=$5
OUTDIR=$6

outsummary=$7

#========================================================
#                        FUNCTIONS
#========================================================
# $1 parameter allows to specify a prefix to identify your tmp folders
function get_tmp_dir(){
  local tmp_dir; tmp_dir=$(mktemp -d -t "$1"_$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXXXX)
  echo $tmp_dir
}

# in debug mode ($1=1), do not delete the temporary directory passed as $2
function clean_tmp_dir(){
  if (( $1==0 )); then
    rm -rf "$2"
  fi
}

function quit_pb_option() {
    printf "\n\nThe config.yml should contain 6 parameters.\n"
    printf "\ntarget_genome: absolute path of fasta file containing the sequence of the genome to annotate"
    printf "\nref_genome: absolute path of fasta file contain the reference genome sequence"
    printf "\nref_gff: absolute path of gff file containing gene annotations of the reference genome"
    printf "\nref_info_locus: absolute path of tab separated file containing for each reference gene to transfer: its identifier, its family, and its class"
    printf "\nlrrome: empty string or absolute path of an already build LRRome folder."
    printf "\nOUTPUTS_DIRNAME: absolute path of output directory that will contain the pipeline results"
    printf "\nyour config file contain incorrect parameters. Please check them."
    printf "\n\nFor further details please check the documentation on pipeline git page: https://github.com/cgottin/GeneModelTransfer/\n\n"
    exit 1
}


function check_in_file_param(){
  local file; file=$2
  local has_problem=0
  if [[ ! -f $file ]]; then
    printf "\nProblem with option $1, File $2 does not exist\n" >&2
    has_problem=1
  fi
  return $has_problem
}


function check_out_dir_param(){
  local dir=$1
  local has_problem=0
  if [ -e  $dir ]; then
    printf "\nWarning: output Directory $1 already  exist\n" >&2
  fi
  return $has_problem
}


function check_mode(){
  local has_problem=0
  if [[ ! "$1" =~ ^(first|best)$ ]];then
    printf "\nProblem with option mode, $1 is not handle, please choose between first and best\n" >&2
    has_problem=1
  fi
  return $has_problem
}


function check_fasta() {
  local has_problem=0
  local nbseq=$(grep -c "^>" $1)
  if (( $nbseq==0 ));then
    has_problem=1
    printf "\nProblem: no sequence found in fasta file $1\n" >&2
  fi
  echo $nbseq
  return $has_problem
}







check_out_dir_param $OUTDIR
check_in_file_param target_genome $target_genome || quit_pb_option
check_in_file_param ref_genome $ref_genome || quit_pb_option
check_in_file_param ref_gff $ref_gff || quit_pb_option
check_in_file_param ref_locus_info $ref_locus_info || quit_pb_option
#check_mode $mode || quit_pb_option



tmpdir=$(get_tmp_dir LRRtransfer_check)


## Control of fasta file
nbseq_target_genome=$(check_fasta $target_genome) || quit_pb_option
nbseq_ref_genome=$(check_fasta $ref_genome) || quit_pb_option

## control of gff
nb_gene=$(cut -f 3 $ref_gff | grep -c "gene")
nb_cds=$(cut -f 3 $ref_gff | grep -c "CDS")
nb_exon=$(cut -f 3 $ref_gff | grep -c "exon")

pb_gff=0
if (( $nb_gene==0 ));then pb_gff=1; printf "no gene feature found in gff file $ref_gff" >&2;fi
if (( $nb_cds==0 ));then pb_gff=1; printf "no CDS feature found in gff file $ref_gff" >&2;fi
if (( $nb_exon==0 ));then pb_gff=1; printf "no exon feature found in gff file $ref_gff" >&2;fi

if (( pb_gff==1 ));then quit_pb_option;fi

## control of chromosome between reference fasta and gff
grep "^>" $ref_genome | cut -d " " -f 1 | sed 's/>//g' | uniq > $tmpdir/chrnames.tmp
cut -f1 $ref_gff | uniq > $tmpdir/chrnamesgff.tmp

pb_chr=$(gawk 'BEGIN{pb=0; not_found=""}{if(NR==FNR){CHR[$1]=1}else{if(! CHR[$1]){if(pb==0){not_found=$1};pb=1}}}END{print(not_found)}'  $tmpdir/chrnames.tmp $tmpdir/chrnamesgff.tmp)
if [ ! -z $pb_chr ];then printf "\nchromosome $pb_chr from the reference gff ($ref_gff) file is missing in the reference fasta file ($ref_genome)\n";quit_pb_option;fi

## check columns of ref_locus_info file
cut -f3,9 $ref_gff | grep "^gene" | cut -f2 | cut -d";" -f1 | sed 's/ID=//g' > $tmpdir/identgenegff.tmp

pb_gene=$(gawk 'BEGIN{pb=0; not_found=""}{if(NR==FNR){gene[$1]=1}else{if(! gene[$1]){if(pb==0){not_found=$1};pb=1}}}END{print(not_found)}' $ref_locus_info $tmpdir/identgenegff.tmp)

if [ ! -z $pb_gene ];then printf "\nwarning: some genes (e.g. $pb_gene) from the reference gff ($ref_gff) are missing from the info locus file ($ref_locus_info)\n" >&2;fi


##log
echo "Input parameters where checked without error." > $outsummary
echo "fasta file $target_genome contains $nbseq_target_genome sequences" >> $outsummary
echo "fasta file $ref_genome contains $nbseq_ref_genome sequences" >> $outsummary
