#!/bin/bash
#========================================================
# PROJET : LRRtransfer
# SCRIPT : candidateLoci.sh
# AUTHOR : Celine Gottin & Thibaud Vicat & Vincent Ranwez
# CREATION : 2020.02.20
#========================================================
# DESCRIPTION : Use of mmseqs to find regions of interest in the
#               target genome from protein sequences contained in LRRome
#               and create query/target pairs returns this LRRome
# ARGUMENTS : o $1 : Path to the Target genome fasta file
#             o $2 : Path to the LRRome directory
#             o $3 : Path to the reference GFF
#             o $4 : Results directory
#             o $5 : Path toward LRR script  directory
# DEPENDENCIES : o python3
#========================================================

# set strict bash mode
set -e -u


#========================================================
#                Environment & variables
#========================================================
TARGET_GENOME=$1
LRRome=$2
REF_GFF=$3
blast_res=$4

CDNA=$LRRome/REF_cDNA
PROTEINS=$LRRome/REF_PEP
EXONS=$LRRome/REF_EXONS

RES_DIR=$5
LRR_SCRIPT=$6
mmseqs="mmseqs"

#========================================================
#                        FUNCTIONS
#========================================================
# $1 parameter allows to specify a prefix to identify your tmp folders
function get_tmp_dir(){
  local tmp_dir;
  tmp_dir=$(mktemp -d -t "$1"_$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXXXX)
  echo $tmp_dir
}

# in debug mode ($1=1), do not delete the temporary directory passed as $2
function clean_tmp_dir(){
  if (( $1==0 )); then
    rm -rf "$2"
  fi
}

function extractSeq {
	##usage :: extractSeq multifasta.file
	##Extracting each sequence from a fasta in separate files
	gawk -F"[;]" '{if($1~/>/){line=$1;gsub(">","");filename=$1;print(line) > filename}else{print > filename}}' $1
}

#========================================================
#                        SCRIPT
#========================================================


tmpdir=$(get_tmp_dir LRRtransfer_candidateLoci)
cd $tmpdir

python ${LRR_SCRIPT}/candidate_loci_VR.py -g $REF_GFF -t ${blast_res} -o filtered_candidatsLRR.gff -l list_query_target.txt


python $LRR_SCRIPT/Extract_sequences_from_genome.py -f $TARGET_GENOME -g filtered_candidatsLRR.gff -o DNA_candidatsLRR.fasta  -t gene

mkdir CANDIDATE_SEQ_DNA
cd CANDIDATE_SEQ_DNA
extractSeq ../DNA_candidatsLRR.fasta
cd ..

#saving files
cp list_query_target.txt $RES_DIR/.
cp filtered_candidatsLRR.gff $RES_DIR/.
cp -r CANDIDATE_SEQ_DNA $RES_DIR/.

echo $tmpdir
clean_tmp_dir 1 $tmpdir
