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


          #------------------------------------------#
          # 1. Find Regions of interest with mmseqs2 #
          #------------------------------------------#

tmpdir=$(get_tmp_dir LRRtransfer_candidateLoci)
cd $tmpdir

filename=$(basename ${TARGET_GENOME%.fasta})

## using mmseqs
#$mmseqs createdb $TARGET_GENOME ${filename}_db -v 0
#$mmseqs createdb $LRRome/REF_proteins.fasta prot_db  -v 0 
#$mmseqs search prot_db ${filename}_db resultDB_aln.m8 tmp -s 8.5 -a -e 0.1 --min-length 10 --merge-query 1 --cov-mode 2 --max-seqs 30000 --sequence-overlap 1000 -v 0
#$mmseqs convertalis prot_db ${filename}_db resultDB_aln.m8 res_candidatsLRR.out --format-output query,target,qlen,alnlen,qstart,qend,tstart,tend,nident,pident,gapopen,evalue,bits  -v 0

## using blast+
#makeblastdb -in $TARGET_GENOME -out ${filename}_db -dbtype nucl
#tblastn -db ${filename}_db -query $LRRome/REF_proteins.fasta -evalue 1 -out res_candidatsLRR.out -outfmt "6 qseqid sseqid qlen length qstart qend sstart send nident pident gapopen evalue bitscore"

cp ${blast_res} res_candidatsLRR.out


          #------------------------------------------#
          # 2. Extract and concat hits               #
          #------------------------------------------#

## Build up global hits by proteins in 2 steps: 
## 1st round with high thresholds to fix "anchor" exons
## 2nd round with lower ths to complete the areas at the ends and define other paralogous type areas
sort -k1,2 -Vk7,7 res_candidatsLRR.out > tmp.tmp
mv tmp.tmp res_candidatsLRR.out
python ${LRR_SCRIPT}/filter_res_align.py -g $REF_GFF -t res_candidatsLRR.out > concat_candidatsLRR.out



          #------------------------------------------#
          # 3. format Candidate regions              #
          #------------------------------------------#

sort -k2,2 -Vk7,7 concat_candidatsLRR.out > tmp.tmp
mv tmp.tmp concat_candidatsLRR.out
python ${LRR_SCRIPT}/create_candidate_from_align.py -t concat_candidatsLRR.out -o filtered_candidatsLRR.gff > list_query_target.txt



          #------------------------------------------#
          # 5. Export files                          #
          #------------------------------------------#

python $LRR_SCRIPT/Extract_sequences_from_genome.py -f $TARGET_GENOME -g filtered_candidatsLRR.gff -o DNA_candidatsLRR.fasta  -t gene 

mkdir CANDIDATE_SEQ_DNA 
cd CANDIDATE_SEQ_DNA
extractSeq ../DNA_candidatsLRR.fasta
cd .. 

#update list_query_target
blastp ../DNA_candidatsLRR.fasta $LRRome/REF_proteins.fasta -outfmt 6 "qseqid sseqid qframe bitscore"

#saving files
cp list_query_target.txt $RES_DIR/.
cp filtered_candidatsLRR.gff $RES_DIR/.
cp -r CANDIDATE_SEQ_DNA $RES_DIR/.

echo $tmpdir
clean_tmp_dir 1 $tmpdir
