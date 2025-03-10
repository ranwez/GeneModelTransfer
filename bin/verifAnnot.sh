#!/bin/bash
#========================================================
# PROJET : lrrtransfer
# SCRIPT : verifAnnot.sh
# AUTHOR : Celine Gottin & Thibaud Vicat & Vincent Ranwez
# CREATION : 2020.02.20
#========================================================
# DESCRIPTION : Check the new gene models
#               annotate for control (presence of start, stop
#               canonical intron, non-overlapping frmaeshift)
# ARGUMENTS : o $1 : Path to a text file with ref locus info
#             o $2 : Target genome
#             o $3 : predicted annotation GFF
#             o $4 : Result directory
#             o $5 : Path toward LRR script  directory
#             o $6 : name of the transfert method used (to differentiate final results) 
#========================================================


infoLocus=$1
dna_seq_file=$2
input_gff=$3

RES_DIR=$4
LRR_SCRIPT=$5
method=$6


source $LRR_SCRIPT/../bin/lib_gff_comment.sh
source $LRR_SCRIPT/../bin/lib_tmp_dir.sh

tmpdir=$(get_tmp_dir LRRtransfer_verifAnnot)
cd $tmpdir

add_origin_info ${input_gff} ${infoLocus} ${input_gff}_with_origin
compute_NC_alerts ${input_gff}_with_origin ${dna_seq} ${input_gff}_NC_alert
add_comment_NC ${input_gff}_with_origin ${input_gff}_NC_alert ${input_gff}_with_origin_comment

## export output files
cat ${input_gff}_with_origin_comment > $RES_DIR/LRRlocus_${method}_predicted.gff
cat ${input_gff}_NC_alert | sort -k1,2 > $RES_DIR/alert_${method}_NC_Locus.txt

clean_tmp_dir 0 $tmpdir
