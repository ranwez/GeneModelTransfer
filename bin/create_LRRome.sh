#!/bin/bash
#========================================================
# PROJET : LRRtransfer
# SCRIPT : create_LRRome.sh
# AUTHOR : Celine Gottin & Thibaud Vicat
# CREATION : 2021.05.07
#========================================================
# DESCRIPTION : Extract fasta files for LRR loci from all
#               given species and then build a LRRome.
#               The process need path of gff and genomic
#               fasta files for each species. 
#               Paths are read from tab separated file with 
#               one line per species.
#               If an LRRome is given as input the process
#               copy data to the working directory.
# ARGUMENTS : o $1 : Path to the reference assembly
#             o $2 : Path to the reference GFF
#             o $3 : Launch directory
#             o $4 : Path to LRRome if one already exist
# DEPENDENCIES : o python3
#========================================================


#========================================================
#                Environment & variables
#========================================================
REF_GENOME=$(readlink -f "$1")
echo $1
echo $REF_GENOME
REF_GFF=$(readlink -f "$2")
echo $2
echo $REF_GFF
RES_DIR=$(readlink -f "$3")
echo $3
echo $RES_DIR
PREBUILT_LRRome=$4
echo $4

echo $PWD
#========================================================
#                        Functions
#========================================================

function extractSeq {
	##usage :: extractSeq multifasta.file
	##Extracting each sequence from a fasta in separate files
	gawk -F"[;]" '{if($1~/>/){line=$1;gsub(">","");filename=$1;print(line) > filename}else{print > filename}}' $1
}

export -f extractSeq


#========================================================
#                Script
#========================================================

mkdir $RES_DIR/LRRome
cd $RES_DIR/LRRome

if [[ $REF_GENOME != 'NULL' ]] && [[ $REF_GFF != 'NULL' ]] && [[ $PREBUILT_LRRome == 'NULL' ]];then

	mkdir -p REF_PEP
	mkdir -p REF_EXONS
	mkdir -p REF_cDNA


	python3 ${LG_SCRIPT}/Extract_sequences_from_genome.py -g ${REF_GFF} -f ${REF_GENOME} -o REF_proteins.fasta -t prot
	python3 ${LG_SCRIPT}/Extract_sequences_from_genome.py -g ${REF_GFF} -f ${REF_GENOME} -o REF_cDNA.fasta -t cdna
	python3 ${LG_SCRIPT}/Extract_sequences_from_genome.py -g ${REF_GFF} -f ${REF_GENOME} -o REF_exons.fasta -t exon

	cd REF_PEP
	extractSeq ../REF_proteins.fasta
	cd ../REF_cDNA
	extractSeq ../REF_cDNA.fasta
	cd ../REF_EXONS
	extractSeq ../REF_exons.fasta
	cd ..

elif [ $PREBUILT_LRRome != 'NULL' ];then
	$PREBUILT_LRRome=$(readlink -f "$4")
	cp -r $PREBUILT_LRRome/* $RES_DIR/LRRome/

fi


