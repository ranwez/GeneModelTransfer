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
# ARGUMENTS : o $1 : Path to the input file (tab separated file)
#             o $2 : Path to LRRome if one already exist
#             o $3 : Launch directory
# DEPENDENCIES : o python3
#========================================================


#========================================================
#                Environment & variables
#========================================================
INFO_FILE=$1
LRRome=$2
LAUNCH_DIR=$3


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

mkdir LRRome
cd LRRome

if [ $INFO_FILE != 'NULL' ] && [ $LRRome == 'NULL' ];then

	mkdir -p REF_PEP
	mkdir -p REF_EXONS
	mkdir -p REF_cDNA

	while read line
	do
		IFS='\t' read code path_gff path_fasta info_locus <<< "$line" ;
		python3 $LG_SCRIPT/Extract_sequences_from_genome.py -g ${path_gff} -f ${path_fasta} -o ${code}_proteins.fasta -t prot
		python3 $LG_SCRIPT/Extract_sequences_from_genome.py -g ${path_gff} -f ${path_fasta} -o ${code}_cDNA.fasta -t cdna
		python3 $LG_SCRIPT/Extract_sequences_from_genome.py -g ${path_gff} -f ${path_fasta} -o ${code}_exons.fasta -t exon
		cd REF_PEP
		extractSeq ../${code}_proteins.fasta
		cd ../REF_cDNA
		extractSeq ../${code}_cDNA.fasta
		cd ../REF_EXONS
		extractSeq ../${code}_exons.fasta
		cd ../
	done < $INFO_FILE

elif [ $LRRome != 'NULL' ];then

	cp -r $LRRome/* ./

fi
