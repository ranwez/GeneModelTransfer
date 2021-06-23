#!/bin/bash
#========================================================
# PROJET : TransfertGeneModel
# SCRIPT : create_LRRome.sh
# AUTHOR : Celine Gottin
# CREATION : 2021.05.07
#========================================================
# DESCRIPTION : Extract fasta files for LRR loci from all
#               given species and then build a LRRome.
#				The process need path of gff
#               and genomic fasta files for each species. 
#               Paths are read from spaced text file with 
#               one line per species.
#				If an LRRome is given as input the process
#               returns this LRRome 
# ARGUMENTS : o $1 : text file with gff and genome fasta path
# DEPENDENCIES : o python3
#========================================================


#========================================================
#                Environment & variables
#========================================================

INFO_FILE=$1
#echo $1
LRRome=$2
#echo $2
LAUNCH_DIR=$3
#echo $3
SCRIPT='/GeneModelTransfer.git/branches/container/SCRIPT/'
#echo "$SCRIPT"
#head $SCRIPT/Extract_sequences_from_genome.py
#========================================================
#                Script
#========================================================
function extractSeq {
	##Extracting each sequence from a fasta in separate files
	gawk -F"[;]" '{if($1~/>/){line=$1;gsub(">","");filename=$1;print(line) > filename}else{print > filename}}' $1
}
if [ $INFO_FILE != 'NULL' ] && [ $LRRome == 'NULL' ]
    then
		mkdir -p LRRome
		cd LRRome
		mkdir -p REF_PEP
		mkdir -p REF_CDS
		mkdir -p REF_cDNA
		while read line
		do
			code=$(echo "${line}" | cut -f1)
			#echo "code"
			#echo "$code"
			mkdir -p $3/Transfert_$code
			path_gff=$(echo "${line}" | cut -f2)
			#echo "-----------------path_gff"
			#head ${path_gff}
			#echo "${path_gff}"
			#echo "----------------------path_fasta"
			path_fasta=$(echo "${line}" | cut -f3)
			#head ${path_fasta}
			#echo "${path_fasta}"

			python3 $SCRIPT/Extract_sequences_from_genome.py -g ${path_gff} -f ${path_fasta} -o ${code}_proteins.fasta -t prot
			cd REF_PEP
			extractSeq ../${code}_proteins.fasta
			cd ../
			python3 $SCRIPT/Extract_sequences_from_genome.py -g ${path_gff} -f ${path_fasta} -o ${code}_cDNA.fasta -t cdna
			cd REF_cDNA
			extractSeq ../${code}_cDNA.fasta
			cd ../
			python3 $SCRIPT/Extract_sequences_from_genome.py -g ${path_gff} -f ${path_fasta} -o ${code}_exons.fasta -t exon
			cd REF_CDS
			extractSeq ../${code}_exons.fasta
			cd ../
		done < $INFO_FILE
elif [ $LRRome != 'NULL' ]
	then
		mkdir -p LRRome
		cd LRRome
		cp -r $LRRome/* ./
		while read line
		do
			code=$(echo "${line}" | cut -f1)
			mkdir -p $3/Transfert_$code
		done < $INFO_FILE
fi
