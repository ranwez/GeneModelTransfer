#!/bin/bash
#========================================================
# PROJET : LRRtransfer
# SCRIPT : create_LRRome.sh
# AUTHOR : Celine Gottin & Thibaud Vicat
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
# ARGUMENTS : o $1 : Path to a text file with 4 columns :
#                    First column contain a code the accession.
#                    Second column contain a path to the reference GFF containing LRR 
#                    Third column contain a path to the referene asembly (fasta format)
#                    Fourth column is not obligatory and should contain a path to a file containing information for LRR (family and class of each location)
#			  o $2 : Path to LRRome if one already exist
#			  o $2 : Launch directory
# DEPENDENCIES : o python3
#========================================================
#========================================================
#                Environment & variables
#========================================================
INFO_FILE=$1
LRRome=$2
LAUNCH_DIR=$3
SCRIPT='/SCRIPT/'
#========================================================
#                Script
#========================================================
function extractSeq {
	##Extracting each sequence from a fasta in separate files
	gawk -F"[;]" '{if($1~/>/){line=$1;gsub(">","");filename=$1;print(line) > filename}else{print > filename}}' $1
}
if [ $INFO_FILE != 'NULL' ] && [ $LRRome == 'NULL' ]
    then
		echo $INFO_FILE
		mkdir -p LRRome
		cd LRRome
		mkdir -p REF_PEP
		mkdir -p REF_CDS
		mkdir -p REF_cDNA
		while read line
		do
			echo $INFO_FILE
			code=$(echo "${line}" | cut -f1)
			echo $code
			mkdir -p $3/Transfert_$code
			path_gff=$(echo "${line}" | cut -f2)
			echo $path_gff
			path_gff=$(realpath "$path_gff")
			echo $path_gff
			path_fasta=$(echo "${line}" | cut -f3)
			echo $path_fasta
			path_fasta=$(realpath $path_fasta)
			echo $path_fasta
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
