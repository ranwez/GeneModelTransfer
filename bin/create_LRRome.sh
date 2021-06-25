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
echo $1
LRRome=$2
echo $2
LAUNCH_DIR=$3
echo $3
SCRIPT='/GeneModelTransfer.git/branches/container/SCRIPT/'
echo "$SCRIPT"
head $SCRIPT/Extract_sequences_from_genome.py
