#!/bin/bash
#========================================================
# PROJET : 
# SCRIPT : 
# AUTHOR : Céline Gottin
# CREATION : 2020.02.20
#========================================================
# DESCRIPTION : Verifier les modeles de gene nouvellement
#               annoter pour controle (presence du start, stop
#               intron canonique, frmaeshift non chevauchant)
# ARGUMENTS : o $1 : species
#             o $2 : GFF issu du transfert d'annot 
#             o $3 : fasta file with genomic sequences
# DEPENDENCIES :
#========================================================


#========================================================
#                Environment & variables
#========================================================

SPECIES=$1
GFF=$2
GENOME=$3
#SCRIPT=$(pwd)/SCRIPT

#========================================================
#                Debut du script
#========================================================
# passage au format table
#gawk 'BEGIN{OFS=";"}{if($3~/gene/){if(line){print(line)};split($9,T,"=");line=T[2]";"$7}else{if($3=="CDS"){line=line";"$4";"$5}}}END{print(line)}' $GFF > geneModel_${SPECIES}.tbl

# Controle des modeles de gene
python ../SCRIPT/Canonical_gene_model_test.py -f $GENOME -t geneModel_${SPECIES}.tbl > ${SPECIES}_alert.txt

#frameshift?
gawk '{if($3=="gene"){end=0}else{if($3=="CDS"){if(end==0){end=$5}else{if($4<(end+25)){split($9,T,/[=:]/);print(T[2])}}}}}' $GFF | sort -u > frameshift.txt

## Canonic/non-canonique

gawk '{if(NR==FNR){F[$1]=1}else{if(F[$2]==1){$5="True";$7="notValid"};print}}' frameshift.txt ${SPECIES}_alert.txt > tmp 

mv tmp ${SPECIES}_alert.txt




