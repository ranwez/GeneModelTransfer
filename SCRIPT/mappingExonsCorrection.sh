#!/bin/bash
SCRIPT=$1
echo $SCRIPT
BLASTDB=$2
echo $BLASTDB
liste_query_target=$3
echo $liste_query_target
SPECIES=$4
echo $SPECIES
##corriger pos gene, lancer correction
gawk -F"\t" 'BEGIN{OFS="\t"}{split($9,T,/[=:;]/);if(NR==FNR){if($3=="gene"){max[T[2]]=$5;min[T[2]]=$5}else{if($5>max[T[2]]){max[T[2]]=$5};if($4<min[T[2]]){min[T[2]]=$4}}}else{if($3=="gene"){$4=min[T[2]];$5=max[T[2]]};print}}' mappingCDS_$SPECIES.gff mappingCDS_$SPECIES.gff > tmp ; mv tmp mappingCDS_$SPECIES.gff
cd ..

python $SCRIPT/Exonerate_correction.py -f $BLASTDB -g $SCRIPT/mappingCDS_$SPECIES.gff > $SCRIPT/mapping_LRRlocus_${SPECIES}.gff
#cat mapping_LRRlocus_${SPECIES}.gff > $SCRIPT/mapping_LRRlocus_${SPECIES}.gff 
##to transfert with cdna2genome
gawk 'BEGIN{OFS="\t"}{if(NR==FNR){if($3=="gene"){split($9,T,";");gsub("ID=","",T[1]);OK[T[1]]=1}}else{if(!OK[$1]==1){print}}}' $SCRIPT/mapping_LRRlocus_$SPECIES.gff $liste_query_target > $SCRIPT/to_transfer_with_cdna.txt

echo "END PART3"