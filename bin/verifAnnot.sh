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
echo ------------
cat $4
info=$(cat $1)
echo info 
echo $info
infoLocus=$(cat "$1" | cut -f4)
echo infoLocus
echo $infoLocus
SPECIES=$(cat "$1" | cut -f1)
echo specie 
echo $specie
GFF=$3/Transfert_$SPECIES/annotation_transfert_${SPECIES}.gff
echo GFF
echo $GFF
GENOME=$2
echo genome
echo $GENOME
SCRIPT='/GeneModelTransfer.git/branches/container/SCRIPT/'

while read line
do
    if [ ${line:0:1} == 'C' ]
    then 
    echo "$line"
    echo "OSJnip_$line" >> output.txt
    else 
        echo "$line"
        echo "$line" >> output.txt
    fi
done < $infoLocus

cat output.txt
#Ajout comment : famille gene Nip, classe gene Nip, +autre
gawk -F"\t" 'BEGIN{OFS="\t"}{
    if(NR==FNR){
        F[$1]=$2;C[$1]=$3}
    else{
        if($3~/gene/){
            split($9,T,/[;/]/);origin=substr(T[2],16);gsub(" ","",origin);$9=$9" / Gene-Fam="F[origin]" / Gene-Class="C[origin]};print}}' output.txt $GFF > LRRlocus_in_${SPECIES}_complet.gff

echo complet
cat LRRlocus_in_${SPECIES}_complet.gff

gawk 'BEGIN{OFS=";"}{if($3~/gene/){if(line){print(line)};split($9,T,";");line=substr(T[1],4)";"$7}else{if($3=="CDS"){line=line";"$4";"$5}}}END{print(line)}' LRRlocus_in_${SPECIES}_complet.gff > geneModel_${SPECIES}.tbl
echo tbl
cat geneModel_${SPECIES}.tbl 
echo "python3 $SCRIPT/Canonical_gene_model_test.py -f $GENOME -t geneModel_${SPECIES}.tbl > alert.txt"
python3 $SCRIPT/Canonical_gene_model_test.py -f $GENOME -t geneModel_${SPECIES}.tbl > alert.txt
echo alert 
cat alert.text
## Couleur des genes bons/pas bons + raison
## Rouge ssi RLP/RLK/NLR et pas D

gawk -F"\t" '{if(NR==FNR){
                if($3=="True"){START[$2]=1;ADD[$2]=ADD[$2]" / noStart"};
                if($4=="True"){STOP[$2]=1;ADD[$2]=ADD[$2]" / noStop"};
                if($5=="True"){OF[$2]=1;ADD[$2]=ADD[$2]" / pbFrameshift"};
                if($3$4$5~/True/){
                  color[$2]=2}
                else{
                  if($6~/True/){color[$2]=10}
                  else{color[$2]=3}}}
              else{
                if($3=="gene"){
                  split($9,T,";");
                  id=substr(T[1],4);
                  if(color[id]==3){
                     print($0";color="color[id])}
                  else{
                     if(color[id]==10 && ($9!~/ident:100/ || $9!~/cov:1/)){color[id]=2};
                         if(($9~/Fam=RLP/ || $9~/Fam=RLK/ || $9~/Fam=NLR/) && $9!~/Class=D/){
                            print($0""ADD[id]";color="color[id])}
                         else{
                            print($0""ADD[id]";color=3")}}}
                  else{print}}}' alert.txt LRRlocus_in_${SPECIES}_complet.gff > tmp
echo tmp 
cat tmp 
#cat alert.txt > /Users/thibaudvicat/pipelinegit/version3old/test/alert.txt
#cat tmp > /Users/thibaudvicat/pipelinegit/version3old/test/tmp
#cat LRRlocus_in_${SPECIES}_complet.gff > /Users/thibaudvicat/pipelinegit/version3old/test/LRRlocus_in_${SPECIES}_complet.gff
                  

gawk 'BEGIN{OFS="\t";p=0}{
  if($3~/CDS/){
    if(p==0){
      line=$0;P4=$4;P5=$5;p=1}
    else{
      if($4<=P5+25 && ($4-P5-1)%3==0){
        $4=P4;line=$0;P4=$4;P5=$5}
      else{print(line);line=$0;P4=$4;P5=$5}
    }
  }else{
    if(p!=0){print(line)};P4=0;P5=0;p=0;print}
}END{if(p!=0){print(line)}}' tmp > LRRlocus_in_${SPECIES}_complet.gff
echo complet2
cat LRRlocus_in_${SPECIES}_complet.gff
#cat LRRlocus_in_${SPECIES}_complet.gff > /Users/thibaudvicat/pipelinegit/version3old/test/LRRlocus_in_${SPECIES}_complet.gff2
#========================================================
#                Debut du script
#========================================================
# passage au format table
echo genetbl
cat geneModel_${SPECIES}.tbl
gawk 'BEGIN{OFS=";"}{if($3~/gene/){if(line){print(line)};split($9,T,"=");line=T[2]";"$7}else{if($3=="CDS"){line=line";"$4";"$5}}}END{print(line)}' LRRlocus_in_${SPECIES}_complet.gff > geneModel_${SPECIES}.tbl
#cat geneModel_${SPECIES}.tbl > /Users/thibaudvicat/pipelinegit/version3old/test/geneModel_${SPECIES}.tbl
# Controle des modeles de gene
echo speciealert
cat ${SPECIES}_alert.txt
echo "python3 $SCRIPT/Canonical_gene_model_test.py -f $GENOME -t geneModel_${SPECIES}.tbl > ${SPECIES}_alert.txt"
python3 $SCRIPT/Canonical_gene_model_test.py -f $GENOME -t geneModel_${SPECIES}.tbl > ${SPECIES}_alert.txt
#cat ${SPECIES}_alert.txt > /Users/thibaudvicat/pipelinegit/version3old/test/${SPECIES}_alert.txt
#frameshift?

gawk '{if($3=="gene"){end=0}else{if($3=="CDS"){if(end==0){end=$5}else{if($4<(end+25)){split($9,T,/[=:]/);print(T[2])}}}}}' LRRlocus_in_${SPECIES}_complet.gff | sort -u > frameshift.txt
echo frameshift
cat frameshift.txt 
## Canonic/non-canonique
gawk '{if(NR==FNR){F[$1]=1}else{if(F[$2]==1){$5="True";$7="notValid"};print}}' frameshift.txt ${SPECIES}_alert.txt > tmp 
mv tmp ${SPECIES}_alert.txt
cat LRRlocus_in_${SPECIES}_complet.gff > $3/Transfert_$SPECIES/LRRlocus_in_${SPECIES}_acurate.gff 
