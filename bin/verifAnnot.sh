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
#echo ok
#echo 1
info=$(cat $1)
#echo $info
#echo 2
#echo $2
#echo 3
#echo $3
#infoLocus=$(echo "${info}" | cut -f4)
infoLocus=$(cat "$1" | cut -f4)
#echo $(cat "$1" | cut -f1)
#echo infoLocus
echo $infoLocus
#infoLocus=$(cat "$infoLocus")
SPECIES=$(cat "$1" | cut -f1)
#echo espece 
#echo $SPECIES
GFF=$3/Transfert_$SPECIES/annotation_transfert_${SPECIES}.gff
#echo le GFF
#echo $GFF
#echo genome
GENOME=$2
#echo $GENOME
SCRIPT=$3/SCRIPT
#echo script
echo $SCRIPT


while read line
do
    echo "OSJnip_$line" >> output.txt
done < $infoLocus

#Ajout comment : famille gene Nip, classe gene Nip, +autre
gawk -F"\t" 'BEGIN{OFS="\t"}{
    if(NR==FNR){
        F[$1]=$2;C[$1]=$3}
    else{
        if($3~/gene/){
            split($9,T,/[;/]/);origin=substr(T[2],16);gsub(" ","",origin);$9=$9" / Gene-Fam="F[origin]" / Gene-Class="C[origin]};print}}' output.txt $GFF > LRRlocus_in_${SPECIES}_complet.gff

echo ----------LRRlocus_in_${SPECIES}_complet.gff
head -n 5 LRRlocus_in_${SPECIES}_complet.gff
gawk 'BEGIN{OFS=";"}{if($3~/gene/){if(line){print(line)};split($9,T,";");line=substr(T[1],4)";"$7}else{if($3=="CDS"){line=line";"$4";"$5}}}END{print(line)}' LRRlocus_in_${SPECIES}_complet.gff > geneModel_${SPECIES}.tbl

echo ---------------------geneModel_${SPECIES}.tbl
head -n 5 geneModel_${SPECIES}.tbl


python3 $SCRIPT/Canonical_gene_model_test.py -f $GENOME -t geneModel_${SPECIES}.tbl > alert.txt

echo ----------------------alert.txt
head -n 5 alert.txt
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
                  
echo --------------tmp
head -n 5 tmp

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

echo -----------LRRlocus_in_${SPECIES}_complet.gff
head -n 5 LRRlocus_in_${SPECIES}_complet.gff
cat LRRlocus_in_${SPECIES}_complet.gff > LRRlocus_complet
cat LRRlocus_in_${SPECIES}_complet.gff  >> $3/Transfert_$SPECIES/lecomplet.gff


#========================================================
#                Debut du script
#========================================================
# passage au format table

gawk 'BEGIN{OFS=";"}{if($3~/gene/){if(line){print(line)};split($9,T,"=");line=T[2]";"$7}else{if($3=="CDS"){line=line";"$4";"$5}}}END{print(line)}' LRRlocus_in_${SPECIES}_complet.gff > geneModel_${SPECIES}.tbl
echo "----------------------geneModel_${SPECIES}.tbl------------------------"
head -n 5 geneModel_${SPECIES}.tbl
# Controle des modeles de gene
python3 $SCRIPT/Canonical_gene_model_test.py -f $GENOME -t geneModel_${SPECIES}.tbl > ${SPECIES}_alert.txt
echo "----------------${SPECIES}_alert.txt-------------------"
head -n 5 ${SPECIES}_alert.txt
#frameshift?
gawk '{if($3=="gene"){end=0}else{if($3=="CDS"){if(end==0){end=$5}else{if($4<(end+25)){split($9,T,/[=:]/);print(T[2])}}}}}' LRRlocus_in_${SPECIES}_complet.gff | sort -u > frameshift.txt
echo "-----------------------frameshift.txt---------------------"
head -n 5 frameshift.txt
## Canonic/non-canonique

gawk '{if(NR==FNR){F[$1]=1}else{if(F[$2]==1){$5="True";$7="notValid"};print}}' frameshift.txt ${SPECIES}_alert.txt > tmp 
echo "----------------tmp--------------"
head -n 5 tmp
mv tmp ${SPECIES}_alert.txt
echo "------------------${SPECIES}_alert.txt-------------------"
head -n 5 ${SPECIES}_alert.txt

cat LRRlocus_in_${SPECIES}_complet.gff 

