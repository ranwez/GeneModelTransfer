#!/bin/bash
#========================================================
# PROJET : 
# SCRIPT : 
# AUTHOR : Céline Gottin & Thibaud Vicat
# CREATION : 2020.02.20
#========================================================
# DESCRIPTION : Check the new gene models
#               annotate for control (presence of start, stop
#               canonical intron, non-overlapping frmaeshift)
# ARGUMENTS : o $1 : Path to a text file with 4 columns :
#                    First column contain a code the accession.
#                    Second column contain a path to the reference GFF containing LRR 
#                    Third column contain a path to the referene asembly (fasta format)
#                    Fourth column is not obligatory and should contain a path to a file containing information for LRR (family and class of each location)
#             o $2 : Target genome
#             o $3 : Launch directory
#========================================================
#                Environment & variables
#========================================================
info=$(cat $1)
infoLocus=$(cat "$1" | cut -f4)
SPECIES=$(cat "$1" | cut -f1)
GFF=$3/Transfert_$SPECIES/annotation_transfert_${SPECIES}.gff
GENOME=$2
SCRIPT='/GeneModelTransfer.git/branches/container/SCRIPT'
#========================================================
#                Beginning of the script
#========================================================
#correction of infoLocus identifiers
while read line
do
    if [ ${line:0:1} == 'C' ]
    then 
    echo "OSJnip_$line" >> output.txt
    else 
        echo "$line" >> output.txt
    fi
done < $infoLocus
#Add comment : Nip gene family, Nip gene class, +other
gawk -F"\t" 'BEGIN{OFS="\t"}{
    if(NR==FNR){
        F[$1]=$2;C[$1]=$3}
    else{
        if($3~/gene/){
            split($9,T,/[;/]/);origin=substr(T[2],16);gsub(" ","",origin);$9=$9" / Gene-Fam="F[origin]" / Gene-Class="C[origin]};print}}' output.txt $GFF > LRRlocus_in_${SPECIES}_complet.gff
gawk 'BEGIN{OFS=";"}{if($3~/gene/){if(line){print(line)};split($9,T,";");line=substr(T[1],4)";"$7}else{if($3=="CDS"){line=line";"$4";"$5}}}END{print(line)}' LRRlocus_in_${SPECIES}_complet.gff > geneModel_${SPECIES}.tbl
python3 $SCRIPT/Canonical_gene_model_test.py -f $GENOME -t geneModel_${SPECIES}.tbl > alert.txt
## Color of genes good/not good + reason
## Red if RLP/RLK/NLR and not D
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
# switch to table format
gawk 'BEGIN{OFS=";"}{if($3~/gene/){if(line){print(line)};split($9,T,"=");line=T[2]";"$7}else{if($3=="CDS"){line=line";"$4";"$5}}}END{print(line)}' LRRlocus_in_${SPECIES}_complet.gff > geneModel_${SPECIES}.tbl
# Control of gene models
python3 $SCRIPT/Canonical_gene_model_test.py -f $GENOME -t geneModel_${SPECIES}.tbl > ${SPECIES}_alert.txt
#frameshift?
gawk '{if($3=="gene"){end=0}else{if($3=="CDS"){if(end==0){end=$5}else{if($4<(end+25)){split($9,T,/[=:]/);print(T[2])}}}}}' LRRlocus_in_${SPECIES}_complet.gff | sort -u > frameshift.txt
## Canonic/non-canonic
gawk '{if(NR==FNR){F[$1]=1}else{if(F[$2]==1){$5="True";$7="notValid"};print}}' frameshift.txt ${SPECIES}_alert.txt > tmp 
mv tmp ${SPECIES}_alert.txt
cat LRRlocus_in_${SPECIES}_complet.gff > $3/Transfert_$SPECIES/LRRlocus_in_${SPECIES}_acurate.gff 
if $4 == 0
then 
rm -r $3/work
fi