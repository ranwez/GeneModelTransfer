#!/bin/bash
#========================================================
# PROJET : lrrtransfer
# SCRIPT : verifAnnot.sh
# AUTHOR : Celine Gottin & Thibaud Vicat
# CREATION : 2020.02.20
#========================================================
# DESCRIPTION : Check the new gene models
#               annotate for control (presence of start, stop
#               canonical intron, non-overlapping frmaeshift)
# ARGUMENTS : o $1 : Path to a text file with ref locus info
#             o $2 : Target genome
#             o $3 : predicted annotation GFF
#             o $4 : Result directory
#========================================================

#========================================================
#                Environment & variables
#========================================================


infoLocus=$1
TARGET_GENOME=$2
GFF=$3

RES_DIR=$4

#========================================================
#                        FUNCTIONS
#========================================================

# $1 parameter allows to specify a prefix to identify your tmp folders
function get_tmp_dir(){
  local tmp_dir; tmp_dir=$(mktemp -d -t "$1"_$(date +%Y-%m-%d-%H-%M-%S)-XXXXXXXXXXXX)
  echo $tmp_dir
}

# in debug mode ($1=1), do not delete the temporary directory passed as $2
function clean_tmp_dir(){
  if (( $1==0 )); then
    rm -rf "$2"
  fi
}

#========================================================
#                        SCRIPT
#========================================================

tmpdir=$(get_tmp_dir LRRtransfer_verifAnnot)
cd $tmpdir


#Add comment : Nip gene family, Nip gene class, +other
gawk -F"\t" 'BEGIN{OFS="\t"}{
    if(NR==FNR){
        F[$1]=$2;C[$1]=$3}
    else{
        if($3~/gene/){
            split($9,T,/[;/]/);origin=substr(T[2],16);gsub(" ","",origin);$9=$9" / Origin-Fam="F[origin]" / Origin-Class="C[origin]};print}}' $infoLocus $GFF > LRRlocus_complet.tmp


## concatenate CDS is less than 25 nucl appart and if in the same frame  
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
}END{if(p!=0){print(line)}}' LRRlocus_complet.tmp > LRRlocus_complet2.tmp

gawk 'BEGIN{OFS=";"}{if($3~/gene/){if(line){print(line)};split($9,T,";");line=substr(T[1],4)";"$7}else{if($3=="CDS"){line=line";"$4";"$5}}}END{print(line)}' LRRlocus_complet2.tmp > geneModel.tbl
python3 ${LRR_SCRIPT}/Canonical_gene_model_test.py -f $TARGET_GENOME -t geneModel.tbl -o alert_NC_Locus.tmp


## Color of genes good/not good + reason
# 2=RED ; 10=Orange ; 3=vert
## Red if RLP/RLK/NLR and Non-canonical
gawk -F"\t" 'BEGIN{OFS="\t"}{if(NR==FNR){
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
                     $9=$9";color="color[id]
                     }
                  else{
                     if(color[id]==10 && ($9!~/ident:100/ || $9!~/cov:1/)){color[id]=2};
                         if(($9~/Fam=RLP/ || $9~/Fam=RLK/ || $9~/Fam=NLR/) && $9!~/Class=Non-canonical/){
                            $9=$9""ADD[id]";color="color[id]
                            }
                         else{
                            $9=$9""ADD[id]";color=3"
                            }}};
                  print}}' alert_NC_Locus.tmp LRRlocus_complet2.tmp > LRRlocus_complet.gff






cat LRRlocus_complet.gff > $RES_DIR/LRRlocus_predicted.gff
cat alert_NC_Locus.tmp | sort -k1,2 > $RES_DIR/alert_NC_Locus.txt

clean_tmp_dir 0 $tmpdir