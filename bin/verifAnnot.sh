#!/bin/bash
#========================================================
# PROJET : lrrtransfer
# SCRIPT : verifAnnot.sh
# AUTHOR : Celine Gottin & Thibaud Vicat & Vincent Ranwez
# CREATION : 2020.02.20
#========================================================
# DESCRIPTION : Check the new gene models
#               annotate for control (presence of start, stop
#               canonical intron, non-overlapping frmaeshift)
# ARGUMENTS : o $1 : Path to a text file with ref locus info
#             o $2 : Target genome
#             o $3 : predicted annotation GFF
#             o $4 : Result directory
#             o $5 : Path toward LRR script  directory
#             o $6 : name of the transfert method used (to differentiate final results) 
#========================================================

#========================================================
#                Environment & variables
#========================================================


infoLocus=$1
TARGET_GENOME=$2
GFF=$3

RES_DIR=$4
LRR_SCRIPT=$5
method=$6
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
            split($9,T,/[;/]/);origin=substr(T[2],16);gsub(" ","",origin);$9=$9" / Origin-Fam:"F[origin]" / Origin-Class:"C[origin]};print}}' $infoLocus $GFF > LRRlocus_complet.tmp


## concatenate CDS if less than 25 nucl appart and if in the same reading frame  
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
# 2=Red ; 10=Orange ; 3=Green
## Red if RLP/RLK/NLR and Non-canonical
gawk -F"\t" 'BEGIN{OFS="\t"}{
              if(NR==FNR){
                COMMENT[$2]="";
                if($3=="True"){COMMENT[$2]=COMMENT[$2]" / noStart"};
                if($4=="True"){COMMENT[$2]=COMMENT[$2]" / noStop"};
                if($5=="True"){COMMENT[$2]=COMMENT[$2]" / pbFrameshift"};
                if($6=="True"){COMMENT[$2]=COMMENT[$2]" / ncIntron"};
                if($7=="True"){COMMENT[$2]=COMMENT[$2]" / stopInFrame"};
                if($8=="True"){COMMENT[$2]=COMMENT[$2]" / pbLength"};
                
                if($3$4$5$6$7$8~/True/){
                    NC[$2]=1;
                    COMMENT[$2]="Gene-Class:Non-canonical"COMMENT[$2]
                }else{
                    COMMENT[$2]="Gene-Class:Canonical"COMMENT[$2]
                };
              }else{
                if($3=="gene"){
                  split($9,infos,";");
                  id=substr(infos[1],4);
                  genecolor=3;
                  if(NC[id]==1){
                    genecolor=2;
                    if($9~/ident:100/ && $9~/cov:1/){genecolor=10};
                  };
                  $9=$9";color="genecolor"; comment="COMMENT[id];
                };
                print}}' alert_NC_Locus.tmp LRRlocus_complet2.tmp > LRRlocus_complet.gff


## sort gff file
grep "gene" LRRlocus_complet.gff | cut -d"=" -f2 | cut -d";" -f1 | sort -g > list_gene.tmp


for gn in $(cat list_gene.tmp)
do
	grep $gn LRRlocus_complet.gff >> LRRlocus_predicted.gff
done


## export output files
cat LRRlocus_predicted.gff > $RES_DIR/LRRlocus_${method}_predicted.gff
cat alert_NC_Locus.tmp | sort -k1,2 > $RES_DIR/alert_${method}_NC_Locus.txt

clean_tmp_dir 0 $tmpdir
