#!/bin/bash

echo $1
BLASTDB=$1
echo $1
LRRome=$2
echo $2
CDNA=$LRRome/REF_cDNA
GFF=$(cat $3| cut -f2)
echo ------------
echo $GFF
echo ---------
PROTEINS=$LRRome/REF_PEP
CDS=$LRR/REF_CDS
SPECIES=$(cat $3| cut -f1)
echo $SPECIES


function extractSeq {
	##Extracting each sequence from a fasta in separate files
	gawk -F"[;]" '{if($1~/>/){line=$1;gsub(">","");filename=$1;print(line) > filename}else{print > filename}}' $1
}

function filter_Blastp {
	##filtering blastp results;remove redonduncies (Hit insides an other hit);concatenate consecutive blast hit modifying stat
	sort -k1,1 -Vk7,7 $1 | gawk -F"\t" 'BEGIN{OFS="\t"}{
		if(FNR==1){
			Q=$1;S=$2;L=$4;Qstart=$5;Qend=$6;Sstart=$7;Send=$8;nident=$9;line=$0;}
		else{
			if($1==Q && $2==S){
				if($7>=Send-30 && $8>Send && $5>Qtart && $5>Qend-30){
					$4=($4+L-(Send-$7));
					$5=Qstart;
					$7=Sstart;
					$9=$9+nident;
					$10=($9/$4)*100
					line=$0}}
			else{
				print(line);Q=$1;S=$2;Qstart=$5;Qend=$6;Sstart=$7;Send=$8;nident=$9;line=$0;}}
	}END{print(line)}' > $2
}


#export -f extractSeq
#export -f filter_Blastp



          #------------------------------------------#
          # 1. Find Regions of interest with tblastn #
          #------------------------------------------#

mkdir  $LRRome/mmseqs
echo la
mmseqs createdb $1 $LRRome/mmseqs/${SPECIES}_genome_db 
echo 1db
#head $LRRome/mmseqs/$1_genome_db  
#ls $LRRome
mmseqs createdb $LRRome/*_proteins.fasta $LRRome/mmseqs/prot_db 
echo 2db 
#head -n 10 $LRRome/mmseqs/prot_db 
mmseqs search $LRRome/mmseqs/prot_db $LRRome/mmseqs/${SPECIES}_genome_db resultDB_aln.m8 tmp -s 8.5 -a -e 0.1 --min-length 10 --merge-query 1 --cov-mode 2 --max-seqs 30000 --sequence-overlap 1000 
mmseqs convertalis $LRRome/mmseqs/prot_db $LRRome/mmseqs/${SPECIES}_genome_db resultDB_aln.m8 $LRRome/mmseqs/res_candidatsLRR_in_$SPECIES.out --format-output query,target,qlen,alnlen,qstart,qend,tstart,tend,nident,pident,gapopen,evalue,bits  
echo 3db 
wc -l resultDB_aln.m8 
#touch res_candidatsLRR_in_$SPECIES.out
cat $LRRome/mmseqs/res_candidatsLRR_in_$SPECIES.out > res_candidatsLRR_in_$SPECIES.out

echo ---------------res_candidatsLRR_in_$SPECIES.out-----------------ok
wc -l res_candidatsLRR_in_$SPECIES.out


#cat res_candidatsLRR_in_$SPECIES.out
#echo ________________________________________________________________
##-----------------------------
## Constituer des hits globaux par proteines en 2 temps : 1er tour avec des seuils hauts pour fixer des exons d"ancrage"
## puis deuxieme tour avec les hsp plus faibles pour completer les zones aux extremites et definir d'autres zones types paralogues


## 1st run : high threshold
##-----------------------------
## filtering and Sorting
gawk 'BEGIN{OFS="\t"}{if($10>=65){print($0)}}' res_candidatsLRR_in_$SPECIES.out | sort -k1,2 -Vk7,7 > sort_65_candidatsLRR_in_$SPECIES.out
#echo -------------sort_65_candidatsLRR_in_$SPECIES.out-------------
#cat sort_65_candidatsLRR_in_$SPECIES.out
#OKKK
## Si un alignement donne la totalite de la sequence -> extraction region
## sinon, est ce que le hit suivant est proche?... oui -> cumul
## si alignement cumule < 60% de la prot --> elim

## Set intron size with max intron in Nip for the prot

gawk 'BEGIN{OFS="\t"}{
		if(NR==FNR){
			if($3=="gene"){gsub("ID=","",$9);new=1;ID=$9}
			else{
				if($3=="CDS"){
					if(new==1){stop=$5;new=0;MAX_INTRON[ID]=0}
					else{
						if($4-stop>MAX_INTRON[ID]){MAX_INTRON[ID]=($4-stop)};stop=$5}}}}
		else{if($7<$8){strand="+"}else{strand="-"};
		limIntron=5000;if(MAX_INTRON[$1]>limIntron){limIntron=MAX_INTRON[$1]+500};
		if(FNR==1){Q=$1;T=$2;P1=$7;P2=$8;old5=$5;old6=$6;S=strand;line=$0;}
		else{
			split(line,tab,"\t");
			if($1==Q && $2==T && strand==S && ((strand=="+" && $6>old6 && old6<$5+100) || (strand=="-" && $5<old5 && $6<old5+100)) && ($(7)-P2<limIntron || $(8)-P1<limIntron)){
				P1=$7;P2=$8;
				$4=$4+tab[4];
				old5=$5;
				old6=$6;
				if($5>tab[5]){$5=tab[5]};
				if($6<tab[6]){$6=tab[6]};
				if(strand=="+"){
					if($7>tab[7]){$7=tab[7]};
					if($8<tab[8]){$8=tab[8]};}
				else{
					if($7<tab[7]){$7=tab[7]};
					if($8>tab[8]){$8=tab[8]};}
				$9=$9+tab[9];
				$11=$11+tab[11];
				$13=$13+tab[13]
				$10=($9/$4)*100;
				line=$0;}
			else{
				print(line);Q=$1;T=$2;P1=$7;P2=$8;old5=$5;old6=$6;S=strand;line=$0}}
}}END{print(line)}' $GFF sort_65_candidatsLRR_in_$SPECIES.out | sort -k2,2 -Vk7,7 > concat_65_candidatsLRR_in_$SPECIES.tmp
echo icimeme
cat concat_65_candidatsLRR_in_$SPECIES.tmp
## 2nd run : lower threshold
##------------------------------
gawk 'BEGIN{OFS="\t"}{if($10>=45 && $10<65){print($0)}}' res_candidatsLRR_in_$SPECIES.out | sort -k1,2 -Vk7,7 > sort_45_candidatsLRR_in_$SPECIES.out

##we discarded hits falling inside already identified regions
cat concat_65_candidatsLRR_in_$SPECIES.tmp sort_45_candidatsLRR_in_$SPECIES.out | sort -k1,2 -Vk7,7 | gawk 'BEGIN{OFS="\t"}{if(NR==1){query=$1;target=$2;p7=$7;p8=$8;print}else{if($1!=query || $2!=target || ($7<$8 && $8>p8) || ($7>$8 && $7>p7)){print;p7=$7;p8=$8;query=$1;target=$2}}}' > concat_candidatsLRR_in_$SPECIES.tmp


gawk 'BEGIN{OFS="\t"}{
		if(NR==FNR){
			if($3=="gene"){gsub("ID=","",$9);new=1;ID=$9}
			else{
				if($3=="CDS"){
					if(new==1){stop=$5;new=0;MAX_INTRON[ID]=0}
					else{
						if($4-stop>MAX_INTRON[ID]){MAX_INTRON[ID]=($4-stop)};stop=$5}}}}
		else{if($7<$8){strand="+"}else{strand="-"};
		limIntron=5000;if(MAX_INTRON[$1]>limIntron){limIntron=MAX_INTRON[$1]+500};
		if(FNR==1){Q=$1;T=$2;P1=$7;P2=$8;old5=$5;old6=$6;S=strand;line=$0;}
		else{
			split(line,tab,"\t");
			if($1==Q && $2==T && strand==S && ((strand=="+" && $6>old6 && old6<$5+100) || (strand=="-" && $5<old5 && $6<old5+100)) && ($(7)-P2<limIntron || $(8)-P1<limIntron)){
				P1=$7;P2=$8;
				$4=$4+tab[4];
				old5=$5;
				old6=$6;
				if($5>tab[5]){$5=tab[5]};
				if($6<tab[6]){$6=tab[6]};
				if(strand=="+"){
					if($7>tab[7]){$7=tab[7]};
					if($8<tab[8]){$8=tab[8]};}
				else{
					if($7<tab[7]){$7=tab[7]};
					if($8>tab[8]){$8=tab[8]};}
				$9=$9+tab[9];
				$11=$11+tab[11];
				$13=$13+tab[13]
				$10=($9/$4)*100;
				line=$0;}
			else{
				if((tab[6]-tab[5]+1)/tab[3]>=0.6 && tab[10]>=50){print(line)};Q=$1;T=$2;P1=$7;P2=$8;old5=$5;old6=$6;S=strand;line=$0}}
}}END{print(line)}' $GFF concat_candidatsLRR_in_$SPECIES.tmp | sort -k2,2 -Vk7,7 > concat_candidatsLRR_in_$SPECIES.out
#concat_candidatsLRR_in_$SPECIES.out > $directory/concat.txt
#echo --------------------$GFF----------------------
#cat $GFF
echo -------------------------concat_candidatsLRR_in_$SPECIES.tmp ---------------
cat concat_candidatsLRR_in_$SPECIES.tmp 
echo ----------------concat_candidatsLRR_in_$SPECIES.out------------------------
cat concat_candidatsLRR_in_$SPECIES.out
## Regions candidates par query --> a filtrer pour enlever les redondances
## Garder pour une zone d'interet, la query donnant la meilleure couverture 

gawk -F"\t" 'BEGIN{OFS="\t"}{
       if(NR==1){T=$2;P1=$7;P2=$8;score=$(13);line=$0}
       else{
          if($2==T && ($(7)<P2 || P1>$(8))){
              if(score<$(13)){T=$2;P1=$7;P2=$8;score=$(13);line=$0} 
          }else{
              print(line);T=$2;P1=$7;P2=$8;line=$0;score=$(13)}}
}END{print(line)}' concat_candidatsLRR_in_$SPECIES.out > filtered_candidatsLRR_in_$SPECIES.out
echo ---------concat_candidatsLRR_in_$SPECIES.out----------------
cat concat_candidatsLRR_in_$SPECIES.out
echo ------------filtered_candidatsLRR_in_$SPECIES.out--------------------
cat filtered_candidatsLRR_in_$SPECIES.out

gawk -F"\t" 'BEGIN{OFS="\t";}{if(NR==FNR){CHR[FNR]=$2;if($7<$8){START[FNR]=$7;STOP[FNR]=$8}else{START[FNR]=$8;STOP[FNR]=$7}}
                            else{if($2!=CHR[FNR-1]){s1="0";current=$2}else{s1=STOP[FNR-1]};
								 if($2!=CHR[FNR+1]){if($7<$8){s2=$8+5000}else{s2=$7+5000}}else{s2=START[FNR+1]};print($0,s1,s2)}}' filtered_candidatsLRR_in_$SPECIES.out filtered_candidatsLRR_in_$SPECIES.out > filtered_candidatsLRR_in_$SPECIES.out2
echo -----------------filtered_candidatsLRR_in_$SPECIES.out------------------------
cat filtered_candidatsLRR_in_$SPECIES.out
echo ----------------------filtered_candidatsLRR_in_$SPECIES.out-----------------------
cat filtered_candidatsLRR_in_$SPECIES.out


# Extract regions of interest + 300bp before and after
# gff format and query/target list
gawk -v sp=$SPECIES 'BEGIN{OFS="\t";}{
           if($7<$8){
               if($5>10){$7=$7-3000}else{$7=$7-300};
               if($7<=$(14)){$7=$(14)+1};
               if($6<$3-10){$8=$8+3000}else{$8=$8+300};
			   if($8>=$(15)){$8=$(15)-1};
               pos=sprintf("%08d",$7);
               print($2,"TransfertAnnot","gene",$7,$8,".","+",".","ID="sp"_"$2"_"pos";Origin="$1);print(sp"_"$2"_"pos,$1,"+")>"liste_query_target.txt"}
           else{
               if($6<$3-10){$8=$8-3000}else{$8=$8-300};
			   if($8<=$(14)){$8=$(14)+1};
               if($5>10){$7=$7+3000}else{$7=$7+300};
               if($7>=$(15)){$7=$(15)-1};
               pos=sprintf("%08d",$7);
               print($2,"TransfertAnnot","gene",$8,$7,".","-",".","ID="sp"_"$2"_"pos";Origin="$1);print(sp"_"$2"_"pos,$1,"-")>"liste_query_target.txt"}
}' filtered_candidatsLRR_in_$SPECIES.out2 > filtered_candidatsLRR_in_$SPECIES.gff
echo ------------------filtered_candidatsLRR_in_$SPECIES.out2------------------
cat filtered_candidatsLRR_in_$SPECIES.out2
echo -----------------filtered_candidatsLRR_in_$SPECIES.gff----------------
cat filtered_candidatsLRR_in_$SPECIES.gff
cat liste_query_target.txt > candidate_loci_to_LRRome
echo -------------------------liste_query_target.txt--------------------------
wc -l liste_query_target.txt
echo liste
cat liste_query_target.txt 
#cat candidate_loci_to_LRRome 
DIR=$(pwd)
cat filtered_candidatsLRR_in_$SPECIES.gff > filtered_candidatsLRR
#cat filtered_candidatsLRR_in_$SPECIES.gff
echo $(wc -l filtered_candidatsLRR_in_$SPECIES.gff)
python3 /Users/thibaudvicat/pipelinegit/version3/SCRIPT/Extract_sequences_from_genome.py -f $BLASTDB -g filtered_candidatsLRR_in_$SPECIES.gff -o ./DNA_candidatsLRR_in_$SPECIES.fasta  -t gene 

cat ./DNA_candidatsLRR_in_$SPECIES.fasta > xDNA_candidatsLRR_in 
#cat xDNA_candidatsLRR_in 
mkdir CANDIDATE_SEQ_DNA ; cd CANDIDATE_SEQ_DNA
#echo $PWD > CANDIDATE_SEQ_DNA
DIR=$(pwd)
extractSeq ../DNA_candidatsLRR_in_$SPECIES.fasta
chmod 777 *
#cat *
cd ..
#echo $DIR
#echo $DIR > CANDIDATE_SEQ_DNA

#echo $DIR 

#rm -rf $4/Transfert_$1/ || true
#mkdir -p $4../toto