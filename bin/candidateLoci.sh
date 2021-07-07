#!/bin/bash

BLASTDB=$1
echo blastdb
echo $1
echo LRRome
LRRome=$2
echo $2
CDNA=$LRRome/REF_cDNA
GFF=$(cat $3| cut -f2)
echo $GFF
PROTEINS=$LRRome/REF_PEP
echo $PROTEINS
CDS=$LRR/REF_CDS
echo $CDS
SPECIES=$(cat $3| cut -f1)
treshold1=$(cat $3| cut -f5)
treshold2=$(cat $3| cut -f6)
echo $SPECIES
SCRIPT='/GeneModelTransfer.git/branches/container/SCRIPT'
echo $SCRIPT
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

          #------------------------------------------#
          # 1. Find Regions of interest with mmseqs2 #
          #------------------------------------------#

mkdir $LRRome/mmseqs
echo "mmseqs createdb $1 $LRRome/mmseqs/${SPECIES}_genome_db -v 0"
mmseqs createdb $1 $LRRome/mmseqs/${SPECIES}_genome_db -v 0
mmseqs createdb $LRRome/*_proteins.fasta $LRRome/mmseqs/prot_db  -v 0
mmseqs search $LRRome/mmseqs/prot_db $LRRome/mmseqs/${SPECIES}_genome_db resultDB_aln.m8 tmp -s 8.5 -a -e 0.1 --min-length 10 --merge-query 1 --cov-mode 2 --max-seqs 30000 --sequence-overlap 1000 -v 0