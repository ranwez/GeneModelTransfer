#!/bin/bash
#========================================================
# PROJET : TransfertGeneModel
# SCRIPT : Transfert.sh
# AUTHOR : Celine Gottin
# CREATION : 2020.01.27
#========================================================
# DESCRIPTION : Transfert des modeles de gene des locus LRR
#               de Nipponbare vers les autres genomes de riz
#               via exonerate.
# ARGUMENTS : o $1 : species
#             o $2 : blast DB of target species
#             o $3 : cDNA sequence from Nip
#             o $4 : gff file of LRR protein from Nip
#             o $5 : protein sequences from Nip
# DEPENDENCIES : o ncbi-blast
#                o exonerate
#========================================================


#========================================================
#                Environment & variables
#========================================================

module load bioinfo/ncbi-blast/2.6.0
module load bioinfo/exonerate/2.4.7

SPECIES=$1
BLASTDB=$2
CDNA=$3
GFF=$4
PROTEINS=$5
CDS=$6

SCRIPT=$(pwd)/SCRIPT
WD=Transfert_$SPECIES
mkdir $WD; 
cd $WD

#========================================================
#                Script
#========================================================

##
## Functions
##

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

##
## Extracting Sequences
##

mkdir REF_CDS ; cd REF_CDS

extractSeq $CDS

cd ..

mkdir REF_PEP ; cd REF_PEP

extractSeq $PROTEINS

cd ..

mkdir REF_cDNA ; cd REF_cDNA

extractSeq $CDNA

cd ..


          #------------------------------------------#
          # 1. Find Regions of interest with tblastn #
          #------------------------------------------#

#tblastn -db $BLASTDB -query $PROTEINS -evalue 1 -out res_candidatsLRR_in_$SPECIES.out -outfmt "6 qseqid sseqid qlen length qstart qend sstart send nident pident gapopen evalue bitscore" #-max_target_seqs 2


##-----------------------------
## Constituer des hits globaux par proteines en 2 temps : 1er tour avec des seuils hauts pour fixer des exons d"ancrage"
## puis deuxieme tour avec les hsp plus faibles pour completer les zones aux extremites et definir d'autres zones types paralogues


## 1st run : high threshold
##-----------------------------
## filtering and Sorting
gawk 'BEGIN{OFS="\t"}{if($10>=65){print($0)}}' res_candidatsLRR_in_$SPECIES.out | sort -k1,2 -Vk7,7 > sort_65_candidatsLRR_in_$SPECIES.out

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

gawk -F"\t" 'BEGIN{OFS="\t";}{if(NR==FNR){CHR[FNR]=$2;if($7<$8){START[FNR]=$7;STOP[FNR]=$8}else{START[FNR]=$8;STOP[FNR]=$7}}
                            else{if($2!=CHR[FNR-1]){s1="0";current=$2}else{s1=STOP[FNR-1]};
								 if($2!=CHR[FNR+1]){if($7<$8){s2=$8+5000}else{s2=$7+5000}}else{s2=START[FNR+1]};print($0,s1,s2)}}' filtered_candidatsLRR_in_$SPECIES.out filtered_candidatsLRR_in_$SPECIES.out > filtered_candidatsLRR_in_$SPECIES.out2

echo "END PART1"


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

python $SCRIPT/Extract_sequences_from_genome.py -f $BLASTDB -g filtered_candidatsLRR_in_$SPECIES.gff -o DNA_candidatsLRR_in_$SPECIES.fasta -t gene 



mkdir TARGET_DNA ; cd TARGET_DNA

extractSeq ../DNA_candidatsLRR_in_$SPECIES.fasta

cd ..

echo "END PART2"

          #------------------------------------------#
          # 2.     Mapping CDS                       #
          #------------------------------------------#
mkdir mapping ; cd mapping

function mapcds {
    # Param 1 : TARGET = fichier sequence genomique d'interet chez la cible
    # Param 2 : QUERY = ID proteine de Nip pour mapping dans la zone
    #echo $1
    #echo $2
    cat ../REF_CDS/${2}* > query.fasta

    blastn -query query.fasta -subject ../TARGET_DNA/$1 -outfmt "6 qseqid sseqid qlen length qstart qend sstart send nident pident gapopen" > blastn.tmp
    
	if [[ -s blastn.tmp ]];then
		cat blastn.tmp >> blastn.save
		## traitement resultats
		# 1. retrait des matchs incoherents par rapport a la query(sort par id cds Nip, taille d'aligenement)
		sort -k1,1 -Vrk4,4 blastn.tmp | gawk 'BEGIN{OFS="\t"}{
				if(NR==1){P5=$5;P6=$6;currentCDS=$1;print}
				else{if(($1==currentCDS && $5>P6-10) || $1!=currentCDS){print;P5=$5;P6=$6;currentCDS=$1}}}' > blastn2.tmp
	
		# 2. sortir GFF ///// /!\ \\\\\ ATTENTION, la determination du chromosome d�eend de la nomenclature de la cible
		sort -Vk7,7 blastn2.tmp | gawk 'BEGIN{OFS="\t";cds=1}{
				if(NR==FNR){
					strand[$1]=$3}
				else{
					split($1,M,":");
					split($2,T,"_");
					if(length(T)==3){chr=T[2]}else{chr=T[2]"_"T[3]};
					pos=T[length(T)];
					if(strand[$2]=="+"){deb=(pos+$7-1);fin=(pos+$8-1)}
					else{deb=(pos-$8+1);fin=(pos-$7+1)};
					if(FNR==1){
						print(chr,"blastCDS","gene","0",fin,".",strand[$2],".","ID="$2";origin="M[1]);
						S=$1;P1=$6;P2=$8;
						print(chr,"blastCDS","CDS",deb,fin,".",strand[$2],".","ID="$2":cds"cds);cds=cds+1}
					else{
						if(($1==S && $7>P2-10 && $5>P1-10) || ($1!=S)){
							print(chr,"blastCDS","CDS",deb,fin,".",strand[$2],".","ID="$2":cds"cds);
							cds=cds+1;S=$1;P1=$6;P2=$8;}}}}' ../liste_query_target.txt - | sed 's/gene/Agene/g' | sort -Vk4,4 | sed 's/Agene/gene/g' > $1.gff


		## verif par blast
		python $SCRIPT/Extract_sequences_from_genome.py -f $BLASTDB -g $1.gff -o $1.fasta -t prot 2>/dev/null

		blastp -query $1.fasta -subject ../REF_PEP/$2 -outfmt "6 qseqid sseqid slen length qstart qend sstart send nident pident gapopen" > blastp.tmp
 
		cat blastp.tmp >> blastp.save
		#si cov > 97% et pid>75% = ok
		if [[ -s blastp.tmp ]];then
			filter_Blastp blastp.tmp blastp2.tmp
			check=$(gawk 'NR==1{if(($8-$7+1)/$3>=0.97 && $10>75){print(1)}else{print(0)}}' blastp2.tmp)
		fi

		if [[ $check -eq 1 ]];then
			# Res blast + ajout GFF global
			gawk -F"\t" 'BEGIN{OFS="\t"}{if(NR==FNR){IDENT=$10;COV=($8-$7+1)/$3}
			else{if($3~/gene/){print($0" / pred:mappingCDS / blast-%-ident:"IDENT" / blast-cov:"COV)}
				else{print}}}' blastp.tmp $1.gff > ${1}_2.gff

			cat ${1}_2.gff >> mappingCDS_$SPECIES.gff
			rm ${1}_2.gff
		fi

		rm $1.gff
		rm $1.fasta
	fi
}

while read line
do
    #echo "$line"
    query=$(echo $line | cut -f1)
    target=$(echo $line | cut -f2)

    mapcds $target $query

done < ../liste_query_target.txt


##corriger pos gene, lancer correction
gawk -F"\t" 'BEGIN{OFS="\t"}{split($9,T,/[=:;]/);if(NR==FNR){if($3=="gene"){max[T[2]]=$5;min[T[2]]=$5}else{if($5>max[T[2]]){max[T[2]]=$5};if($4<min[T[2]]){min[T[2]]=$4}}}else{if($3=="gene"){$4=min[T[2]];$5=max[T[2]]};print}}' mappingCDS_$SPECIES.gff mappingCDS_$SPECIES.gff > tmp ; mv tmp mappingCDS_$SPECIES.gff

cd ..

python $SCRIPT/Exonerate_correction.py -f $BLASTDB -g mapping/mappingCDS_$SPECIES.gff > mapping_LRRlocus_${SPECIES}.gff

##to transfert with cdna2genome
gawk 'BEGIN{OFS="\t"}{if(NR==FNR){if($3=="gene"){split($9,T,";");gsub("ID=","",T[1]);OK[T[1]]=1}}else{if(!OK[$1]==1){print}}}' mapping_LRRlocus_$SPECIES.gff liste_query_target.txt > to_transfer_with_cdna.txt

echo "END PART3"

          #------------------------------------------#
          # 3.     Run exonerate cdna2genome         #
          #------------------------------------------#



mkdir exonerate ; cd exonerate

# exonerate lines
gawk -F"\t" -v species=$SPECIES '{target=$1;query=$2;print("exonerate -m cdna2genome --bestn 1 --showalignment no --showvulgar no --showtargetgff yes --annotation",query".an --query ../REF_cDNA/"query,"--target ../TARGET_DNA/"target,">> LRRlocus_in_"species"_cdna.out")}' ../to_transfer_with_cdna.txt > exe ##cdna

## annotation files
gawk -F"\t" 'BEGIN{OFS="\t"}{if($3=="gene"){start=1;split($9,T,";");id=substr(T[1],4);filename=id".an"}else{if($3=="CDS"){len=$5-$4+1;print(id,"+",start,len)>>filename;start=start+len}}}' $GFF 

chmod +x exe 
./exe

rm *.an

cd ..

echo "END PART4"

function parseExonerate {
    # with $1 type of exonerate model --> cdna or prot
    gawk -F"\t" 'BEGIN{OFS="\t"}{if($7=="+" && ($3=="gene" || $3=="similarity")){print}}' exonerate/LRRlocus_in_${SPECIES}_$1.out > exonerate/LRRlocus_in_${SPECIES}_$1.tmp ##prot,cdna

    ## Reconstruct gff from align section
    gawk -F"\t" 'BEGIN{OFS="\t"}{if($3=="similarity"){split($9,T,";");for(i=3;i<=length(T);i++){split(T[i],M," ");start=M[2];end=start+M[4]-1;print(chr,proj,"CDS",start,end,".",strand,".",$9)}}else{print;proj=$2;chr=$1;strand=$7}}' exonerate/LRRlocus_in_${SPECIES}_$1.tmp > exonerate/LRRlocus_in_${SPECIES}_$1.gff ##cdna,prot

   ## definir ID et Parent et strand
    gawk -F"\t" 'BEGIN{OFS="\t"}{
               if(NR==FNR){
                  split($9,M,/[=;]/);strand[M[2]]=$7} 
               else{
                  split($1,T,"_");$7=strand[$1];
                  if(length(T)==3){pos=T[3];name=T[2]}else{pos=T[4];name=T[2]"_"T[3]}
                  if(strand[$1]=="+"){
                     $4=pos+$4-1;$5=pos+$5-1}
                  else{
                     o4=$4;o5=$5;$5=pos-o4+1;$4=pos-o5+1};
                  if($3=="gene"){$9="ID="$1};
                  if($3=="CDS"){$3="CDS";$9="Parent="$1};
                  $1=name;print}
}' filtered_candidatsLRR_in_$SPECIES.gff exonerate/LRRlocus_in_${SPECIES}_$1.gff > filtered_LRRlocus_in_${SPECIES}_$1.gff

    ## Eliminer redondance des genes

    gawk -F"\t" 'BEGIN{OFS="\t"}{if(NR==FNR){if($3=="gene"){if(!START[$9] || $4<START[$9]){START[$9]=$4};if($5>STOP[$9]){STOP[$9]=$5}}}else{if($3=="CDS"){T[$1,$4,$5]++;if(T[$1,$4,$5]==1){print}}else{M[$9]++;if(M[$9]==1){$4=START[$9];$5=STOP[$9];print}}}}' filtered_LRRlocus_in_${SPECIES}_$1.gff filtered_LRRlocus_in_${SPECIES}_$1.gff > filtered2_LRRlocus_in_${SPECIES}_$1.gff

    sed 's/gene/Agene/g' filtered2_LRRlocus_in_${SPECIES}_$1.gff | sort -k1,1 -Vk4,4 -k3,3 | sed 's/Agene/gene/g' > filtered3_LRRlocus_in_${SPECIES}_$1.gff

    ## eliminer redondance et chevauchement des CDS si sur la même phase ( on check la phase par $4 si strand + et par $5 si strand -)
    # 1. retrait cds inclus dans autres cds et cds colle
    gawk -F"\t" 'BEGIN{OFS="\t"}{if($3~/gene/){print;lim=0}else{if($5>lim){print;lim=$5}}}' filtered3_LRRlocus_in_${SPECIES}_$1.gff | gawk -F"\t" 'BEGIN{OFS="\t";line=""}{if($3=="gene"){if(line!=""){print(line)};print;lim=0;line=""}else{if($4==(lim+1)){$4=old4;line=$0;lim=$5}else{if(line!=""){print(line)};old4=$4;line=$0;lim=$5}}}' | gawk -F"\t" 'BEGIN{OFS="\t"}{if($5-$4>3){print}}' > filtered4_LRRlocus_in_${SPECIES}_$1.gff


   # 2. retrait chevauchement et intron de moins de 15 bases si même phase 
   ## checker le changement de phase 
    gawk -F"\t" 'BEGIN{OFS="\t";line=""}{
                if($3~/gene/){if(line!=""){print(line)};print;p=0}
                else{
                   if(p==0){
                       p=1;line=$0;start=$4;stop=$5;mod=($(5)+1)%3}
                   else{
                       if($5<stop+15 || $4<start+15){$4=start;stop=$5;line=$0;mod=($(5)+1)%3}
                       else{
                          if(($4>stop+15) || $4%3!=mod){print(line);line=$0;start=$4;stop=$5;mod=($(5)+1)%3}
                          else{$4=start;stop=$5;line=$0;mod=($(5)+1)%3}
 }}}}END{print(line)}' filtered4_LRRlocus_in_${SPECIES}_$1.gff > filtered5_LRRlocus_in_${SPECIES}_$1.tmp

}

parseExonerate cdna

echo "END PART5"

#Correct PROT
python $SCRIPT/Exonerate_correction.py -f $BLASTDB -g filtered5_LRRlocus_in_${SPECIES}_cdna.tmp > filtered5_LRRlocus_in_${SPECIES}_cdna.gff


## verif prot
# extraction, alignement des prot
python $SCRIPT/Extract_sequences_from_genome.py -f $BLASTDB -g filtered5_LRRlocus_in_${SPECIES}_cdna.gff -o PROT_predicted_from_cdna_in_$SPECIES.fasta -t prot 

# BLAST
mkdir Blast
cd Blast

extractSeq ../PROT_predicted_from_cdna_in_$SPECIES.fasta

# generer les lignes executables
gawk -F"\t" -v species=$SPECIES 'BEGIN{OFS=""}{query=$1;subject=$2;print("blastp -query ",query," -subject ../REF_PEP/"subject," -outfmt \"6 qseqid sseqid slen length qstart qend sstart send nident pident gapopen\" >> res_predicted_from_cdna_in_",species,".out")}' ../to_transfer_with_cdna.txt > exe

chmod +x exe ; ./exe

filter_Blastp res_predicted_from_cdna_in_$SPECIES.out res_predicted_from_cdna_in_$SPECIES.out2

cd ..

# extraction seq cible a repredir en protein2genome
gawk -F"\t" '{if(NR==FNR){if($10>=70 && ($8-$7+1)/$3>=0.97){OK[$1]=1}}else{if(!OK[$1]){print}}}' Blast/res_predicted_from_cdna_in_$SPECIES.out2 to_transfer_with_cdna.txt > to_transfer_with_prot.txt
 
# gff avec info origin + blast dans section comment
gawk -F"\t" 'BEGIN{OFS="\t"}{if(NR==FNR){Nip[$1]=$2;ID[$1]=$10;COV[$1]=($8-$7+1)/$3}else{if($3~/gene/){split($9,T,";");locname=substr(T[1],4);gsub("comment=","",T[2]);$9=T[1];print($0";comment=Origin:"Nip[locname]" / pred:cdna2genome / blast-%-ident:"ID[locname]" / blast-cov:"COV[locname]" / "T[2])}else{print}}}' Blast/res_predicted_from_cdna_in_$SPECIES.out2 filtered5_LRRlocus_in_${SPECIES}_cdna.gff > filtered6_LRRlocus_in_${SPECIES}_cdna.gff

echo "END PART6"

          #------------------------------------------#
          # 4.     Run exonerate prot2genome         #
          #------------------------------------------#


# 4. Transfert avec exonerate protein2genome
cd exonerate ; rm exe

gawk -v species=$SPECIES '{target=$1;query=$2;print("exonerate -m protein2genome --showalignment no --showvulgar no --showtargetgff yes -q ../REF_PEP/"query,"-t ../TARGET_DNA/"target,">> LRRlocus_in_"species"_prot.out")}' ../to_transfer_with_prot.txt > exe 

chmod +x exe; ./exe

cd ..

echo "END PART 7"

parseExonerate prot

echo "END PART8"

#Correct PROT
python $SCRIPT/Exonerate_correction.py -f $BLASTDB -g filtered5_LRRlocus_in_${SPECIES}_prot.tmp > filtered6_LRRlocus_in_${SPECIES}_prot.gff

echo "END PART9"
####  BLAST + ajouter res blast au gff dans section comment + method=prot2genome

# extraction, alignement des prot
python $SCRIPT/Extract_sequences_from_genome.py -f $BLASTDB -g filtered6_LRRlocus_in_${SPECIES}_prot.gff -o PROT_predicted_from_prot_in_$SPECIES.fasta -t prot 

### blast
cd Blast
rm ${SPECIES}_*

extractSeq ../PROT_predicted_from_prot_in_$SPECIES.fasta

# generer les lignes executables
rm exe
gawk -v species=$SPECIES 'BEGIN{OFS=""}{query=$1;subject=$2;print("blastp -query ",query" -subject ../REF_PEP/"subject" -outfmt \"6 qseqid sseqid slen length qstart qend sstart send nident pident gapopen\" >> res_predicted_from_prot_in_",species,".out")}' ../to_transfer_with_prot.txt > exe

# lancer les blasts
chmod +x exe ; ./exe

filter_Blastp res_predicted_from_prot_in_$SPECIES.out res_predicted_from_prot_in_$SPECIES.out2

cd ..
 
# gff avec info origin + blast dans section comment
gawk -F"\t" 'BEGIN{OFS="\t"}{
       if(NR==FNR){
          Nip[$1]=$2;ID[$1]=$10;COV[$1]=($8-$7+1)/$3}
       else{if($3~/gene/){  
              gsub("comment=","");
              split($9,T,";");
              locname=substr(T[1],4);
              $9=T[1];print($0";comment=Origin:"Nip[locname]" / pred:prot2genome / blast-%-ident:"ID[locname]" / blast-cov:"COV[locname]" / "T[2])}else{print}}}' Blast/res_predicted_from_prot_in_$SPECIES.out2 filtered6_LRRlocus_in_${SPECIES}_prot.gff > filtered7_LRRlocus_in_${SPECIES}_prot.gff



          #------------------------------------------#
          # 5.     Concatenate all results           #
          #------------------------------------------#



echo "END PART10"
#### CONCATENATION GFF cdna2genome ET protein2genome AVEC INFO LOCUS NIP

# liste locus prot2genome a enlever des resultats cdna2genome
cut -f1 to_transfer_with_prot.txt > list_prot2genome.txt

# retrait des sequences predites en prot2genome des resultats cdna2genome
grep -v -f list_prot2genome.txt filtered6_LRRlocus_in_${SPECIES}_cdna.gff > filtered7_LRRlocus_in_${SPECIES}_cdna.gff

# Concatenation avec prot2genome
cat filtered7_LRRlocus_in_${SPECIES}_cdna.gff filtered7_LRRlocus_in_${SPECIES}_prot.gff | sed 's/CDS/zCDS/g' | sort -k1,1 -Vk4,4 -k3,3 | sed 's/zCDS/CDS/g' > LRRlocus_in_${SPECIES}_complet.tmp

## Concatener avec res mapping
#cat ${SPECIES}_liftover_clean.gff LRRlocus_in_${SPECIES}_complet.tmp | sed 's/CDS/zCDS/g' | sort -k1,1 -Vk4,4 -k3,3 | sed 's/zCDS/CDS/g' > LRRlocus_in_${SPECIES}_complet2.tmp
cat mapping_LRRlocus_$SPECIES.gff LRRlocus_in_${SPECIES}_complet.tmp | sed 's/CDS/zCDS/g' | sort -k1,1 -Vk4,4 -k3,3 | sed 's/zCDS/CDS/g' > LRRlocus_in_${SPECIES}_complet2.tmp


# Ajout comment : famille gene Nip, classe gene Nip, +autre
gawk -F"\t" 'BEGIN{OFS="\t"}{if(NR==FNR){F[$1]=$2;C[$1]=$3}else{if($3~/gene/){split($9,T,/[;/]/);origin=substr(T[2],16);gsub(" ","",origin);$9=$9" / Gene-Fam="F[origin]" / Gene-Class="C[origin]};print}}' ../Info_locus_Nipponbare.txt LRRlocus_in_${SPECIES}_complet2.tmp > LRRlocus_in_${SPECIES}_complet.gff

## verif annot
gawk 'BEGIN{OFS=";"}{if($3~/gene/){if(line){print(line)};split($9,T,";");line=substr(T[1],4)";"$7}else{if($3=="CDS"){line=line";"$4";"$5}}}END{print(line)}' LRRlocus_in_${SPECIES}_complet.gff > geneModel_${SPECIES}.tbl

python $SCRIPT/Canonical_gene_model_test.py -f $BLASTDB -t geneModel_${SPECIES}.tbl > alert.txt


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

#mv tmp LRRlocus_in_${SPECIES}_complet.gff


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


# rm *.tmp


## Res liftover
#gawk -F"\t" 'BEGIN{OFS="\t";p=0}{if($3=="CDS" && p==1){print}else{if($3=="gene"){line=$0}else{if($9~/color=3/ || $9~/color=9/){print(line);print;p=1}else{p=0}}}}' ../${SPECIES}_liftover/* > tmp

#gawk -F"\t" 'BEGIN{OFS="\t"}{if(NR==FNR){if($3=="gene"){name=$9;min[$9]=$4;max[$9]=$5}else{if($4<min[name]){min[name]=$4};if($5>max[name]){max[name]=$5}}}else{if($3=="gene"){name=$9;$4=min[name];$5=max[name]}else{if($3=="mRNA"){$4=min[name];$5=max[name]}};print}}' tmp tmp > ${SPECIES}_liftover_clean.tmp

#sed 's/mRNA/AmRNA/g' ${SPECIES}_liftover_clean.tmp | sort -k3,3 | gawk -F"\t" 'BEGIN{OFS="\t"}{if($3=="AmRNA"){split($9,T,";");feat[substr(T[2],8)]=T[4]";"T[5];print($1,$2,$3,$4,$5,$6,$7,$8,T[1]";"T[2]";"T[3])}else{if($3=="gene"){split($9,T,";");prot=substr(T[1],4);print($0";"feat[prot])}else{print}}}' | sed 's/AmRNA/mRNA/g' | sed 's/CDS/zCDS/g' | sort -k1,1 -Vk4,4 -k3,3 | sed 's/zCDS/CDS/g' > ${SPECIES}_liftover.gff

##gawk -F"\t" 'BEGIN{OFS="\t"}{if($3=="gene"){split($9,T,";");gsub("ID=","",T[1]);id=T[1];split(T[2],M," / ");gsub("comment=origin:","",M[1]);gsub("comment=Origin:","",M[1]);ori=M[1];print(id,ori)}}' Transfert_new_kitaake/LRRlocus_in_kitaake_complet.gff > liste_gene.txt
##gawk 'BEGIN{OFS="\t"}{if(NR==FNR){F[$1]=$2;C[$1]=$3}else{print($0,F[$2],C[$2])}}' Info_locus_Nipponbare.txt liste_gene.txt > tmp
##gawk 'BEGIN{OFS="\t"}{if(NR==FNR){F[$1]=$3;C[$1]=$4}else{print($0,F[$2],C[$2])}}' liste_gene.txt alert.txt > alert2.txt