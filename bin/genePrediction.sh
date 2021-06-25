#!/bin/bash
          #------------------------------------------#
          # 2.     Mapping CDS                       #
          #------------------------------------------#
export line=$(echo | cat $1)
echo $line > file
echo "line"
echo $line 
echo "file"
cat file
export TARGET_DNA=$2
export BLASTDB=$3
export SPECIES=$(cat $8 | cut -f1)
export mode=$4
export filtered_candidatsLRR=$5
export resDir=$6/Transfert_$SPECIES
export LRRome=$7
export SCRIPT='/GeneModelTransfer.git/branches/container/SCRIPT/'
export REF_PEP=$LRRome/REF_PEP
export REF_CDS=$LRRome/REF_CDS
export REF_cDNA=$LRRome/REF_cDNA
export GFF=$(cat $8 | cut  -f2)
echo "$line" > to_transfer_with_cdna.txt
echo "$line" > to_transfer_with_prot.txt
export infoLocus=$(cat $8 | cut  -f4)
mkdir mapping ; cd mapping
function extractSeq {
	##Extracting each sequence from a fasta in separate files
	gawk -F"[;]" '{if($1~/>/){line=$1;gsub(">","");filename=$1;print(line) > filename}else{print > filename}}' $1
}
echo $line >> $resDir/line
cat file >> $resDir/file
export -f extractSeq
function mapcds {
   # Param 1 : TARGET = fichier sequence genomique d'interet chez la cible
   # Param 2 : QUERY = ID proteine de Nip pour mapping dans la zone
   cat $REF_CDS/$1* > query.fasta
   blastn -query query.fasta -subject $TARGET_DNA/$2 -outfmt "6 qseqid sseqid qlen length qstart qend sstart send nident pident gapopen" > blastn.tmp
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
							cds=cds+1;S=$1;P1=$6;P2=$8;}}}}' ../file - | sed 's/gene/Agene/g' | sort -Vk4,4 | sed 's/Agene/gene/g' > $2.gff
		#cat $2.gff >> $resDir/blastCDS.gff
		#cat $2.gff > test.gff
		## verif par blast
		python3 $SCRIPT/Extract_sequences_from_genome.py -f $BLASTDB -g $2.gff -o $2.fasta -t prot 2>/dev/null

		blastp -query $2.fasta -subject $REF_PEP/$1 -outfmt "6 qseqid sseqid slen length qstart qend sstart send nident pident gapopen" > blastp.tmp
		cat blastp.tmp >> blastp.save
		#si cov > 97% et pid>75% = ok
		if [[ -s blastp.tmp ]];then
			sh $SCRIPT/filter_Blastp.sh blastp.tmp blastp2.tmp
			check=$(gawk 'NR==1{if(($8-$7+1)/$3>=0 && $10>0){print(1)}else{print(0)}}' blastp2.tmp)
		fi
		if [[ $check -eq 1 ]];then
			# Res blast + ajout GFF global
			gawk -F"\t" 'BEGIN{OFS="\t"}{if(NR==FNR){IDENT=$10;COV=($8-$7+1)/$3}
			else{if($3~/gene/){print($0" / pred:mappingCDS / blast-%-ident:"IDENT" / blast-cov:"COV)}
				else{print}}}' blastp.tmp $2.gff > $2_2.gff
			cat $2_2.gff >> mappingCDS_$SPECIES.gff
		fi
		rm $2.gff
		rm $2.fasta
	fi
}
export query=$(echo $line | cut -d ' ' -f1)
export target=$(echo $line | cut -d ' ' -f2)
#mapcds $target $query
