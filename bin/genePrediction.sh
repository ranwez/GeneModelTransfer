#!/bin/bash
#========================================================
# PROJET : LRRtransfer
# SCRIPT : genePrediction.sh
# AUTHOR : Celine Gottin & Thibaud Vicat
# CREATION : 2020.02.20
#========================================================
# DESCRIPTION : Use of blast and exonerate to predict gene models from 
#               query/target pairs and build an annotation based on selected mode
# ARGUMENTS : o $1 : Query/target couple
#             o $2 : Directory with extracted genomic regions of interest
#             o $3 : Target genome
#             o $4 : Filtered_candidatsLRR
#             o $5 : Path to LRRome
#             o $6 : Path to the query GFF
#             o $7 : Path to te infofile
#             o $8 : Results directory
#             o $9 : Path to outfile
#             o $10 : Selected mode

#========================================================
#                Environment & variables
#========================================================


pairID=$1
TARGET_DNA=$2
TARGET_GENOME=$3
filtered_candidatsLRR=$4
LRRome=$5
GFF=$6
infoLocus=$7
RES_DIR=$8

outfile=$9
mode=${10}

REF_PEP=$LRRome/REF_PEP
REF_EXONS=$LRRome/REF_EXONS
REF_cDNA=$LRRome/REF_cDNA



target=$(cat $pairID | cut -f1)
query=$(cat $pairID | cut -f2)


#========================================================
#                        Functions
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

function filter_Blastp {
	##usage :: filter_Blastp blastp.res blastp_filter.res
	##filtering blastp results;remove redonduncies (Hit insides an other hit);
	##concatenate consecutive blast hit modifying stat
	sort -k1,1 -Vk5,5 $1 | gawk -F"\t" 'BEGIN{OFS="\t"}{
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




function parseExonerate {
	# with $1 type of exonerate model --> cdna or prot
	gawk -F"\t" 'BEGIN{OFS="\t"}{if($7=="+" && ($3=="gene" || $3=="similarity")){print}}' LRRlocus_$1.out > LRRlocus_$1.tmp ##prot,cdna

	## Reconstruct gff from align section
	gawk -F"\t" 'BEGIN{OFS="\t"}{
      if($3=="similarity"){
          split($9,T,";");
          for(i=3;i<=length(T);i++){
              split(T[i],M," ");
              start=M[2];
              end=start+M[4]-1;
              print(chr,proj,"CDS",start,end,".",strand,".",$9)}}
      else{
          print;
          proj=$2;
          chr=$1;
          strand=$7}
	}' LRRlocus_$1.tmp > LRRlocus_$1.gff ##cdna,prot

	## define ID and Parent and strand
	gawk -F"\t" 'BEGIN{OFS="\t"}{
				if(NR==FNR){
					split($9,M,/[=;]/);strand[M[2]]=$7} 
				else{
					split($1,T,"_");$7=strand[$1];
					if(length(T)==3){pos=T[3];name=T[2]}else{pos=T[2];name=T[1]}
					if(strand[$1]=="+"){
						$4=pos+$4-1;$5=pos+$5-1}
					else{
						o4=$4;o5=$5;$5=pos-o4+1;$4=pos-o5+1};
					if($3=="gene"){$9="ID="$1};
					if($3=="CDS"){$9="Parent="$1};
					$1=name;
					print}
	}' $filtered_candidatsLRR LRRlocus_$1.gff > filtered_LRRlocus_$1.gff

	## Eliminate gene redundancy
	gawk -F"\t" 'BEGIN{OFS="\t"}{
          if(NR==FNR){
              if($3=="gene"){
                  if(!START[$9] || $4<START[$9]){
                      START[$9]=$4};
                  if($5>STOP[$9]){
                      STOP[$9]=$5}}}
          else{
              if($3=="CDS"){
                  T[$1,$4,$5]++;
                      if(T[$1,$4,$5]==1){
                          print}}
              else{
                  M[$9]++;
                  if(M[$9]==1){
                      $4=START[$9];
                      $5=STOP[$9];
                      print}}}
	}' filtered_LRRlocus_$1.gff filtered_LRRlocus_$1.gff > filtered2_LRRlocus_$1.gff

	sed 's/gene/Agene/g' filtered2_LRRlocus_$1.gff | sort -k1,1 -Vk4,4 -k3,3 | sed 's/Agene/gene/g' > filtered3_LRRlocus_$1.gff

	## eliminate redundancy and overlap of CDS if on the same phase (we check the phase by $4 if strand + and by $5 if strand -)
	# 1. remove cds included in other cds and glued cds

	#cat filtered3_LRRlocus_$1.gff > filtered4_LRRlocus_$1.gff
	gawk -F"\t" 'BEGIN{OFS="\t"}{
		if($3~/gene/){
			print;
			lim=0}
		else{
			if($5>lim){
				print;
				lim=$5}}}' filtered3_LRRlocus_$1.gff | gawk 'BEGIN{OFS="\t"; line=""}{
		if($3=="gene" ){
			if(line!=""){
				print(line)};
			print;
			lim=0;
			line=""}
		else{
			if($4==(lim+1)){
				$4=old4;
				line=$0;
				lim=$5}
			else{
				if(line!=""){
					print(line)};
				old4=$4;
				line=$0;
				lim=$5}}}END{print(line)}' | gawk -F"\t" 'BEGIN{OFS="\t"}{
		if($5-$4>3){
			print }}' > filtered4_LRRlocus_$1.gff

	# 2. removal of overlap and intron of less than 15 bases if same phase 
	## check the phase change 
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
	}}}}END{print(line)}' filtered4_LRRlocus_$1.gff > filtered5_LRRlocus_$1.tmp
}



#========================================================
#                SCRIPT
#========================================================

tmpdir=$(get_tmp_dir LRRtransfer)
cd $tmpdir

          #------------------------------------------#
          # 1.     Mapping CDS                       #
          #------------------------------------------#

mkdir mapping 
cd mapping

cat $REF_EXONS/$query* > query.fasta
blastn -query query.fasta -subject $TARGET_DNA/$target -outfmt "6 qseqid sseqid qlen length qstart qend sstart send nident pident gapopen" > blastn.tmp

if [[ -s blastn.tmp ]];then
	## processing results
	# 1. remove inconsistent matches from the query (sort by id cds Nip, size of aligenement)
	sort -k1,1 -Vrk4,4 blastn.tmp | gawk 'BEGIN{OFS="\t"}{
			if(NR==1){P5=$5;P6=$6;currentCDS=$1;print}
			else{if(($1==currentCDS && $5>P6-10) || $1!=currentCDS){print;P5=$5;P6=$6;currentCDS=$1}}}' > blastn2.tmp

	# 2. output GFF ///// /!\ \\\\\ ATTENTION, the determination of the chromosome depends on the nomenclature of the target
	sort -Vk7,7 blastn2.tmp | gawk 'BEGIN{OFS="\t";cds=1}{
			if(NR==FNR){
				strand[$1]=$3}
			else{
				split($1,Qid,":");
				split($2,Tid,"_");
				if(length(Tid)==3){chr=Tid[2]}else{chr=Tid[1]};
				pos=Tid[length(Tid)];
				if(strand[$2]=="+"){deb=(pos+$7-1);fin=(pos+$8-1)}
				else{deb=(pos-$8+1);fin=(pos-$7+1)};
				if(FNR==1){
					print(chr,"blastCDS","gene","0",fin,".",strand[$2],".","ID="$2";origin="Qid[1]);
					S=$1;P1=$6;P2=$8;
          print(chr,"blastCDS","CDS",deb,fin,".",strand[$2],".","ID="$2":cds"cds);
          cds=cds+1}
				else{
					if(($1==S && $7>P2-10 && $5>P1-10) || ($1!=S)){
            print(chr,"blastCDS","CDS",deb,fin,".",strand[$2],".","ID="$2":cds"cds);
						cds=cds+1;S=$1;P1=$6;P2=$8;}}}}' $pairID - | sed 's/gene/Agene/g' | sort -Vk4,4 | sed 's/Agene/gene/g' > $target.gff

	gawk -F"\t" 'BEGIN{OFS="\t"}{split($9,T,/[=:;]/);if(NR==FNR){if($3=="gene"){max[T[2]]=$5;min[T[2]]=$5}else{if($5>max[T[2]]){max[T[2]]=$5};if($4<min[T[2]]){min[T[2]]=$4}}}else{if($3=="gene"){$4=min[T[2]];$5=max[T[2]]};print}}' $target.gff $target.gff > tmp ; mv tmp $target.gff

	python3 $LRR_SCRIPT/Exonerate_correction.py -f $TARGET_GENOME -g $target.gff > mapping_LRRlocus.tmp


	## 3. blast verification
	python3 $LRR_SCRIPT/format_GFF.py -g mapping_LRRlocus.tmp -o mapping_LRRlocus.tmp2

	python3 $LRR_SCRIPT/Extract_sequences_from_genome.py -f $TARGET_GENOME -g mapping_LRRlocus.tmp2 -o $target.fasta -t cdna 2>/dev/null
	blastx -query $target.fasta -subject $REF_PEP/$query -outfmt "6 qseqid sseqid slen length qstart qend sstart send nident pident gapopen" > blastx.tmp

	if [[ -s blastx.tmp ]];then
		filter_Blastp blastx.tmp blastx2.tmp
    # gff with origin info + blast in comment section 
    gawk -F"\t" 'BEGIN{OFS="\t"}{if(NR==FNR){Nip[$1]=$2;ID[$1]=$10;COV[$1]=($8-$7+1)/$3}else{if($3~/gene/){split($9,T,";");locname=substr(T[1],4);gsub("comment=","",T[2]);$9=T[1];print($0";comment=Origin:"Nip[locname]" / pred:mapping / blast-%-ident:"ID[locname]" / blast-cov:"COV[locname]" / "T[2])}else{print}}}' blastx2.tmp mapping_LRRlocus.tmp2 > mapping_LRRlocus.gff
		idMapping=$(gawk 'NR==1{print($10)}' blastx2.tmp)
		covMapping=$(gawk 'NR==1{print(100*($8-$7+1)/$3)}' blastx2.tmp)
		ScoreMapping=$(echo "0.6*${idMapping}+0.4*${covMapping}" | bc -l)
  else
  ## coller gff mapping_LRRlocus.gff dans le dossier parent
  cp mapping_LRRlocus.tmp2 mapping_LRRlocus.gff
	fi

fi


cd ..


          #------------------------------------------#
          # 2.     Run exonerate cdna2genome         #
          #------------------------------------------#

mkdir exonerateCDNA
cd exonerateCDNA

## annotation files
gawk -F"\t" 'BEGIN{OFS="\t"}{if($3=="gene"){start=1;split($9,T,";");id=substr(T[1],4);filename=id".an"}else{if($3=="CDS"){len=$5-$4+1;print(id,"+",start,len)>>filename;start=start+len}}}' $GFF
chmod +x $query.an

exonerate -m cdna2genome --bestn 1 --showalignment no --showvulgar no --showtargetgff yes --annotation $query.an --query $REF_cDNA/$query --target $TARGET_DNA/$target > LRRlocus_cdna.out

parseExonerate cdna

#Correcting gene model
python3 $LRR_SCRIPT/Exonerate_correction.py -f $TARGET_GENOME -g filtered5_LRRlocus_cdna.tmp > filtered5_LRRlocus_cdna.gff


# extraction, alignement of nucl on query prot
python3 $LRR_SCRIPT/format_GFF.py -g filtered5_LRRlocus_cdna.gff -o filtered6_LRRlocus_cdna.gff
python3 $LRR_SCRIPT/Extract_sequences_from_genome.py -f $TARGET_GENOME -g filtered6_LRRlocus_cdna.gff -o CDNA_predicted_from_cdna.fasta -t cdna 

blastx -query CDNA_predicted_from_cdna.fasta -subject $REF_PEP/$query -outfmt "6 qseqid sseqid slen length qstart qend sstart send nident pident gapopen" > res_predicted_from_cdna.out
filter_Blastp res_predicted_from_cdna.out res_predicted_from_cdna.out2

idCdna2genome=$(gawk 'NR==1{print($10)}' res_predicted_from_cdna.out2)
covCdna2genome=$(gawk 'NR==1{print(100*($8-$7+1)/$3)}' res_predicted_from_cdna.out2)
ScoreCdna2genome=$(echo "0.6*${idCdna2genome}+0.4*${covCdna2genome}" | bc -l)

# gff with origin info + blast in comment section
gawk -F"\t" 'BEGIN{OFS="\t"}{if(NR==FNR){Nip[$1]=$2;ID[$1]=$10;COV[$1]=($8-$7+1)/$3}else{if($3~/gene/){split($9,T,";");locname=substr(T[1],4);gsub("comment=","",T[2]);$9=T[1];print($0";comment=Origin:"Nip[locname]" / pred:cdna2genome / blast-%-ident:"ID[locname]" / blast-cov:"COV[locname]" / "T[2])}else{print}}}' res_predicted_from_cdna.out2 filtered6_LRRlocus_cdna.gff > cdna2genome_LRRlocus.gff

cd ..


          #------------------------------------------#
          # 3.     Run exonerate prot2genome         #
          #------------------------------------------#

mkdir exoneratePROT 
cd exoneratePROT

exonerate -m protein2genome --showalignment no --showvulgar no --showtargetgff yes --query $REF_PEP/$query --target $TARGET_DNA/$target > LRRlocus_prot.out

parseExonerate prot

#Correct PROT
python3 $LRR_SCRIPT/Exonerate_correction.py -f $TARGET_GENOME -g filtered5_LRRlocus_prot.tmp > filtered5_LRRlocus_prot.gff

#### BLAST + add res blast to gff in comment section + method=prot2genome
# extraction, alignment
python3 $LRR_SCRIPT/format_GFF.py -g filtered5_LRRlocus_prot.gff -o filtered6_LRRlocus_prot.gff
python3 $LRR_SCRIPT/Extract_sequences_from_genome.py -f $TARGET_GENOME -g filtered6_LRRlocus_prot.gff -o PROT_predicted_from_prot.fasta -t cdna 

blastx -query PROT_predicted_from_prot.fasta -subject $REF_PEP/$query -outfmt "6 qseqid sseqid slen length qstart qend sstart send nident pident gapopen" > res_predicted_from_prot.out
filter_Blastp res_predicted_from_prot.out res_predicted_from_prot.out2

idProt2genome=$(gawk 'NR==1{print($10)}' res_predicted_from_prot.out2)
covProt2genome=$(gawk 'NR==1{print(100*($8-$7+1)/$3)}'  res_predicted_from_prot.out2)
ScoreProt2genome=$(echo "0.6*${idProt2genome}+0.4*${covProt2genome}" | bc -l)

# gff with origin info + blast in comment section
gawk -F"\t" 'BEGIN{OFS="\t"}{
       if(NR==FNR){
          Nip[$1]=$2;ID[$1]=$10;COV[$1]=($8-$7+1)/$3}
       else{if($3~/gene/){  
              gsub("comment=","");
              split($9,T,";");
              locname=substr(T[1],4);
              $9=T[1];print($0";comment=Origin:"Nip[locname]" / pred:prot2genome / blast-%-ident:"ID[locname]" / blast-cov:"COV[locname]" / "T[2])}else{print}}}' res_predicted_from_prot.out2 filtered6_LRRlocus_prot.gff > prot2genome_LRRlocus.gff


cd ..

          #------------------------------------------#
          # 4.     Parse All Predictions             #
          #------------------------------------------#

if [ $mode == "first" ];then

	if [ -s mapping_LRRlocus.gff ];then 
	  cat mapping/mapping_LRRlocus.gff >> $RES_DIR/annotation_transfer.gff
	  cat mapping/mapping_LRRlocus.gff > $outfile
	elif [ -s filtered6_LRRlocus_cdna.gff ];then 
	  cat exonerateCDNA/cdna2genome_LRRlocus.gff >> $RES_DIR/annotation_transfer.gff
	  cat exonerateCDNA/cdna2genome_LRRlocus.gff > $outfile
	elif [ -s filtered6_LRRlocus_prot.gff ];then 
	  cat exoneratePROT/prot2genome_LRRlocus.gff  >> $RES_DIR/annotation_transfer.gff
	  cat exoneratePROT/prot2genome_LRRlocus.gff > $outfile
	fi
 
elif [ $mode == "best" ];then

	if (( $(echo "$ScoreMapping > $ScoreCdna2genome" |bc -l) )) && (( $(echo "$ScoreMapping > $ScoreProt2genome" |bc -l) ));then 
	  cat mapping/mapping_LRRlocus.gff >> $RES_DIR/annotation_transfer.gff
	  cat mapping/mapping_LRRlocus.gff > $outfile
	elif (( $(echo "$ScoreMapping < $ScoreCdna2genome" |bc -l) )) && (( $(echo "$ScoreProt2genome < $ScoreCdna2genome" |bc -l) ));then 
	  cat exonerateCDNA/cdna2genome_LRRlocus.gff >> $RES_DIR/annotation_transfer.gff
	  cat exonerateCDNA/cdna2genome_LRRlocus.gff > $outfile
	else
	  cat exoneratePROT/prot2genome_LRRlocus.gff >> $RES_DIR/annotation_transfer.gff
	  cat exoneratePROT/prot2genome_LRRlocus.gff > $outfile
	fi

fi

clean_tmp_dir 0 $tmpdir
