#!/bin/bash
#========================================================
# PROJET : LRRtransfer
# SCRIPT : genePrediction.sh
# AUTHOR : Celine Gottin & Thibaud Vicat & Vincent Ranwez
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
#             o $11 : Path toward LRR script  directory

#========================================================
#                Environment & variables
#========================================================

#set -euo pipefail
set -eu
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
LRR_SCRIPT=${11}

REF_PEP=$LRRome/REF_PEP
REF_EXONS=$LRRome/REF_EXONS
REF_cDNA=$LRRome/REF_cDNA



target=$(cat $pairID | cut -f1)
query=$(cat $pairID | cut -f2)
pairStrand=$(cat $pairID | cut -f3)

mmseqs="mmseqs"

#========================================================
#                        Functions
#========================================================

source $LRR_SCRIPT/../bin/lib_gff_comment.sh
source $LRR_SCRIPT/../bin/lib_tmp_dir.sh

function exonerate_GFF_from_similarity {
	local input_exonerate_res=$1
	local output_exonerate_gff=$2
	gawk -F"\t" 'BEGIN{OFS=FS}{if($7=="+" && ($3=="gene" || $3=="similarity")){print}}' ${input_exonerate_res} > ${input_exonerate_res}.tmp ##prot,cdna
	gawk -F"\t" 'BEGIN{OFS=FS}{
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
	}' ${input_exonerate_res}.tmp > ${output_exonerate_gff}
}

function exonerate_GFF_from_exon {
	local input_exonerate_res=$1
	local output_exonerate_gff=$2
	gawk -F"\t" 'BEGIN{OFS=FS}{if($7=="+"){ if ($3=="gene") {print} ; if ($3=="exon"){ $3="CDS";print} } }' ${input_exonerate_res} > ${output_exonerate_gff}
}

function parseExonerate {
	local exoneRate_input=$1
	local gff_output=$2
	local gff_from=$3 # "similarity" or "exon"
	
	if [[ ${gff_from} == "similarity" ]]; then
		exonerate_GFF_from_similarity ${exoneRate_input} ${exoneRate_input}.gff
	else
		exonerate_GFF_from_exon ${exoneRate_input} ${exoneRate_input}.gff
	fi

	## define ID and Parent and strand
	gawk -F"\t" 'BEGIN{OFS=FS}{
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
	}' $filtered_candidatsLRR ${exoneRate_input}.gff > ${exoneRate_input}_filtered.gff

	## Eliminate gene redundancy
	gawk -F"\t" 'BEGIN{OFS=FS}{
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
	}' ${exoneRate_input}_filtered.gff ${exoneRate_input}_filtered.gff > ${exoneRate_input}_filtered2.gff

	sed 's/gene/Agene/g' ${exoneRate_input}_filtered2.gff | sort -k1,1 -Vk4,4 -k3,3 | sed 's/Agene/gene/g' > ${exoneRate_input}_filtered3.gff

	## eliminate redundancy and overlap of CDS if on the same phase (we check the phase by $4 if strand + and by $5 if strand -)
	# 1. remove cds included in other cds and glued cds

	#cat filtered3_LRRlocus_$1.gff > filtered4_LRRlocus_$1.gff
	gawk -F"\t" 'BEGIN{OFS=FS}{
		if($3~/gene/){
			print;
			lim=0}
		else{
			if($5>lim){
				print;
				lim=$5}}}' ${exoneRate_input}_filtered3.gff | gawk 'BEGIN{OFS="\t"; line=""}{
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
			print }}' > ${exoneRate_input}_filtered4.gff

	# 2. removal of overlap and intron of less than 15 bases if same phase 
	## check the phase change 
	gawk -F"\t" 'BEGIN{OFS=FS;line=""}{
                if($3~/gene/){if(line!=""){print(line)};print;p=0}
                else{
                   if(p==0){
                       p=1;line=$0;start=$4;stop=$5;mod=($(5)+1)%3}
                   else{
                       if($5<stop+15 || $4<start+15){$4=start;stop=$5;line=$0;mod=($(5)+1)%3}
                       else{
                          if(($4>stop+15) || $4%3!=mod){print(line);line=$0;start=$4;stop=$5;mod=($(5)+1)%3}
                          else{$4=start;stop=$5;line=$0;mod=($(5)+1)%3}
	}}}}END{if (NR > 0) {print(line)}}' ${exoneRate_input}_filtered4.gff > ${gff_output}
}

function parse_blast_to_gff {
	local blast_input=$1
	local gff_output=$2
	# 1. remove inconsistent matches from the query (sort by id cds Nip, size of aligenement)
	sort -k1,1 -Vrk4,4 ${blast_input} | gawk -F"\t" 'BEGIN{OFS=FS}{
			if(NR==1){P5=$5;P6=$6;currentCDS=$1;print}
			else{if(($1==currentCDS && $5>P6-10) || $1!=currentCDS){print;P5=$5;P6=$6;currentCDS=$1}}}' > blastn2.tmp

	## VR modif gene coordinates are min max of CDS coordinates 
	# 2. output GFF ///// /!\ \\\\\ ATTENTION, the determination of the chromosome depends on the nomenclature of the target
	sort -Vk7,7 blastn2.tmp | gawk -F"\t" 'BEGIN{OFS=FS;cds=1}{
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
					geneDeb=deb;
					geneFin=fin;
					geneStrand=strand[$2];
					geneId=$2;
					geneOrigin=Qid[1];
					S=$1;P1=$6;P2=$8;
		 			print(chr,"blastCDS","CDS",deb,fin,".",strand[$2],".","ID="$2":cds"cds);
		  			cds=cds+1
				}
				else{
					if(($1==S && $7>P2-10 && $5>P1-10) || ($1!=S)){
						if( deb < geneDeb) {geneDeb=deb};
						if(fin > geneFin){geneFin=fin};
						print(chr,"blastCDS","CDS",deb,fin,".",strand[$2],".","ID="$2":cds"cds);
						cds=cds+1;S=$1;P1=$6;P2=$8;}
					}
				}
			} 
			END{print(chr,"blastCDS","gene",geneDeb,geneFin,".",geneStrand,".","ID="geneId";origin="geneOrigin);}' $pairID - | sed 's/gene/Agene/g' | sort -Vk4,4 | sed 's/Agene/gene/g' > ${target}_1.gff


	gawk -F"\t" 'BEGIN{OFS=FS}{if($4>$5){max=$4;$4=$5;$5=max};print}' ${target}_1.gff > ${gff_output} 
}

function try_merging_CDS {
	local input_draft_gff_m=$1
	local dna_seq_file=$2
	local output_draft_mergedCDS_gff=$3

	# merge CDS that are nearby and in the same RF
	gawk -v dmax=25 'BEGIN{OFS="\t";p=0}{
		if($3~/CDS/){
			if(p==0){
			line=$0;P4=$4;P5=$5;p=1}
			else{
			if($4<=P5+dmax && ($4-P5-1)%3==0){
				$4=P4;line=$0;P4=$4;P5=$5}
			else{print(line);line=$0;P4=$4;P5=$5}
			}
		}else{
			if(p!=0){print(line)};P4=0;P5=0;p=0;print}
		}END{if(p!=0){print(line)}}' ${input_draft_gff_m} > ${output_draft_mergedCDS_gff}

	nbCDS_origin=$(grep -w "CDS" ${input_draft_gff_m} | grep -c ".");
	nbCDS_merged=$(grep -w "CDS" ${output_draft_mergedCDS_gff} | grep -c ".");

	merged=0;
	if (( $nbCDS_merged < $nbCDS_origin )); then
		seqAA_origin=$(python3 $LRR_SCRIPT/Extract_sequences_from_genome.py -f ${dna_seq_file} -g ${input_draft_gff_m} -t FSprot 2>/dev/null)
		seqAA_merged=$(python3 $LRR_SCRIPT/Extract_sequences_from_genome.py -f ${dna_seq_file} -g ${output_draft_mergedCDS_gff} -t FSprot 2>/dev/null)
		nb_stop_origin=$(echo "$seqAA_origin" | grep -o '\*' | grep -c ".");
		nb_stop_merged=$(echo "$seqAA_merged" | grep -o '\*' | grep -c ".");
		if (( $nb_stop_merged <= $nb_stop_origin )); then
			merged=$(echo "$nbCDS_origin - $nbCDS_merged" | bc)
			sed -i "/\bgene\b/s/$/;merge_cds=${merged}/" ${output_draft_mergedCDS_gff}
		fi
	fi
	
	if (( merged==0 )); then
		cp ${input_draft_gff_m} ${output_draft_mergedCDS_gff};
	fi
}

function improve_annot {
	local input_draft_gff=$1
	local dna_seq_file=$2
	local output_improved_gff=$3
	
	if [[ -s ${input_draft_gff} ]];then
		### using target region is faster than using full genome
		gff_genome_to_target ${input_draft_gff} ${input_draft_gff}_onTarget
		try_merging_CDS  ${input_draft_gff}_onTarget  ${dna_seq_file} ${input_draft_gff}_cdsmerged_onTarget
		python3 ${LRR_SCRIPT}/Exonerate_correction.py -f ${dna_seq_file} -g ${input_draft_gff}_cdsmerged_onTarget > ${input_draft_gff}_tmp1_onTarget.gff
		gff_target_to_genome ${input_draft_gff}_tmp1_onTarget.gff ${input_draft_gff}_tmp1.gff
		###
		gawk -F"\t" 'BEGIN{OFS=FS}{if($4>$5){max=$4;$4=$5;$5=max};print}' ${input_draft_gff}_tmp1.gff > ${input_draft_gff}_tmp2.gff
		python3 $LRR_SCRIPT/format_GFF.py -g ${input_draft_gff}_tmp2.gff -o ${input_draft_gff}_tmp3.gff
		gawk -F"\t" 'BEGIN{OFS=FS}{if($4>$5){max=$4;$4=$5;$5=max};print}' ${input_draft_gff}_tmp3.gff > ${output_improved_gff}
	else
		touch ${output_improved_gff}
	fi
}

function gff_genome_to_target {
	local input_gff_on_genome=$1
	local output_gff_on_target=$2

	target_start=$( echo ${target}| sed -e 's/.*_[0]*//');
	if [[ $pairStrand == '+' ]]; then
		gawk -F"\t" -v tstart=${target_start} -v chrT=${target} 'BEGIN{OFS=FS} {$1=chrT; $4=$4 - tstart +1 ; $5=$5 - tstart +1 ; print}' ${input_gff_on_genome} > ${output_gff_on_target}
	else
		gawk  -F"\t" -v tstart=${target_start} -v chrT=${target} 'BEGIN{OFS=FS} {
			bound1=tstart-$4 +1; bound2= tstart -$5 +1; 
			$1=chrT; $4=bound2; $5=bound1;
			if ($7 == "-") $7 = "+"; else $7 = "-";
			print}' ${input_gff_on_genome} | sort -s -n -k4,4 > ${output_gff_on_target}
	fi
}

function gff_target_to_genome {
	local input_gff_on_genome=$1
	local output_gff_on_target=$2

	target_start=$( echo ${target}| sed -e 's/.*_[0]*//');
	chr_genome=$( echo ${target}| sed -e 's/_.*//');
	if [[ $pairStrand == '+' ]]; then
		gawk -F"\t" -v tstart=${target_start} -v chrG=${chr_genome} 'BEGIN{OFS=FS} {$1=chrG; $4=tstart +$4 -1 ; $5= tstart +$5 -1 ; print}' ${input_gff_on_genome} > ${output_gff_on_target}
	else
		gawk -F"\t" -v tstart=${target_start} -v chrG=${chr_genome} 'BEGIN{OFS=FS} {
			bound1=tstart-$4 +1; bound2= tstart -$5 +1; 
			$1=chrG; $4=bound2; $5=bound1;
			if ($7 == "-") $7 = "+"; else $7 = "-";
			print}' ${input_gff_on_genome} | sort -s -n -k4,4  > ${output_gff_on_target}
	fi
}

function non_canonical_penalty {
	local input_gff=$1
	local dna_seq=$2
	local output_alert_NC_info=$3

	compute_NC_alerts ${input_gff} ${dna_seq} ${output_alert_NC_info}
	nb_issues=$(tail -1 ${output_alert_NC_info} | grep -w -o 'True' | grep -c '.' )
	penalty=$((nb_issues * 2))
	echo $penalty

}

function evaluate_annotation {
	local input_gff=$1
	local output_alert_NC_info=$2
	local cov_denom=$3
	# return identity coverage and combined (with and without NC penalty) score of the newly annotated prot w.r.t the model prot
	res="0 0 0 0" 
	if [[ -s ${input_gff} ]];then
		# extract the predicted protein sequence corresponding to the input gff
		gff_genome_to_target ${input_gff}  ${input_gff}_onTarget
		python3 $LRR_SCRIPT/Extract_sequences_from_genome.py -f ${TARGET_DNA}/$target -g ${input_gff}_onTarget -o ${input_gff}_prot.fasta -t FSprot 2>/dev/null
	
		# evaluate the similarity between the newly predicted protein and the reference one
		bestHit=$(blastp -query $REF_PEP/$query -subject ${input_gff}_prot.fasta -outfmt "6 length qlen slen pident bitscore" | sort -n -k 5,5 | tail -1) 
		if [[ -n "$bestHit" ]];then
			penalty=$(non_canonical_penalty ${input_gff}_onTarget ${TARGET_DNA}/$target ${output_alert_NC_info})
			res=$(echo "$bestHit" | gawk -F"\t" -v penalty=$penalty -v covDenom=${cov_denom} '{
				maxL = ($2>$3 ? $2 : $3);covFull=(100*$1)/maxL;
				ident=$4; cov=(100*$1)/(covDenom); score=0.6*ident+0.4*cov; scoreNC=0.6*(ident-penalty)+0.4*cov;
				print ident,covFull,score,scoreNC}')
		fi
	fi
	echo  $res
}

function set_gff_comments {
	local input_gff=$1
	local infoLocus=$2
	local cov_denom=$3
	local method=$4
	local updated_gff=$5


	read ident cov score scoreNC< <(evaluate_annotation ${input_gff} ${input_gff}_NC_alert.tsv ${cov_denom})
	
	if [[ -s ${input_gff} ]]; then
		# add scoring comments
		gawk -F"\t" -v query_id=$query -v target_id=$target -v ident=$ident -v cov=$cov -v method=$method -v score=$score -v scoreNC=$scoreNC 'BEGIN{OFS=FS}{
				if($3~/gene/){
					split($9,T,";");
					gsub("comment=","",T[2]);
					$9="ID="target_id";comment=Origin:"query_id" / pred:"method" / prot-%-ident:"ident" / prot-%-cov:"cov" / score:"score" / scoreNC:"scoreNC" / "T[2]};
				print}' ${input_gff} > ${input_gff}_w_scoring
		
		# add origin details and NC comments 
		add_origin_info ${input_gff}_w_scoring ${infoLocus} ${input_gff}_w_scoring_origin
		add_comment_NC ${input_gff}_w_scoring_origin ${input_gff}_NC_alert.tsv ${updated_gff}
	else
		touch ${updated_gff}
	fi
	echo $score
}

function mrna_length {
    local input_gff=$1
    grep -w "CDS" "$input_gff" | gawk -F"\t" 'BEGIN{L=0} {L += ($5>$4 ? $5-$4 : $4-$5)} END{print (L+3)/3}' | cut -f1 -d"."
}

#========================================================
#                SCRIPT
#========================================================

tmpdir=$(get_tmp_dir LRRtransfer)
cd $tmpdir

          #---------------------------------------------------------#
          #      Build draft gff estimate for each method           #
          #---------------------------------------------------------#

methods="mapping cdna2genome cdna2genomeExon prot2genome prot2genomeExon"
for method in $(echo $methods); do mkdir $method; done

cd mapping
#TODO various separators are use to separate cds numbers, should be improved
cat $REF_EXONS/${query}[:_-]* > query.fasta
blastn -query query.fasta -subject $TARGET_DNA/$target -outfmt "6 qseqid sseqid qlen length qstart qend sstart send nident pident gapopen" > blastn.tmp
if [[ -s blastn.tmp ]]; then
	parse_blast_to_gff blastn.tmp ${target}_draft.gff
fi

cd ../cdna2genome
grep $query $GFF | gawk -F"\t" 'BEGIN{OFS=FS}{if($3=="gene"){start=1;split($9,T,";");id=substr(T[1],4);filename=id".an"}else{if($3=="CDS"){len=$5-$4+1;print(id,"+",start,len)>>filename;start=start+len}}}'
chmod +x $query.an
exonerate -m cdna2genome --bestn 1 --showalignment no --showvulgar no --showtargetgff yes --annotation $query.an --query $REF_cDNA/$query --target $TARGET_DNA/$target > LRRlocus_cdna.out
if [[ -s LRRlocus_cdna.out ]]; then
	parseExonerate LRRlocus_cdna.out ${target}_draft.gff "similarity"
	parseExonerate LRRlocus_cdna.out ../cdna2genomeExon/${target}_draft.gff "exon"
fi

cd ../prot2genome
exonerate -m protein2genome --showalignment no --showvulgar no --showtargetgff yes --query $REF_PEP/$query --target $TARGET_DNA/$target > LRRlocus_prot.out
if [[ -s LRRlocus_prot.out ]]; then
	parseExonerate LRRlocus_prot.out ${target}_draft.gff "similarity"
	parseExonerate LRRlocus_prot.out ../prot2genomeExon/${target}_draft.gff "exon"
fi

cd ..

          #--------------------------------------------------------------#
          #    Build gff for each method and search the best             #
          #--------------------------------------------------------------#



# use draft gff to get prot length estimate and better coverage estimation for scoring
LG_REF=0
for method in $(echo $methods);do 
	if [[ -s ${method}/${target}_draft.gff ]];then
		lg=$(mrna_length ${method}/${target}_draft.gff)
		if (( $lg > $LG_REF )); then LG_REF=$lg ;fi
	fi
done

bestScore=0; bestGff="";
for method in $(echo $methods);do 
	cd $method
	if [[ -s ${target}_draft.gff ]];then
		improve_annot  ${target}_draft.gff ${TARGET_DNA}/$target ${target}.gff 
		scoreMethod=$(set_gff_comments ${target}.gff $infoLocus $LG_REF "$method" ${outfile}_${method}.gff )
		if (( $(echo "$scoreMethod > $bestScore" | bc -l) )); then 
			bestGff=$(realpath ${outfile}_${method}.gff)
			bestScore=$scoreMethod
		fi
	else : # or do nothing to raise the error ?
		touch  ${outfile}_${method}.gff 
	fi
	cd ..
done

if [ $mode == "best" ];then
	if [[ -s $bestGff ]]; then
		cp $bestGff ${outfile}_best.gff
	else : # or do nothing to raise the error ?
		touch  ${outfile}_best.gff
	fi
fi

echo $tmpdir
clean_tmp_dir 0 $tmpdir
