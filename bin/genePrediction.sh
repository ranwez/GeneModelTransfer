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


ALLOWED_EXONERATE_ERRORS=( # /!\ Special characters such as "*" or "[" "]" need to be escaped
  "^\*\* FATAL ERROR \*\*: Initial HSP score \[-1\] less than zero"
)

#========================================================
#                        Functions
#========================================================

source $LRR_SCRIPT/../bin/lib_gff_comment.sh
source $LRR_SCRIPT/../bin/lib_tmp_dir.sh
source $LRR_SCRIPT/../bin/lib_gff_utils.sh

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

function exonerate_GFF_from_cds {
	local input_exonerate_res=$1
	local output_exonerate_gff=$2
	gawk -F"\t" 'BEGIN{OFS=FS}{if($7=="+"){ if ($3=="gene") {print} ; if ($3=="cds"){ $3="CDS";print} } }' ${input_exonerate_res} > ${output_exonerate_gff}
}

function parseExonerate {
	local exoneRate_input=$1
	local gff_output=$2
	local gff_from=$3 # "similarity" or "exon" or "cds"

	if [[ ${gff_from} == "similarity" ]]; then
		exonerate_GFF_from_similarity ${exoneRate_input} ${exoneRate_input}.gff
	else
		if [[ ${gff_from} == "exon" ]]; then
			exonerate_GFF_from_exon ${exoneRate_input} ${exoneRate_input}.gff
		else
			exonerate_GFF_from_cds ${exoneRate_input} ${exoneRate_input}.gff
		fi
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

	local nbCDS_origin=$(grep -w "CDS" ${input_draft_gff_m} | grep -c ".");
	local nbCDS_merged=$(grep -w "CDS" ${output_draft_mergedCDS_gff} | grep -c ".");

	local merged=0;
	if (( $nbCDS_merged < $nbCDS_origin )); then
		local seqAA_origin=$(python3 $LRR_SCRIPT/Extract_sequences_from_genome.py -f ${dna_seq_file} -g ${input_draft_gff_m} -t FSprot &>> log_Extract_sequences_from_genome.txt)
		local seqAA_merged=$(python3 $LRR_SCRIPT/Extract_sequences_from_genome.py -f ${dna_seq_file} -g ${output_draft_mergedCDS_gff} -t FSprot &>> log_Extract_sequences_from_genome.txt)
		local nb_stop_origin=$(echo ${seqAA_origin} | grep -o '\*' | grep -c ".");
		local nb_stop_merged=$(echo ${seqAA_merged} | grep -o '\*' | grep -c ".");
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
	local lg_max=$4
	# return identity coverage and combined (with and without NC penalty) score of the newly annotated prot w.r.t the model prot
	res="0 0 0 0"
	if [[ -s ${input_gff} ]];then
		# extract the predicted protein sequence corresponding to the input gff
		gff_genome_to_target ${input_gff}  ${input_gff}_onTarget
		python3 $LRR_SCRIPT/Extract_sequences_from_genome.py -f ${TARGET_DNA}/$target -g ${input_gff}_onTarget -o ${input_gff}_prot.fasta -t FSprot &>> log_Extract_sequences_from_genome.txt
		# detect issues
		penalty=$(non_canonical_penalty ${input_gff}_onTarget ${TARGET_DNA}/$target ${output_alert_NC_info})
		# evaluate the similarity between the newly predicted protein and the reference one
		#bestHit=$(blastp -query $REF_PEP/$query -subject ${input_gff}_prot.fasta -outfmt "6 length qlen slen pident positive bitscore" | sort -n -k 6,6 | tail -1)
		#if [[ -n "$bestHit" ]];then
		#	res=$(echo "$bestHit" | gawk -F"\t" -v penalty=$penalty -v covDenom=${cov_denom} '{
		#		maxL = ($2>$3 ? $2 : $3);covFull=(100*$1)/maxL;
		#		ident=$4; positive=$5; score=positive/covDenom; scoreNC=score-(0.01*penalty);
		#		print ident,covFull,score,scoreNC}')
		#fi
		read nbPositives nbIdentity RawNbPositives pcHomology< <(python3 $LRR_SCRIPT/VR/prot_prediction_scoring.py ${input_gff}_prot.fasta $REF_PEP/$query)
		lgQuery=$(grep -v ">" $REF_PEP/$query | sed 's,\n,,' | wc -c)
		if (( $RawNbPositives > 0 )) ; then
			res=$(awk -v nbPos=${nbPositives} -v covDenom=${cov_denom} -v nbIdent=${nbIdentity} -v lgmax=$lg_max -v lgQuery=$lgQuery 'BEGIN{
				score=nbPos/covDenom;
				pident=nbIdent/covDenom;
				pidentQ=nbIdent/lgQuery;
				cov=nbPos/lgmax;
				covQ=nbPos/lgQuery;
				scoreNC=score-(0.01*penalty);
				print pident,pidentQ,cov,covQ,score,scoreNC}')
		fi
	fi
	echo  $res
}

function set_gff_comments {
	local input_gff=$1
	local infoLocus=$2
	local cov_denom=$3
	local lg_max=$4
	local method=$5
	local updated_gff=$6

	local score=0
	if [[ -s ${input_gff} ]]; then
		read ident identQ cov covQ score scoreNC< <(evaluate_annotation ${input_gff} ${input_gff}_NC_alert.tsv ${cov_denom} ${lg_max})
		# add scoring comments
		gawk -F"\t" -v query_id=$query -v target_id=$target -v ident=$ident -v identQ=$identQ -v cov=$cov -v covQ=$covQ -v method=$method -v score=$score -v scoreNC=$scoreNC 'BEGIN{OFS=FS}{
				if($3~/gene/){
					split($9,T,";");
					gsub("comment=","",T[2]);
					$9="ID="target_id";comment=Origin:"query_id" / pred:"method" / prot-%-ident:"ident" / prot-%-identQ:"identQ" / prot-%-cov:"cov" / prot-%-covQ:"covQ" / score:"score" / scoreNC:"scoreNC" / "T[2]};
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

function get_new_template {
	local input_fasta=$1
	bestHit=$(blastp -query ${input_fasta} -subject $REF_PEP/../REF_proteins.fasta  -outfmt "6 length qlen slen pident positive bitscore sseqid" | sort -n -k 6,6 | tail -1)
	best_template=$(echo $bestHit | awk '{print $7}')
	echo ${best_template}
}

function run_exonerate() {
	trap "set -eu" RETURN
	local out_file="$1"
    shift
    local cmd="$*"
    local error_message
	set +eu
    error_message=$(${cmd} 2>&1 > "$out_file")

    # If exonerate failed, check if the error is allowed
    if [[ $? -ne 0 ]]; then
        for allowed in "${ALLOWED_EXONERATE_ERRORS[@]}"; do
            if echo "$error_message" | grep -E "$allowed" > /dev/null; then
                echo "Warning: Exonerate encountered known error:"
                echo "$error_message" | head -1
				echo "When running:"
				echo ${cmd}
                echo "Continuing..."
                return 0
            fi
        done
        echo "Exonerate failed with an unexpected error:"
        echo "$error_message"
        exit 1
    fi
}


#========================================================
#                SCRIPT
#========================================================



methods="mapping cdna2genome cdna2genomeExon cds2genome cds2genomeExon prot2genome prot2genomeExon"
lg=$(grep -v ">" $TARGET_DNA/$target | wc -c)
if (( $lg <= 1 )); then
	for method in $(echo $methods);do
		touch ${outfile}_${method}.gff
	done
	touch ${outfile}_best.gff
	echo " WARNING empty target file $TARGET_DNA/$target"
	exit 0
fi

tmpdir=$(get_tmp_dir LRRtransfer)
echo $tmpdir
cd $tmpdir
          #---------------------------------------------------------#
          #      Build draft gff estimate for each method           #
          #---------------------------------------------------------#


for method in $(echo $methods); do mkdir $method; done

cp $REF_cDNA/$query query_cDNA.fasta
cp $REF_PEP/$query query_PEP.fasta

cd mapping
#TODO various separators are use to separate cds numbers, should be improved
cat $REF_EXONS/${query}[:_-]* > query.fasta
blastn -query query.fasta -subject $TARGET_DNA/$target -outfmt "6 qseqid sseqid qlen length qstart qend sstart send nident pident gapopen" > blastn.tmp
if [[ -s blastn.tmp ]]; then
	parse_blast_to_gff blastn.tmp ${target}_draft.gff
fi


cd ../cdna2genome
extract_gene_from_sortedGFF $query $GFF | gawk -F"\t" 'BEGIN{OFS=FS}{if($3=="gene"){start=1;split($9,T,";");id=substr(T[1],4)}else{if($3=="CDS"){len=$5-$4+1;print(id,"+",start,len);start=start+len}}}' > query.an
chmod +x query.an
run_exonerate LRRlocus_cdna.out exonerate -m cdna2genome --bestn 1 --showalignment no --showvulgar no --showtargetgff yes --annotation query.an --query ../query_cDNA.fasta --target $TARGET_DNA/$target

if [[ -s LRRlocus_cdna.out ]]; then
	parseExonerate LRRlocus_cdna.out ${target}_draft.gff "similarity"
	parseExonerate LRRlocus_cdna.out ../cdna2genomeExon/${target}_draft.gff "exon"
fi

cd ../cds2genome
# in REF_cDNA we only got the coding fragment (concatenation of CDS in the + direction) so the query.an is simply:
# query_name + 1 query_length
query_lg=$(sed 's/[[:space:]]//g' ../query_cDNA.fasta  | sed '/^>/d' | wc -c)
echo -e "$query\t+\t1\t${query_lg}" > query.an
chmod +x query.an
run_exonerate LRRlocus_cds.out exonerate -m coding2genome --bestn 1 --showalignment no --showvulgar no --showtargetgff yes --annotation query.an --query ../query_cDNA.fasta --target $TARGET_DNA/$target --refine full
if [[ -s LRRlocus_cds.out ]]; then
	parseExonerate LRRlocus_cds.out ${target}_draft.gff "similarity"
	parseExonerate LRRlocus_cds.out ../cds2genomeExon/${target}_draft.gff "cds"
fi

cd ../prot2genome
run_exonerate LRRlocus_prot.out exonerate -m protein2genome --showalignment no --showvulgar no --showtargetgff yes --query ../query_PEP.fasta --target $TARGET_DNA/$target
if [[ -s LRRlocus_prot.out ]]; then
	parseExonerate LRRlocus_prot.out ${target}_draft.gff "similarity"
	parseExonerate LRRlocus_prot.out ../prot2genomeExon/${target}_draft.gff "cds"
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

# to get coverage use max of template and predicted protein length
lg_max=$(sed 's/[[:space:]]//g' $REF_PEP/$query  | sed '/^>/d' | wc -c)
if(( LG_REF > $lg_max )); then lg_max=$LG_REF; fi

bestScore=0; bestGff=""; bestProtFasta=""
for method in $(echo $methods);do
	cd $method
	if [[ -s ${target}_draft.gff ]];then
		improve_annot  ${target}_draft.gff ${TARGET_DNA}/$target ${target}.gff
		scoreMethod=$(set_gff_comments ${target}.gff $infoLocus $LG_REF $lg_max "$method" ${outfile}_${method}.gff )
		#if (( $(echo "$scoreMethod > $bestScore" | bc -l) )); then
		# comparison of numbers in scientific notation does not work with bc -l so we use awk instead:
		isBetter=$(echo -e "$scoreMethod\t$bestScore" | awk '{if ($1 > $2){print 1} else {print 0}}')
		if (( $isBetter == 1 )); then
			bestGff=$(realpath ${outfile}_${method}.gff)
			bestScore=$scoreMethod
			bestProtFasta=$(realpath ${target}.gff_prot.fasta)
		fi
	else : # or do nothing to raise the error ?
		touch  ${outfile}_${method}.gff
	fi
	cd ..
done

#echo $bestScore
#cat $bestProtFasta
#echo ""
if [ $mode == "best2rounds" ];then
	cat $REF_PEP/$query | sed -e 's/>/>temp1_/' > ${outfile}_protinfo.fasta
	cat ${bestProtFasta}| sed -e 's/>/>pred1_/' >> ${outfile}_protinfo.fasta
	if [[ -s $bestGff ]]; then
		cp $bestGff ${outfile}_best1.gff
		new_template=$(get_new_template $bestProtFasta)
		# if we find a better template prot use it in best mode
		if [ $new_template != $query ]; then
			for method in $(echo $methods);do
				rm ${outfile}_${method}.gff
			done
			echo -e "$target\t${new_template}\t${pairStrand}" > ${outfile}_pairID
			$0 ${outfile}_pairID $2 $3 $4 $5 $6 $7 $8 $9 best ${LRR_SCRIPT}
		# else directly switch back to best mode
		else :
			mode="best"
		fi
	else :
		touch ${outfile}_best1.gff
		mode="best"
	fi
fi

if [ $mode == "best" ];then
	if [[ -s $bestGff ]]; then
		cp $bestGff ${outfile}_best.gff
		if [[ -s ${outfile}_pairID ]]; then
			cat $REF_PEP/$query | sed -e 's/>/>temp2_/'>> ${outfile}_protinfo.fasta
			cat ${bestProtFasta} | sed -e 's/>/>pred2_/' >> ${outfile}_protinfo.fasta
		fi
	else : # or do nothing to raise the error ?
		touch  ${outfile}_best.gff
	fi
fi

echo $tmpdir
#clean_tmp_dir 0 $tmpdir
