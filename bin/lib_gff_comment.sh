################################################
## useful functions to add comments to GFF file
################################################

#Add comments regarding model protein at the origin of the annotation
function add_origin_info {
  local input_gff=$1
  local input_info_locus_orig=$2
  local output_gff_with_origin=$3
  gawk -F"\t" 'BEGIN{OFS="\t"}{
    if(NR==FNR){
        F[$1]=$2;C[$1]=$3}
    else{
        if($3~/gene/){
            split($9,T,/[;/]/);origin=substr(T[2],16);gsub(" ","",origin);$9=$9" / Origin-Fam:"F[origin]" / Origin-Class:"C[origin]};print}}' ${input_info_locus_orig} ${input_gff} > ${output_gff_with_origin}
}

# compute a tsv file with true/false for each non canonical events (missing start, frameshift etc.)
function compute_NC_alerts {
  local input_gff=$1
  local dna_seq=$2
  local output_alert_NC_info=$3

  gawk 'BEGIN{OFS=";"}{if($3~/gene/){if(line){print(line)};split($9,T,";");line=substr(T[1],4)";"$7}else{if($3=="CDS"){line=line";"$4";"$5}}}END{print(line)}' ${input_gff}> ${input_gff}_cds_bounds.tbl
  python3 ${LRR_SCRIPT}/Canonical_gene_model_test.py -f ${dna_seq} -t ${input_gff}_cds_bounds.tbl -o ${output_alert_NC_info}

}

# add comment regarding non canonical events and colour information
function add_comment_NC {
  local input_gff=$1
  local alert_NC_info=$2
  local output_commented_gff=$3
    gawk -F"\t" 'BEGIN{OFS="\t"}{
              if(NR==FNR){
                COMMENT[$2]="";
                if($3=="True"){COMMENT[$2]=COMMENT[$2]" / noStart"};
                if($4=="True"){COMMENT[$2]=COMMENT[$2]" / noStop"};
                if($5=="True"){COMMENT[$2]=COMMENT[$2]" / pbFrameshift"};
                if($6=="True"){COMMENT[$2]=COMMENT[$2]" / unexpectedSplicingSite"};
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
                if ($9 ~ /comment=/) { sub(/(comment=[^;]*)/, "& " COMMENT[id], $9);} 
                else { $9 = $9 "; comment=" COMMENT[id];}
                $9=$9";color="genecolor; 
                };
                print}}' ${alert_NC_info} ${input_gff} > ${output_commented_gff}

}
