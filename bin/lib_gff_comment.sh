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

  #gawk 'BEGIN{OFS=";"}{if($3~/gene/){if(line){print(line)};split($9,T,";");line=substr(T[1],4)";"$7}else{if($3=="CDS"){line=line";"$4";"$5}}}END{print(line)}' ${input_gff}> ${input_gff}_cds_bounds.tbl
  gawk 'BEGIN{OFS=";"}{if($3~/gene/){if(line){print(line)};split($9,T,";");line=substr(T[1],4)"@"$1";"$7}else{if($3=="CDS"){line=line";"$4";"$5}}}END{print(line)}' ${input_gff}> ${input_gff}_cds_bounds.tbl
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

#Add comments regarding LRR family type
function add_family_info {
  local input_gff=$1
  local input_LRR_profiler_classification=$2
  local output_gff_with_LRR_classification=$3

  gawk -F";" 'BEGIN{OFS="\t"} {$1=$1;print $0}' ${input_LRR_profiler_classification} > __LRR_family.tmp

  gawk -F"\t" 'BEGIN{OFS="\t"}{
              if(NR==FNR){
              FAMILY[$1]=$2
              }
              else
              {
                if($3=="gene"){
                  split($9,infos,";");
                  id=substr(infos[1],4);
                  if (id in FAMILY){
                    $9=$9" / Fam="FAMILY[id]
                  }
                }
                print
              }
            }' __LRR_family.tmp ${input_gff} | sed -e 's/Fam=other/Fam=UC/' -e 's/Fam=RLK/Fam=LRR-RLK/' -e 's/Fam=RLP/Fam=LRR-RLP/' -e 's/Fam=NLR/Fam=NBS-LRR/'> ${output_gff_with_LRR_classification}
}

function basic_NC_family_stat {
  local input_gff=$1
  nbG=$(grep -cw gene ${input_gff});
  echo "$nbG genes"
  for m in noStart noStop pbFrameshift unexpectedSplicingSite stopInFrame pbLength; do echo -n "$m :" ;grep -w gene ${input_gff} | grep -c $m; done
  for fam in LRR-RLK LRR-RLP NBS-LRR UC ; do echo -n "$fam :" ;grep -w gene ${input_gff} | grep -c "Fam=$fam"; done
  echo -n "not LRR :"; grep -w gene ${input_gff} | grep -vc "Fam=" 
}

function remove_genes_from_gff {
  local input_gff=$1
  local input_gene_list_file=$2 # 1 id per line
  local output_gff=$3

  awk '
  NR==FNR {
    genes_to_remove[$1]
  }
  {
    # Update variables for lines of type "gene"
    if ($3 == "gene") {
      match($9, /ID=([^;]+)/, arr)
      current_gene_id = arr[1]
      if (current_gene_id in genes_to_remove) {
        skip = 1
      } else {
        skip = 0
      }
    }

    # Extract the ID of the current feature
    match($9, /ID=([^;]+)/, arr)
    feature_id = arr[1]

    # Check if the feature ID matches the current gene ID
    if (feature_id !~ "^"current_gene_id"(_|$)") {
      print "Error: Feature ID " feature_id " does not match current gene ID " current_gene_id > "/dev/stderr"
    } else if (skip == 0) {
      print $0
    }
  }
  ' "$input_gene_list_file" "$input_gff" > "$output_gff"
}