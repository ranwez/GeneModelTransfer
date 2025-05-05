## extract all features of a gene from a sorted gff
extract_gene_from_sortedGFF(){
  local gene_id=$1
  local sorted_gff=$2
  awk -v ID="ID=${gene_id}" -F '[\t;]' '
    BEGIN { OFS="\t"; within_gene=0 }
    {
        if ($3 == "gene") {
            if ($9 == ID) {
                within_gene = 1
                print $0
            } else if (within_gene == 1) {
                exit  # Stop processing when encountering another "gene"
            }
        } 
        else if (within_gene == 1) {
            print $0
        }
    }' "$sorted_gff"    
}


