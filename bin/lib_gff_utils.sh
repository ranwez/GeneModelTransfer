## extract all features of a gene from a sorted gff
extract_gene_from_sortedGFF() {
  local gene_id=$1
  local sorted_gff=$2

  awk -v gene_id="$gene_id" -F '\t' '
    BEGIN {
      within_gene = 0
    }

    $3 == "gene" {
      current_id = ""
      n = split($9, attributes, ";")

      for (i = 1; i <= n; i++) {
        if (attributes[i] ~ /^ID=/) {
          current_id = substr(attributes[i], 4)
          break
        }
      }

      if (current_id == gene_id) {
        within_gene = 1
        print
        next
      }

      if (within_gene) {
        exit
      }
    }

    within_gene {
      print
    }
  ' "$sorted_gff"
}
