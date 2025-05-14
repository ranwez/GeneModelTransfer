## Example of a 3-gene cluster that was annotated as a single big gene
# - the expected gff is from the 012025 annotation, which did well at separating the 3 genes
# - the (problematic) LRRt output gff is from the 022025 annotation (after the candidateLoci refacto)

GMT=/storage/replicated/cirad_users/girodollej/2024_LRR/03_scripts/LRRtransfer/GeneModelTransfer
GMT_SINGULARITY_IMAGE=/storage/replicated/cirad/projects/GE2POP/2023_LRR/LRRtransfer_image/LRRtransfer.sif

chr=Chr6B
start=732108000
end=732134000
zone_suffix="_07321XXXXX"
LRRt_suffix="_202502"

gff_LRRt=/storage/replicated/cirad/projects/GE2POP/2023_LRR/SVEVO3_LRR_ANNOT_2025_02/02_OUTPUT_2025_02_06/OUTPUTS_partial/annotate_best_CLVR_06022025_chr.gff
gff_input=/storage/replicated/cirad/projects/GE2POP/2023_LRR/SVEVO3_LRR_ANNOT_2025_01/01_INPUT/04_final_GFF/IRGSP_DWSvevo3January_LRR.gff
exp_gff=/storage/replicated/cirad/projects/GE2POP/2023_LRR/SVEVO3_LRR_ANNOT_2025_01/02_OUTPUT_2025_01_20/OUTPUTS/annot_best.gff
blast=/storage/replicated/cirad/projects/GE2POP/2023_LRR/SVEVO3_LRR_ANNOT_2025_02/02_OUTPUT_2025_02_06/OUTPUTS/blast_refProt.tsv
input_gene_list=${chr}${zone_suffix}_input_genes.list # list of genes used as template for annotating of the zone

module load singularity/3.6.3

extract_gff() {
  local gff=$1
  local chr=$2
  local start=$3
  local end=$4

  awk -v chr=$chr -v start=$start -v end=$end '{
    if ($1 == chr) {
      if ((start <= $4 && $4 <= end) || (start <= $5 && $5 <= end)) {
        print $0
      }
    }
  }' $gff
}

extract_blast() {
  local blast=$1
  local chr=$2
  local start=$3
  local end=$4

  awk -v chr=$chr -v start=$start -v end=$end '{
    if ($2 == chr) {
      if ((start <= $7 && $7 <= end) || (start <= $8 && $8 <= end)) {
        print $0
      }
    }
  }' $blast
}

# Extract LRRtransfer output gff
extract_gff $gff_LRRt $chr $start $end >raw/${chr}${zone_suffix}_LRRt${LRRt_suffix}.gff

# Extract LRRtransfer input gff
singularity exec $GMT_SINGULARITY_IMAGE python $GMT/SCRIPT/CANDIDATE_LOCI/filter_gff_by_geneID.py -g $gff_input -l ${input_gene_list} -o raw/${chr}${zone_suffix}_input.gff

# Extract expected gff
extract_gff $exp_gff $chr $start $end >raw/${chr}${zone_suffix}_exp.gff

# Sort gffs
for gff in raw/*.gff; do
  sorted_gff=$(basename $gff)
  singularity exec $GMT_SINGULARITY_IMAGE python $GMT/SCRIPT/VR/gff_cleaner.py -g ${gff} -o ${sorted_gff} -a
done

# Extract blast
extract_blast $blast $chr $start $end >${chr}${zone_suffix}_tblastn.tsv
