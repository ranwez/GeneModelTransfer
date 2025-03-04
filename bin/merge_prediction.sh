set -euo pipefail

inDir=$1
LRR_SCRIPT=$2
method=$3
prefix=$4
outDir=$5

cat ${inDir}/annotate_one_*_${method}.gff > ${outDir}/${method}_tmp
python3 ${LRR_SCRIPT}/VR/gff_cleaner.py -a -p ${prefix} -g ${outDir}/${method}_tmp -o ${outDir}/annot_${method}.gff > ${outDir}/annot_${method}_cleaning.log
awk '{if ($1 != "contig"){print $0}}' ${outDir}/annot_${method}.gff > ${outDir}/annot_${method}_chr.gff
rm ${outDir}/${method}_tmp
