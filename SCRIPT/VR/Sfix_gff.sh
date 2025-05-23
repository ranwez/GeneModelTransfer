# wk -F'\t' 'BEGIN{OFS=FS}{$1="chr2A"; print $0}' EXP_chr2A_LRRlocus_predicted_best_VR_NC_20240209.gff > EXP_chr2A_LRRlocus_curated_20240209.gff
# ./Sfix_gff.sh EXP_chr2A_LRRlocus_curated_20240209.gff curated_chr2A_LRRlocus_20240209_cleaned.gff
input_gff=$1
prefix=$2
cleaned_gff=$3
encode=$( file $input_gff | cut -f2 -d" "); echo $encode
if [[ $encode -eq "ISO-8859" ]]; then
	iconv -f ISO-8859-1 -t UTF-8 ${input_gff} > ${input_gff}.utf8
else cp ${input_gff} ${input_gff}.utf8
fi

SCRIPTDIR=$(realpath $(dirname "${BASH_SOURCE[0]}") )
python3 ${SCRIPTDIR}/gff_cleaner.py -a -p ${prefix} -g ${input_gff}.utf8 -o ${cleaned_gff} > ${cleaned_gff}.log
rm ${input_gff}.utf8
# launch with ./Sfix_gff.sh EXP_chr2A_LRRlocus_curated_20240209.gff DWSvevo1  curated_chr2A_LRRlocus_20240209_cleaned.gff
