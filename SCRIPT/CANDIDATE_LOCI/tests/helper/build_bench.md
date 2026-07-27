extract blast hits and gff features from a region to create a unit test 
```bash
awk -F'\t' '$2=="Chr4"{s=($7<$8?$7:$8); e=($7>$8?$7:$8); if(s<=20000000 && e>=19790000) print}' blast.tsv
awk -F'\t' '$1=="Chr4" && $4<=20000000 && $5>=19790000' annotation.gff
```