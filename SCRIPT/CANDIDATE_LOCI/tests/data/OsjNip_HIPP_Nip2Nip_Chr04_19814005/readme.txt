When using softmasking tblastn
```
python ../../SCRIPT/candidate_loci_VR.py -g OsjNip_HIPPHPP_Chr04_19814005/ref.gff -t OsjNip_HIPPHPP_Chr04_19814005/tblastn_softmasking.tsv -o _test_chr4.gff -l _test_list_chr4.txt
````
Part of the gene locus is missing
```
../helper/check_gene_coverage.sh   --ref_gff OsjNip_HIPPHPP_Chr04_19814005/ref.gff   --candidate_gff _test_chr4.gff   --out _all_chr_coverage.tsv
more _all_chr_coverage.tsv
chr     start   end     gene_id gene_length     covered_bases   coverage_percent        fully_covered
Chr4    19796632        19797845        OsjNip_Chr04_19797845   1214    1214    100.00  YES
Chr4    19812328        19814005        OsjNip_Chr04_19814005   1678    918     54.71   NO
Chr4    19819617        19820436        OsjNip_Chr04_19820436   820     820     100.00  YES
Chr4    19822056        19822949        OsjNip_Chr04_19822949   894     894     100.00  YES
Chr4    19826726        19827860        OsjNip_Chr04_19827860   1135    1135    100.00  YES
Chr4    19835501        19836437        OsjNip_Chr04_19836437   937     937     100.00  YES
Chr4    19841797        19842670        OsjNip_Chr04_19842670   874     874     100.00  YES
```
