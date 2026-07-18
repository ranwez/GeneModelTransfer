# task make a plan to add a new annotation transfert to the pipeline

## Objective
No coding, just a plan to add a new annotation transfer to the pipeline.

## current state of the pipeline
Lrrtransfer is a complex gene family annotation transfer pipeline designed for plant LRR (Leucine-Rich Repeat) containing receptors. 

- The pipeline use a reference annotation of LRR gene families to annotate LRR of a closely related genome.  It procedes in several steps @SNAKEMAKE/LRR_transfer.smk (blast to identify candidate loci, for each pair candidate locus-modelProtein, diffeent annotation transfert strategies are used based on exonerate or blast; improved and compare to get the best possible annotation for this locus).
- the core of the pipeline is the annotation of a locus based on the associated protein (annotation transfer) and is done in @bin/genePrediction.sh.
- the current version already extract the genomic region of the model and target locus using Extract_sequences_from_genome.py  and handle feature coordinates shifting and strand orientation when moving from the full region to the extracted region handled by the script (functions gff_genome_to_target and gff_target_to_genome in @bin/genePrediction.sh).

Current pipeline miss short exon and introduce errors regions that are highly similar at nucleotide level. The idea is to add a new annotation transfer strategy that will be able to handle these cases and improve the overall quality of the annotation.

## the new annotation transfer strategy; locus-alignment based annotation transfer
1) get the locus corresponding to the mode protein (from the reference species) and the target locus (from the target species).
2) blast them at the nucleotide level, to detect if the strategy is relevant (i.e. if the two loci are highly similar at nucleotide level). Criteria: the best blast hit should cover at least 85% of the model locus length and have at least 95% identity. (we allow small trucations, but the aligned region should be highly similar).
3) if the blast criteria are met, then we will use a new basic strategy to transfer the annotation from the model locus to the target locus. 
a) align the two loci using MAFFT.
b) transfer the annotation from the model locus (single mRNA only non overlapping CDS) to the target locus based on the alignment, order the CDS features then traverse the alignment at each position keep positions in the model/target locus (i.e. unchanged if there is a gap; +1 otherwise), if current coordinate in the model locus is a CDS start/end feature, then add the corresponding coordinate in the target locus as a start/end CDS feature.
c) add this has a new annotatition transfer method in @bin/genePrediction.sh "methods="mapping cdna2genome cdna2genomeExon cds2genome cds2genomeExon prot2genome prot2genomeExon"; so that the rest is automatically handled by the pipeline (improving the splicing site positions, searching missed start etc; and then comparison with other alignment methods).

## pipeline integration 

1) the new transfert should first be implemented as separate scripts python/bash/awk and the sole modification in the current pipeline file will be to add a new method in the list of methods in @bin/genePrediction.sh. 
2) LLRome currently do not heve the locus/gff per model protein. A folder could be added to the LRRome folder with these extra information (locus/gff per model protein) to be used by the new annotation transfer strategy.
3) the current locus extraction for pre-exploration of this feature has been done from full genome using
```bash 
SCRIPT/CANDIDATE_LOCI/extract_loci.sh OsjNip_HIPPHPP_202606v1.gff Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fasta REF_LOCI
SCRIPT/CANDIDATE_LOCI/extract_loci.sh OsjNip_HIPPHPP_202606v1.gff Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fasta REF_LOCI 
```
This locus extraction is much faster than the current python script Extract_sequences_from_genome.py and once locus of reference genes will be available in the LRRome folder, the genePrediction.sh script will be able to  utilize these pre-extracted loci (reducing memory and CPU usage, this is for a next pipeline improvement iteration).
4) The new annotation could first be tested alone using a set of model/target loci pairs to validate the method and then integrated in the pipeline. These tested should relies on data provides in @HIPP/HIPP_NIP_KIT and include:
- pairs of loci from same species (to validate the method on identical loci), the raw predicted annotation should be identical to the model locus annotation eg. OsjNip_Chr04_02891360 and OsjNip_Chr04_02891360.
- pairs of loci similar, on same and opposite strand, to validate the method on loci that are similar but not identical, e.g. annotating OSJkit_chr01_0023883423 using OsjNip_Chr01_23324465 (same strand) and annotating OSJkit_chr04_0020834228 using OsjNip_Chr04_19814005 (different strand, detail in @HIPP/HIPP_NIP_KIT/blast_refProt.tsv | grep OsjNip_Chr04_19814005 )
- pairs of loci that diverge where the method should fail and the blast criteria should not be met, e.g. OSJkit_chr04_0003866905 and OsjNip_Chr04_02891360

5) in a second step the Extract_sequences_from_genome.py called from @bin/genePrediction.sh should be modified to use the locus.fasta files now available in the LRRome folder. This will be done latter to avoid breaking the current pipeline and to allow testing of the new annotation transfer strategy independently, but any relevant discover during this stage should be documented in the nex task specification file TASKS/planning_speedUpLocus_extraction_byLRRome.md.