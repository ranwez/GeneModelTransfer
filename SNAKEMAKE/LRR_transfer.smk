#!/usr/bin/env python

####################               SINGULARITY CONTAINER              ####################

#singularity:"library://cgottin/default/lrrtransfer:2.1"
#singularity:"../lrrtransfer_2.1.sif"
####################   DEFINE CONFIG VARIABLES BASED ON CONFIG FILE   ####################

localrules: transfer_stats

target_genome = config["target_genome"]
ref_genome = config["ref_genome"]
ref_gff = config["ref_gff"]
ref_locus_info = config["ref_locus_info"]
mode = "best2rounds"
prefix = "DWSvevo3"
gP_methods = ["best", "best1", "mapping", "cdna2genome", "cdna2genomeExon", "cds2genome", "cds2genomeExon", "prot2genome", "prot2genomeExon"]
preBuildLRRomeDir = config["lrrome"]
outDir = config["OUTPUTS_DIRNAME"]

if (len(preBuildLRRomeDir) == 0):
    preBuildLRRomeDir='NULL'
outLRRomeDir = outDir+"/LRRome"

### Define paths
path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
LRR_SCRIPT = snakefile_dir+"/../SCRIPT"
LRR_BIN = snakefile_dir+"/../bin"
working_directory = os.getcwd()


####################                  RUNNING PIPELINE                ####################


rule All:
    input:
        outDir+"/stats/GFFstats.txt",
        outDir+"/stats/jobsStats_out.txt",
        outDir+"/stats/jobsStats_err.txt"
        #outDir+"/refProts"
        #outDir+"/annot_best.gff"
        #outDir+"/LRRlocus_predicted_best.gff",
        #outDir+"/LRRlocus_predicted_mapping.gff",
        #outDir+"/LRRlocus_predicted_cdna2genome.gff",
        #outDir+"/LRRlocus_predicted_prot2genome.gff"


 # ------------------------------------------------------------------------------------ #

rule checkFiles:
    input:
        target_genome,
        ref_genome,
        ref_gff,
        ref_locus_info
    output:
        outDir+"/input_summary.log"
    conda:
        "./conda_tools.yml"
    shell:
        "{LRR_BIN}/check_files.sh {input} {outLRRomeDir} {outDir} {output};"

 # ------------------------------------------------------------------------------------ #

rule buildLRROme:
    input:
        ref_genome=ref_genome,
        ref_gff=ref_gff,
        log_file=outDir+"/input_summary.log"
    output:
        directory(outLRRomeDir),
        outLRRomeDir+"/REF_proteins.fasta"
    conda:
        "./conda_tools.yml"
    shell:
        "{LRR_BIN}/create_LRRome.sh {input.ref_genome} {input.ref_gff} {outDir} {preBuildLRRomeDir} {LRR_SCRIPT}"

# ------------------------------------------------------------------------------------ #
# TODO do not lsmix output from makeblastdb with input genome fasta file
rule makeBlastdb:
    input:
        target_genome=target_genome,
        ref_prots=outDir+"/refProts"
    params:
        outDir=outDir,
    output:
        target_genome+".nal"
    conda:
        "./conda_tools.yml"
    shell:
        "cd {input.ref_prots};"
        "echo makeblastdb -in {input.target_genome} -dbtype nucl -out {input.target_genome};"
        "makeblastdb -in {input.target_genome} -dbtype nucl -out {input.target_genome};"
        "cd -"


checkpoint split_blast:
    input:
        outLRRomeDir+"/REF_proteins.fasta"
    output:
        refProts=directory(outDir+"/refProts")
    conda:
        "./conda_tools.yml"
    shell:
        "mkdir {outDir}/refProts; cd {outDir}/refProts; split -a 5 -d -l 20 {input} REF_proteins_split."

def aggregate_blast(wildcards):
    checkpoint_output = checkpoints.split_blast.get(**wildcards).output[0]
    return expand(outDir+"/refProts/blast_split_{id}_res.tsv",
        id=glob_wildcards(os.path.join(checkpoint_output, "REF_proteins_split.{id}")).id)


rule blastProt:
    input:
        ref_prots=outDir+"/refProts/REF_proteins_split.{id}",
        target_genome=target_genome,
        blast_db=target_genome+".nal",
    params:
        outDir=outDir,
        resFile="blast_split_{id}_res.tsv"
    output:
        outDir+"/refProts/blast_split_{id}_res.tsv"
    conda:
        "./conda_tools.yml"
    #threads:
    #    6
    shell:
        ### WARNING TRICK TO NOT RECOMPUTE BLAST
        #"cp /lustre/ranwezv/RUN_LRROME/LRR_TRANSFERT_OUTPUTS_BUG/refProts/{params.resFile} {output}"
        "tblastn -db {input.target_genome} -query {input.ref_prots} -evalue 1 -out {output} -outfmt '6 qseqid sseqid qlen length qstart qend sstart send nident pident gapopen evalue bitscore positive' "
        #"touch {output}"

rule merge_blast:
    input:
        aggregate_blast
    params:
        outDir=outDir,
    output:
        outDir+"/blast_refProt.tsv"
    shell:
        "cat {outDir}/refProts/blast_split_*_res.tsv > {output}"
        #"cp {outDir}/blast_refProt_save.tsv {output}"

rule candidateLoci:
    input:
        target_genome=target_genome,
        ref_gff=ref_gff,
        blast_res=outDir+"/blast_refProt.tsv"
    output:
        outDir+"/list_query_target.txt",
        outDir+"/filtered_candidatsLRR.gff",
        directory(outDir+"/CANDIDATE_SEQ_DNA")
    conda:
        "./conda_tools.yml"
    envmodules:
        "bedtools/2.30.0"
    params:
        min_sim = config["CL_min_sim"]
    shell:
        """
        ## amelio : split par chromosome de target_genome et parallélisation
        #"{LRR_BIN}/candidateLoci.sh {input} {outDir} {LRR_SCRIPT}"
        python {LRR_SCRIPT}/candidate_loci_VR.py -g {input.ref_gff} -t {input.blast_res} -o {outDir}/filtered_candidatsLRR.gff -l {outDir}/list_query_target.txt -s {params.min_sim}
        {LRR_SCRIPT}/CANDIDATE_LOCI/extract_loci.sh {outDir}/filtered_candidatsLRR.gff {input.target_genome} {outDir}/CANDIDATE_SEQ_DNA
        """


 # ------------------------------------------------------------------------------------ #

checkpoint split_candidates:
    input:
        outDir+"/list_query_target.txt"
    output:
        queryTargets=directory(outDir+"/queryTargets")
    conda:
        "./conda_tools.yml"
    shell:
        "mkdir {outDir}/queryTargets; cd {outDir}/queryTargets; split -a 5 -d -l 1 {input} list_query_target_split."
 # ------------------------------------------------------------------------------------ #

def aggregate_best(wildcards):
    checkpoint_output = checkpoints.split_candidates.get(**wildcards).output[0]
    return expand(outDir+"/annotate_one/annotate_one_{id}_best.gff",
           id=glob_wildcards(os.path.join(checkpoint_output, "list_query_target_split.{id}")).id)


rule sortGFF:
    input:
        ref_gff
    output:
        outDir+"/ref_sorted.gff"
    conda:
        "./conda_tools.yml"
    shell:
        "python3 {LRR_SCRIPT}/sort_gff.py -g {input} -o {output}"


rule genePrediction:
    input:
        outDir+"/queryTargets/list_query_target_split.{split_id}",
        outDir+"/CANDIDATE_SEQ_DNA",
        target_genome,
        outDir+"/filtered_candidatsLRR.gff",
        outLRRomeDir,
        #ref_gff,
        outDir+"/ref_sorted.gff",
        ref_locus_info,
    params:
        outDir=outDir,
        mode=mode
    output:
        best=outDir+"/annotate_one/annotate_one_{split_id}_best.gff",
        best1=outDir+"/annotate_one/annotate_one_{split_id}_best1.gff",
        mapping=outDir+"/annotate_one/annotate_one_{split_id}_mapping.gff",
        cdna=outDir+"/annotate_one/annotate_one_{split_id}_cdna2genome.gff",
        cds=outDir+"/annotate_one/annotate_one_{split_id}_cds2genome.gff",
        prot=outDir+"/annotate_one/annotate_one_{split_id}_prot2genome.gff",
        cdnaExon=outDir+"/annotate_one/annotate_one_{split_id}_cdna2genomeExon.gff",
        cdsExon=outDir+"/annotate_one/annotate_one_{split_id}_cds2genomeExon.gff",
        protExon=outDir+"/annotate_one/annotate_one_{split_id}_prot2genomeExon.gff"
    conda:
        "./conda_tools.yml"
    shell:
        "{LRR_BIN}/genePrediction.sh {input} {params.outDir} {outDir}/annotate_one/annotate_one_{wildcards.split_id} {params.mode} {LRR_SCRIPT}"

 # ------------------------------------------------------------------------------------ #

rule merge_prediction:
    input:
        aggregate_best
    output:
        outDir+"/annot_best.gff",
        outDir+"/annot_mapping.gff",
        outDir+"/annot_cdna2genome.gff",
        outDir+"/annot_cds2genome.gff",
        outDir+"/annot_prot2genome.gff",
        outDir+"/annot_cdna2genomeExon.gff",
        outDir+"/annot_cds2genomeExon.gff",
        outDir+"/annot_prot2genomeExon.gff",
        #outDir+"/stats/GFFstats.txt",
        #outDir+"/stats/jobsStats_out.txt",
        #outDir+"/stats/jobsStats_err.txt"
    conda:
        "./conda_tools.yml"
    shell:
        """
        for method in {gP_methods}; do
            {LRR_BIN}/merge_prediction.sh {outDir}/annotate_one {LRR_SCRIPT} ${{method}} {prefix} {outDir}
        done

        #{LRR_SCRIPT}/STATS_OUTPUTS/stats_transfer.sh {outDir} {outDir}/.. {outDir}/stats
        """

rule transfer_stats:
    input:
        outDir+"/annot_best.gff",
        outDir+"/annot_mapping.gff",
        outDir+"/annot_cdna2genome.gff",
        outDir+"/annot_cds2genome.gff",
        outDir+"/annot_prot2genome.gff",
        outDir+"/annot_cdna2genomeExon.gff",
        outDir+"/annot_cds2genomeExon.gff",
        outDir+"/annot_prot2genomeExon.gff"
    output:
        outDir+"/stats/GFFstats.txt",
        outDir+"/stats/jobsStats_out.txt",
        outDir+"/stats/jobsStats_err.txt"
    conda:
        "./conda_tools.yml"
    shell:
        """
        {LRR_SCRIPT}/STATS_OUTPUTS/stats_transfer.sh {outDir} {outDir}/.. {outDir}/stats
        """

 # ------------------------------------------------------------------------------------ #
