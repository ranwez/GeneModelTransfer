#!/usr/bin/env python

####################               SINGULARITY CONTAINER              ####################

#singularity:"library://cgottin/default/lrrtransfer:2.1"
#singularity:"../lrrtransfer_2.1.sif"
####################   DEFINE CONFIG VARIABLES BASED ON CONFIG FILE   ####################

target_genome = config["target_genome"]
ref_genome = config["ref_genome"]
ref_gff = config["ref_gff"]
ref_locus_info = config["ref_locus_info"]
mode = "best"
lrrome = config["lrrome"]
outDir = config["OUTPUTS_DIRNAME"]

if (len(lrrome) == 0):
    lrrome = "NULL"

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
        outDir+"/LRRlocus_predicted_best.gff",
        outDir+"/LRRlocus_predicted_mapping.gff",
        outDir+"/LRRlocus_predicted_cdna2genome.gff",
        outDir+"/LRRlocus_predicted_prot2genome.gff"


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
        "{LRR_BIN}/check_files.sh {input} {lrrome} {outDir} {output};"

 # ------------------------------------------------------------------------------------ #

rule buildLRROme:
    input:
        ref_genome=ref_genome,
        ref_gff=ref_gff,
        log_file=outDir+"/input_summary.log"
    output:
        directory(outLRRomeDir)
    conda:
        "./conda_tools.yml"
    shell:
        "{LRR_BIN}/create_LRRome.sh {input.ref_genome} {input.ref_gff} {outDir} {lrrome} {LRR_SCRIPT}"
# ------------------------------------------------------------------------------------ #

rule candidateLoci:
    input:
        target_genome,
        outLRRomeDir,
        ref_gff
    output:
        outDir+"/list_query_target.txt",
        temp(outDir+"/filtered_candidatsLRR.gff"),
        temp(directory(outDir+"/CANDIDATE_SEQ_DNA"))
    conda:
        "./conda_tools.yml"
    shell:
        ## amelio : split par chromosome de target_genome et parallélisation
        "{LRR_BIN}/candidateLoci.sh {input} {outDir} {LRR_SCRIPT}"

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
    return expand(outDir+"/annotate_one_{id}_best.gff",
           id=glob_wildcards(os.path.join(checkpoint_output, "list_query_target_split.{id}")).id)

rule genePrediction:
    input:
        outDir+"/queryTargets/list_query_target_split.{split_id}",
        outDir+"/CANDIDATE_SEQ_DNA",
        target_genome,
        outDir+"/filtered_candidatsLRR.gff",
        outLRRomeDir,
        ref_gff,
        ref_locus_info,
    params:
        outDir=outDir,
        mode=mode
    output:
        best=outDir+"/annotate_one_{split_id}_best.gff",
        mapping=outDir+"/annotate_one_{split_id}_mapping.gff",
        cdna=outDir+"/annotate_one_{split_id}_cdna2genome.gff",
        prot=outDir+"/annotate_one_{split_id}_prot2genome.gff"
    conda:
        "./conda_tools.yml"
    shell:
        "{LRR_BIN}/genePrediction.sh {input} {params.outDir} {outDir}/annotate_one_{wildcards.split_id} {params.mode} {LRR_SCRIPT}"
 # ------------------------------------------------------------------------------------ #

rule merge_prediction:
    input:
        aggregate_best
    output:
        best=temp(outDir+"/annot_best.gff"),
        mapping=temp(outDir+"/annot_mapping.gff"),
        cdna=temp(outDir+"/annot_cdna2genome.gff"),
        prot=temp(outDir+"/annot_prot2genome.gff")
    conda:
        "./conda_tools.yml"
    shell:
        "cat {input} > {output.best};"
        "cat {outDir}/annotate_one_*_mapping.gff > {output.mapping};"
        "cat {outDir}/annotate_one_*_cdna2genome.gff > {output.cdna};"
        "cat {outDir}/annotate_one_*_prot2genome.gff > {output.prot};"
        #"rm {outDir}/annotate_one_*.gff;"

 # ------------------------------------------------------------------------------------ #

rule verif_annotation:
    input:
        ref_locus_info,
        target_genome,
        outDir+"/annot_{method}.gff",

    output:
        outDir+"/LRRlocus_predicted_{method}.gff",
        outDir+"/alert_NC_Locus_{method}.txt"
    params:
        outDir
    conda:
        "./conda_tools.yml"
    shell:
        "{LRR_BIN}/verifAnnot.sh {input} {params} {LRR_SCRIPT} {wildcards.method}"
#modif wildcard