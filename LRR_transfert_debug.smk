#!/usr/bin/env python

import pandas as pd
import os,sys
from itertools import compress

singularity:"library://cgottin/default/lrrtransfer:2.0"

####################   DEFINE CONFIG VARIABLES BASED ON CONFIG FILE   ####################

### Variables from config file

target_genome:     "example/Punctata_chr1.fasta"
ref_genome:        "example/Nipponbare_Chr1.fasta"
ref_gff:           "example/Nipponbare_LRR-CR_Chr1.gff"
ref_locus_info:    "example/Info_locus_Nipponbare.txt"
# Choosing the annotation (mostly for testing purpose): "best"
mode: "best"
# optional folder containing a prebuild LRROME
lrrome:            ""
OUTPUTS_DIRNAME: "LRR_TRANSFERT_OUTPUTS"

target_genome = config["target_genome"]
ref_genome = config["ref_genome"]
ref_gff = config["ref_gff"]
ref_locus_info = config["ref_locus_info"]
mode = config["mode"]
lrrome = config["lrrome"]
outDir = config["OUTPUTS_DIRNAME"]

if (len(lrrome) == 0):
    lrrome = "NULL"

### Define paths
path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
#scripts_dir = snakefile_dir+"/SCRIPTS"
### Define outputs subfolders
outLRRomeDir = outDir+"/LRRome"



### PIPELINE ###
#rule order: genePred_one > genePrediction

rule FinalTargets:
    input:
        #outLRRomeDir+"/REF_proteins.fasta"
        #outDir+"/filtered_candidatsLRR.gff"
        #outDir+"/annot.gff"
        #outDir+"/annotate_one_00000.gff"
        outDir+"/LRRlocus_predicted.gff"
 # ----------------------------------------------------------------------------------------------- #


rule buildLRROme:
    input:
        ref_genome=ref_genome,
        ref_gff=ref_gff
    output:
        #outLRRomeDir+"/REF_cDNA.fasta",
        #outLRRomeDir+"/REF_exons.fasta",
        #outLRRomeDir+"/REF_proteins.fasta"
        directory(outLRRomeDir)
#    singularity:
#        "library://cgottin/default/lrrtransfer:2.0"
    shell:
        "${{LRR_BIN}}/create_LRRome.sh {input.ref_genome} {input.ref_gff} {outDir} {lrrome}"

rule candidateLoci:
    input:
    	target_genome,
    	outLRRomeDir,
    	ref_gff
    output:
    	outDir+"/list_query_target.txt",
    	outDir+"/filtered_candidatsLRR.gff",
        directory(outDir+"/CANDIDATE_SEQ_DNA")
 #   singularity:
 #       "library://cgottin/default/lrrtransfer:2.0"
    shell:
    	"${{LRR_BIN}}/candidateLoci.sh {input} {outDir}"
        #"/work_home/cgottin/SNAKEMAKE/candidateLoci_test.sh {input} {outDir}"

# rule genePred_one:
#     input:
#     	outDir+"/list_query_target_split.00000",
#     	outDir+"/CANDIDATE_SEQ_DNA",
#     	target_genome,
#     	outDir+"/filtered_candidatsLRR.gff",
#     	outLRRomeDir,
#     	ref_gff,
#     	ref_locus_info,
#     params:
#         outDir=outDir,
#         mode=mode
#     output:
#     	outDir+"/annotate_one_00000.gff"
#     shell:
#         #"${{LRR_BIN}}/genePrediction.sh {input} {output} {params}"
#         "bash -xv /work_home/cgottin/SNAKEMAKE/genePrediction_test.sh {input} {params.outDir} {output} {params.mode}"



rule split_candidates:
	input:
		outDir+"/list_query_target.txt"
	output:
		dynamic(outDir+"/list_query_target_split.{split_id}")
	shell:
		"cd {outDir}; head {input} > {input}.head ; split -a 5 -d -l 1 {input}.head list_query_target_split."

rule genePrediction:
    input:
    	outDir+"/list_query_target_split.{split_id}",
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
    	outDir+"/annotate_one_{split_id}.gff"
    shell:
        "${{LRR_BIN}}/genePrediction.sh {input} {params.outDir} {output} {params.mode}"
        #"bash -xv /work_home/cgottin/SNAKEMAKE/genePrediction_test.sh {input} {params.outDir} {output} {params.mode}"

rule merge_prediction:
	input:
		dynamic(outDir+"/annotate_one_{split_id}.gff")
	output:
		outDir+"/annot.gff"
	shell:
		"cat {input}>>{output}"

rule verif_annotation:
    input:
        ref_locus_info,
        target_genome,
        outDir,
        outDir+"/annot.gff"
    output:
        outDir+"/LRRlocus_predicted.gff",
        outDir+"/alert_NC_Locus.txt"
    shell:
        #"bash -xv /work_home/cgottin/SNAKEMAKE/verifAnnot_test.sh {input}"
        "${{LRR_BIN}}/verifAnnot.sh {input}"
