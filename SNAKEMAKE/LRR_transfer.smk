#!/usr/bin/env python

####################               SINGULARITY CONTAINER              ####################

singularity:"library://cgottin/default/lrrtransfer:2.1"

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


####################                  RUNNING PIPELINE                ####################


rule All:
    input:
        outDir+"/LRRlocus_predicted.gff"


 # ------------------------------------------------------------------------------------ #

rule checkFiles:
    input:
        target_genome,
        ref_genome,
        ref_gff,
        ref_locus_info
    output:
        outDir+"/input_summary.log"
    shell:
        "${{LRR_BIN}}/check_files.sh {input} {lrrome} {outDir} {output};"

 # ------------------------------------------------------------------------------------ #

rule buildLRROme:
    input:
        ref_genome=ref_genome,
        ref_gff=ref_gff,
        log_file=outDir+"/input_summary.log"
    output:
        directory(outLRRomeDir)
    shell:
        "${{LRR_BIN}}/create_LRRome.sh {input.ref_genome} {input.ref_gff} {outDir} {lrrome}"

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
    shell:
    	"${{LRR_BIN}}/candidateLoci.sh {input} {outDir}"

 # ------------------------------------------------------------------------------------ #

rule split_candidates:
	input:
		outDir+"/list_query_target.txt"
	output:
		dynamic(temp(outDir+"/list_query_target_split.{split_id}"))
	shell:
		"cd {outDir}; split -a 5 -d -l 1 {input} list_query_target_split."

 # ------------------------------------------------------------------------------------ #

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
    	best=temp(outDir+"/annotate_one_{split_id}_best.gff"),
        mapping=temp(outDir+"/annotate_one_{split_id}_mapping.gff"),
        cdna=temp(outDir+"/annotate_one_{split_id}_cdna2genome.gff"),
        prot=temp(outDir+"/annotate_one_{split_id}_prot2genome.gff")
    shell:
        "${{LRR_BIN}}/genePrediction.sh {input} {params.outDir} {outDir}/annotate_one_{wildcards.split_id} {params.mode}"

 # ------------------------------------------------------------------------------------ #

rule merge_prediction:
	input:
		best=dynamic(outDir+"/annotate_one_{split_id}_best.gff"),
        mapping=dynamic(outDir+"/annotate_one_{split_id}_mapping.gff"),
        cdna=dynamic(outDir+"/annotate_one_{split_id}_cdna2genome.gff"),
        prot=dynamic(outDir+"/annotate_one_{split_id}_prot2genome.gff")
	output:
		best=temp(outDir+"/annot_best.gff"),
        mapping=temp(outDir+"/annot_mapping.gff"),
        cdna=temp(outDir+"/annot_cdna2genome.gff"),
        prot=temp(outDir+"/annot_prot2genome.gff")

	shell:
		"cat {input.best}>>{output.best};"
        "cat {input.mapping}>>{output.mapping};"
        "cat {input.cdna}>>{output.cdna};"
        "cat {input.prot}>>{output.prot};"

 # ------------------------------------------------------------------------------------ #

rule verif_annotation:
    input:
        ref_locus_info,
        target_genome,
        outDir+"/annot_{method}.gff"
    output:
        outDir+"/LRRlocus_predicted_{method}.gff",
        outDir+"/alert_NC_Locus_{method}.txt"
    params:
        outDir
    shell:
        "${{LRR_BIN}}/verifAnnot.sh {input} {params}"
