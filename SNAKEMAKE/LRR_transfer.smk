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
        #outDir+"/refProts"
        outDir+"/annot_best.gff"
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
        target_genome+".nsq"
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
        #"mkdir {outDir}/refProts; cd {outDir}/refProts; split -a 5 -d -l 4 {input} REF_proteins_split."
        "mkdir {outDir}/refProts; cd {outDir}/refProts; split -a 5 -d -l 20 {input} REF_proteins_split."

def aggregate_blast(wildcards):
    checkpoint_output = checkpoints.split_blast.get(**wildcards).output[0]
    return expand(outDir+"/refProts/blast_split_{id}_res.tsv",
        id=glob_wildcards(os.path.join(checkpoint_output, "REF_proteins_split.{id}")).id)


rule blastProt:
    input:
        ref_prots=outDir+"/refProts/REF_proteins_split.{id}",
        target_genome=target_genome,
        blast_db=target_genome+".nsq",
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
        "cp /lustre/ranwezv/RUN_LRROME/LRR_TRANSFERT_OUTPUTS_BUG/refProts/{params.resFile} {output}"
        #"tblastn -db {input.target_genome} -query {input.ref_prots} -evalue 1 -out {output} -outfmt '6 qseqid sseqid qlen length qstart qend sstart send nident pident gapopen evalue bitscore' "
        #"tblastn -num_threads 4 -db {input.target_genome} -query {input.ref_prots} -evalue 1 -out {output} -outfmt '6 qseqid sseqid qlen length qstart qend sstart send nident pident gapopen evalue bitscore' "
 
rule merge_blast:
    input:
        aggregate_blast
    params:
        outDir=outDir,
    output:
        outDir+"/blast_refProt.tsv"
    shell:
        #"find . -name {outDir}/refProts/blast_split_*_res.tsv -print0 | xargs -0 cat > {output}"
        "cat {outDir}/refProts/blast_split_*_res.tsv > {output}"

rule candidateLoci:
    input:
        target_genome,
        outLRRomeDir,
        ref_gff,
        outDir+"/blast_refProt.tsv"
    output:
        outDir+"/list_query_target.txt",
        outDir+"/filtered_candidatsLRR.gff",
        directory(outDir+"/CANDIDATE_SEQ_DNA")
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
        prot=outDir+"/annotate_one_{split_id}_prot2genome.gff",
        cdnaExon=outDir+"/annotate_one_{split_id}_cdna2genomeExon.gff",
        protExon=outDir+"/annotate_one_{split_id}_prot2genomeExon.gff"
    conda:
        "./conda_tools.yml"
    shell:
        "{LRR_BIN}/genePrediction.sh {input} {params.outDir} {outDir}/annotate_one_{wildcards.split_id} {params.mode} {LRR_SCRIPT}"
 # ------------------------------------------------------------------------------------ #

rule merge_prediction:
    input:
        aggregate_best
    output:
        best=outDir+"/annot_best.gff",
        mapping=outDir+"/annot_mapping.gff",
        cdna=outDir+"/annot_cdna2genome.gff",
        prot=outDir+"/annot_prot2genome.gff",
        cdnaExon=outDir+"/annot_cdna2genomeExon.gff",
        protExon=outDir+"/annot_prot2genomeExon.gff"
    conda:
        "./conda_tools.yml"
    shell:
        """
        #cat {input} > {output.best}_tmp;
        
        for method in best mapping cdna2genome cdna2genomeExon prot2genome prot2genomeExon; do
            cat {outDir}/annotate_one_*_${{method}}.gff > {outDir}/${{method}}_tmp;
            {LRR_SCRIPT}/VR/Sfix_gff.sh {outDir}/${{method}}_tmp DWSvevo3 {outDir}/annot_${{method}}.gff;
            rm  {outDir}/${{method}}_tmp;
        done

        #"rm {outDir}/annotate_one_*.gff;"
        """

 # ------------------------------------------------------------------------------------ #
