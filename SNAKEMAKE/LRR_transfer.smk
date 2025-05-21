#!/usr/bin/env python
import os 

####################               SINGULARITY CONTAINER              ####################

#singularity:"library://cgottin/default/lrrtransfer:2.1"
#singularity:"/storage/replicated/cirad/projects/GE2POP/2023_LRR/LRRtransfer_image/LRRtransfer.sif"
singularity_image = config["singularity_image"]
singularity: singularity_image

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
outDir = os.path.abspath(config["OUTPUTS_DIRNAME"])
ignore_exonerate_errors = str(config["IGNORE_EXONERATE_ERRORS"]).lower() 

if (len(preBuildLRRomeDir) == 0):
    preBuildLRRomeDir='NULL'
outLRRomeDir = outDir+"/LRRome"

## Functions
def remove_file_ext(file_name, nb_ext=1):
    parts = file_name.split(".")
    if len(parts) <= nb_ext:
        return parts[0]
    return ".".join(parts[:-nb_ext])

def get_config_file():
    if '--configfile' in sys.argv:
        i = sys.argv.index('--configfile')
    elif '--configfiles' in sys.argv:
        i = sys.argv.index('--configfiles')
    config_file = sys.argv[i+1]
    return config_file


### Define paths
path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
LRR_SCRIPT = snakefile_dir+"/../SCRIPT"
LRR_BIN = snakefile_dir+"/../bin"
working_directory = os.getcwd()
config_file = get_config_file()

# db directory
target_genome_file_name=os.path.basename(target_genome)
target_genome_basename=remove_file_ext(file_name=target_genome_file_name)
target_genome_dir=os.path.dirname(target_genome)
target_genome_db_dir=f"{target_genome_dir}/{target_genome_basename}_db"

####################                  RUNNING PIPELINE                ####################


rule All:
    input:
        outDir+"/stats/GFFstats.txt",
        outDir+"/stats/jobsStats_out.txt",
        outDir+"/stats/jobsStats_err.txt"



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
        "{LRR_BIN}/check_files.sh {input} {outLRRomeDir} {outDir} {output};"

# ------------------------------------------------------------------------------------ #

rule writeWorkflowLog:
    input:
        config_file
    output:
        temp(outDir+"/log.sentinel")
    params:
        log_file = outDir+"/LRRtransfer.log",
    shell:
        "{LRR_BIN}/write_workflow_log.sh {snakefile_dir} {input} {params.log_file} {output}"

 # ------------------------------------------------------------------------------------ #

rule buildLRROme:
    input:
        ref_genome=ref_genome,
        ref_gff=ref_gff,
        log_file=outDir+"/input_summary.log",
        log_sentinel=outDir+"/log.sentinel"
    output:
        directory(outLRRomeDir),
        outLRRomeDir+"/REF_proteins.fasta"
    shell:
        "{LRR_BIN}/create_LRRome.sh {input.ref_genome} {input.ref_gff} {outDir} {preBuildLRRomeDir} {LRR_SCRIPT}"

# ------------------------------------------------------------------------------------ #
# TODO do not lsmix output from makeblastdb with input genome fasta file
rule makeBlastdb:
    input:
        target_genome=target_genome,
        ref_prots=outDir+"/refProts"
    output:
        blast_db_dir=directory(target_genome_db_dir)
    shell:
        "makeblastdb -in {input.target_genome} -dbtype nucl -out {output.blast_db_dir}/{target_genome_basename};"


checkpoint split_blast:
    input:
        outLRRomeDir+"/REF_proteins.fasta"
    output:
        refProts=directory(outDir+"/refProts")
    shell:
        "mkdir {outDir}/refProts; cd {outDir}/refProts; split -a 5 -d -l 20 {input} REF_proteins_split."

def aggregate_blast(wildcards):
    checkpoint_output = checkpoints.split_blast.get(**wildcards).output[0]
    return expand(outDir+"/refProts/blast_split_{id}_res.tsv",
        id=glob_wildcards(os.path.join(checkpoint_output, "REF_proteins_split.{id}")).id)


rule blastProt:
    input:
        ref_prots=outDir+"/refProts/REF_proteins_split.{id}",
        blast_db_dir=rules.makeBlastdb.output.blast_db_dir
    params:
        outDir=outDir,
        resFile="blast_split_{id}_res.tsv"
    output:
        outDir+"/refProts/blast_split_{id}_res.tsv"
    shell:
        ### WARNING TRICK TO NOT RECOMPUTE BLAST
        #"cp /lustre/ranwezv/RUN_LRROME/LRR_TRANSFERT_OUTPUTS_BUG/refProts/{params.resFile} {output}"
        "tblastn -db {input.blast_db_dir}/{target_genome_basename} -query {input.ref_prots} -evalue 1 -out {output} -outfmt '6 qseqid sseqid qlen length qstart qend sstart send nident pident gapopen evalue bitscore positive' "
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
    shell:
        "python3 {LRR_SCRIPT}/sort_gff.py -g {input} -o {output}"


rule genePrediction:
    input:
        outDir+"/queryTargets/list_query_target_split.{split_id}",
        outDir+"/CANDIDATE_SEQ_DNA",
        target_genome,
        outDir+"/filtered_candidatsLRR.gff",
        outLRRomeDir,
        outDir+"/ref_sorted.gff",
        ref_locus_info,
    params:
        outDir=outDir,
        mode=mode,
        ignore_exonerate_errors=ignore_exonerate_errors
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
    shell:
        "{LRR_BIN}/genePrediction.sh {input} {params.outDir} {outDir}/annotate_one/annotate_one_{wildcards.split_id} {params.mode} {LRR_SCRIPT} {params.ignore_exonerate_errors}"
        
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
    shell:
        """
        {LRR_SCRIPT}/STATS_OUTPUTS/stats_transfer.sh {outDir} {outDir}/.. {outDir}/stats
        """

 # ------------------------------------------------------------------------------------ #
