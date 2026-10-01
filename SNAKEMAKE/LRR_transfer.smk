#!/usr/bin/env python

import os
import sys
from typing import Any, Optional

from snakemake.utils import validate

snakefile_dir = os.path.dirname(os.path.abspath(workflow.snakefile))
repo_root = os.path.abspath(os.path.join(snakefile_dir, ".."))
lrr_script = os.path.join(repo_root, "SCRIPT")
lrr_bin = os.path.join(repo_root, "bin")

sys.path.insert(0, lrr_script)

from run_info import create_run_id, finish_run_record, start_run_record

####################                  CONFIGURATION                  ####################

validate(config, "schemas/config.schema.yaml")

singularity_image = os.path.abspath(config["singularity_image"])
singularity: singularity_image

def optional_abspath(path: Optional[str]) -> Optional[str]:
    return os.path.abspath(path) if path is not None else None


target_genome = os.path.abspath(config["target_genome"])
ref_genome = optional_abspath(config["ref_genome"])
ref_gff = os.path.abspath(config["ref_gff"])
ref_locus_info = optional_abspath(config["ref_locus_info"])

provided_lrrome = optional_abspath(config["lrrome"])
provided_tblastn = optional_abspath(config["tblastn_results"])
provided_blastn = optional_abspath(config["blastn_results"])

out_dir = os.path.abspath(config["out_dirname"])
prefix = config["out_feature_id_prefix"]

mode = "best2rounds"

prediction_methods = [
    "best",
    "best1",
    "mapping",
    "locusAlignment",
    "cdna2genome",
    "cdna2genomeExon",
    "cds2genome",
    "cds2genomeExon",
    "prot2genome",
    "prot2genomeExon",
]


####################                     PATHS                      ####################

# Fasta indexes
target_genome_fai = f"{target_genome}.fai"
ref_genome_fai = f"{ref_genome}.fai" if ref_genome else None

# LRRome
generated_lrrome = os.path.join(out_dir, "LRRome")
lrrome = provided_lrrome or generated_lrrome

ref_proteins = os.path.join(lrrome, "REF_proteins.fasta")
ref_loci_fasta = os.path.join(lrrome, "REF_loci.fasta")
ref_pep_dir = os.path.join(lrrome, "REF_PEP")
ref_exons_dir = os.path.join(lrrome, "REF_EXONS")
ref_cdna_dir = os.path.join(lrrome, "REF_cDNA")
ref_loci_dir = os.path.join(lrrome, "REF_LOCI")
ref_loci_gff_dir = os.path.join(lrrome, "REF_LOCI_GFF")

# BLAST results
generated_tblastn = os.path.join(out_dir, "tblastn_refProt.tsv")
generated_blastn = os.path.join(out_dir, "blastn_refProt.tsv")

tblastn_results = provided_tblastn or generated_tblastn
blastn_results = provided_blastn or generated_blastn

# Persistent target BLAST database
target_genome_filename = os.path.basename(target_genome)
target_genome_basename = os.path.splitext(target_genome_filename)[0]
target_genome_dir = os.path.dirname(target_genome)
target_blast_db_dir = os.path.join(
    target_genome_dir,
    f"{target_genome_basename}_db",
)

# BLAST intermediate directories
tblastn_dir = os.path.join(out_dir, "refProts", "tblastn")
blastn_dir = os.path.join(out_dir, "refProts", "blastn")

# Candidate loci
candidate_pairs = os.path.join(out_dir, "list_query_target.txt")
candidate_gff = os.path.join(out_dir, "filtered_candidatsLRR.gff")
candidate_sequences = os.path.join(out_dir, "CANDIDATE_SEQ_DNA")
candidate_chunks = os.path.join(out_dir, "queryTargets")

# Reference GFF
sorted_ref_gff = os.path.join(out_dir, "ref_sorted.gff")

# Predictions
predictions_dir = os.path.join(out_dir, "annotate_one")
prediction_prefix = os.path.join(predictions_dir,"annotate_one_{split_id}")
prediction_outputs = [
    f"{prediction_prefix}_{method}.gff"
    for method in prediction_methods
]

# Statistics
stats_dir = os.path.join(out_dir, "stats")

# Workflow log
workflow_log = os.path.join(out_dir, "LRRtransfer.log")
run_history_dir = os.path.join(out_dir, "run_history")
latest_run_info = os.path.join(out_dir, "run_info.yaml")

run_id = create_run_id()
run_history_file = os.path.join(
    run_history_dir,
    f"{run_id}.yaml",
)

# Validation
input_check = os.path.join(out_dir, "input_check.log")
validation_dir = os.path.join(out_dir, "validation")
ref_gff_validation = os.path.join(validation_dir, "ref_gff.tsv")
lrrome_validation = os.path.join(validation_dir,"lrrome.tsv")
blast_results_validation = os.path.join(validation_dir,"blast_results.tsv")
inputs_validation = os.path.join(validation_dir,"inputs_validated.tsv")


####################                 LOCAL RULES                     ####################

localrules: transfer_stats


####################                   HELPERS                       ####################

def get_split_ids(checkpoint_output: str, pattern: str) -> list[str]:
    """Return sorted split IDs from a checkpoint output directory."""
    return sorted(
        glob_wildcards(
            os.path.join(checkpoint_output, pattern)
        ).id
    )


def aggregate_tblastn(wildcards: Any) -> list[str]:
    """Return all TBLASTN result chunks produced from the checkpoint."""
    checkpoint_output = checkpoints.split_tblastn.get(**wildcards).output.chunks
    split_ids = get_split_ids(
        checkpoint_output,
        "REF_proteins_split.{id}",
    )

    return expand(
        os.path.join(tblastn_dir, "blast_split_{id}_res.tsv"),
        id=split_ids,
    )


def aggregate_blastn(wildcards: Any) -> list[str]:
    """Return all BLASTN result chunks produced from the checkpoint."""
    checkpoint_output = checkpoints.split_blastn.get(**wildcards).output.chunks
    split_ids = get_split_ids(
        checkpoint_output,
        "REF_loci_split.{id}",
    )

    return expand(
        os.path.join(blastn_dir, "blast_split_{id}_res.tsv"),
        id=split_ids,
    )


def aggregate_predictions(wildcards: Any) -> list[str]:
    """Return every prediction output produced for candidate pairs."""
    checkpoint_output = checkpoints.split_candidates.get(**wildcards).output.chunks
    split_ids = get_split_ids(
        checkpoint_output,
        "list_query_target_split.{id}",
    )

    return expand(
        prediction_outputs,
        split_id=split_ids,
    )


####################                    WORKFLOW                     ####################

rule all:
    input:
        os.path.join(stats_dir, "GFFstats.txt"),
        os.path.join(stats_dir, "jobsStats_out.txt"),
        os.path.join(stats_dir, "jobsStats_err.txt")


# -------------------------------------------------------------------------------------- #

rule check_files:
    input:
        target_genome=target_genome,
        ref_gff=ref_gff,
        ref_genome=[ref_genome] if ref_genome else [],
        ref_locus_info=[ref_locus_info] if ref_locus_info else [],
        lrrome=[provided_lrrome] if provided_lrrome else [],
        tblastn=[provided_tblastn] if provided_tblastn else [],
        blastn=[provided_blastn] if provided_blastn else []
    output:
        input_check
    shell:
        """
        "{lrr_bin}/check_files.sh" \
            --target-genome "{input.target_genome}" \
            --ref-gff "{input.ref_gff}" \
            --ref-genome "{input.ref_genome}" \
            --ref-locus-info "{input.ref_locus_info}" \
            --lrrome "{input.lrrome}" \
            --tblastn-results "{input.tblastn}" \
            --blastn-results "{input.blastn}" \
            --output "{output}"
        """


# -------------------------------------------------------------------------------------- #


rule index_target_fasta:
    input:
        fasta=target_genome,
        checked=input_check
    output:
        fai=target_genome_fai
    shell:
        """
        samtools faidx "{input.fasta}"
        """

if ref_genome is not None and ref_genome != target_genome:

    rule index_ref_fasta:
        input:
            fasta=ref_genome,
            checked=input_check
        output:
            fai=ref_genome_fai
        shell:
            """
            samtools faidx "{input.fasta}"
            """


# -------------------------------------------------------------------------------------- #

rule validate_ref_gff:
    input:
        checked=input_check,
        gff=ref_gff,
        ref_fai=[ref_genome_fai] if ref_genome_fai else [],
        ref_locus_info=[ref_locus_info] if ref_locus_info else []
    output:
        report=ref_gff_validation
    shell:
        """
        python3 "{lrr_script}/INPUT_VALIDATION/validate_ref_gff.py" \
            --gff "{input.gff}" \
            --ref-fai "{input.ref_fai}" \
            --ref-locus-info "{input.ref_locus_info}" \
            --output "{output.report}"
        """

if provided_lrrome is not None:

    rule validate_lrrome:
        input:
            lrrome=provided_lrrome,
            gff=ref_gff,
            gff_validation=ref_gff_validation
        output:
            report=lrrome_validation
        shell:
            """
            python3 "{lrr_script}/INPUT_VALIDATION/validate_lrrome.py" \
                --lrrome "{input.lrrome}" \
                --gff "{input.gff}" \
                --output "{output.report}"
            """

if provided_tblastn is not None or provided_blastn is not None:

    rule validate_blast_results:
        input:
            gff=ref_gff,
            gff_validation=ref_gff_validation,
            target_fai=target_genome_fai,
            tblastn=[provided_tblastn] if provided_tblastn else [],
            blastn=[provided_blastn] if provided_blastn else []
        output:
            report=blast_results_validation
        shell:
            """
            python3 "{lrr_script}/INPUT_VALIDATION/validate_blast_results.py" \
                --gff "{input.gff}" \
                --target-fai "{input.target_fai}" \
                --tblastn "{input.tblastn}" \
                --blastn "{input.blastn}" \
                --output "{output.report}"
            """


rule validate_inputs:
    input:
        (
            [ref_gff_validation]
            + ([lrrome_validation] if provided_lrrome else [])
            + (
                [blast_results_validation]
                if provided_tblastn or provided_blastn
                else []
            )
        )
    output:
        inputs_validation
    shell:
        """
        printf "status\tOK\n" > "{output}"
        """
        
# -------------------------------------------------------------------------------------- #

if ref_genome is not None:

    rule build_lrrome:
        input:
            validated=inputs_validation,
            ref_genome=ref_genome,
            ref_genome_fai=ref_genome_fai,
            ref_gff=ref_gff
        output:
            ref_proteins=os.path.join(generated_lrrome, "REF_proteins.fasta"),
            ref_loci_fasta=os.path.join(generated_lrrome, "REF_loci.fasta"),
            ref_cdna_fasta=os.path.join(generated_lrrome, "REF_cDNA.fasta"),
            ref_exons_fasta=os.path.join(generated_lrrome, "REF_exons.fasta"),
            ref_pep=directory(os.path.join(generated_lrrome, "REF_PEP")),
            ref_exons=directory(os.path.join(generated_lrrome, "REF_EXONS")),
            ref_cdna=directory(os.path.join(generated_lrrome, "REF_cDNA")),
            ref_loci=directory(os.path.join(generated_lrrome, "REF_LOCI")),
            ref_loci_gff=directory(os.path.join(generated_lrrome, "REF_LOCI_GFF")),
            provenance=os.path.join(generated_lrrome, "REF_LOCI_PROVENANCE.tsv")
        shell:
            """
            "{lrr_bin}/create_LRRome.sh" "{input.ref_genome}" "{input.ref_gff}" "{generated_lrrome}" "{lrr_script}"
            """


# -------------------------------------------------------------------------------------- #

rule make_blastdb:
    input:
        validated=inputs_validation,
        target_genome=target_genome
    output:
        blast_db=directory(target_blast_db_dir)
    shell:
        """
        makeblastdb -in "{input.target_genome}" -dbtype nucl -out "{output.blast_db}/{target_genome_basename}"
        """


# -------------------------------------------------------------------------------------- #

checkpoint split_tblastn:
    input:
        validated=inputs_validation,
        ref_proteins=ref_proteins
    output:
        chunks=directory(tblastn_dir)
    shell:
        """
        mkdir -p "{output.chunks}"
        split -a 5 -d -l 10 "{input.ref_proteins}" "{output.chunks}/REF_proteins_split."
        """


rule tblastn:
    input:
        ref_proteins=os.path.join(tblastn_dir, "REF_proteins_split.{id}"),
        blast_db=rules.make_blastdb.output.blast_db
    output:
        os.path.join(tblastn_dir, "blast_split_{id}_res.tsv")
    shell:
        """
        tblastn \
            -soft_masking true \
            -db "{input.blast_db}/{target_genome_basename}" \
            -query "{input.ref_proteins}" \
            -evalue 1 \
            -out "{output}" \
            -outfmt '6 qseqid sseqid qlen length qstart qend sstart send nident pident gapopen evalue bitscore positive'
        """


rule merge_tblastn:
    input:
        aggregate_tblastn
    output:
        generated_tblastn
    shell:
        """
        cat {input:q} > {output:q}
        """


# -------------------------------------------------------------------------------------- #

checkpoint split_blastn:
    input:
        validated=inputs_validation,
        ref_loci=ref_loci_fasta
    output:
        chunks=directory(blastn_dir)
    shell:
        """
        mkdir -p "{output.chunks}"
        split -a 5 -d -l 20 "{input.ref_loci}" "{output.chunks}/REF_loci_split."
        """


rule blastn:
    input:
        ref_loci=os.path.join(blastn_dir, "REF_loci_split.{id}"),
        blast_db=rules.make_blastdb.output.blast_db
    output:
        os.path.join(blastn_dir, "blast_split_{id}_res.tsv")
    shell:
        """
        blastn \
            -db "{input.blast_db}/{target_genome_basename}" \
            -query "{input.ref_loci}" \
            -evalue 1 \
            -out "{output}" \
            -outfmt "6 qseqid sseqid qlen length qstart qend sstart send nident pident gapopen evalue bitscore positive"
        """


rule merge_blastn:
    input:
        aggregate_blastn
    output:
        generated_blastn
    shell:
        """
        cat {input:q} > {output:q}
        """


# -------------------------------------------------------------------------------------- #

rule candidate_loci:
    input:
        validated=inputs_validation,
        ref_gff=ref_gff,
        tblastn=tblastn_results,
        blastn=blastn_results
    output:
        pairs=candidate_pairs,
        gff=candidate_gff
    params:
        min_similarity=config["CL_min_sim"]
    shell:
        """
        python3 "{lrr_script}/candidate_loci_VR.py" \
            --gff_file "{input.ref_gff}" \
            --table "{input.tblastn}" \
            --blastn_table "{input.blastn}" \
            --output_gff "{output.gff}" \
            --output_list "{output.pairs}" \
            --min_sim {params.min_similarity}
        """


rule extract_candidate_sequences:
    input:
        candidate_gff=candidate_gff,
        target_genome=target_genome,
        target_genome_fai=target_genome_fai
    output:
        sequences=directory(candidate_sequences)
    shell:
        """
        "{lrr_script}/CANDIDATE_LOCI/extract_loci.sh" \
            "{input.candidate_gff}" \
            "{input.target_genome}" \
            "{output.sequences}"
        """


# -------------------------------------------------------------------------------------- #

checkpoint split_candidates:
    input:
        pairs=candidate_pairs
    output:
        chunks=directory(candidate_chunks)
    shell:
        """
        mkdir -p "{output.chunks}"
        split -a 5 -d -l 1 "{input.pairs}" "{output.chunks}/list_query_target_split."
        """


# -------------------------------------------------------------------------------------- #

rule sort_reference_gff:
    input:
        validated=inputs_validation,
        ref_gff=ref_gff
    output:
        sorted_gff=sorted_ref_gff
    shell:
        """
        python3 "{lrr_script}/sort_gff.py" --gff "{input.ref_gff}" --output "{output.sorted_gff}"
        """


# -------------------------------------------------------------------------------------- #

rule gene_prediction:
    input:
        pair=os.path.join(candidate_chunks,"list_query_target_split.{split_id}"),
        target_loci=candidate_sequences,
        ref_gff=sorted_ref_gff,
        ref_locus_info=[ref_locus_info] if ref_locus_info else [],
        ref_proteins=ref_proteins,
        ref_pep=ref_pep_dir,
        ref_exons=ref_exons_dir,
        ref_cdna=ref_cdna_dir,
        ref_loci=ref_loci_dir,
        ref_loci_gff=ref_loci_gff_dir
    output:
        prediction_outputs
    params:
        lrrome=lrrome,
        mode=mode,
        ignore_exonerate_errors=str(config["ignore_exonerate_errors"]).lower(),
        output_prefix=prediction_prefix
    shell:
        """
        "{lrr_bin}/genePrediction.sh" \
            --pair-file "{input.pair}" \
            --target-loci-dir "{input.target_loci}" \
            --lrrome "{params.lrrome}" \
            --ref-gff "{input.ref_gff}" \
            --ref-locus-info "{input.ref_locus_info}" \
            --output-prefix "{params.output_prefix}" \
            --mode "{params.mode}" \
            --script-dir "{lrr_script}" \
            --ignore-exonerate-errors "{params.ignore_exonerate_errors}"
        """


# -------------------------------------------------------------------------------------- #

rule merge_predictions:
    input:
        aggregate_predictions
    output:
        expand(
            os.path.join(out_dir, "annot_{method}{suffix}"),
            method=prediction_methods,
            suffix=[".gff", "_chr.gff", "_cleaning.log"],
        )
    params:
        methods=" ".join(prediction_methods)
    shell:
        """
        for method in {params.methods}; do
            "{lrr_bin}/merge_prediction.sh" \
                "{predictions_dir}" \
                "{lrr_script}" \
                "$method" \
                "{prefix}" \
                "{out_dir}"
        done
        """


# -------------------------------------------------------------------------------------- #

rule transfer_stats:
    input:
        os.path.join(out_dir, "annot_best_chr.gff")
    output:
        gff_stats=os.path.join(stats_dir, "GFFstats.txt"),
        jobs_out=os.path.join(stats_dir, "jobsStats_out.txt"),
        jobs_err=os.path.join(stats_dir, "jobsStats_err.txt")
    shell:
        """
        "{lrr_script}/STATS_OUTPUTS/stats_transfer.sh" "{out_dir}" "{out_dir}/.." "{stats_dir}"
        """



# -------------------------------------------------------------------------------------- #
# Workflow provenance
# -------------------------------------------------------------------------------------- #

onstart:
    start_run_record(
        run_id=run_id,
        history_path=run_history_file,
        latest_path=latest_run_info,
        log_path=workflow_log,
        config=config,
        repo_root=repo_root,
        argv=sys.argv,
        working_directory=os.getcwd(),
    )


onsuccess:
    finish_run_record(
        history_path=run_history_file,
        latest_path=latest_run_info,
        log_path=workflow_log,
        status="success",
    )


onerror:
    finish_run_record(
        history_path=run_history_file,
        latest_path=latest_run_info,
        log_path=workflow_log,
        status="failed",
    )
