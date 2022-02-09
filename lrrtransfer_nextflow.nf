#!/usr/bin/env nextflow
help=false
params.lrrome = "NULL"
LAUNCH_DIR="$workflow.launchDir"
WD="${LAUNCH_DIR}/LRRtransfer_output"
LG_SCRIPT="/home/ubuntu/LRRtransfer/SCRIPT"
LG_BIN="/home/ubuntu/LRRtransfer/bin"

if(!params.lrrome || !params.ref_gff) {
    println """please specify one of the following parameters
                --lrrome: the path to a LRRome directory if you already have it
                --ref_gff: path to the gff containing LRR annotations for the reference genome"""}

if(!params.target_genome) {
    println """ please specify a path of genome to annotate  using --genome"""
    help=true;
}

if(!params.mode) {
    println """ --mode is not specified, default is <best>"""
}

if( help == true){
  println """\
    Usage:
        ======================================
        nextflow run lrrtransfer_v0_nextflow.nf --ref_genome <reference_genome_fasta> --ref_gff <reference_gff> --ref_locus_info <reference_locus_info.txt> --target_genome <genome_to_annotate.fasta> [ --lrrome <PATH_to_precompute_LRRome> --mode <first, best>]
        ======================================
    """
    exit 1
}else{
  println """\
LRRtransfer pipeline runnning ...
        ===================================
        target_genome:     ${params.target_genome}
        ref_genome:        ${params.ref_genome}
        ref_gff:           ${params.ref_gff}
        ref_locus_info:    ${params.ref_locus_info}
        mode:              ${params.mode}
        lrrome:            ${params.lrrome}
        """
}


def helpMessage() 
{
    log.info"""
  Usage: 
    nextflow run lrrtransfer_v0_nextflow.nf --ref_genome <reference_genome_fasta> --ref_gff <reference_gff> --target_genome <genome_to_annotate.fasta> [ --lrrome <PATH_to_precompute_LRRome> --mode <first, best>]
  Mandatory arguments:
    --target_genome      Path to the fasta genome that need to be annotate 
    --ref_genome         Path to the fasta genome use as reference
    --ref_gff            Path to the LRR locus annotations from the reference species
    --ref_locus_info     Path to the text file with reference LRR locus information (gene family and class)
  Optional arguments:
    --lrrome             Path to a precompute LRRome folder
    --mode               mode [first, best] default is <best>
"""
}


/*The following process builds an LRRome only if an
  LRRome is not given as input and a results directory
  is built.*/

process buildLRROme { 
    input:
      val refgenome from file(params.ref_genome)
    output:
      path LRRome into LRRome_dirch
    script:
    """
    mkdir -p $WD
    ${LG_BIN}/create_LRRome.sh ${refgenome} ${params.ref_gff} ${WD} ${params.lrrome}
    """
}


/*The following process find regions of interest in the 
  target genome*/


process candidateLoci  { 
    input:
      val LRRome from LRRome_dirch
      val targetgenome from file(params.target_genome)
      val refgff from file(params.ref_gff)
    output:
      path CANDIDATE_SEQ_DNA into CANDIDATE_SEQ_DNAch
      file 'list_query_target.txt' into list_query_targetch
      file 'filtered_candidatsLRR.gff' into filtered_candidatsLRRch
    script:
    """
    ${LG_BIN}/candidateLoci.sh ${targetgenome} ${LRRome} ${refgff} ${WD}
    """
} 


/*Individual recuperation of all "query target" couples in order
  to parallelize the genePrediction process for each couple.*/

list_query_targetch.splitText(file: true).set{ one_candidatech }


/*The following process produce a gene model for all regions of interest (GFF file)*/

process genePrediction {
    errorStrategy 'ignore'
    input:
      file one_candidate from one_candidatech
      val CANDIDATE_SEQ_DNA from CANDIDATE_SEQ_DNAch
      val targetgenome from file(params.target_genome)
      val filtered_candidatsLRR from filtered_candidatsLRRch
      val LRRome from LRRome_dirch
      val refgff from file(params.ref_gff)
      val refinfo from file(params.ref_locus_info)
    output:
      path one_candidate_gff into one_candidate_gffch
    script:
    """
    ${LG_BIN}/genePrediction.sh ${one_candidate} ${CANDIDATE_SEQ_DNA} ${targetgenome} ${params.mode}  ${filtered_candidatsLRR} ${LRRome} ${WD} ${refgff} ${refinfo}
    """
}

one_candidate_gffch.collectFile(name: 'annot.gff').set{ genePredictionch }

/*process foo {
  input:
    file x from genePredictionch
  output:
    file 'annot' into completeGFF_ch
  script:
    """
    < $x zcat > $annot
    """
}

completeGFF_ch.collectFile()*/


/*The following process produce a currated GFF file*/

process verifAnnot {
    errorStrategy 'ignore'
    input:
      val refinfo from file(params.ref_locus_info)
      val targetgenome from file(params.target_genome)
      file 'annot' from genePredictionch
    script:
    """
      ${LG_BIN}/verifAnnot.sh ${refinfo} ${targetgenome} ${WD} ${annot}
    """
}