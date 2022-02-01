#!/usr/bin/env nextflow
help=false
params.lrrome = "NULL"
LAUNCH_DIR="$workflow.launchDir"


if(!params.lrrome || !params.input) {
    println """please specify one of the following parameters
				--lrrome: the path to a LRRome directory if you already have it
				--input:text file with 4 columns.
				First column contain a code the accession.
				Second column contain a path to the reference GFF containing LRR
				Third column contain a path to the referene asembly (fasta format)
				Fourth column is not mandatory and should contain a path to a file containing information for LRR (family and class of each locus)"""
}

if(!params.genome) {
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
          nextflow run pipeline.nf --genome <genome_to_annotate.fasta> --input <tab-delimited_file.txt>  [--mode <first, best>]
          ======================================
    """
    exit 1
}else{
  println """\
LRRtransfer pipeline runnning ...
         ===================================
         genome:            ${params.genome}
         input:             ${params.input}
         mode:              ${params.mode}
         lrrome:            ${params.lrrome}
         debug:             ${params.debug}
         """
}

def helpMessage() 
{
    log.info"""
  Usage: 
  nextflow run pipeline.nf --genome <genome_to_annotate.fasta> --input <tab-delimited_file.txt>  [--mode <mode_type>]
  Mandatory arguments:
    --genome      File path to genome that need to be annotate 
    --input       A tab separated file with 4 columns: 1 species code, 2 path to GFF (only LRR), 3 path to genome fasta file, 4 path to Info_LRR file (family and class for each locus)
  Other options:
    --mode        mode [first, best]
"""
}



//inputch = file(params.input)

/*The following process builds an LRRome only if an
  LRRome is not given as input and a results directory
  is built.*/

process buildLRROme { 
    input:
      val input_file from file(params.input)
    output:
      path LRRome into LRRome_dirch
    script:
    """
    mkdir -p ${LAUNCH_DIR}/LRRtransfer_$(date +"%Y%m%d")
    ${LG_BIN}/create_LRRome.sh ${input_file} ${params.lrrome} ${LAUNCH_DIR}
    """
}


/*The following process find regions of interest in the 
  target genome*/
  
process candidateLoci  { 
    input:
      val LRRome from LRRome_dirch
    output:
      path CANDIDATE_SEQ_DNA into CANDIDATE_SEQ_DNAch
      path candidate_loci_to_LRRome into candidate_loci_to_LRRomech
      path filtered_candidatsLRR into filtered_candidatsLRRch
    script:
    """
    ${LG_BIN}/candidateLoci.sh ${params.genome} ${LRRome} ${params.input} ${LAUNCH_DIR}
    """
} 

/*Individual recuperation of all "query target" couples in order
  to parallelize the genePrediction process for each couple.*/

candidate_loci_to_LRRomech.splitText(file: true).set{ candidate_locich }


/*The following process produce a gene model 
  for all regions of interest (GFF file)*/

process genePrediction {
    errorStrategy 'ignore'
    input:
      val filtered_candidatsLRR from filtered_candidatsLRRch
      val LRRome from LRRome_dirch
      val CANDIDATE_SEQ_DNA from CANDIDATE_SEQ_DNAch
      file one_candidate from candidate_locich
    output:
      path one_candidate_gff into one_candidate_gffch
    script:
    """
    ${LG_BIN}/genePrediction.sh ${one_candidate} ${CANDIDATE_SEQ_DNA} ${params.genome} ${params.mode}  ${filtered_candidatsLRR} ${LAUNCH_DIR} ${LRRome} ${params.input}
    """
}

one_candidate_gffch.collect().set{ genePredictionch }

/*The following process produce a currated GFF file*/

process verifAnnot {
  errorStrategy 'ignore'
  input:
    val one_prediction_gff from genePredictionch
  script:
  """
   ${LG_BIN}/verifAnnot.sh ${params.input} ${params.genome} ${LAUNCH_DIR}
  """
}