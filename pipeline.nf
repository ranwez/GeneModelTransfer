#!/usr/bin/env nextflow
help=false
params.lrrome = "NULL"
params.debug = 0
LAUNCH_DIR="$workflow.launchDir"

if(!params.lrrome) {
    println """ please specify the path to a LRRome directory if you already have it using  --lrrome"""
    help=true;
}

if(!params.input) {
    println """ please specify a path to a text file with 4 columns.
    First column contain a code the accession.
    Second column contain a path to the reference GFF containing LRR 
    Third column contain a path to the referene asembly (fasta format)
    Fourth column is not obligatory and should contain a path to a file containing information for LRR (family and class of each location) using --input"""
    help=true;
}

if(!params.genome) {
    println """ please specify a path of genome to annotate  using --genome"""
    help=true;
}

if(!params.mode) {
    println """ please specify the mode of execution (first, best or consensus) using --mode"""
    help=true;
}

if( help == true)
{
  println """\

    Usage:
          ======================================
          nextflow run pipeline.nf --genome_to_annotate.fasta --input tab-delimited_file.txt  --mode chosen_mode
          ======================================
    """
    exit 1
}
else{
  println """\
         
GeneModelTransfer pipeline runnning ...
         ===================================
         genome:            ${params.genome}
         input:             ${params.input}
         mode:           ${params.mode}
         lrrome:             ${params.lrrome}
         debug:              ${params.debug}
         """
}

def helpMessage() 
{
	log.info"""

  Usage: 
  nextflow run pipeline.nf --genome_to_annotate.fasta --input tab-delimited_file.txt  --mode chosen_mode
  Mandatory arguments:
    --genome      File path to genome that need to be annotate 
    --input // Un fichier texte à 4 colonnes : 1 code pour l'accesion et trois chemin d'accès : GFF(only LRR), Assemblage.fasta, Info_LRR(famille et classe de chaque locus)
    --mode mode 
  Other options:
    *--treshold            treshold    treshold      
"""
}

//The following process builds an LRRome only if an LRRome is not given as input and a results directory is built.
process buildLRROme { 
   //echo true
    input:
    val input_file from file(params.input)
    output:
    path LRRome into LRRome_dirch
    script:
    """
    /GeneModelTransfer.git/branches/container/bin/create_LRRome.sh $input_file ${params.lrrome} $LAUNCH_DIR 
    """
}

//The following process find regions of interest in target genome
process candidateLoci  { 
    //echo true
    input:
    val LRRome from LRRome_dirch
    output:
    path CANDIDATE_SEQ_DNA into CANDIDATE_SEQ_DNAch
    path candidate_loci_to_LRRome into candidate_loci_to_LRRomech
    path filtered_candidatsLRR into filtered_candidatsLRRch
    script:
    """
    /GeneModelTransfer.git/branches/container/bin/candidateLoci.sh ${params.genome} $LRRome ${params.input} $LAUNCH_DIR 
    """
} 

//Individual recuperation of all "query target" couples in order to parallelize the genePrediction process for each couple
candidate_loci_to_LRRomech.splitText().set{ candidate_locich }


//The following process produce a gene prediction for all regions of interest (GFF file)
process genePrediction {
    errorStrategy 'ignore'
    //echo true
    input:
    val filtered_candidatsLRR from filtered_candidatsLRRch
    val LRRome from LRRome_dirch
    val CANDIDATE_SEQ_DNA from CANDIDATE_SEQ_DNAch
    file one_candidate from candidate_locich
    output:
    path one_candidate_gff into one_candidate_gffch
    script:
    """
    /GeneModelTransfer.git/branches/container/bin/genePrediction.sh $one_candidate $CANDIDATE_SEQ_DNA ${params.genome} ${params.mode}  $filtered_candidatsLRR $LAUNCH_DIR $LRRome ${params.input}
    """
}

one_candidate_gffch.collect().set{ genePredictionch }
//The following process produce a currated GFF file 
process verifAnnot {
  errorStrategy 'ignore'
  echo true
  input:
  val one_prediction_gff from genePredictionch
  script:
  """
   /GeneModelTransfer.git/branches/container/bin/verifAnnot.sh ${params.input} ${params.genome} $LAUNCH_DIR ${params.debug}
  """
}
