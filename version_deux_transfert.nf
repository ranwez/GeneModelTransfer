#!/usr/bin/env nextflow
def helpMessage() 
{
	log.info"""
  Usage: 
  nextflow run transfer.nf --reads1 R1 --reads2 R2 --prefix prefix --outdir path [options] 
  
  Mandatory arguments:
    --specie   SN  Name of the specie for the analysis naming
    --genome    genome   file path to genome that need to be annotate
    *--cdna  cdna    file path to the cdna of reference
    *--gff   gff      file path to the gff containing the LRR gene model of reference
    *--pep   pep file path to the fasta file containing the peptidique sequence of LRR gene of reference
    *--exon  exon   file path to file containing the exon of LRR gene of reference 
    !
    --info   info file path containing info locus of nipponbare 
  
  Other options:
    --treshold            treshold    treshold      
"""
}
wd = file("$baseDir/database")
wd.mkdirs()
//Specie 
params.specie = "UnknowSpecie"
specie = (params.specie)
process wichsepcie {
    output:
    stdout speciech
    """
    printf "Specie : "'${params.specie}'
    """
}
params.mode = "first"
mode = (params.mode)
process wichmode {
    output:
    stdout modech
    """
    printf "mode : "'${params.mode}'
    """
}
//speciech.view()
//Genome
params.genome= "$baseDir/../../pipeline/data/ZS.fasta"
genome = file(params.genome)
process wichgenome {
    output:
    stdout genomech
    """
    printf "Genome : "'${params.genome}'
    """
}
//genomech.view()
//cdna
params.cdna="$baseDir/../../pipeline/data/3cdna.fasta"
cdna = file(params.cdna)
process wichcdna {
    output:
    stdout cdnach
    """
    printf "cdna : "'${params.cdna}'
    """
}
//cdnach.view()
//gff
params.gff="$baseDir/../../pipeline/data/4gff.gff"
gff = file(params.gff)
process wichgff {
    output:
    stdout gffch
    """
    printf "gff : "'${params.gff}'
    """
}
//gffch.view()
//pep
params.pep="$baseDir/../../pipeline/data/5pep.fasta"
pep = file(params.pep)
process wichpep {
    output:
    stdout pepch
    """
    printf "pep : "'${params.pep}'
    """
}
//pepch.view()
//exon
params.exon="$baseDir/../../pipeline/data/exon.fasta"
exon = file(params.exon)
process wichexon {
    output:
    stdout exonch
    """
    printf "exon : "'${params.exon}'
    """
}
params.infoLocus="$baseDir/../../pipeline/data/Info_locus_Nipponbare.txt"
infoLocus = file(params.infoLocus)
process wichInfoLocus {
    output:
    stdout infoLocusch
    """
    printf "infoLocus : "'${params.infoLocus}'
    """
}

//exonch.view()

//setupworkingdirectory
resDir=file("$baseDir/Transfert_${params.specie}")
resDir.mkdirs()


scriptdir="$workflow.launchDir/SCRIPT"



//vegetable_datasets = Channel.fromPath(params.genome).splitFasta(by: 1, file: true)

process part1 {
    //echo true
    //input:
    //file chr from vegetable_datasets
    output:
    path liste_query_target into liste_query_targetch
    //file xDNA_candidatsLRR_in into number
    path TARGET_DNA into dirtargetdnach
    path REF_PEP into refpepdir
    path REF_CDS into refcdsdir
    path REF_cDNA into refcdnadir
    path filtered_candidatsLRR into filtered_candidatsLRRch
    script:
    """
    part1.sh ${params.specie} ${params.genome} ${params.cdna} ${params.gff} ${params.pep} ${params.exon} '$scriptdir' $wd 
    """
}

//filtered_candidatsLRRch.view()

//liste_query_targetch.into {
 // liste_query_targetch1
//  liste_query_targetch2
//}



liste_query_targetch.splitText().set{ toto_ch}
process mappingcds {
  echo true
  input:
  val filtered_candidatsLRR from filtered_candidatsLRRch
  val listqt from liste_query_targetch
  val refcds from refcdsdir
  val refpep from refpepdir
  val refcdna from refcdnadir
  val targetdnadir from dirtargetdnach
  file x from toto_ch
  output:
  path LRRlocus_complet into LRRlocus_completch
  script:
  """
  part3.sh $x $scriptdir $targetdnadir ${params.genome} ${refpep} ${params.specie} $refcds $listqt ${params.mode} ${params.gff} $refcdna $filtered_candidatsLRR ${params.infoLocus} ${params.mode} $resDir
  """
}

LRRlocus_completch.collect().set{ totoch }

process concat {
  //echo true
  input:
  each x from totoch
  script:
  """
  cat $x >> $baseDir/final.gff
  echo ----------------------------------------------------------------- >> $baseDir/final.gff
  """
}
/*    
process checker {
    echo true
    input: 
    env y from line1_ch 
    script:
    """
    echo $y
    """
}
/*
str = liste_query_targetch.readLines()
expl1 = Channel.value(str)
process printEnv {
    echo true
    input:
    env HELLO from str
    script :
    '''
    echo $HELLO
    '''
}
*/

/*
process printlerefpepdir{
    echo true
    input:
    path liste_query_target from liste_query_targetch
    path toto from refpepdir
    exec : 
    //liste_query_targetch.println()
    toto2=refpepdir.getParent()
    script:
    """
    echo '$toto2'
    echo '$liste_query_target'
    """
}*/