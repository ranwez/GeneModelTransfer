# GeneModelTransfer
Lrrtransfer is a complex gene family annotation transfer pipeline designed for LRR (Leucine-Rich Repeat) containing plant receptors. 
The pipeline allows to annotate the LRRs of a new genome from a quality annotation of this gene family of a nearby genome  

This pipeline is mostly bash and python scripts encapsulated within Singularity containers and combined into NextFlow workflows.

The pipeline use a Singularity container and Nextflow for execution.  
The pipeline works with nextflow 19.10.0 and higher.


## Singularity (.sif) container

The file LRRtransfer.def provides the recipe used to build this singularity container.

The container can be downloaded from [the sylabs singularity repository](https://sylabs.io/) using the following command: 
```
$ singularity pull library://thiabud/default/lrrtransfer:v1
```

Alternatively, you can specify the address of this container in your nextflow config file (see next section) by adding the following line:
```
'library://thiabud/default/lrrtransfer:v1'    
```


## Nextflow 
lrrtransfer.nf requires the customization of the nextflow.config file according to your execution environment. Some *sample configuration files are provided in this repository*. The key parameters to specify are:   
-The job manager e.g."executor='sge"; "executor='slurm" or nothing  nothing for a local execution.  
-The number of CPU to be used (e.g. CPU=100) .  
-The container to be used to run the pipeline (i.e. the address or path to the llrtransfer.sif container). 
See [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more details.

## Running LRRtransfer
### The program can be run with the command line :
```
nextflow run lrrtransfer.nf --genome <target_genome> --mode <choosen_mode> --input <tab-separated file>
```
 < target_genome >  is the absolute path to the genome to annotate in fasta format.   
 < choosen_mode >  can be 'first' or 'best'. The 'best' method allows to obtain gene models with a higher identity to the reference than the 'first' method.  
 < tab-separated file >  is a path to a text file with 4 columns.  
	First column contains a code the accession.
 	Second column contains a path to the reference GFF containing LRR 
 	Third column contains a path to the referene asembly (fasta format)
	Fourth column is not obligatory and should contains a path to a file containing information for LRR (available in this repository for Nipponbare: 'Info_locus_Nipponbare.txt)
### Using the example files :   
Execute the following lines to get the files needed for the test.  
You can copy and paste the following instructions into your terminal.
```
#Download the Info_locus file available the lrrtransfer repository :
wget https://github.com/cgottin/GeneModelTransfer/raw/container/Info_locus_Nipponbare.txt

#Download the input.txt file availablethe the lrrtransfer repository:
wget https://github.com/cgottin/GeneModelTransfer/raw/container/input.txt

#Download the chromosome1_punctata.fasta available on lrrtransfer repository :
wget https://github.com/cgottin/GeneModelTransfer/raw/container/chromosome1_punctata.fasta

#Download the reference genome : 
wget https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_genome.fasta.gz;
gunzip IRGSP-1.0_genome.fasta.gz

#Download the reference genome annotation : 
wget https://github.com/cgottin/GeneModelTransfer/raw/container/Oryza_Nipponbare_IRGSP-1.0_LRR-CR__20210715.gff

#Download the lrrtransfer.nf
wget https://github.com/cgottin/GeneModelTransfer/raw/container/lrrtransfer.nf

#Download the nextflow.config file 
wget https://github.com/cgottin/GeneModelTransfer/raw/container/nextflow.config

#Replace the name of the files contained in input.txt by their absolute paths 
realpath $(cat input.txt) > input_tmp.txt ; awk -vRS="\n" -vORS="\t" '1' input_tmp.txt > input.txt ; rm input_tmp.txt ; cat input.txt |   sed 's/.*T/T/' > input_tmp.txt ; cat input_tmp.txt > input.txt ; rm input_tmp.txt ; echo "" >> input.txt

#Format the name of chromosomes
sed -i -e 's/>c/>C/g' IRGSP-1.0_genome.fasta  ; sed -i -e 's/r0/r/g' IRGSP-1.0_genome.fasta
```


Then execute the following command to run the pipeline with the test files
```
nextflow run lrrtransfer.nf --genome $PWD/chromosome1_punctata.fasta --mode best --input $PWD/input.txt
```

## Singularity overview

A singularity container [[Kurtzer, 2017]](#Kurtzer_2017) contains everything that is needed to execute a specific task. The person building the container has to handle dependencies and environment configuration so that the end-user do not need to bother. The file specifying the construction of the container is a simple text file called a recipe (we provide the recipe of our container as well as the containers). As our scripts/pipelines often relies on several other scripts and external tools (e.g. MAFFT) singularity container is very handy as the end user just need to install singularity and download the container without having to care for installing dependencies or setting environment variables.

A brief introduction to singularity is available [here](https://bioweb.supagro.inra.fr/macse/index.php?menu=pipelines).

If you got an error message stating that your input file does not exist it is probably related to the fact that the folder containing them is not visible from the singularity container. A solution found by one user is to use the [SINGULARITY_BINDPATH variable](https://sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html):   
```
export SINGULARITY_BINDPATH="/path/to/fasta"
```

## Nextflow overview

Nextflow [[Di Tommaso, 2017]](#Di_Tommaso_2017) enables scalable and reproducible scientific workflowsusing software containers allowing the adaptation of pipelines written in the most commonscripting languages.

Nextflow separates the workflow itself from the directive regarding the correct way to execute it in the environment. One key advantage of Nextflow is that, by changing slightly the “nextflow.config” file, the same workflow will be parallelized and launched to exploit the full resources of a high performance computing (HPC) cluster.

## References

Github repository   
	https://github.com/cgottin/LRRprofiler/

Initial works  
	A New Comprehensive Annotation of Leucine-Rich Repeat-Containing Receptors in Rice. Gottin C., Diévart A., Summo M., Droc G., Périn C., Ranwez V. and Chantret N. preprint on bioRxiv. doi: https://doi.org/10.1101/2021.01.29.428842

MMseqs2  
	Steinegger, M., & Söding, J. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature biotechnology, 35(11), 1026-1028.

exonerate  
	Slater, G. S. C., & Birney, E. (2005). Automated generation of heuristics for biological sequence comparison. BMC bioinformatics, 6(1), 1-11.

BLAST+  
	Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421
	
 Singularity   
 <a id="Di_Tommaso_2017"></a> Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., and Notredame, C.(2017). Nextflow enables reproducible computational workflows. Nature Biotechnology,35(4):316–319. [Nextflow web site](https://www.nextflow.io/)
Nextflow   
<a id="Kurtzer_2017"></a> Kurtzer, G. M., Sochat, V., and Bauer, M. W. (2017). Singularity: Scientific containers formobility of compute. PloS One, 12(5):e0177459. [singularity web site](https://sylabs.io/)
