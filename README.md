# LRRtransfer
Lrrtransfer is a complex gene family annotation transfer pipeline designed for plant LRR (Leucine-Rich Repeat) containing receptors. 
The pipeline use a reference annotation of LRR gene families to annotate LRR of a closely related genome.  

This pipeline is mostly bash and python scripts encapsulated within Singularity containers and combined into NextFlow workflows.

The pipeline use a Singularity container and Nextflow for execution.  
The pipeline works with nextflow 19.10.0 and higher.


## Singularity (.sif) container

The file lrrtransfer_singularity.def provides the recipe used to build this singularity container.

The container will be load automatically if the address of this container is specified in your nextflow config file (see next section):
```
container = 'library://cgottin/default/lrrtransfer:2.0'    
```

If needed, you can download it from [the sylabs singularity repository](https://sylabs.io/) using the following command: 
```
$ singularity pull library://cgottin/default/lrrtransfer:2.0
```




## Nextflow 
lrrtransfer_nexflow.nf requires the customization of the nextflow.config file according to your execution environment. Some sample configuration files are provided in this repository. The key parameters to specify are:   
-The job manager "executor" e.g.'sge', 'slurm' or 'local'.  
-The number of CPU to be used (e.g. CPU=100) .  
-The container to be used to run the pipeline (i.e. the address or path to the lrrtransfer_v2.sif container). 
See [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more details.

## Running LRRtransfer
### Clone the repository
```
git clone https://github.com/cgottin/GeneModelTransfer.git
```
### The program can be run with the command line :
```
nextflow run lrrtransfer_nextflow.nf --ref_genome <fasta> --ref_gff <gff> --ref_locus_info <txt> --target_genome <fasta> --mode <[first, best]>
```
### Using the example files :  
Example files are provided in the "example" folder of this git repository.
You have to modify the nextflow.config file according to your execution environment. 
```
cd GeneModelTransfer/
nextflow run lrrtransfer_nextflow.nf --ref_genome example/Nipponbare_Chr1.fasta --ref_gff example/Nipponbare_LRR-CR_Chr1.gff --ref_locus_info example/Info_locus_Nipponbare.txt --target_genome example/Punctata_chr1.fasta --mode best
```

The gff file resulting from the test ("LRRlocus_predicted.gff") is located in a new directory named 'LRRtransfer_output'.

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
	
 Nextflow    
 <a id="Di_Tommaso_2017"></a> Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., and Notredame, C.(2017). Nextflow enables reproducible computational workflows. Nature Biotechnology,35(4):316–319. [Nextflow web site](https://www.nextflow.io/) 
 
Singularity        
<a id="Kurtzer_2017"></a> Kurtzer, G. M., Sochat, V., and Bauer, M. W. (2017). Singularity: Scientific containers formobility of compute. PloS One, 12(5):e0177459. [singularity web site](https://sylabs.io/)
