# LRRtransfer
Lrrtransfer is a complex gene family annotation transfer pipeline designed for plant LRR (Leucine-Rich Repeat) containing receptors. 
The pipeline use a reference annotation of LRR gene families to annotate LRR of a closely related genome.  

This pipeline is mostly bash and python scripts encapsulated within Singularity containers and combined into Snakemake workflows.

The pipeline use a Singularity container and Snakemake for execution.  
The pipeline works with snakemake 6.9.1 and higher.



## Running LRRtransfer
The LRRtransfer pipeline can be run with:

```
profile="./SNAKEMAKE/PROFILES/<your_profile_folder>"
cluster_config="./SNAKEMAKE/<your_cluster_config.yaml>"
config_file="./SNAKEMAKE/<your_config_file.yaml>"
smk_Pipeline="./SNAKEMAKE/LRR_transfer.smk"
snakemake  --profile $profile --cluster-config $cluster_config --snakefile $smk_Pipeline --configfile $config_file
```

Examples for each configuration file are provided in the folder SNAKEMAKE/.



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

## Singularity overview

A singularity container [[Kurtzer, 2017]](#Kurtzer_2017) contains everything that is needed to execute a specific task. The person building the container has to handle dependencies and environment configuration so that the end-user do not need to bother. The file specifying the construction of the container is a simple text file called a recipe (we provide the recipe of our container as well as the containers). As our scripts/pipelines often relies on several other scripts and external tools (e.g. MAFFT) singularity container is very handy as the end user just need to install singularity and download the container without having to care for installing dependencies or setting environment variables.

A brief introduction to singularity is available [here](https://bioweb.supagro.inra.fr/macse/index.php?menu=pipelines).

If you got an error message stating that your input file does not exist it is probably related to the fact that the folder containing them is not visible from the singularity container. A solution found by one user is to use the [SINGULARITY_BINDPATH variable](https://sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html):   
```
export SINGULARITY_BINDPATH="/path/to/fasta"
```




## Clone the repository
```
git clone https://github.com/cgottin/GeneModelTransfer.git
```

## References

Github repository   
	https://github.com/cgottin/LRRprofiler/

Initial work describing the gene model transfer pipeline:
	A New Comprehensive Annotation of Leucine-Rich Repeat-Containing Receptors in Rice. Gottin C., Diévart A., Summo M., Droc G., Périn C., Ranwez V. and Chantret N. preprint on bioRxiv. doi: https://doi.org/10.1101/2021.01.29.428842

MMseqs2  
	Steinegger, M., & Söding, J. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature biotechnology, 35(11), 1026-1028.

exonerate  
	Slater, G. S. C., & Birney, E. (2005). Automated generation of heuristics for biological sequence comparison. BMC bioinformatics, 6(1), 1-11.

BLAST+  
	Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421
	
Snakemake
	Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.
Singularity        
<a id="Kurtzer_2017"></a> Kurtzer, G. M., Sochat, V., and Bauer, M. W. (2017). Singularity: Scientific containers formobility of compute. PloS One, 12(5):e0177459. [singularity web site](https://sylabs.io/)
