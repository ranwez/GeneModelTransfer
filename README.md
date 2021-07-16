# GeneModelTransfer
Lrrtransfert is a complex gene family annotation transfer pipeline designed for LRR (Leucine-Rich Repeat) containing plant receptors. 
The pipeline allows to annotate the LRRs of a new genome from a quality annotation of this gene family of a nearby genome  

These pipelines are mostly bash and python scripts encapsulated within Singularity containers and combined into NextFlow workflows.

The program use a Singularity container and Nextflow for execution.


## Singularity (.sif) container

The file LRRtransfert.def provides the recipe used to build this singularity container.

The container can be downloaded from [the sylabs singularity repository](https://sylabs.io/) using the following command: 
```
$ singularity pull --arch amd64 library://thiabud/default/lrrtransfert:v1
```

Alternatively, you can specify the address of this container in your nextflow config file (see next section) by adding the following line:
```
'library://thiabud/default/lrrtransfert:v1'    
```


## Nextflow 
LRRtransfert.nf requires the customization of the nextflow.config file according to your execution environment. Some *sample configuration files are provided in this repository*. The key parameters to specify are:   
-The job manager e.g."executor='sge"; "executor='slurm" or nothing  nothing for a local execution.  
-The number of CPU to be used (e.g. CPU=100) .  
-The container to be used to run the pipeline (i.e. the address or path to the llrtransfer.sif container). 
See [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more details.

## Running LRRtransfer
### The program can be run with the command line :
```
nextflow run lrrtransfer.nf --genome <target_genome> --mode <choosen_mode> --input <tab-separated file>
```
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
```
Then execute the following command to run the pipeline
```
nextflow run lrrtransfer.nf --genome chromosome1_punctata.fasta --mode best --input input.txt
```
## References
