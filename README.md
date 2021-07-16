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
## References
