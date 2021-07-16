# GeneModelTransfer
Lrrtransfert is a complex gene family annotation transfer pipeline designed for LRR (Leucine-Rich Repeat) containing plant receptors. 
The pipeline allows to annotate the LRRs of a new genome from a quality annotation of this gene family of a nearby genome  

These pipelines are mostly bash and python scripts encapsulated within Singularity containers and combined into NextFlow workflows.

The program use a Singularity container and Nextflow for execution.


## Singularity (.sif) container

The Singularity container can be downloaded to the machine to be used by nextflow in one of the following ways 
### Pull with Singularity
$ singularity pull --arch amd64 library://thiabud/default/lrrtransfert:v1
### Pull by unique ID (reproducible even if tags change)
$ singularity pull library://thiabud/default/lrrtransfert:sha256.6ae4453dd18a36367800eba83c62d829845b1c8b7a3214d4534adfbd52450e59

The Singularity container
The container singularity address can also be specified in the Nextflow configuration file.
Please add 'library://thiabud/default/lrrtransfert:v1'
Or library://thiabud/default/lrrtransfert:sha256.6ae4453dd18a36367800eba83c62d829845b1c8b7a3214d4534adfbd52450e59

The file LRRtransfert.def provide the corresponding singularity recipe.

## Nextflow 
LRRtransfert.nf requires the customization of the nextflow.config file according to your execution environment and the resources to allocate to the run.
-The field "executor='sge" should be replaced by your environment executor for example 'slurm' or nothing for a local execution
-The field CPU allows you to choose the number of cpus to use
-The fields container must contain the address or path to the llrtransfer.sif container
## References
