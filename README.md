# LRRtransfer
Lrrtransfer is a complex gene family annotation transfer pipeline designed for plant LRR (Leucine-Rich Repeat) containing receptors. 
The pipeline use a reference annotation of LRR gene families to annotate LRR of a closely related genome.  

This pipeline is mostly bash and python scripts orchestrated via Snakemake workflows to leverage HPC.

The pipeline use Conda for package management. The version of blast, mmseq and exoneate used by the pipeline are listed in SNAKEMAKE/conda_tools.yml.

The pipeline works with snakemake 6.9.1 and higher.

##  Installing LRRtransfer
The easiest way to install LRRtransfer is to clone this git repository: 
```
git clone https://github.com/cgottin/GeneModelTransfer.git
```
To actually run LRRtransfer snakemake pipeline you need to have:
- snakemake and
- conda installed. 

Conda and snakemake will take care of all other dependencies.

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

Note that, the first time you run LRRtransfer, a conda environment will be created. This could take some times, but it is done only once. 


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
