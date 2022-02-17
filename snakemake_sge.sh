#!/bin/bash


snakePipe=$1
cluster_config=$2
profile=$3

mkdir jobinfo/

snakemake  --profile $profile --cluster-config $cluster_config --snakefile $snakePipe --configfile config.yml

mv LRRtransfer*.err jobinfo/
mv LRRtransfer*.out jobinfo/

## run exemple
#./snakemake.sh LRR_transfert.smk /work_home/cgottin/SNAKEMAKE/cluster_config.yml /work_home/cgottin/SNAKEMAKE/PROFILES
