#!/bin/bash
#
#SBATCH -J LRRtransfer
#SBATCH -o LRRtransfer."%j".out
#SBATCH -e LRRtransfer."%j".err

#SBATCH -p agap_normal

module purge
module load snakemake/5.13.0
module load singularity/3.6.3

snakePipe=$1
cluster_config=$2
profile=$3
input_config=$4

mkdir jobinfo/


snakemake  --profile $profile --cluster-config $cluster_config --snakefile $snakePipe --configfile $input_config 


mv LRRtransfer*.err jobinfo/
mv LRRtransfer*.out jobinfo/
