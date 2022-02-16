#!/bin/bash
#

snakePipe=$1



snakemake  --use-singularity --latency-wait 40 --snakefile $snakePipe --configfile config.yml --cluster "qsub -V -cwd -b y" --jobs 200 --printshellcmds --use-envmodules


