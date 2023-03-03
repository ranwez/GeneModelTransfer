rm -rf .snakemake/log .snakemake/locks
module load snakemake/7.15.1-conda
profile="./SNAKEMAKE/PROFILES/slurm_muse/"
cluster_config="./SNAKEMAKE/slurm_config.yaml"
config_file="./SNAKEMAKE/config.yaml"
smk_Pipeline="./SNAKEMAKE/LRR_transfer.smk"
snakemake  --profile $profile --cluster-config $cluster_config --snakefile $smk_Pipeline --configfile $config_file
#sbatch --partition=agap_long --wrap "./launch_transfert.sh"
