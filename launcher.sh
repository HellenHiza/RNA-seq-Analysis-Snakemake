#!/bin/bash
#SBATCH --job-name tbtz
#SBATCH --qos short
#SBATCH --time 1-00:00:00
#SBATCH --cpus-per-task 16
#SBATCH --mem 32GB
#SBATCH --output slurm-%x-%j.out
#SBATCH --mail-type ALL
#SBATCH --mail-user hehi@uv.es

#conda activate snakemake9
snakemake -c 16 --use-conda --rerun-incomplete  # --use-conda only if u use conda: "...yaml"

# $ sbatch launcher.sh
