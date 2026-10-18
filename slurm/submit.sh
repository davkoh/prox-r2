#!/bin/bash
#SBATCH --job-name=prox-r2
#SBATCH --time=25:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --output=snakemake-%j.out

# Make R, CmdStan and Snakemake available here, for example with "module load".
snakemake --profile slurm
