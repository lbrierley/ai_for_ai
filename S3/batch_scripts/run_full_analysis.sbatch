#!/bin/sh
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=2
#SBATCH --time=12:00:00
#SBATCH --partition=walkerpass
#SBATCH --job-name=zoonosis-risk-ai-modelling
#SBATCH --error=S3/batch_scripts/full_analysis.%J.stderr
#SBATCH --output=S3/batch_scripts/full_analysis.%J.stdout


source ~/miniconda3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate zoonosis-risk-ml 


Rscript S3/scripts/runall.R
