#!/bin/sh
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=3-00:00:00
#SBATCH --partition=walkerpass
#SBATCH --job-name=zoonosis-risk-ai-modelling
#SBATCH --error=batch_scripts/outputs/full_analysis.%J.stderr
#SBATCH --output=batch_scripts/outputs/full_analysis.%J.stdout


eval "$(conda shell.bash hook)"
conda activate zoonosis-risk-ml 


Rscript S3/scripts/runall.R
