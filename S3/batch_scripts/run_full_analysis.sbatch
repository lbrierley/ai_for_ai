#!/bin/sh
#SBATCH --ntasks-per-node=20
#SBATCH --nodes=1
#SBATCH --time=3-00:00:00
#SBATCH --partition=walkerpass
#SBATCH --job-name=zoonosis-risk-ai-modelling
#SBATCH --error=S3/batch_scripts/full_analysis.%J.stderr
#SBATCH --output=S3/batch_scripts/full_analysis.%J.stdout


eval "$(conda shell.bash hook)"
conda activate zoonosis-risk-ml 


Rscript S3/scripts/runall.R
