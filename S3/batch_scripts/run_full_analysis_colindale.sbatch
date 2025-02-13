#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:30:00
#SBATCH --partition=standard-wp
#SBATCH --job-name=zoonosis-risk-ai-modelling
#SBATCH --error=S3/batch_scripts/full_analysis.%J.stderr
#SBATCH --output=S3/batch_scripts/full_analysis.%J.stdout


source /hpscol02/tenant1/zoonosis-risk-ai/miniconda3/etc/profile.d/conda.sh
conda activate zoonosis-risk-ml

Rscript S3/scripts/runall.R
