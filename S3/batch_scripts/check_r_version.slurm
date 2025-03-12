#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=10-00:00:00
#SBATCH --partition=standard-wp
#SBATCH --job-name=zoonosis-risk-ai-modelling
#SBATCH --error=S3/batch_scripts/full_analysis.%J.stderr
#SBATCH --output=S3/batch_scripts/full_analysis.%J.stdout


source /hpscol02/tenant1/zoonosis-risk-ai/miniconda3/etc/profile.d/conda.sh
conda activate zoonosis-risk-ml

Rscript --version