#!/bin/sh
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task 32
#SBATCH --time=7-00:00:00
#SBATCH --partition=walkerpass
#SBATCH --job-name=zoonosis-risk-ai-modelling
#SBATCH --error=S3/batch_scripts/full_analysis.%J.stderr
#SBATCH --output=S3/batch_scripts/full_analysis.%J.stdout


source /data/projects/zoonosis-risk-ai/miniconda3/etc/profile.d/conda.sh
conda activate zoonosis-risk-ml 


Rscript S3/scripts/runall.R