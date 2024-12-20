#!/bin/sh
#SBATCH -n 96
#SBATCH --ntasks-per-node=20
#SBATCH --time=7-00:00:00
#SBATCH --partition=walkerpass
#SBATCH --job-name=zoonosis-risk-ai-modelling
#SBATCH --error=S3/batch_scripts/full_analysis.%J.stderr
#SBATCH --output=S3/batch_scripts/full_analysis.%J.stdout


source /hpspor01/tenant1/projects/zoonosis-risk-ai/miniconda3/etc/profile.d/conda.sh
conda activate zoonosis-risk-ml


mpirun -np 1 Rscript -e "source('S3/scripts/runall.R')"