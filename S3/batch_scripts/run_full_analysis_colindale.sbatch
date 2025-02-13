#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=7-00:00:00
#SBATCH --partition=standard-wp
#SBATCH --job-name=zoonosis-risk-ai-modelling
#SBATCH --error=S3/batch_scripts/full_analysis.%J.stderr
#SBATCH --output=S3/batch_scripts/full_analysis.%J.stdout


source /hpspor01/tenant1/projects/zoonosis-risk-ai/miniconda3/etc/profile.d/conda.sh
conda activate zoonosis-risk-ml


echo =========================================================   
echo SLURM job: submitted  date = `date`      
date_start=`date +%s`

Rscript S3/scripts/runall.R

echo   
echo ---------------                                           
echo Job output ends                                           
date_end=`date +%s`
seconds=$((date_end-date_start))
minutes=$((seconds/60))
seconds=$((seconds-60*minutes))
hours=$((minutes/60))
minutes=$((minutes-60*hours))
echo =========================================================   
echo SLURM job: finished   date = `date`   
echo Total run time : $hours Hours $minutes Minutes $seconds Seconds
echo =========================================================  
exit