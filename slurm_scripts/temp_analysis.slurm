#!/bin/bash
#SBATCH --job-name=60N_temp_analysis
#SBATCH --partition=cpu-invest
#SBATCH --qos=job_icpu-stocker
#SBATCH --account=invest
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=160G
#SBATCH --time=08:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/temp_analysis_%A_%a.out
#SBATCH --error=logs/temp_analysis_%A_%a.err

mkdir -p logs

module load R/4.4.2-gfbf-2024a
module load GDAL/3.10.0-foss-2024a

SSPS=(126 245 370 585)
SSP=${SSPS[$SLURM_ARRAY_TASK_ID]}

echo "Running temperature analysis for SSP ${SSP}"

Rscript temp_analysis.R ${SSP}

echo "Finished SSP ${SSP}"