#!/bin/bash
#SBATCH --job-name=vert_thawing
#SBATCH --output=vert_thaw_%A_%a.out
#SBATCH --error=vert_thaw_%A_%a.err
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --time=9:00:00
#SBATCH --mem=140G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu-invest
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker

module purge
module load R/4.4.2-gfbf-2024a
module load GDAL/3.10.0-foss-2024a

SSPS=(126 245 370 585)
SSP=${SSPS[$SLURM_ARRAY_TASK_ID-1]}

echo "Running SSP ${SSP}"

Rscript cum_thawed_N_by_depth.R "${SSP}"