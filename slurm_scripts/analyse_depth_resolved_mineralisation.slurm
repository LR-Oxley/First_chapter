#!/bin/bash
#SBATCH --job-name=arctic_min_depth
#SBATCH --output=arctic_min_depth_%A_%a.out
#SBATCH --error=arctic_min_depth_%A_%a.err
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=254G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu-invest

# Investment account
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker

module load R/4.4.2-gfbf-2024a
module load GDAL/3.10.0-foss-2024a

# Map array index to SSP
SSPS=(126 245 370 585)
SSP=${SSPS[$SLURM_ARRAY_TASK_ID-1]}

echo "Running SSP $SSP (task $SLURM_ARRAY_TASK_ID)"

Rscript analyse_depth_resolved_mineralisation.R $SSP