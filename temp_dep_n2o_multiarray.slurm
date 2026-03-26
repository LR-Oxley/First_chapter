#!/bin/bash
#SBATCH --job-name=n2o_emissions
#SBATCH --output=logs/n2o_emissions_%A_%a.out
#SBATCH --error=logs/n2o_emissions_%A_%a.err
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=512G

#SBATCH --partition=cpu-invest
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker

# ---- Job array: 4 SSPs ----
#SBATCH --array=0-3

module load R/4.4.2-gfbf-2024a
module load GDAL/3.10.0-foss-2024a

# Map array index -> SSP
SSPS=(126 245 370 585)
SSP=${SSPS[$SLURM_ARRAY_TASK_ID]}

Rscript temp_dep_n2o_multiarray.R ${SLURM_ARRAY_TASK_ID}
