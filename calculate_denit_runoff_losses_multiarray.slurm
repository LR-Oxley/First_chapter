#!/bin/bash
#SBATCH --job-name=compute_losses_denit_runoff_loss
#SBATCH --output=compute_losses_denit_runoff_loss_%A_%a.out
#SBATCH --error=compute_losses_denit_runoff_loss_%A_%a.err
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=512G
#SBATCH --partition=cpu-invest

# Investment account
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker

# ---- Job array: 4 SSPs ----
#SBATCH --array=0-3

module load R/4.4.2-gfbf-2024a
module load GDAL/3.10.0-foss-2024a

# Map array index -> SSP
SSPS=(126 245 370 585)
SSP=${SSPS[$SLURM_ARRAY_TASK_ID]}

echo "Running compute_losses_denit_runoff_loss for SSP${SSP}"

# Run R script with SSP argument
Rscript compute_losses_denit_runoff_multiarray_3m_lim.R ${SSP}