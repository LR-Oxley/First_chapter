#!/bin/bash
#SBATCH --job-name=compute_losses_denit_runoff_loss
#SBATCH --output=compute_losses_denit_runoff_loss_%j.out
#SBATCH --nodes=1
#SBATCH --error=compute_losses_denit_runoff_loss_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=512G
#SBATCH --partition=cpu-invest

# Investment account
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker

module load R/4.4.2-gfbf-2024a
module load GDAL/3.10.0-foss-2024a

Rscript compute_losses_denit_runoff_loss.R
