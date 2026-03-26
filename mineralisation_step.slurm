#!/bin/bash
#SBATCH --job-name=n_mineralisation
#SBATCH --output=n_mineralisation_%j.out
#SBATCH --nodes=1
#SBATCH --error=n_mineralisation_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --partition=cpu-invest

# Investment account
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker

module load R/4.4.2-gfbf-2024a
module load GDAL/3.10.0-foss-2024a

Rscript mineralisation_step.R
