#!/bin/bash
#SBATCH --job-name=n_mineralisation
#SBATCH --output=n_mineralisation_%j.out
#SBATCH --array=1-4
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --error=n_mineralisation_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --partition=cpu-invest

# Investment account
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker

module load R/4.4.2-gfbf-2024a
module load GDAL/3.10.0-foss-2024a

# Define SSPs
SSPS=(126 245 370 585)

# Map array task ID to SSP
SSP=${SSPS[$SLURM_ARRAY_TASK_ID-1]}

echo "Running SSP = $SSP (task $SLURM_ARRAY_TASK_ID)"

# Run R script, pass SSP as argument
Rscript mineralisation_new.R $SSP