#!/bin/bash
#SBATCH --job-name=temp_sens_analysis
#SBATCH --output=temp_sens_analysis_%j.out
#SBATCH --nodes=1
#SBATCH --error=temp_sens_analysis_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
# Investment account
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker

module load R/4.4.2-gfbf-2024a
module load GDAL/3.10.0-foss-2024a

Rscript temp_sensitivity_analysis.R 

