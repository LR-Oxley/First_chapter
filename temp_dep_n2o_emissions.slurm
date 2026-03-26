#!/bin/bash
#SBATCH --job-name=n2o_emissions
#SBATCH --output=n2o_emissions_%j.out
#SBATCH --nodes=1
#SBATCH --error=n2o_emissions_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=512G
#SBATCH --partition=cpu-invest

# Investment account
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker

module load R/4.4.2-gfbf-2024a
module load GDAL/3.10.0-foss-2024a

Rscript temp_dependent_n2o.R
