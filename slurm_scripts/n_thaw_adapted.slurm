#!/bin/bash
#SBATCH --job-name=arctic_N_thaw
#SBATCH --output=N_thaw_%A_%a.out
#SBATCH --error=N_thaw_%A_%a.err
#SBATCH --array=0-11%4
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=85G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu-invest

# Investment account
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker

module load R/4.4.2-gfbf-2024a
module load GDAL/3.10.0-foss-2024a

# Map array index to SSP and run: 4 SSPs x 3 runs = 12 tasks
#   0 = 126 mean, 1 = 126 plus, 2 = 126 minus, 3 = 245 mean, ...
SSPS=(126 245 370 585)
RUNS=(mean plus minus)

SSP=${SSPS[$(( SLURM_ARRAY_TASK_ID / 3 ))]}
RUN=${RUNS[$(( SLURM_ARRAY_TASK_ID % 3 ))]}

echo "Running N thaw for SSP $SSP, run $RUN (task $SLURM_ARRAY_TASK_ID)"

Rscript n_thaw_adapted.R $SSP $RUN