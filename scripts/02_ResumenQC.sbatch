#!/bin/bash
#SBATCH --job-name=QC_summary
#SBATCH --output=QC_summary_%j.out
#SBATCH --error=QC_summary_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1

set -euo pipefail

module purge
module load R

Rscript scripts/R/summary_QC.R
