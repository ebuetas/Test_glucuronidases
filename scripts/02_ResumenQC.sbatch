#!/bin/bash

#SBATCH --job-name=QC_summary
#SBATCH --output=/data/home/molivares/elena/glucuronidases/logs/QC_summary_%j.out
#SBATCH --error=/data/home/molivares/elena/glucuronidases/logs/QC_summary_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1

set -euo pipefail



Rscript /data/home/molivares/elena/glucuronidases/scripts/R/summary_QC.R
