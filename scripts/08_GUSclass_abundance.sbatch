#!/bin/bash
#SBATCH --job-name=GUS_relative
#SBATCH --output=/data/home/molivares/elena/glucuronidases/logs/RELATIVE_%j.out
#SBATCH --error=/data/home/molivares/elena/glucuronidases/logs/RELATIVE_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

PROJECT_DIR="/data/home/molivares/elena/glucuronidases"

COUNTS="${PROJECT_DIR}/results/GUS_class_x_sample_counts.tsv"

QC="${PROJECT_DIR}/results/QC_summary.tsv"

# Ajusta esta ruta si clinical_metadata.tsv está en otro directorio
CLINICAL="${PROJECT_DIR}/data/clinical_metadata.tsv"

SCRIPT="${PROJECT_DIR}/scripts/R/relative_abundance_GUS.R"

OUTPUT_DIR="${PROJECT_DIR}/statistical_analysis/relative_abundance"

mkdir -p "${OUTPUT_DIR}"

module purge
module load R

echo "Starting GUS relative abundance analysis"

if [[ ! -s "${COUNTS}" ]]; then
    echo "ERROR: Missing counts matrix: ${COUNTS}"
    exit 1
fi

if [[ ! -s "${QC}" ]]; then
    echo "ERROR: Missing QC summary: ${QC}"
    exit 1
fi

if [[ ! -s "${CLINICAL}" ]]; then
    echo "ERROR: Missing clinical metadata: ${CLINICAL}"
    exit 1
fi

if [[ ! -s "${SCRIPT}" ]]; then
    echo "ERROR: Missing R script: ${SCRIPT}"
    exit 1
fi

Rscript "${SCRIPT}" \
    "${COUNTS}" \
    "${QC}" \
    "${CLINICAL}" \
    "${OUTPUT_DIR}"

echo "GUS relative abundance analysis completed"
