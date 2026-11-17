#!/bin/bash
#SBATCH --job-name=GUS_counts
#SBATCH --output=/data/home/molivares/elena/glucuronidases/logs/COUNTS_%A.out
#SBATCH --error=/data/home/molivares/elena/glucuronidases/logs/COUNTS_%A.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

PROJECT_DIR="/data/home/molivares/elena/glucuronidases"

METADATA="${PROJECT_DIR}/data/samples.tsv"
RESULTS="${PROJECT_DIR}/results/GUS_diamond"
OUTPUT_DIR="${PROJECT_DIR}/results"

SCRIPT="${PROJECT_DIR}/scripts/R/GUS_matriz.R"

mkdir -p "${OUTPUT_DIR}"

module purge
module load R

if [[ ! -s "${SCRIPT}" ]]; then
    echo "ERROR: No existe el script R: ${SCRIPT}"
    exit 1
fi

echo "Iniciando recuento de GUS y clases..."

Rscript "${SCRIPT}" \
    "${METADATA}" \
    "${RESULTS}" \
    "${OUTPUT_DIR}"

echo "Recuento completado."
