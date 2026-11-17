#!/bin/bash
#SBATCH --job-name=GUS_annot
#SBATCH --output=/data/home/molivares/elena/glucuronidases/logs/ANNOT_%A_%a.out
#SBATCH --error=/data/home/molivares/elena/glucuronidases/logs/ANNOT_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-48%4

set -euo pipefail

# --------------------------------------------------
# DIRECTORIES
# --------------------------------------------------

PROJECT_DIR="/data/home/molivares/elena/glucuronidases"

METADATA="${PROJECT_DIR}/data/samples.tsv"

RESULTS="${PROJECT_DIR}/results/GUS_diamond"

CLASSES="${PROJECT_DIR}/data/GUS_metadata.tsv"

SCRIPT="${PROJECT_DIR}/scripts/R/add_GUS_annot.R"

LOG_DIR="${PROJECT_DIR}/logs"


# --------------------------------------------------
# LOAD R
# --------------------------------------------------

module purge
module load R

echo "=========================================="
echo "GUS ANNOTATION"
echo "Task: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "=========================================="

# --------------------------------------------------
# READ SAMPLE
# --------------------------------------------------

SAMPLE_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${METADATA}")

if [[ -z "${SAMPLE_LINE}" ]]; then
    echo "ERROR: Sample not found in samples.tsv"
    exit 1
fi

IFS=$'\t' read -r SAMPLE ORIGINAL_FASTQ <<< "${SAMPLE_LINE}"

SAMPLE="${SAMPLE//$'\r'/}"

echo "Sample: ${SAMPLE}"

# --------------------------------------------------
# CHECK INPUT FILES
# --------------------------------------------------

INPUT_TSV="${RESULTS}/${SAMPLE}.GUS.tsv"

if [[ ! -s "${INPUT_TSV}" ]]; then
    echo "ERROR: DIAMOND result not found: ${INPUT_TSV}"
    exit 1
fi

if [[ ! -s "${CLASSES}" ]]; then
    echo "ERROR: GUS classes file not found: ${CLASSES}"
    exit 1
fi

if [[ ! -s "${SCRIPT}" ]]; then
    echo "ERROR: R script not found: ${SCRIPT}"
    exit 1
fi

# --------------------------------------------------
# RUN ANNOTATION
# --------------------------------------------------

echo "Annotating ${SAMPLE}..."

Rscript "${SCRIPT}" \
    "${INPUT_TSV}" \
    "${CLASSES}"

echo "=========================================="
echo "Finished sample: ${SAMPLE}"
echo "Date: $(date)"
echo "=========================================="
