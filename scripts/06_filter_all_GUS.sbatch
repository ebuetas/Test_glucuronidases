#!/bin/bash
#SBATCH --job-name=GUS_filter
#SBATCH --output=/data/home/molivares/elena/glucuronidases/logs/FILTER_%A_%a.out
#SBATCH --error=/data/home/molivares/elena/glucuronidases/logs/FILTER_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-48%4

set -euo pipefail

PROJECT_DIR="/data/home/molivares/elena/glucuronidases"
METADATA="${PROJECT_DIR}/data/samples.tsv"
RESULTS="${PROJECT_DIR}/results/GUS_diamond"
SCRIPT="${PROJECT_DIR}/scripts/R/filter_GUS.R"

# Cargar R
module purge
module load R

# Comprobar archivos necesarios
if [[ ! -s "${METADATA}" ]]; then
    echo "ERROR: No existe el metadata: ${METADATA}"
    exit 1
fi

if [[ ! -s "${SCRIPT}" ]]; then
    echo "ERROR: No existe el script R: ${SCRIPT}"
    exit 1
fi

# Obtener identificador de muestra según el índice del array
SAMPLE_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${METADATA}")

if [[ -z "${SAMPLE_LINE}" ]]; then
    echo "ERROR: No se ha encontrado la muestra ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

IFS=$'\t' read -r SAMPLE ORIGINAL_FASTQ <<< "${SAMPLE_LINE}"
SAMPLE="${SAMPLE//$'\r'/}"

INPUT_TSV="${RESULTS}/${SAMPLE}.GUS.tsv"
OUTPUT_TSV="${RESULTS}/${SAMPLE}.GUS_filtered_best.tsv"

echo "Muestra: ${SAMPLE}"
echo "Entrada: ${INPUT_TSV}"
echo "Salida esperada: ${OUTPUT_TSV}"

if [[ ! -s "${INPUT_TSV}" ]]; then
    echo "ERROR: No existe el TSV de entrada o está vacío: ${INPUT_TSV}"
    exit 1
fi

Rscript "${SCRIPT}" "${INPUT_TSV}"

if [[ ! -f "${OUTPUT_TSV}" ]]; then
    echo "ERROR: No se ha generado el archivo esperado: ${OUTPUT_TSV}"
    exit 1
fi

echo "Filtrado completado para ${SAMPLE}"
