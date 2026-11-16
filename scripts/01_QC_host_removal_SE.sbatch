#!/bin/bash
#SBATCH --job-name=QC_HostRemoval
#SBATCH --output=/data/home/molivares/elena/glucuronidases/logs/QC_%A_%a.out
#SBATCH --error=/data/home/molivares/elena/glucuronidases/logs/QC_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8

#SBATCH --array=1-48%4

###############################################################################
# Control de calidad y eliminación de ADN humano
# FASTQ SINGLE-END
#
# Herramientas:
#   FastQC
#   Trimmomatic
#   Bowtie2
#
###############################################################################

set -euo pipefail
module purge

############################
# 1. CONFIGURACIÓN
############################

PROJECT_DIR="/data/home/molivares/elena/glucuronidases/"

METADATA="${PROJECT_DIR}/data/samples.tsv"

# Índice humano construido en el trabajo anterior
HUMAN_INDEX="/data/home/molivares/references/human/GRCh38_Ensembl116/bowtie2_index/GRCh38"

RESULTS="${PROJECT_DIR}/data/processed"

QC_RAW="${RESULTS}/fastqc_raw"
TRIMMED="${RESULTS}/trimmed"
QC_TRIMMED="${RESULTS}/fastqc_trimmed"
NON_HUMAN="${RESULTS}/non_human"
QC_FINAL="${RESULTS}/fastqc_nonhuman"
REPORTS="${RESULTS}/reports"
LOG_DIR="${PROJECT_DIR}/logs"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

############################
# 2. DIRECTORIOS
############################

mkdir -p "$QC_RAW" "$TRIMMED" "$QC_TRIMMED" \
         "$NON_HUMAN" "$QC_FINAL" "$REPORTS"

############################
# 3. CARGAR MÓDULOS
############################


# Adaptar a los nombres y versiones de tu HPC
#module avail 2>&1 | grep -Ei 'fastp|trimmomatic|bbmap|bbtools|cutadapt|kneaddata|sickle'
module load FastQC
module load Trimmomatic/0.40-Java-17
module load Bowtie2/2.5.4-GCC-13.2.0

# ============================================================
# LOCALIZAR TRIMMOMATIC
# ============================================================
TRIMMOMATIC_DIR="${EBROOTTRIMMOMATIC}"

TRIMMOMATIC_JAR=$(find "${TRIMMOMATIC_DIR}" \
    -maxdepth 2 \
    -type f \
    -name "*.jar" \
    | head -n 1)

if [[ -z "${TRIMMOMATIC_JAR}" || ! -f "${TRIMMOMATIC_JAR}" ]]; then
    echo "ERROR: No se encuentra el JAR de Trimmomatic"
    exit 1
fi

echo "Trimmomatic JAR: ${TRIMMOMATIC_JAR}"

# ============================================================
# LEER INFORMACION DE LA MUESTRA
# ============================================================

SAMPLE_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${METADATA}")

if [[ -z "${SAMPLE_LINE}" ]]; then
    echo "ERROR: No existe la muestra ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

IFS=$'\t' read -r SAMPLE FASTQ <<< "${SAMPLE_LINE}"

if [[ ! -f "${FASTQ}" ]]; then
    echo "ERROR: No encuentro el FASTQ: ${FASTQ}"
    exit 1
fi

echo "========================================"
echo "Muestra: ${SAMPLE}"
echo "FASTQ: ${FASTQ}"
echo "Fecha: $(date)"
echo "========================================"

# ============================================================
# ARCHIVOS DE SALIDA
# ============================================================

TRIMMED1="${TRIMMED}/${SAMPLE}.trimmed.fastq.gz"
NONHUMAN1="${NON_HUMAN}/${SAMPLE}.nonhuman.fastq.gz"

TRIM_LOG="${LOG_DIR}/${SAMPLE}_trimmomatic.log"
BOWTIE_LOG="${LOG_DIR}/${SAMPLE}_bowtie2.log"

# ============================================================
# 1. FASTQC INICIAL
# ============================================================

echo "1. FastQC inicial"

fastqc \
    --threads 2 \
    --outdir "${QC_RAW}" \
    "${FASTQ}"

# ============================================================
# 2. TRIMMOMATIC
# ============================================================

echo "2. Trimmomatic"

java -jar "${TRIMMOMATIC_JAR}" SE \
    -threads "${SLURM_CPUS_PER_TASK}" \
    -phred33 \
    "${FASTQ}" \
    "${TRIMMED1}" \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:50 \
    2>&1 | tee "${TRIM_LOG}"

# ============================================================
# 3. FASTQC DESPUES DEL FILTRADO
# ============================================================

echo "3. FastQC posterior a Trimmomatic"

fastqc \
    --threads 2 \
    --outdir "${QC_TRIMMED}" \
    "${TRIMMED1}"

# ============================================================
# 4. ELIMINACION DE LECTURAS HUMANAS
# ============================================================

echo "4. Eliminacion de lecturas humanas con Bowtie2"

bowtie2 \
    --very-sensitive \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -x "${HUMAN_INDEX}" \
    -U "${TRIMMED1}" \
    --un-gz "${NONHUMAN1}" \
    -S /dev/null \
    2> "${BOWTIE_LOG}"

# ============================================================
# 5. FASTQC DE LECTURAS NO HUMANAS
# ============================================================

echo "5. FastQC de lecturas no humanas"

fastqc \
    --threads 2 \
    --outdir "${QC_FINAL}" \
    "${NONHUMAN1}"

# ============================================================
# FINAL
# ============================================================

echo "========================================"
echo "Muestra completada: ${SAMPLE}"
echo "FASTQ filtrado: ${TRIMMED1}"
echo "FASTQ no humano: ${NONHUMAN1}"
echo "========================================"
