#!/bin/bash
#SBATCH --job-name=GUS_makedb
#SBATCH --output=/data/home/molivares/elena/glucuronidases/logs/GUS_makedb_%j.out
#SBATCH --error=/data/home/molivares/elena/glucuronidases/logs/GUS_makedb_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2


set -euo pipefail

PROJECT_DIR="/data/home/molivares/elena/glucuronidases"
GUS_FASTA="${PROJECT_DIR}/data/HMGC279.fasta"
DB_DIR="${PROJECT_DIR}/data/references/diamond_db"
DB_PREFIX="${DB_DIR}/HMGC279"

mkdir -p "${DB_DIR}"

module purge
module load DIAMOND/2.1.11-GCC-13.3.0

echo "DIAMOND version:"
diamond version

if [[ ! -s "${GUS_FASTA}" ]]; then
    echo "ERROR: No encuentro el FASTA de proteínas: ${GUS_FASTA}"
    exit 1
fi

echo "Building GUS protein database..."

diamond makedb \
    --in "${GUS_FASTA}" \
    --db "${DB_PREFIX}" \
    --threads "${SLURM_CPUS_PER_TASK}"

if [[ ! -s "${DB_PREFIX}.dmnd" ]]; then
    echo "ERROR: No se ha generado la base de datos."
    exit 1
fi

echo "Database successfully created:"
ls -lh "${DB_PREFIX}.dmnd"

echo "Checking database information..."
diamond dbinfo --db "${DB_PREFIX}"

echo "Finished."
