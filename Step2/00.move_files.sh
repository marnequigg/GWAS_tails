#!/bin/bash
#SBATCH --job-name=move_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=32G
#SBATCH -A ACF-UTK0011
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --output=logs/cp_fqgz_%j.out
#SBATCH --error=logs/cp_fqgz_%j.err
#SBATCH --time=2:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mquigg1@vols.utk.edu

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# ─── USER CONFIGURATION ──────────────────────────────────────────────────────
SOURCE_DIR="/lustre/isaac24/proj/UTK0032/TSIP_ash/NRS_data_master"
DEST_PATH="/lustre/isaac24/scratch/mquigg1/tails_GWAS/00.input_data"
PARALLEL_JOBS=15
# ─────────────────────────────────────────────────────────────────────────────

# Activate conda environment properly in a non-interactive shell
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate basics

# Create log dir if it doesn't exist
mkdir -p logs

echo "Job started: $(date)"
echo "Searching for .fq.gz files under: ${SOURCE_DIR}"

# Find all .fq.gz files recursively (following symlinks with -L),
# then copy each file in parallel, dereferencing symlinks with rsync -L
find -L "${SOURCE_DIR}" -type f -name "*.fq.gz" | \
    parallel --jobs "${PARALLEL_JOBS}" \
             --plain \
             --joblog logs/parallel_transfer_log_${SLURM_JOB_ID}.txt \
             --resume \
    "rsync -aL {} ${DEST_PATH}/"

EXIT_CODE=$?

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "All copies completed successfully."
else
    echo "One or more copies failed. Check logs/parallel_transfer_log_${SLURM_JOB_ID}.txt for details."
fi

echo "Job finished: $(date)"
exit ${EXIT_CODE}
