#!/bin/bash
#SBATCH -J SNP_call_attempt2
#SBATCH -A ACF-UTK0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=6-00:00:00
#SBATCH --error=logs/job.SNP_call.e%J
#SBATCH --output=logs/job.SNP_call.o%J
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mquigg1@vols.utk.edu

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate vary_cool

export NXF_OPTS="-Xms500M -Xmx2G"
export NXF_ANSI_LOG=false

nextflow /lustre/isaac24/scratch/mquigg1/tails_GWAS/00.vary_cool/vary_cool/main.nf \
    --publish_dir /lustre/isaac24/scratch/mquigg1/tails_GWAS/01.vary_cool_output \
    --input /lustre/isaac24/scratch/mquigg1/tails_GWAS/metadata.csv \
    --genome /lustre/isaac24/scratch/mquigg1/tails_GWAS/00.references/pe57_v-T-B.H.C.C.A.A.FINAL.hap1.fasta \
    --skip_qc false \
    --skip_trim false \
    --skip_mark_dupe false \
    --aligner bwa_mem \
    --ploidy 2 \
    --chunks 12 \
    --caller bcftools \
    --bp_intervals 10000000 \
    -profile slurm,custom \
    -resume
