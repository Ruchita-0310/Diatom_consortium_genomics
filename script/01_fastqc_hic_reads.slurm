#!/bin/bash
#SBATCH --job-name=hic_fastqc
#SBATCH --output=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/01_qc/hic_fastqc_%j.out
#SBATCH --error=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/01_qc/hic_fastqc_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G

set -euo pipefail

# Logic: load conda so the hic_diatom environment can be activated inside the SLURM job.
source ~/miniforge3/etc/profile.d/conda.sh

# Logic: activate the environment containing FastQC and MultiQC.
conda activate hic_diatom

# Logic: define the Hi-C working directory and QC output folder.
HIC_DIR=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms
OUT=${HIC_DIR}/01_qc

# Logic: create the output directory if it does not already exist.
mkdir -p ${OUT}

# Logic: run FastQC on both Hi-C read files to assess raw read quality.
fastqc -t 8 \
    ${HIC_DIR}/1574499_S6_L001_R1_001.fastq.gz \
    ${HIC_DIR}/1574499_S6_L001_R2_001.fastq.gz \
    -o ${OUT}

# Logic: combine FastQC reports into one MultiQC summary.
multiqc ${OUT} -o ${OUT}
