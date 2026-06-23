#!/bin/bash
#SBATCH --job-name=hic_yahs_whole
#SBATCH --output=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/03_yahs_scaffolding/yahs_whole_%j.out
#SBATCH --error=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/03_yahs_scaffolding/yahs_whole_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G

set -euo pipefail

# Logic: load conda so the hic_diatom environment can be activated inside the SLURM job.
source ~/miniforge3/etc/profile.d/conda.sh

# Logic: activate the environment containing YaHS and seqkit.
conda activate hic_diatom

# Logic: define the local assembly link, name-sorted Hi-C BAM, and YaHS output directory.
ASM=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/02_map_to_whole_assembly/pypolca_corrected.hic_input.fasta
BAM=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/02_map_to_whole_assembly/hic_to_whole_assembly.name_sorted.bam
OUT=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/03_yahs_scaffolding

# Logic: create and enter the YaHS output directory.
mkdir -p ${OUT}
cd ${OUT}

# Logic: run exploratory whole-assembly Hi-C scaffolding.
yahs -o DL_diatom_whole_hic_yahs ${ASM} ${BAM}

# Logic: compare input and YaHS output assembly statistics.
seqkit stats \
    ${ASM} \
    DL_diatom_whole_hic_yahs_scaffolds_final.fa \
    > DL_diatom_whole_hic_yahs.seqkit_stats.txt
