#!/bin/bash
#SBATCH --job-name=hic_map_whole
#SBATCH --output=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/02_map_to_whole_assembly/hic_map_whole_%j.out
#SBATCH --error=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/02_map_to_whole_assembly/hic_map_whole_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G

set -euo pipefail

# Logic: load conda so the hic_diatom environment can be activated inside the SLURM job.
source ~/miniforge3/etc/profile.d/conda.sh

# Logic: activate the environment containing BWA, samtools, and seqkit.
conda activate hic_diatom

# Logic: define the original assembly, paired Hi-C reads, and mapping output directory.
ASM_ORIG=/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_output/pypolca_corrected.fasta
R1=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/1574499_S6_L001_R1_001.fastq.gz
R2=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/1574499_S6_L001_R2_001.fastq.gz
OUT=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/02_map_to_whole_assembly

# Logic: create and enter the mapping directory where all index and BAM files will be written.
mkdir -p ${OUT}
cd ${OUT}

# Logic: create a local symlink to avoid writing BWA and samtools index files into the original assembly directory.
ASM_LOCAL=${OUT}/pypolca_corrected.hic_input.fasta
ln -sfn ${ASM_ORIG} ${ASM_LOCAL}

# Logic: remove old partial outputs so reruns start cleanly.
rm -f hic_to_whole_assembly.raw.bam \
      hic_to_whole_assembly.name_sorted.bam \
      hic_to_whole_assembly.coord_sorted.bam \
      hic_to_whole_assembly.coord_sorted.bam.bai \
      hic_to_whole_assembly.flagstat.txt \
      hic_to_whole_assembly.idxstats.txt \
      hic_to_whole_assembly.samtools_stats.txt

# Logic: index the assembly for BWA read mapping.
bwa index ${ASM_LOCAL}

# Logic: index the assembly with samtools to generate contig lengths.
samtools faidx ${ASM_LOCAL}

# Logic: save a two-column contig-size table for downstream contact-map or plotting tools.
cut -f1,2 ${ASM_LOCAL}.fai > pypolca_corrected.chrom.sizes

# Logic: map Hi-C paired-end reads to the assembly and convert the SAM stream to BAM.
bwa mem -SP5M -t 24 ${ASM_LOCAL} ${R1} ${R2} \
    | samtools view -@ 8 -bS - \
    > hic_to_whole_assembly.raw.bam

# Logic: make a name-sorted BAM for Hi-C scaffolding tools such as YaHS.
samtools sort -@ 24 -n \
    -o hic_to_whole_assembly.name_sorted.bam \
    hic_to_whole_assembly.raw.bam

# Logic: make a coordinate-sorted BAM for mapping statistics and contact extraction.
samtools sort -@ 24 \
    -o hic_to_whole_assembly.coord_sorted.bam \
    hic_to_whole_assembly.raw.bam

# Logic: index the coordinate-sorted BAM for fast access by samtools.
samtools index hic_to_whole_assembly.coord_sorted.bam

# Logic: summarize read-level mapping statistics.
samtools flagstat hic_to_whole_assembly.coord_sorted.bam \
    > hic_to_whole_assembly.flagstat.txt

# Logic: summarize mapped and unmapped reads per contig.
samtools idxstats hic_to_whole_assembly.coord_sorted.bam \
    > hic_to_whole_assembly.idxstats.txt

# Logic: save a detailed BAM statistics report.
samtools stats hic_to_whole_assembly.coord_sorted.bam \
    > hic_to_whole_assembly.samtools_stats.txt
