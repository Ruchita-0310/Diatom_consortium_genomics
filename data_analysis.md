# 1. Assembly
Basecalled using Guppy
```
flye --nano-raw pass_trim.fastq.gz --meta -g 50m --min-overlap 5000 --out-dir flye_out_new -i 3 --threads 8
```
# 2. Mapping
## Mapping Short Reads (SR) to the assembled reads
```
minimap2 -ax sr /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/guppy_flye_assembly.fasta Diatoms_merged.fastq.gz > sr_alignment.sam
samtools view -S -b sr_alignment.sam > alignment.bam
samtools sort alignment.bam -o alignment_sorted.bam
samtools index alignment_sorted.bam
samtools flagstat alignment_sorted.bam > mapping_stats.txt
samtools idxstats alignment_sorted.bam | sort -k3,3rn > sr_all_nanopore_hits.tsv
samtools depth alignment_sorted.bam > sr_depth.txt
awk '{sum[$1]+=$3; count[$1]++} END {for (c in sum) print c, sum[c]/count[c]}' sr_depth.txt | sort -k2,2nr > sr_mean_depth.tsv
```
# 3. Polishing
## Medaka - LR polising
```
medaka_consensus \
  -i pass_trim.fastq.gz \
  -d guppy_flye_assembly.fasta \
  -o medaka_euk_polished \
  -t 12
```
# Polishing - SR
```
## Mapping ###
bwa mem -t 16 -a /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/medaka_euk_polished/consensus.fasta /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz > alignments_1.sam
bwa mem -t 16 -a /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/medaka_euk_polished/consensus.fasta /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz > alignments_2.sam

### Polypolish filter ###
polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam

### Polypolish filter ###
polypolish polish \
  /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/medaka_euk_polished/consensus.fasta \
  filtered_1.sam filtered_2.sam \
  > sr_poly.fasta

### Pypolca ###
pypolca run -a sr_poly.fasta \
  -1 /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz \
  -2 /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz \
  -t 12 -o sr_pypolca_output --careful
```
Run BUSCO on pyloca_corrected.fasta
```
busco -i pyloca_corrected.fasta -l busco_downloads/lineages/eukaryota_odb10 -o busco_report -m genome
```
# 4. Metabat2
```
conda create -n metabat2_v2 -c conda-forge -c bioconda metabat2 libdeflate=1.10

minimap2 -ax map-ont -t 16 1_sr_pypolca_output/pypolca_corrected.fasta pass_trim.fastq.gz | \
samtools view -@ 16 -bS - | \
samtools sort -@ 16 -m 10G -o aligned_reads.sorted.bam

samtools index -@ 16 aligned_reads.sorted.bam

samtools index -@ 16 aligned_reads.sorted.bam

# Summarize Depth
jgi_summarize_bam_contig_depths --outputDepth depth.txt --percentIdentity 85 aligned_reads.sorted.bam

# Binning
mkdir -p 2_metabat2_bins
metabat2 -i 1_sr_pypolca_output/pypolca_corrected.fasta -a depth.txt -o 2_metabat2_bins/bin -m 1500 -t 16 --unbinned
```
# 5. CheckM2
```
checkm2 predict --threads 16 --input 2_metabat2_bins/ --output_directory 3_checkm2_results
```
# 6. GTDB classification
```
# Assign taxonomy using the Genome Taxonomy Database
gtdbtk classify_wf --genome_dir 2_metabat2_bins/ --out_dir 4_gtdbtk_output --cpus 16 -x fa
```
# 7. CoverM
```
conda create -n coverm -c bioconda -c conda-forge coverm
conda activate coverm

coverm genome \
    --genome-fasta-directory 2_metabat2_bins/bac_bins \
    --genome-fasta-extension fa \
    -1 /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz \
    -2 /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz \
    --mapper bwa-mem \
    -m mean relative_abundance covered_fraction \
    --threads 8 \
    --min-read-percent-identity 95 \
    -o bac_output_coverm.tsv
```
# 8. MetaEUK combined with polished assembly
```
# 1. Load the module
module load blast/2.10.0

# 2. Run TBLASTN (Query: Protein | DB: DNA)
# This uses the DB you already made with makeblastdb
tblastn -query /work/ebg_lab/eb/diatom_consortia/sr_diatoms/metaeuk_result/sr_contigs_metaeuk_output.fna \
-db /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/sr_pypolca_corrected.fasta \
-outfmt 6 -max_target_seqs 1 -num_threads 16 -out bridge_results.txt

# 3. Filter for High-Confidence Hits (E-value < 1e-10)
awk '$5 < 1e-10 {print $2}' bridge_results.txt | sort | uniq > euk_contig_names_final.txt

# 4. Extract the Final Eukaryotic Bin
# Ensure seqtk is available (usually in your miniconda path)
seqtk subseq /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/sr_pypolca_corrected.fasta euk_contig_names_final.txt > Diatom_Euk_Bin.fasta
```
