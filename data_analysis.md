# 1. Assembly
Basecalled using Guppy
```
flye --nano-raw pass_trim.fastq.gz --meta -g 50m --min-overlap 5000 --out-dir flye_out_new -i 3 --threads 8
```
# 2. Mapping
## Mapping Hi-C reads to the assembled reads
```
minimap2 -ax sr /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/guppy_flye_assembly.fasta Diatoms_hi-c_merged.fastq.gz > hi-c_alignment.sam
samtools view -S -b hi-c_alignment.sam > hi-c_alignment.bam
samtools sort hi-c_alignment.bam -o hi-c_alignment_sorted.bam
samtools index hi-c_alignment_sorted.bam
samtools flagstat hi-c_alignment_sorted.bam > mapping_stats.txt
samtools idxstats hi-c_alignment_sorted.bam | sort -k3,3rn > hi-c_all_nanopore_hits.tsv
samtools depth hi-c_alignment_sorted.bam > depth.txt
awk '{sum[$1]+=$3; count[$1]++} END {for (c in sum) print c, sum[c]/count[c]}' depth.txt > mean_depth.tsv
```
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
# Polishing
## Medaka - LR polising
```
medaka_consensus \
  -i pass_trim.fastq.gz \
  -d guppy_flye_assembly.fasta \
  -o medaka_euk_polished \
  -t 12
```
## Polypolish - SR polishing
For SR and Hi-C reads - do it seperately        
```
# --- 1. Variables ---
THREADS=16
DRAFT="medaka_euk_polished/consensus.fasta"
R1="/work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz"
R2="/work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz"

# --- 2. Polypolish Step --- do this step only once as its already indexed
echo "Indexing draft..."
bwa index "$DRAFT"

echo "Aligning reads..."
# Aligning separately as required by Polypolish filter
bwa mem -t $THREADS -a "$DRAFT" "$R1" > sr_alignments_1.sam
bwa mem -t $THREADS -a "$DRAFT" "$R2" > sr_alignments_2.sam

echo "Filtering alignments..."
polypolish filter --in1 sr_alignments_1.sam --in2 sr_alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam

echo "Polishing..."
# Ensure the output filename is correctly captured
polypolish polish "$DRAFT" filtered_1.sam filtered_2.sam > sr_assembly.polypolish.fasta

# CHECK: Did the file actually get created?
if [ ! -s sr_assembly.polypolish.fasta ]; then
    echo "ERROR: Polypolish failed to create sr_assembly.polypolish.fasta"
    exit 1
fi

rm sr_alignments_1.sam sr_alignments_2.sam filtered_1.sam filtered_2.sam
```
## Pypolca - SR polishing
For SR and Hi-C reads - do it seperately        
```
# --- 3. pypolca Step ---
echo "Starting pypolca..."
pypolca run \
    -a sr_assembly.polypolish.fasta \
    -1 "$R1" \
    -2 "$R2" \
    -t $THREADS \
    --memory_limit 16G \
    -o sr_pypolca_output
```
# BUSCO
On Hi-C and SR 
```
busco -i pypolca_corrected.fasta -o busco_results -m genome -l busco_downloads/lineages/eukaryota_odb10 --cpu 8
```
# SALSA2
Only for Hi-C reads
```
bedtools bamtobed -i /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/hi-c_alignment.bam > /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/hic_alignments.bed
awk '/^>/ {if (seqlen){print seqlen}; printf substr($0,2) "\t"; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/hi-c_corrected_results/hi-c_pypolca_output/pypolca_corrected.fasta  > /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/hi-c_corrected_results/hi-c_pypolca_output/contig_lengths.txt
python ~/SALSA/run_pipeline.py \
  -a /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/hi-c_corrected_results/hi-c_pypolca_output/pypolca_corrected.fasta \
  -l /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/hi-c_corrected_results/hi-c_pypolca_output/contig_lengths.txt \
  -b /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/hic_alignments.bed \
  -o /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/hi-c_corrected_results/salsa2_output
```
Output you need: 
1. assembly.cleaned.fasta
2. scaffolds_FINAL.agp
3. scaffolds_FINAL.fasta
Run BUSCO again 
```
busco -i scaffolds_FINAL.fasta -l busco_downloads/lineages/eukaryota_odb10 -o busco_report -m genome
```
