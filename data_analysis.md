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
echo "Summarizing depth..."
jgi_summarize_bam_contig_depths --outputDepth depth.txt --percentIdentity 85 aligned_reads.sorted.bam

# Binning
echo "Starting MetaBAT2..."
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
# Diatom way
## Map coverage
```
bwa index pypolca_corrected.fasta
bwa mem -t 32 pypolca_corrected.fasta R1.fastq.gz R2.fastq.gz | \
samtools sort -o illumina.bam
samtools index illumina.bam
```
### BlobToolKit
```
conda create -n blobtoolkit -c conda-forge -c bioconda blobtoolkit
conda activate blobtoolkit

blobtools create \
  --fasta pypolca_corrected.fasta \
  diatom_blob

blobtools add \
  --cov illumina.bam \
  diatom_blob

# Extract contig IDs, Extract lengths, Extract GC, Extract coverage (Illumina), 
jq -r '.values[]' diatom_blob/identifiers.json > ids.txt
jq -r '.values[]' diatom_blob/length.json > lengths.txt
jq -r '.values[]' diatom_blob/gc.json > gc.txt
jq -r '.values[]' diatom_blob/illumina_cov.json > cov.txt

# Combine into one TSV: ID, Length, GC, Coverage
paste ids.txt lengths.txt gc.txt cov.txt > diatom_blob_view.tsv

# Nuclear diatom contigs
awk '$3>=0.45 && $3<=0.52 && $4>=20 {print $1}' diatom_blob_view.tsv > nuclear_contigs.txt

# Plastid contigs (high coverage, low GC example)
awk '$2>=100000 {print $1}' diatom_blob_view.tsv > plastid_large_contigs.txt

# Mitochondrial contigs (medium coverage, smaller size)
awk '$2<=100000 && $3>=0.42 && $3<=0.44 {print $1}' diatom_blob_view.tsv > mito_contigs.txt

seqtk subseq pypolca_corrected.fasta nuclear_contigs.txt > diatom_nuclear.fasta
seqtk subseq pypolca_corrected.fasta plastid_contigs.txt > diatom_plastid.fasta
seqtk subseq pypolca_corrected.fasta mito_contigs.txt > diatom_mito.fasta

stats.sh in=diatom_nuclear.fasta out=diatom_nuclear_stats.txt
stats.sh in=diatom_plastid.fasta out=plastid_stats.txt
stats.sh in=diatom_mito.fasta out=mito_stats.txt
```
### HMM way
```
#!/bin/bash
#SBATCH --job-name=diatom_manual_blob
#SBATCH --output=diatom_manual_blob.%j.out
#SBATCH --error=diatom_manual_blob.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=120:00:00
#SBATCH --mem=150G
###SBATCH --partition=cpu2023

# -----------------------------------------
# Load modules / activate conda env
# -----------------------------------------
module load miniconda3
module load bwa/0.7.17
module load minimap2/2.24
module load samtools/1.17
module load bbmap/38.84
module load diamond/2.0.6
module load repeatmodeler/2.0.1
module load repeatmasker/4.1.1
module load ebg_java/11.0.1
module load hmmer/3.3.2

conda activate diatom_env   # env with polypolish, busco, python tools

# -----------------------------------------
# User variables
# -----------------------------------------
ASSEMBLY="/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_out/pypolca_corrected.fasta"
ONT_READS="/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/pass_trim.fastq.gz"
ILLUMINA_R1="/work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz"
ILLUMINA_R2="/work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz"
THREADS=32

NUCLEAR_FASTA="/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_out/diatom_nuclear.fasta"
POLISHED_FASTA="/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_out/diatom_nuclear_polished.fasta"
BUSCO_LINEAGE="/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/busco_downloads/lineages/stramenopiles_odb10"

MITO_HMM="/work/ebg_lab/eb/diatom_consortia/markers/mito_markers.hmm"
PLASTID_HMM="/work/ebg_lab/eb/diatom_consortia/markers/plastid_markers.hmm"

# -----------------------------------------
# Step 1: Map reads and calculate coverage
# -----------------------------------------
echo "[1] Mapping Illumina reads and calculating coverage..."
bwa index $ASSEMBLY
bwa mem -t $THREADS $ASSEMBLY $ILLUMINA_R1 $ILLUMINA_R2 | samtools sort -@ $THREADS -o illumina.bam
samtools index illumina.bam

samtools depth -aa illumina.bam | \
awk '{cov[$1]+=$3; len[$1]++} END {for(c in cov) print c, cov[c]/len[c]}' > illumina_cov.tsv

echo "[2] Mapping ONT reads and calculating coverage..."
minimap2 -ax map-ont -t $THREADS $ASSEMBLY $ONT_READS | samtools sort -@ $THREADS -o ont.bam
samtools index ont.bam

samtools depth -aa ont.bam | \
awk '{cov[$1]+=$3; len[$1]++} END {for(c in cov) print c, cov[c]/len[c]}' > ont_cov.tsv

# -----------------------------------------
# Step 2: Calculate GC content
# -----------------------------------------
echo "[3] Calculating GC content per contig..."
stats.sh in=$ASSEMBLY out=assembly_stats.txt

# Extract contig name + GC fraction for downstream filtering
awk '/^>/{name=$1; getline seq; gsub(/[^GCgc]/,"",seq); gc=length(seq)/length($0); print name, gc}' $ASSEMBLY > gc_content.tsv

# -----------------------------------------
# Step 3: Search for organellar marker genes
# -----------------------------------------
echo "[4] Searching for mitochondrial and plastid contigs..."
hmmpress $MITO_HMM
hmmpress $PLASTID_HMM

hmmsearch --tblout mito_hits.tbl $MITO_HMM $ASSEMBLY
hmmsearch --tblout plastid_hits.tbl $PLASTID_HMM $ASSEMBLY

awk '{if($1!~/^#/){print $1}}' mito_hits.tbl | sort | uniq > mito_contigs.txt
awk '{if($1!~/^#/){print $1}}' plastid_hits.tbl | sort | uniq > plastid_contigs.txt

# -----------------------------------------
# Step 4: Assign nuclear contigs using coverage + GC + marker hits
# -----------------------------------------
echo "[5] Defining nuclear contigs..."
# Start from all contigs not matching mitochondrial/plastid markers
awk 'NR==FNR{org[$1]=1; next} !($1 in org)' mito_contigs.txt $ASSEMBLY | \
awk 'NR==FNR{org[$1]=1; next} !($1 in org)' plastid_contigs.txt > nuclear_contigs_raw.txt

# Filter by GC content (example range for diatom nuclear genome: 0.35–0.45)
awk 'NR==FNR{gc[$1]=$2; next} ($1 in gc) && gc[$1]>=0.35 && gc[$1]<=0.45 {print $1}' gc_content.tsv nuclear_contigs_raw.txt > nuclear_contigs.txt

# Extract nuclear contigs fasta
awk 'NR==FNR{c[$1]=1; next} /^>/{f=($1 in c)?1:0} f' nuclear_contigs.txt $ASSEMBLY > $NUCLEAR_FASTA

# -----------------------------------------
# Step 5: Polypolish nuclear contigs
# -----------------------------------------
echo "[6] Mapping Illumina reads to nuclear contigs..."
bwa index $NUCLEAR_FASTA
bwa mem -t $THREADS $NUCLEAR_FASTA $ILLUMINA_R1 $ILLUMINA_R2 | samtools sort -@ $THREADS -o illumina_nuclear.bam
samtools index illumina_nuclear.bam

echo "[7] Running Polypolish..."
polypolish $NUCLEAR_FASTA illumina_nuclear.bam > $POLISHED_FASTA

# Clean up BAM/SAM files
rm -f *.bam *.bam.bai *.sam

# -----------------------------------------
# Step 6: BUSCO completeness check
# -----------------------------------------
echo "[8] Running BUSCO..."
busco -i $POLISHED_FASTA -l $BUSCO_LINEAGE -m genome -c $THREADS -o busco_nuclear

# -----------------------------------------
# Step 7: Repeat identification and masking
# -----------------------------------------
echo "[9] Building repeat library and masking..."
BuildDatabase -name nuclear_db $POLISHED_FASTA
RepeatModeler -database nuclear_db -pa $THREADS -LTRStruct
RepeatMasker -pa $THREADS -lib nuclear_db-families.fa $POLISHED_FASTA

echo "[10] Pipeline complete!"
```

