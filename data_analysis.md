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

# -----------------------------------------
# Load modules
# -----------------------------------------
module load bbmap/38.84              
module load bwa/0.7.17               
module load miniconda3/samtools      
module load hmmer/v3.3               
module load prodigal/v2.6.3          
module load repeatmodeler/2.0.1      
module load repeatmasker/4.1.1

conda activate diatom_env 

# -----------------------------------------
# Variables
# -----------------------------------------
ASSEMBLY="/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_out/pypolca_corrected.fasta"
ILLUMINA_R1="/work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz"
ILLUMINA_R2="/work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz"
THREADS=32

MITO_HMM="/work/ebg_lab/eb/diatom_consortia/markers/mito_markers.hmm"
PLASTID_HMM="/work/ebg_lab/eb/diatom_consortia/markers/plastid_markers.hmm"
BUSCO_LINEAGE="/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/busco_downloads/lineages/stramenopiles_odb10"

# -----------------------------------------
# Step 1: Mapping and Coverage
# -----------------------------------------
echo "[1] Calculating Illumina coverage..."
bwa index $ASSEMBLY
bwa mem -t $THREADS $ASSEMBLY $ILLUMINA_R1 $ILLUMINA_R2 | samtools sort -@ $THREADS -o illumina_full.bam
samtools index illumina_full.bam

# Output: ContigID <tab> AverageDepth
samtools depth -aa illumina_full.bam | awk '{cov[$1]+=$3; len[$1]++} END {for(c in cov) print c"\t"cov[c]/len[c]}' > illumina_cov.tsv

# -----------------------------------------
# Step 2: Marker Search (Organelle Detection)
# -----------------------------------------
echo "[2] Predicting proteins and searching for markers..."
prodigal -i $ASSEMBLY -a proteins.faa -p meta -q

hmmsearch --tblout mito.tbl --cpu $THREADS $MITO_HMM proteins.faa
hmmsearch --tblout plastid.tbl --cpu $THREADS $PLASTID_HMM proteins.faa

# Create a combined list of organelle contig IDs
awk '!/^#/ {print $1}' mito.tbl plastid.tbl | sed 's/_[0-9]*$//' | sort -u > organelle_ids.txt

# -----------------------------------------
# Step 3: Filtering (AWK-based Master Table)
# -----------------------------------------
echo "[3] Filtering for nuclear contigs using AWK..."
# Get Length and GC via BBTools (format 6: name length gc ...)
stats.sh in=$ASSEMBLY format=6 > assembly_stats.tsv

# Join stats with coverage and filter
# Logic: GC 0.30-0.50, Depth > 5, Not in organelle_ids.txt
awk 'BEGIN {FS="\t"; OFS="\t"} 
    # Load organelle IDs into an array
    NR==FNR {org[$1]=1; next} 
    # Load coverage into an array (from illumina_cov.tsv)
    FILENAME=="illumina_cov.tsv" {cov[$1]=$2; next}
    # Process assembly_stats.tsv
    FILENAME=="assembly_stats.tsv" {
        name=$1; gc=$3; 
        if (name in cov && !(name in org) && gc >= 0.30 && gc <= 0.50 && cov[name] > 5) {
            print name
        }
    }' organelle_ids.txt illumina_cov.tsv assembly_stats.tsv > nuclear_ids.txt

# Extract Nuclear FASTA using AWK
awk 'NR==FNR {a[$1]; next} /^>/ {f=0; id=$1; sub(/^>/,"",id); if (id in a) f=1} f' nuclear_ids.txt $ASSEMBLY > diatom_nuclear.fasta

# -----------------------------------------
# Step 4: Polypolish
# -----------------------------------------
echo "[4] Polishing nuclear genome..."
bwa index diatom_nuclear.fasta
bwa mem -t $THREADS -a diatom_nuclear.fasta $ILLUMINA_R1 > aln1.sam
bwa mem -t $THREADS -a diatom_nuclear.fasta $ILLUMINA_R2 > aln2.sam
polypolish diatom_nuclear.fasta aln1.sam aln2.sam > diatom_nuclear_polished.fasta

# -----------------------------------------
# Step 5: BUSCO & Repeats
# -----------------------------------------
echo "[5] Final QC and Masking..."
busco -i diatom_nuclear_polished.fasta -l $BUSCO_LINEAGE -m genome -c $THREADS -o busco_out

BuildDatabase -name nuc_db diatom_nuclear_polished.fasta
RepeatModeler -database nuc_db -pa $THREADS -LTRStruct
RepeatMasker -pa $THREADS -lib nuc_db-families.fa diatom_nuclear_polished.fasta

echo "Pipeline complete."
```

