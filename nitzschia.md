#!/bin/bash
#SBATCH --job-name=diatom_manual_blob
#SBATCH --output=diatom_manual_blob.%j.out
#SBATCH --error=diatom_manual_blob.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=120:00:00
#SBATCH --mem=150G
#SBATCH --exclusive

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

conda activate diatom_env   # polypolish, busco, python tools

# -----------------------------------------
# User variables
# -----------------------------------------
ASSEMBLY="/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/sr_pypolca_out/pypolca_corrected.fasta"
ONT_READS="/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/pass_trim.fastq.gz"
ILLUMINA_R1="/work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz"
ILLUMINA_R2="/work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz"
THREADS=32

NUCLEAR_FASTA="/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/sr_pypolca_out/diatom_nuclear.fasta"
POLISHED_FASTA="/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/sr_pypolca_out/diatom_nuclear_polished.fasta"
BUSCO_LINEAGE="/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/busco_downloads/lineages/bacillariophyta_odb10"

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
# Step 2: Calculate GC content (multi-line FASTA)
# -----------------------------------------
echo "[3] Calculating GC content per contig..."
awk '
  /^>/{if(seq!=""){gsub(/[^GCgc]/,"",seq); print name, length(seq); seq=""} name=$1; seq=""; next} 
  {seq=seq $0} 
  END{gsub(/[^GCgc]/,"",seq); print name, length(seq)}
' $ASSEMBLY > contig_lengths.tsv

awk '
  NR==FNR{len[$1]=$2; next} 
  {seqs[$1]=seqs[$1] $0} 
  END{for(c in seqs){gsub(/[^GCgc]/,"",seqs[c]); gc=length(seq)/(len[c]); print c, gc}}
' $ASSEMBLY contig_lengths.tsv > gc_content.tsv

# -----------------------------------------
# Step 3: Search for organellar marker genes
# -----------------------------------------
echo "[4] Searching for mitochondrial and plastid contigs..."
# hmmpress only if .h3* files don't exist
[ ! -f "${MITO_HMM}.h3f" ] && hmmpress $MITO_HMM
[ ! -f "${PLASTID_HMM}.h3f" ] && hmmpress $PLASTID_HMM

hmmsearch --cpu $THREADS -E 1e-5 --tblout mito_hits.tbl $MITO_HMM $ASSEMBLY
hmmsearch --cpu $THREADS -E 1e-5 --tblout plastid_hits.tbl $PLASTID_HMM $ASSEMBLY

awk '{if($1!~/^#/){print $1}}' mito_hits.tbl | sort | uniq > mito_contigs.txt
awk '{if($1!~/^#/){print $1}}' plastid_hits.tbl | sort | uniq > plastid_contigs.txt

# -----------------------------------------
# Step 4: Assign nuclear contigs (GC + coverage + marker filtering)
# -----------------------------------------
echo "[5] Defining nuclear contigs..."
# Remove organelle contigs
awk 'NR==FNR{org[$1]=1; next} !($1 in org)' mito_contigs.txt $ASSEMBLY | \
awk 'NR==FNR{org[$1]=1; next} !($1 in org)' plastid_contigs.txt > nuclear_contigs_raw.txt

# Filter by GC (Nitzschia nuclear GC ~0.38–0.46)
awk 'NR==FNR{gc[$1]=$2; next} ($1 in gc) && gc[$1]>=0.38 && gc[$1]<=0.46 {print $1}' gc_content.tsv nuclear_contigs_raw.txt > nuclear_gc_filtered.txt

# Filter by Illumina coverage (optional: remove very low/high coverage)
awk 'NR==FNR{cov[$1]=$2; next} ($1 in cov) && cov[$1]>10 && cov[$1]<50 {print $1}' illumina_cov.tsv nuclear_gc_filtered.txt > nuclear_final.txt

# Extract nuclear contigs fasta
awk 'NR==FNR{c[$1]=1; next} /^>/{f=($1 in c)?1:0} f' nuclear_final.txt $ASSEMBLY > $NUCLEAR_FASTA

# -----------------------------------------
# Step 5: Polypolish nuclear contigs
# -----------------------------------------
echo "[6] Mapping Illumina reads to nuclear contigs..."
bwa index $NUCLEAR_FASTA
bwa mem -t $THREADS $NUCLEAR_FASTA $ILLUMINA_R1 $ILLUMINA_R2 | samtools sort -@ $THREADS -o illumina_nuclear.bam
samtools index illumina_nuclear.bam

echo "[7] Running Polypolish..."
polypolish $NUCLEAR_FASTA illumina_nuclear.bam > $POLISHED_FASTA

# Clean up intermediate BAM/SAM
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
# Optionally filter small contigs (<5 kb) for repeat modeling
awk 'NR==FNR{len[$1]=$2; next} ($1 in len) && len[$1]>=5000 {print $1}' contig_lengths.tsv nuclear_final.txt > nuclear_long.txt
awk 'NR==FNR{c[$1]=1; next} /^>/{f=($1 in c)?1:0} f' nuclear_long.txt $NUCLEAR_FASTA > nuclear_long.fasta

BuildDatabase -name nuclear_db nuclear_long.fasta
RepeatModeler -database nuclear_db -pa $THREADS -LTRStruct
RepeatMasker -pa $THREADS -lib nuclear_db-families.fa $POLISHED_FASTA

echo "[10] Pipeline complete!"
