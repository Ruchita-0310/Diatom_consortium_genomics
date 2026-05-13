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
## Polishing - SR
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
# 8. Phylogenetic tree
```
clustalo -i 18S_new.fasta -o 18S_aligned.fasta
trimal -in 18S_aligned.fasta -out 18S_trimmed.fasta -automated1
/home/ruchita.solanki/iqtree-2.2.2.7-Linux/bin/iqtree2 -s 18S_trimmed.fasta -B 1000 -T 4
```
Once you have the tree file, run this to get names next to the accession IDs
```
python3 -c '
import re
import glob

# Using your specific tree filename
tree_filename = "18S_trimmed.fasta.treefile"
name_map = {}

# 1. Loop through all .fasta files (seqdump, Nitzschia_18S_full, and individual IDs)
fasta_files = glob.glob("*.fasta")
print(f"Reading names from: {len(fasta_files)} files...")

for fasta in fasta_files:
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                parts = line.strip().split()
                acc = parts[0].replace(">", "")
                if len(parts) >= 3:
                    # Extracts Genus_species and cleans special characters
                    species = re.sub(r"[^a-zA-Z0-9_]", "", f"{parts[1]}_{parts[2]}")
                    name_map[acc] = f"{acc}_{species}"

# 2. Perform the swap on the tree
try:
    with open(tree_filename, "r") as t:
        content = t.read()
    
    # Sort keys by length (longest first) to prevent partial matching bugs
    for acc in sorted(name_map.keys(), key=len, reverse=True):
        content = content.replace(acc, name_map[acc])
    
    output_name = "final_species_tree.tre"
    with open(output_name, "w") as out:
        out.write(content)
    
    print(f"\nDone! Mapped {len(name_map)} unique IDs.")
    print(f"Your fixed tree is: {output_name}")

except FileNotFoundError:
    print(f"\nError: Could not find {tree_filename}")
'
```
You can next remove all the "Uncultured", to clean up the tree:
```
from Bio import Phylo

# Load the tree
tree = Phylo.read("new18_species_tree.tre", "newick")

# Identify all terminals (leaves) containing "Uncultured"
to_prune = [leaf for leaf in tree.get_terminals() if "Uncultured" in leaf.name]

# Remove them
for leaf in to_prune:
    tree.prune(leaf)

# Save the cleaned tree
Phylo.write(tree, "cleaned_species_tree.tre", "newick")

print(f"Removed {len(to_prune)} nodes.")
```
# 9. Transcriptome analysis 
[Nf core metadenovo](https://github.com/nf-core/metatdenovo)
```
# --- 1. JAVA SETUP ---
module purge
module load java/openjdk-23.0.1
export JAVA_HOME=$(dirname $(dirname $(readlink -f $(which java))))
export PATH=$JAVA_HOME/bin:$PATH

# --- 2. EXECUTION ---
~/nextflow run nf-core/metatdenovo \
    -profile singularity \
    --input samplesheet.csv \
    --outdir /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/new_results \
    -w /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/work \
    --assembler spades \
    --orf_caller transdecoder \
    --eggnog_dbpath /home/ruchita.solanki/eggnog_db \
    --skip_kofam true \
    --hmmfiles /home/ruchita.solanki/Pfam-A.hmm \
    --eukulele_dbpath /home/ruchita.solanki/eukulele_db \
    --eukulele_db mmetsp \
    -resume \
    -with-report report_skipK.html \
    -with-timeline timeline_skipK.html
```

# 10. Identifying rRNA genes from transcriptome
```
barrnap --kingdom euk --threads 4 spades.transcripts.fa --outseq euk_transcript_rRNA.fna > diatom_euk_rRNA.gff
barrnap --kingdom bac spades.transcripts.fa --outseq bac_transcript_rRNA.fna > diatom_bac_rRNA.gff
barrnap --kingdom mito spades.transcripts.fa --outseq mito_transcript_rRNA.fna > diatom_mito_rRNA.gff
```
# 11. Identifying organelle genome
Download mitogenome - MT742552 & chloroplast genome - MT742551
```
conda create -n quast_env quast
metaquast.py 8_diatom.fasta -R /work/ebg_lab/eb/diatom_consortia/organelle/ref/ -o ./8_metaquast_output
metaquast.py /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_output/pypolca_corrected.fasta -R /work/ebg_lab/eb/diatom_consortia/organelle/ref/ -o ./whole_metaquast_output
```
# 12. Diatom Genome Annotation Pipeline
## 1. Software Environment and Dependencies
Structural genome annotation was performed using the BRAKER4 framework, which integrates RNA-seq evidence with ab initio gene prediction. Transcript-supported refinement of gene models was performed using TSEBRA. Repeat identification and masking were conducted prior to annotation to minimize false-positive gene predictions arising from repetitive genomic regions.
```
# Create conda environment
conda create -n braker_env python=3.9 -y
conda activate braker_env
```
```
# Install required tools
conda install -c bioconda \
    repeatmodeler \
    repeatmasker \
    star \
    samtools \
    stringtie \
    augustus \
    tsebra \
    -y
```
GeneMark-ETP was installed separately following license registration.
```
module load ebg_perl/5.32.0 ebg_perl_modules/5.32.0 miniconda3/4.8.3
module load GeneMark/GeneMark-ES/v4
```
## 2. Repeat Identification and Genome Masking
Repetitive elements were identified de novo using RepeatModeler and subsequently soft-masked using RepeatMasker. Soft masking converts repetitive regions to lowercase sequence while preserving nucleotide information, thereby reducing spurious gene predictions during ab initio annotation.
### 2.1 Repeat Library Construction
```
module load ebg_perl/5.32.0 ebg_perl_modules/5.32.0 recon/1.08 miniconda3/h5py repeatmasker/4.1.1 miniconda3/4.8.3 repeatscout/1.0.6 rmblast/2.10.0 trf/4.09
module load repeatmodeler/2.0.1

BuildDatabase \
    -name genomeDB \
    18_diatom.fasta

RepeatModeler \
    -database genomeDB \
    -pa 32
```
This generated a custom repeat library: ```consensi.fa.classified```
### 2.2 Soft Masking of the Genome
```
RepeatMasker \
    -pa 32 \
    -lib consensi.fa.classified \
    -xsmall \
    18_diatom.fasta
```
The resulting soft-masked genome: ```18_diatom.fasta.masked``` was used for all downstream analyses.

## 3. RNA-seq Alignment to the Soft-Masked Genome
RNA-seq triplicates were pooled and aligned to the soft-masked genome using the splice-aware aligner STAR. RNA-seq evidence improves prediction of exon–intron boundaries and transcript structures.
### 3.1 Genome Index Generation
```
STAR \
    --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir genome_index \
    --genomeFastaFiles 18_diatom.fasta.masked \
    --genomeSAindexNbases 10
```

### 3.2 RNA-seq Alignment
```
BASE="/work/ebg_lab/eb/diatom_consortia/metatranscriptomics"

STAR \
    --runThreadN 12 \
    --genomeDir genome_index \
    --readFilesCommand zcat \
    --readFilesIn \
        $BASE/Li57991-Diatoms-1-4C_S1_R1_001.fastq.gz,\
$BASE/Li57992-Diatoms-2-4C_S2_R1_001.fastq.gz,\
$BASE/Li57993-Diatoms-3-4C_S3_R1_001.fastq.gz \
        $BASE/Li57991-Diatoms-1-4C_S1_R2_001.fastq.gz,\
$BASE/Li57992-Diatoms-2-4C_S2_R2_001.fastq.gz,\
$BASE/Li57993-Diatoms-3-4C_S3_R2_001.fastq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix Diatoms_Combined_
```
The resulting coordinate-sorted BAM file was used as transcriptomic evidence for downstream annotation.

## 4. Generation of RNA-seq Hints
Intron hints were extracted from RNA-seq alignments using bam2hints. These hints provide extrinsic splice-site evidence during AUGUSTUS prediction and TSEBRA refinement.
```
bam2hints \
    --intronsonly \
    --in=Diatoms_Combined_Aligned.sortedByCoord.out.bam \
    --out=rna_hints.gff
```
The RNA-seq evidence tag (src=E) was appended for compatibility with TSEBRA scoring.
```
sed -i 's/$/src=E;/' rna_hints.gff
```
## 5. Genome Annotation with BRAKER4
BRAKER4 was used to generate evidence-supported structural gene predictions from the soft-masked genome and RNA-seq alignments.
### 5.1 BRAKER4 Setup
```
git clone https://github.com/Gaius-Augustus/BRAKER4.git
cd BRAKER4
```
### 5.2 Sample Configuration
A samples.csv configuration file was prepared specifying the genome assembly and transcriptomic evidence.
```
echo "Diatom_18,/work/ebg_lab/eb/metatranscriptomics/18_diatom.fasta.masked,/work/ebg_lab/eb/metatranscriptomics/Diatoms_Combined_Aligned.sortedByCoord.out.bam," > samples.csv
```
### 5.3 BRAKER4 Execution
```
# Install snakemake in braker env
conda install -c conda-forge -c bioconda snakemake

snakemake \
    --use-singularity \
    --singularity-args "--bind /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index:/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index" \
    --cores 24 \
    --config sif_image=braker3.sif
```
BRAKER4 internally performs:
- GeneMark self-training
- AUGUSTUS-based ab initio prediction
- Incorporation of RNA-seq splice evidence
- Evidence-guided transcript selection

## 6. Transcript-Supported Refinement Using TSEBRA
To further refine gene models, TSEBRA was used to integrate transcript-supported predictions with ab initio AUGUSTUS predictions.
### 6.1 Transcript Assembly and Candidate Prediction
Transcript-supported gene models were generated using StringTie.
```
stringtie \
    Diatoms_Combined_Aligned.sortedByCoord.out.bam \
    -o stringtie_preds.gtf
```
### 6.2 Ab Initio Prediction
Additional ab initio predictions were generated using AUGUSTUS with Phaeodactylum tricornutum as the closest available reference species.
```
augustus \
    --species=phaeodactylum_tricornutum \
    18_diatom.fasta.masked \
    > augustus_preds.gtf
```
### 6.3 Standardization of GTF Files
```
fix_gtf_ids.py \
    --gtf stringtie_preds.gtf \
    --out set1.gtf

fix_gtf_ids.py \
    --gtf augustus_preds.gtf \
    --out set2.gtf
```
### 6.4 TSEBRA Integration
TSEBRA was used to select optimal gene models based on transcriptomic support.
```
tsebra \
    -g set1.gtf,set2.gtf \
    -e rna_hints.gff \
    -c default.cfg \
    -o final_diatom_annotation.gtf
```
### 6.5 Final Annotation Formatting
```
rename_gtf.py \
    --gtf final_diatom_annotation.gtf \
    --prefix Diatom_Consortia \
    --out 18_diatom_final.gtf
```
The final output consisted of a transcript-supported structural annotation of the soft-masked diatom genome in GTF format.
