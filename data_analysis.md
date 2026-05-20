# Diatom Consortia: Metagenomic & Metatranscriptomic Pipeline
This repository contains the end-to-end workflow for the assembly, polishing, binning, and annotation of Diatom-associated microbial consortia.                            
🛠 Prerequisites & Environment
The following software environments are required. It is recommended to manage these via Conda or Singularity as noted:
1. Assembly/Polishing: Flye, Medaka, Polypolish, Pypolca
2. Quality & Validation: BUSCO, CheckM2, QUAST/MetaQUAST
3. Binning & Taxonomy: MetaBAT2, GTDB-Tk, CoverM
4. Phylogenetics: ClustalO, TrimAl, IQ-TREE 2, Biopython
5. Annotation: RepeatModeler2, RepeatMasker, STAR, BRAKER4, TSEBRA, StringTie
6. Metatranscriptomics: Nextflow, nf-core/metatdenovo

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
cat *.fasta > 18S_new.fasta
clustalo -i 18S_new.fasta -o 18S_aligned.fasta
trimal -in 18S_aligned.fasta -out 18S_trimmed.fasta -automated1
/home/ruchita.solanki/iqtree-2.2.2.7-Linux/bin/iqtree2 -s 18S_trimmed.fasta -m MFP -bb 1000 -alrt 1000 -nt AUTO
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
## Pipeline Logic
```
Genome
   ↓
Repeat masking
   ↓
RNA-seq alignment to genome
   ↓
BAM files
   ↓
BRAKER4
   ↓
Gene predictions
   ↓
TSEBRA refinement
   ↓
Functional annotation
```
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
You are building a species-specific repeat database from your diatom genome.
```
module load ebg_perl/5.32.0 ebg_perl_modules/5.32.0 recon/1.08 miniconda3/h5py repeatmasker/4.1.1 miniconda3/4.8.3 repeatscout/1.0.6 rmblast/2.10.0 trf/4.09
module load repeatmodeler/2.0.1
```
This:
- converts your FASTA genome into a searchable database
- creates index files needed by RepeatModeler
- does not identify repeats yet
```
BuildDatabase \
    -name genomeDB \
    18_diatom.fasta
```
This is the actual repeat discovery step. RepeatModeler scans the genome and tries to find:
- transposable elements (TEs)
- tandem repeats
- low-complexity regions
- repetitive fragments
- novel repeats specific to your organism
```
RepeatModeler \
    -database genomeDB \
    -pa 32
```
This generated a custom repeat library: ```consensi.fa.classified```
### 2.2 Soft Masking of the Genome
Now you use the repeat library to locate repeats across the genome.
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
Why STAR is “splice-aware”                     
Eukaryotic genes contain introns.
RNA-seq reads often span exon junctions: ```Exon1 ---- intron ---- Exon2```               
A read may align like: ```[Exon1][Exon2]```                   
Normal aligners fail because part of the read is missing from genomic sequence continuity.                    
STAR detects splice junctions and aligns across introns correctly.                  
That is why STAR is ideal for:
- eukaryotic transcriptomes
- BRAKER
- exon prediction
### 3.1 Genome Index Generation
This step prepares the genome for rapid RNA-seq alignment.
```
STAR \
    --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir genome_index \
    --genomeFastaFiles 18_diatom.fasta.masked \
    --genomeSAindexNbases 10
```
### 3.2 RNA-seq Alignment
Now STAR maps RNA reads back onto the genome.
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
This BAM file contains:
- where every RNA read aligned
- splice junctions
- transcript evidence

## 4. Generation of RNA-seq Hints
Intron hints were extracted from RNA-seq alignments using bam2hints. These hints provide extrinsic splice-site evidence during AUGUSTUS prediction and TSEBRA refinement.
```
bam2hints \
    --intronsonly \
    --in=Diatoms_Combined_Aligned.sortedByCoord.out.bam \
    --out=rna_hints.gff
```
Create BAM index
```
samtools index genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam
```
The RNA-seq evidence tag (src=E) was appended for compatibility with TSEBRA scoring.
```
sed -i 's/$/src=E;/' rna_hints.gff
```

## 5. Genome Annotation with BRAKER4
Structural gene annotation of the soft-masked diatom genome assembly was performed using the BRAKER4 workflow (Gaius-Augustus/BRAKER4). Gene prediction was based on a combination of RNA-seq splice evidence and cross-species protein homology under an evidence-guided training framework. The genome assembly was soft-masked prior to annotation using repeat annotations generated from RepeatModeler and RepeatMasker, as described in Section 2. Masking reduced spurious ab initio gene predictions in repetitive regions while preserving nucleotide-level information required for downstream evidence mapping.
### 5.1 Installation and Environment Setup
BRAKER4 was executed within a Singularity container to ensure reproducibility and dependency isolation.
#### 5.1.1 Clone BRAKER4 and Pull Container
```
git clone https://github.com/Gaius-Augustus/BRAKER4.git
singularity pull braker3.sif docker://teambraker/braker3:latest
```
#### 5.1.2 Conda Environment and Snakemake Installation
A dedicated conda environment was used for workflow execution:
```
source ~/miniforge3/etc/profile.d/conda.sh
conda create -n braker_env -c conda-forge -c bioconda snakemake
conda activate braker_env
```
5.1.3 GeneMark License Configuration
GeneMark-ETP requires a valid license key:
```
export GENEMARK_KEY=/home/ruchita.solanki/.gm_key
```
### 5.2 Input Preparation
#### 5.2.1 Genome, RNA-seq, and Protein Evidence
The following inputs were used:
Soft-masked genome assembly: 18_diatom.fasta.masked (linked as data/genome.fa)
RNA-seq evidence: Coordinate-sorted STAR alignment BAM derived from pooled Illumina libraries: Diatoms_Combined_Aligned.sortedByCoord.out.bam
Protein evidence: A curated multi-species diatom and stramenopile UniProt FASTA dataset (data/proteins.fa) containing homologous proteins from multiple Bacillariophyta and related taxa. This dataset was used to enable protein-guided gene structure refinement under ETP mode.
```
wget -O diatoms.faa \
"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=taxonomy_id:2836"

seqkit seq -m 50 diatoms.faa > diatoms.clean.faa
seqkit rmdup -s diatoms.clean.faa > diatoms.nr.faa
```
A data/ directory was created within the BRAKER workflow, and symbolic links were used for standardized input handling:
```
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4
mkdir -p data

ln -s /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta.masked data/genome.fa
ln -s /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam data/rnaseq1.bam
ln -s /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/data/diatoms_etp.clean.faa data/proteins.fa
```
#### 5.2.2 RNA-seq BAM Validation
The RNA-seq BAM file was confirmed to be coordinate-sorted prior to annotation:
```
samtools view -H data/rnaseq1.bam | grep SO
```
Expected: SO:coordinate
This ensures compatibility with GeneMark-ETP intron inference and BRAKER splice hint generation.
#### 5.2.3 Protein Evidence for ETP Mode
Protein evidence was derived from a curated diatom-stramenopile FASTA dataset. The dataset includes homologous proteins from multiple genera (e.g., Pseudo-nitzschia, Nitzschia, Seminavis, Ditylum, Grammatophora), providing phylogenetically informed support for gene boundary refinement. Protein evidence is used exclusively to guide gene model training and does not impose hard constraints on predicted gene structures.
#### 5.2.4 Input Constraint (Single BAM Requirement)
Due to Snakemake rule-level input validation in the BRAKER4 workflow, RNA-seq BAM inputs must be provided as a single coordinate-sorted file per sample. Initial attempts to supply multiple BAMs as semicolon-separated entries resulted in Snakemake parsing failures during BAM validation steps (check_bam_sorted). Therefore, RNA-seq evidence was consolidated into a single pooled BAM alignment for consistent execution.
### 5.3 BRAKER4 Configuration
BRAKER4 was executed in ETP mode, enabling integrated use of:
- RNA-seq splice evidence
- protein homology evidence
- ab initio gene prediction (GeneMark-ETP + AUGUSTUS)
#### 5.3.1 samples.csv
```
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage
example_species,data/genome.fa,,data/proteins.fa,data/rnaseq1.bam,,,,,,,,eukaryota_odb12
```
#### 5.3.2 config.ini
```
[paths]
braker_container = /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/braker3.sif
genemark_key = /home/ruchita.solanki/.gm_key

[DATA]
samples = samples.csv

[PARAMS]
fungus = False
min_contig = 1000
run_red = False
species = diatom_v1
mode = etp

[SLURM_ARGS]
cpus_per_task = 32
mem_of_node = 350
max_runtime = 120
```
#### 5.3.3 Annotation Mode Justification
ETP mode was selected because both RNA-seq and protein homology evidence were available. This allows GeneMark-ETP to integrate splice junction evidence with cross-species protein alignments, improving exon boundary resolution relative to ET-only workflows.
### 5.4 BRAKER4 Execution
```
snakemake \
    --use-singularity \
    --singularity-args "-B /home/ruchita.solanki:/home/ruchita.solanki \
    -B /work/ebg_lab/eb/diatom_consortia/metatranscriptomics:/work/ebg_lab/eb/diatom_consortia/metatranscriptomics" \
    --cores 32 \
    --latency-wait 60 \
    --rerun-incomplete \
    --printshellcmds
```
5.5 BRAKER4 Annotation Workflow
Within BRAKER4:
- RNA-seq BAM was validated and used to generate splice hints
- GeneMark-ETP performed evidence-guided training using RNA-seq and protein homology data
- AUGUSTUS generated ab initio gene predictions informed by trained parameters
- Final gene models were exported in GTF and GFF3 formats
Final annotation outputs were generated in: ```output/example_species/results/```
## 6. Post-processing and Sequence Extraction
Protein and CDS sequences were extracted from the final gene models:
```
gffread -y diatom_proteins.faa -x diatom_cds.fna -g genome.fa braker.gtf
```
The resulting protein FASTA (```diatom_proteins.faa```) represents the final gene set used for downstream functional annotation and BUSCO evaluation.
