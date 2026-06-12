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
# 12. Diatom genome annotation pipeline
Pipeline overview: 
Gene models were generated using a repeat-masked genome, RNA-seq alignments, and protein homology evidence within the BRAKER4 framework (ETP mode). Transcript-supported refinement was optionally performed using TSEBRA.
```
Genome (soft-masked)
   ↓
Repeat annotation (RepeatModeler2 + RepeatMasker)
   ↓
RNA-seq alignment (STAR → coordinate BAM)
   ↓
Protein evidence (UniProt Bacillariophyta + Stramenopiles)
   ↓
BRAKER4 (GeneMark-ETP + AUGUSTUS)
   ↓
Gene models (GTF/GFF3)
   ↓
TSEBRA (optional refinement with StringTie)
   ↓
CDS + protein extraction (gffread)
```
## 1. Software environment and dependencies
Genome annotation was performed using BRAKER4, which integrates RNA-seq evidence and protein homology for gene prediction. Repeat masking was applied prior to annotation to reduce spurious gene calls in repetitive regions.
A conda environment was used to manage dependencies:
```
conda create -n braker_env python=3.9 -y
conda activate braker_env

conda install -c bioconda \
    repeatmodeler repeatmasker star samtools stringtie augustus tsebra -y
```
GeneMark-ETP was installed separately due to licensing constraints and configured via environment variable:
```
export GENEMARK_KEY=/home/ruchita.solanki/.gm_key
```
## 2. Repeat identification and genome masking
Repetitive elements were identified de novo and masked prior to gene prediction to reduce false-positive gene models arising from transposable elements and low-complexity regions.
### 2.1 Repeat library construction
A species-specific repeat library was generated:
```
module load repeatmodeler/2.0.1
BuildDatabase -name genomeDB 18_diatom.fasta
RepeatModeler -database genomeDB -pa 32
```
RepeatModeler performed unsupervised discovery of:
- transposable elements (LINEs, SINEs, LTRs)
- tandem repeats
- low-complexity genomic regions
- lineage-specific repetitive elements                       
The output ```consensi.fa.classified``` served as the repeat library.
### 2.2 Genome soft masking
RepeatMasker was applied using the custom library:
```
RepeatMasker -pa 32 \
    -lib consensi.fa.classified \
    -xsmall \
    18_diatom.fasta
```
The resulting genome ```18_diatom.fasta.masked``` retained nucleotide sequence while masking repeats in lowercase. This masked assembly was used for all downstream steps.
## 3. RNA-seq alignment and evidence generation
RNA-seq reads from multiple biological replicates were pooled prior to alignment.
### 3.1 STAR genome indexing
```
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir genome_index \
     --genomeFastaFiles 18_diatom.fasta.masked \
     --genomeSAindexNbases 10
```
### 3.2 STAR alignment (pooled RNA-seq)
```
STAR --runThreadN 12 \
     --genomeDir genome_index \
     --readFilesCommand zcat \
     --readFilesIn R1_rep1.gz,R1_rep2.gz,R1_rep3.gz \
                    R2_rep1.gz,R2_rep2.gz,R2_rep3.gz \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix Diatoms_Combined_
```
The resulting BAM file provided:
- exon–intron junction evidence
- splice-aware read placement
- transcript coverage profiles across gene loci
### 3.3 BAM consolidation constraint
Although multiple RNA-seq libraries were available, they were merged into a single coordinate-sorted BAM prior to BRAKER execution due to Snakemake input validation constraints in the BRAKER4 workflow (check_bam_sorted rule requires explicit BAM listing per sample without ambiguous multi-file expansion).
## 4. Protein evidence (ETP input construction)
Protein evidence was constructed from UniProtKB reference proteomes using a taxonomy-constrained query targeting Bacillariophyta (TaxID: 2836).
```
wget -O diatoms.faa \
"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=taxonomy_id:2836"
```
### 4.1 Protein filtering and curation
The raw dataset was processed prior to BRAKER input:
```
seqkit seq -m 50 diatoms.faa > diatoms.clean.faa
seqkit rmdup -s diatoms.clean.faa > diatoms.nr.faa
```
This step:
- removed short sequences (<50 aa)
- collapsed exact duplicates
- reduced redundancy in homolog search space
### 4.2 Composition of protein evidence set
The resulting dataset contained homologous proteins from:
- Pseudo-nitzschia
- Nitzschia
- Seminavis
- Ditylum
- Grammatophora
- additional stramenopile lineages

This dataset provided phylogenetically distributed protein homology signals rather than a single-reference proteome, improving sensitivity across divergent gene families.
## 5. Genome Annotation with BRAKER4

BRAKER4 was used to generate evidence-supported structural gene predictions using the genome assembly, RNA-seq alignments, and protein evidence. Because the genome assembly was not soft-masked, repeat masking was performed internally by BRAKER4 using RED. RNA-seq and protein evidence were provided through the `samples.csv` file, which triggered the GeneMark-ETP and AUGUSTUS workflow.
### 5.1 BRAKER4 Setup
BRAKER4 was cloned from GitHub, and the BRAKER3 Singularity image was downloaded for containerized execution.
```
git clone https://github.com/Gaius-Augustus/BRAKER4.git
cd BRAKER4

singularity pull braker3.sif docker://teambraker/braker3:latest
```
Snakemake was installed in the working environment:
```
conda install -c conda-forge -c bioconda snakemake
```
### 5.2 Protein Evidence Preparation
Protein evidence was prepared from diatom and stramenopile protein databases. The final combined protein FASTA was filtered to remove short sequences and duplicate protein sequences.
```
seqkit seq -m 30 diatom_protein_evidence.raw.faa \
    | seqkit rmdup -s \
    > diatom_protein_evidence.clean.faa
```
The final protein evidence file was checked before use:
```
seqkit stats diatom_protein_evidence.clean.faa
grep -n "\*" diatom_protein_evidence.clean.faa | head
grep -n -v -E '^>|^[A-Z]+$' diatom_protein_evidence.clean.faa | head
```
No stop codons or malformed sequence lines were detected. The cleaned protein evidence file was used for BRAKER4:
```
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/proteins/diatom_protein_evidence.clean.faa
```
### 5.3 Sample Configuration
A `samples.csv` file was prepared to specify the genome assembly, protein evidence, RNA-seq alignment file, and BUSCO lineage. Both `protein_fasta` and `bam_files` were provided, which triggers ETP mode in BRAKER4.
```
nano samples.csv
```
```
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage
DL_diatom,/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta,,/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/proteins/diatom_protein_evidence.clean.faa,/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam,,,,,,,,stramenopiles_odb12
```
The `genome_masked` column was left empty because the genome was not soft-masked:
```
grep -v "^>" /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta \
    | grep -q '[a-z]' && echo "soft-masked" || echo "not soft-masked"
```
The `samples.csv` file was checked to confirm that the expected number of columns was present:
```
awk -F',' '{print NR, NF}' samples.csv
```
Expected output:
```
1 13
2 13
```
### 5.4 BRAKER4 Configuration
The `config.ini` file was edited to specify the Singularity image, GeneMark license key, sample file, repeat masking option, and run parameters.
```
nano config.ini
```
```
[paths]
braker_container = /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/braker3.sif
genemark_key = /home/ruchita.solanki/.gm_key

[DATA]
samples = samples.csv

[PARAMS]
fungus = False
min_contig = 1000
run_red = True
species = diatom_v1
mode = etp

[SLURM_ARGS]
cpus_per_task = 32
mem_of_node = 350000
max_runtime = 7200
```
The key setting was:
```
run_red = True
```
This was used because the genome assembly was not soft-masked.
### 5.5 Dry Run
A Snakemake dry run was performed before launching the full analysis:
```
snakemake \
    -s Snakefile \
    --use-singularity \
    --singularity-args "--bind /work/ebg_lab/eb/diatom_consortia/metatranscriptomics" \
    --cores 24 \
    --latency-wait 120 \
    --printshellcmds \
    -n
```
The dry run successfully built the workflow DAG and included the expected BRAKER4 rules, including:
```
run_masking
run_genemark_etp
run_tsebra
run_augustus_hints
busco_genome
busco_proteins
collect_results
```
The presence of `run_genemark_etp` confirmed that BRAKER4 recognized the RNA-seq and protein evidence and configured the run in ETP mode.
### 5.6 BRAKER4 Execution
After the dry run completed successfully, BRAKER4 was executed using the official BRAKER4 `Snakefile`.
```bash
snakemake \
    -s Snakefile \
    --use-singularity \
    --singularity-args "--bind /work/ebg_lab/eb/diatom_consortia/metatranscriptomics" \
    --cores 24 \
    --latency-wait 120 \
    --printshellcmds
```
If the GeneMark key is not visible inside the container, the home directory can also be bound:
```bash
snakemake \
    -s Snakefile \
    --use-singularity \
    --singularity-args "--bind /work/ebg_lab/eb/diatom_consortia/metatranscriptomics,/home/ruchita.solanki" \
    --cores 24 \
    --latency-wait 120 \
    --printshellcmds
```
BRAKER4 internally performs repeat masking, GeneMark-ETP training, AUGUSTUS training and prediction, evidence integration, TSEBRA refinement, BUSCO/compleasm assessment, and final result collection. The final annotation files generated by BRAKER4 include GTF/GFF3 gene models, predicted coding sequences, and predicted protein sequences for downstream functional annotation and BUSCO assessment.
