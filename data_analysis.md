# Diatom Consortia: Metagenomic and Metatranscriptomic Pipeline
This repository contains the workflow used for assembly, polishing, binning, taxonomic classification, organelle identification, transcriptome analysis, and genome annotation of a diatom-associated microbial consortium.
## Prerequisites and Software Environment
The following tools were used across different stages of the analysis. Software was installed using Conda, Singularity, or local module systems depending on availability on the HPC cluster.
1. **Assembly and polishing:** Flye, Medaka, Polypolish, Pypolca
2. **Assembly quality and validation:** BUSCO, CheckM2, QUAST/MetaQUAST
3. **Binning and taxonomy:** MetaBAT2, GTDB-Tk, CoverM
4. **Phylogenetics:** Clustal Omega, TrimAl, IQ-TREE 2, Biopython
5. **Genome annotation:** STAR, BRAKER4, GeneMark-ETP, AUGUSTUS, TSEBRA, BUSCO/compleasm
6. **Metatranscriptomics:** Nextflow, nf-core/metatdenovo, TransDecoder, eggNOG-mapper

---
# 1. Genome Assembly
Long-read assembly was performed using Nanopore reads basecalled with Guppy. The assembly was generated using Flye in metagenome mode.
```bash
flye \
    --nano-raw pass_trim.fastq.gz \
    --meta \
    -g 50m \
    --min-overlap 5000 \
    --out-dir flye_out_new \
    -i 3 \
    --threads 8
```

---
# 2. Read Mapping
## 2.1 Mapping Short Reads to the Assembly
Short reads were mapped to the Nanopore assembly to assess read support, coverage, and assembly quality.
```bash
minimap2 -ax sr \
    /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/guppy_flye_assembly.fasta \
    Diatoms_merged.fastq.gz \
    > sr_alignment.sam

samtools view -S -b sr_alignment.sam > alignment.bam

samtools sort alignment.bam -o alignment_sorted.bam

samtools index alignment_sorted.bam

samtools flagstat alignment_sorted.bam > mapping_stats.txt

samtools idxstats alignment_sorted.bam \
    | sort -k3,3rn \
    > sr_all_nanopore_hits.tsv

samtools depth alignment_sorted.bam > sr_depth.txt

awk '{sum[$1]+=$3; count[$1]++} END {for (c in sum) print c, sum[c]/count[c]}' sr_depth.txt \
    | sort -k2,2nr \
    > sr_mean_depth.tsv
```
The resulting files were used to evaluate mapping rate, contig-level coverage, and short-read support across the assembly.

---
# 3. Assembly Polishing
Assembly polishing was performed using long-read polishing with Medaka followed by short-read polishing with Polypolish and Pypolca.
## 3.1 Long-read Polishing with Medaka
```bash
medaka_consensus \
    -i pass_trim.fastq.gz \
    -d guppy_flye_assembly.fasta \
    -o medaka_euk_polished \
    -t 12
```
The Medaka-polished assembly was used as input for short-read polishing.
## 3.2 Short-read Mapping for Polypolish
Short reads were aligned to the Medaka-polished assembly using BWA-MEM.
```bash
bwa mem -t 16 -a \
    /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/medaka_euk_polished/consensus.fasta \
    /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz \
    > alignments_1.sam

bwa mem -t 16 -a \
    /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/medaka_euk_polished/consensus.fasta \
    /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz \
    > alignments_2.sam
```
## 3.3 Polypolish Filtering
```bash
polypolish filter \
    --in1 alignments_1.sam \
    --in2 alignments_2.sam \
    --out1 filtered_1.sam \
    --out2 filtered_2.sam
```
## 3.4 Polypolish Polishing
```bash
polypolish polish \
    /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/medaka_euk_polished/consensus.fasta \
    filtered_1.sam \
    filtered_2.sam \
    > sr_poly.fasta
```
## 3.5 Pypolca Polishing
```bash
pypolca run \
    -a sr_poly.fasta \
    -1 /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz \
    -2 /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz \
    -t 12 \
    -o sr_pypolca_output \
    --careful
```
The final corrected assembly was used for downstream assembly assessment, binning, organelle identification, and genome annotation.
## 3.6 BUSCO Assessment of the Polished Assembly
BUSCO was used to assess eukaryotic gene completeness in the polished genome assembly.
```bash
busco \
    -i pypolca_corrected.fasta \
    -l busco_downloads/lineages/stramenopiles_odb10 \
    -o busco_report \
    -m genome
```
---
# 4. Metagenomic Binning with MetaBAT2
MetaBAT2 was used to recover genome bins from the polished assembly using Nanopore read coverage.
## 4.1 Install MetaBAT2
```bash
conda create -n metabat2_v2 -c conda-forge -c bioconda metabat2 libdeflate=1.10
conda activate metabat2_v2
```
## 4.2 Map Nanopore Reads to the Polished Assembly
```bash
minimap2 -ax map-ont -t 16 \
    1_sr_pypolca_output/pypolca_corrected.fasta \
    pass_trim.fastq.gz \
    | samtools view -@ 16 -bS - \
    | samtools sort -@ 16 -m 10G -o aligned_reads.sorted.bam

samtools index -@ 16 aligned_reads.sorted.bam
```
## 4.3 Generate Contig Depth File
```bash
jgi_summarize_bam_contig_depths \
    --outputDepth depth.txt \
    --percentIdentity 85 \
    aligned_reads.sorted.bam
```
## 4.4 Run MetaBAT2
```bash
mkdir -p 2_metabat2_bins

metabat2 \
    -i 1_sr_pypolca_output/pypolca_corrected.fasta \
    -a depth.txt \
    -o 2_metabat2_bins/bin \
    -m 1500 \
    -t 16 \
    --unbinned
```
The resulting bins were used for genome quality assessment and taxonomic classification.

---
# 5. Bin Quality Assessment with CheckM2
CheckM2 was used to estimate completeness and contamination of recovered genome bins.
```bash
checkm2 predict \
    --threads 16 \
    --input 2_metabat2_bins/ \
    --output_directory 3_checkm2_results
```
---
# 6. Taxonomic Classification with GTDB-Tk
GTDB-Tk was used to assign taxonomy to recovered bins using the Genome Taxonomy Database.
```bash
gtdbtk classify_wf \
    --genome_dir 2_metabat2_bins/ \
    --out_dir 4_gtdbtk_output \
    --cpus 16 \
    -x fa
```
---
# 7. Genome Coverage and Relative Abundance with CoverM
CoverM was used to calculate coverage and relative abundance of bacterial genome bins using paired-end short reads.
```bash
conda create -n coverm -c bioconda -c conda-forge coverm
conda activate coverm
```
```bash
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
The output table contained mean coverage, relative abundance, and covered fraction for each bacterial bin.

---
# 8. 18S rRNA Phylogenetic Tree
A phylogenetic tree was generated from 18S rRNA sequences using Clustal Omega, TrimAl, and IQ-TREE 2.

```bash
cat *.fasta > 18S_new.fasta

clustalo \
    -i 18S_new.fasta \
    -o 18S_aligned.fasta

trimal \
    -in 18S_aligned.fasta \
    -out 18S_trimmed.fasta \
    -automated1

/home/ruchita.solanki/iqtree-2.2.2.7-Linux/bin/iqtree2 \
    -s 18S_trimmed.fasta \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -nt AUTO
```
The best-fit model was selected by IQ-TREE using ModelFinder, and branch support was estimated using ultrafast bootstrap and SH-aLRT support values.

---
# 9. Transcriptome Analysis
Transcriptome assembly and annotation were performed using the nf-core/metatdenovo workflow.
## 9.1 Java Setup
```bash
module purge
module load java/openjdk-23.0.1

export JAVA_HOME=$(dirname $(dirname $(readlink -f $(which java))))
export PATH=$JAVA_HOME/bin:$PATH
```
## 9.2 nf-core/metatdenovo Execution
```bash
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
The workflow generated transcript assemblies, predicted ORFs, functional annotations, and taxonomic classifications.

---
# 10. rRNA Gene Identification from Transcriptome
Barrnap was used to identify rRNA genes from the assembled transcriptome.

```bash
barrnap \
    --kingdom euk \
    --threads 4 \
    spades.transcripts.fa \
    --outseq euk_transcript_rRNA.fna \
    > diatom_euk_rRNA.gff

barrnap \
    --kingdom bac \
    spades.transcripts.fa \
    --outseq bac_transcript_rRNA.fna \
    > diatom_bac_rRNA.gff

barrnap \
    --kingdom mito \
    spades.transcripts.fa \
    --outseq mito_transcript_rRNA.fna \
    > diatom_mito_rRNA.gff
```
The resulting rRNA FASTA and GFF files were used to identify eukaryotic, bacterial, and mitochondrial rRNA transcripts.

---
# 11. Organelle Genome Identification
Organelle contigs were identified by comparing the polished assembly against reference mitochondrial and chloroplast genomes. The reference mitogenome and chloroplast genome used were:
```text
Mitogenome: MT742552
Chloroplast genome: MT742551
```
MetaQUAST was used to compare assembly contigs against organelle references.
```bash
conda create -n quast_env quast
conda activate quast_env
```
```bash
metaquast.py \
    8_diatom.fasta \
    -R /work/ebg_lab/eb/diatom_consortia/organelle/ref/ \
    -o ./8_metaquast_output
```
```bash
metaquast.py \
    /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_output/pypolca_corrected.fasta \
    -R /work/ebg_lab/eb/diatom_consortia/organelle/ref/ \
    -o ./whole_metaquast_output
```
The MetaQUAST output was used to identify candidate chloroplast and mitochondrial contigs for downstream organelle genome refinement and annotation.

---
# 12. Diatom Genome Annotation with BRAKER4
Gene models were generated using BRAKER4 with the genome assembly, RNA-seq alignments, and protein homology evidence. The genome assembly was not soft-masked prior to annotation; therefore, repeat masking was performed internally by the BRAKER4 workflow using RepeatModeler, RepeatMasker, and TRF. RNA-seq and protein evidence were supplied through `samples.csv`, which triggered BRAKER4 ETP mode using GeneMark-ETP, AUGUSTUS, and TSEBRA.
```text
Genome FASTA
   ↓
STAR genome indexing
   ↓
STAR RNA-seq alignment → coordinate-sorted BAM
   ↓
BRAKER4 internal repeat masking: RepeatModeler + RepeatMasker + TRF
   ↓
Protein homology evidence: OrthoDB Stramenopiles + Phaeodactylum tricornutum
   ↓
BRAKER4 ETP mode: GeneMark-ETP + AUGUSTUS + TSEBRA
   ↓
Gene models: GTF/GFF3
   ↓
Predicted proteins + CDS
   ↓
BUSCO/compleasm assessment
   ↓
Functional annotation and expression quantification
```
## 12.1 BRAKER4 Setup
BRAKER4 was cloned from GitHub, and the BRAKER3 Singularity image was downloaded for containerized execution.
```bash
git clone https://github.com/Gaius-Augustus/BRAKER4.git
cd BRAKER4
singularity pull braker3.sif docker://teambraker/braker3:latest
```
Snakemake was installed in the working environment:
```bash
conda install -c conda-forge -c bioconda snakemake
```
GeneMark requires a license key. The GeneMark key was stored at:
```bash
/home/ruchita.solanki/.gm_key
```
The GeneMark key path was specified in `config.ini`.
## 12.2 RNA-seq Alignment Evidence for BRAKER4
RNA-seq reads were aligned to the diatom genome assembly using STAR before running BRAKER4. BRAKER4 was therefore supplied with a precomputed coordinate-sorted BAM file rather than raw RNA-seq FASTQ files.
The genome assembly used for STAR indexing was:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
```
A STAR genome index was generated using the diatom genome assembly:
```bash
STAR \
    --runThreadN 24 \
    --runMode genomeGenerate \
    --genomeDir /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index \
    --genomeFastaFiles /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta \
    --genomeSAindexNbases 10
```
RNA-seq reads from the diatom consortium were aligned to the genome using STAR. The output was written directly as a coordinate-sorted BAM file:
```bash
STAR \
    --runThreadN 24 \
    --genomeDir /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index \
    --readFilesCommand zcat \
    --readFilesIn R1_rep1.fastq.gz,R1_rep2.fastq.gz,R1_rep3.fastq.gz \
                  R2_rep1.fastq.gz,R2_rep2.fastq.gz,R2_rep3.fastq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_
```
This generated the coordinate-sorted BAM file used as RNA-seq evidence for BRAKER4:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam
```
The BAM file was indexed:
```bash
samtools index /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam
```
The corresponding BAM index was:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam.bai
```
The BAM and BAM index were checked before BRAKER4 execution:
```bash
ls -lh /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam
ls -lh /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam.bai
```
The STAR BAM provides RNA-seq evidence for splice junctions, transcript coverage, and exon-intron boundaries during BRAKER4 gene prediction.
## 12.3 Protein Evidence Preparation
Protein evidence was prepared from broad stramenopile and diatom-related protein databases rather than from a single reference species. The final evidence set included OrthoDB Stramenopiles proteins and *Phaeodactylum tricornutum* proteins.
```bash
mkdir -p /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/proteins
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/proteins
```
### 12.3.1 Download OrthoDB Stramenopiles Proteins
```bash
wget -c https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Stramenopiles.fa.gz

gzip -t Stramenopiles.fa.gz

gunzip -c Stramenopiles.fa.gz > OrthoDB12_Stramenopiles.fa
```
The OrthoDB Stramenopiles FASTA file was checked:
```bash
grep -c "^>" OrthoDB12_Stramenopiles.fa
seqkit stats OrthoDB12_Stramenopiles.fa
```
### 12.3.2 Download *Phaeodactylum tricornutum* Proteins

```bash
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-63/fasta/phaeodactylum_tricornutum/pep/Phaeodactylum_tricornutum.ASM15095v2.pep.all.fa.gz

gzip -t Phaeodactylum_tricornutum.ASM15095v2.pep.all.fa.gz

gunzip -c Phaeodactylum_tricornutum.ASM15095v2.pep.all.fa.gz \
    > Phaeodactylum_tricornutum.ASM15095v2.pep.all.fa
```
The *Phaeodactylum* protein FASTA file was checked:
```bash
grep -c "^>" Phaeodactylum_tricornutum.ASM15095v2.pep.all.fa
seqkit stats Phaeodactylum_tricornutum.ASM15095v2.pep.all.fa
```
### 12.3.3 Combine and Clean Protein Evidence
The protein files were combined:
```bash
cat \
    OrthoDB12_Stramenopiles.fa \
    Phaeodactylum_tricornutum.ASM15095v2.pep.all.fa \
    > diatom_protein_evidence.raw.faa
```
The combined protein FASTA was filtered to remove short sequences and exact duplicate protein sequences:
```bash
seqkit seq -m 30 diatom_protein_evidence.raw.faa \
    | seqkit rmdup -s \
    > diatom_protein_evidence.clean.faa
```
The final protein evidence file was checked before use:
```bash
seqkit stats diatom_protein_evidence.raw.faa diatom_protein_evidence.clean.faa

grep -n "\*" diatom_protein_evidence.clean.faa | head

grep -n -v -E '^>|^[A-Z]+$' diatom_protein_evidence.clean.faa | head
```
No stop codons or malformed sequence lines were detected. The cleaned protein evidence file used for BRAKER4 was:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/proteins/diatom_protein_evidence.clean.faa
```
## 12.4 Genome Input and Masking Strategy
The genome assembly used for BRAKER4 was:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
```
The genome was checked for soft masking:
```bash
grep -v "^>" /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta \
    | grep -q '[a-z]' && echo "soft-masked" || echo "not soft-masked"
```
The output was:
```text
not soft-masked
```
Therefore, the `genome_masked` column in `samples.csv` was left empty, and internal repeat masking was enabled in `config.ini`. During the BRAKER4 run, the masking rule used RepeatModeler to identify repeat families, RepeatMasker to soft-mask repeat regions, and TRF to mask tandem repeats.
The masking step generated:
```bash
output/DL_diatom/preprocessing/genome.fa.masked
```
External pre-masking with RepeatModeler and RepeatMasker was not performed before BRAKER4; masking was handled inside the BRAKER4 workflow.
## 12.5 BRAKER4 Sample Configuration
A `samples.csv` file was prepared to specify the genome assembly, protein evidence, RNA-seq BAM file, and BUSCO lineage. The STAR alignment step was performed before BRAKER4, so the precomputed coordinate-sorted BAM file was supplied in the `bam_files` column. Providing both `protein_fasta` and `bam_files` triggered ETP mode in BRAKER4.
```bash
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4
nano samples.csv
```
```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage
DL_diatom,/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta,,/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/proteins/diatom_protein_evidence.clean.faa,/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam,,,,,,,,stramenopiles_odb12
```
The `samples.csv` file was checked to confirm the expected number of columns:
```bash
awk -F',' '{print NR, NF}' samples.csv
```
Expected output:
```text
1 13
2 13
```
## 12.6 BRAKER4 Configuration
The `config.ini` file was edited to specify the Singularity image, GeneMark license key, sample file, repeat masking option, and run parameters.
```bash
nano config.ini
```
```ini
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
The key settings were:
```ini
run_red = True
mode = etp
```
These were used because the genome assembly was not soft-masked and because the annotation was run with both RNA-seq and protein evidence.

## 12.7 Snakemake Dry Run
A Snakemake dry run was performed before launching the full BRAKER4 analysis:
```bash
snakemake \
    -s Snakefile \
    --use-singularity \
    --singularity-args "--bind /work/ebg_lab/eb/diatom_consortia/metatranscriptomics" \
    --cores 24 \
    --latency-wait 120 \
    --printshellcmds \
    -n
```
The dry run successfully built the workflow DAG and included the expected BRAKER4 rules:
```text
run_masking
run_genemark_etp
run_tsebra
run_augustus_hints
busco_genome
busco_proteins
collect_results
```
The presence of `run_genemark_etp` confirmed that BRAKER4 recognized the RNA-seq and protein evidence and configured the analysis in ETP mode. The presence of `run_tsebra` confirmed that TSEBRA refinement would be performed internally by BRAKER4. Therefore, TSEBRA does not need to be run separately after BRAKER4 completion.
## 12.8 BRAKER4 Execution
After the dry run completed successfully, BRAKER4 was executed using the official BRAKER4 `Snakefile`.
```bash
snakemake \
    -s Snakefile \
    --use-singularity \
    --singularity-args "--bind /work/ebg_lab/eb/diatom_consortia/metatranscriptomics" \
    --cores 24 \
    --latency-wait 120 \
    --printshellcmds \
    --rerun-incomplete
```
If the GeneMark key is not visible inside the container, the home directory can also be bound:
```bash
snakemake \
    -s Snakefile \
    --use-singularity \
    --singularity-args "--bind /work/ebg_lab/eb/diatom_consortia/metatranscriptomics,/home/ruchita.solanki" \
    --cores 24 \
    --latency-wait 120 \
    --printshellcmds \
    --rerun-incomplete
```
The following optional check can be used to confirm whether the GeneMark key is visible inside the container:
```bash
singularity exec \
    --bind /work/ebg_lab/eb/diatom_consortia/metatranscriptomics \
    braker3.sif \
    ls -lh /home/ruchita.solanki/.gm_key
```
BRAKER4 internally performs repeat masking, GeneMark-ETP training, AUGUSTUS training and prediction, evidence integration, TSEBRA refinement, BUSCO/compleasm assessment, and final result collection.
## 12.9 Extraction of Final Annotation Files
After BRAKER4 completed successfully, the final annotation files were extracted from the BRAKER4 results directory.
Output directory:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/output/DL_diatom/results
```
The compressed final files were extracted as follows:
```bash
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/output/DL_diatom/results

gunzip -c braker.gtf.gz > DL_diatom.braker4.gtf
gunzip -c braker.gff3.gz > DL_diatom.braker4.gff3
gunzip -c braker.aa.gz > DL_diatom.braker4.proteins.faa
gunzip -c braker.codingseq.gz > DL_diatom.braker4.cds.fna
```
The final files are:
```text
DL_diatom.braker4.gtf
DL_diatom.braker4.gff3
DL_diatom.braker4.proteins.faa
DL_diatom.braker4.cds.fna
```
No additional TSEBRA run was performed because the BRAKER4 workflow had already included TSEBRA refinement internally.
## 12.10 Basic Annotation Statistics
Basic annotation statistics were generated from the final GFF3, protein FASTA, and CDS FASTA files.
```bash
grep -c $'\tgene\t' DL_diatom.braker4.gff3
grep -c $'\tmRNA\t' DL_diatom.braker4.gff3
grep -c $'\tCDS\t' DL_diatom.braker4.gff3

grep -c "^>" DL_diatom.braker4.proteins.faa
grep -c "^>" DL_diatom.braker4.cds.fna

seqkit stats DL_diatom.braker4.proteins.faa
seqkit stats DL_diatom.braker4.cds.fna
```
Predicted proteins were checked for stop codons:
```bash
grep -n "\*" DL_diatom.braker4.proteins.faa | head
```
Very short predicted proteins were inspected:
```bash
seqkit fx2tab -n -l DL_diatom.braker4.proteins.faa \
    | awk '$2 < 50' \
    | head
```
## 12.11 BRAKER4 Reports and Evidence Support
BRAKER4 report files and quality-control summaries were located using:
```bash
find . -type f | grep -Ei "report|summary|busco|compleasm|statistics|support"
```
Relevant outputs include:
```text
braker_report.html
gene_support.tsv
BUSCO/compleasm summaries
gene set statistics
```
These files were used to assess annotation completeness, evidence support, and overall gene model quality.
## 12.12 Independent BUSCO Assessment
An independent BUSCO assessment was performed on the predicted protein set.
```bash
busco \
    -i DL_diatom.braker4.proteins.faa \
    -l stramenopiles_odb12 \
    -m proteins \
    -o busco_DL_diatom_braker4_proteins \
    -c 24
```
The BUSCO output was used to evaluate the completeness of the predicted protein set.
## 12.13 Longest Isoform Extraction
If BRAKER4 generated longest-isoform files, these were extracted for downstream functional annotation:
```bash
find . -type f | grep -Ei "longest|aa|codingseq"
```
If present, the longest-isoform protein and CDS files were extracted:
```bash
gunzip -c braker.longest.aa.gz > DL_diatom.braker4.longest.proteins.faa
gunzip -c braker.longest.codingseq.gz > DL_diatom.braker4.longest.cds.fna
```
The longest protein isoform file was used preferentially for downstream functional annotation:
```text
DL_diatom.braker4.longest.proteins.faa
```
If longest-isoform files were not available, the full BRAKER4 protein set was used:
```text
DL_diatom.braker4.proteins.faa
```
## 12.14 Functional Annotation
The final protein set was used for downstream functional annotation. Annotation tools may include:
```text
eggNOG-mapper
InterProScan
DIAMOND/BLASTP against UniProt or Swiss-Prot
KEGG/KO annotation
```
Functional annotation focused on pathways relevant to the diatom-dominated consortium, including photosynthesis, carbon-concentrating mechanisms, silica and frustule formation, nitrogen assimilation, lipid metabolism, vitamin and cofactor metabolism, stress responses, motility, and extracellular polymeric substance production.
## 12.15 Expression Quantification Against the Final Gene Models
After the final BRAKER4 annotation was accepted, RNA-seq expression was quantified against the final GTF file.
For gene-level counts, `featureCounts` can be used:
```bash
featureCounts \
    -T 24 \
    -p \
    -t exon \
    -g gene_id \
    -a DL_diatom.braker4.gtf \
    -o DL_diatom.braker4.featureCounts.txt \
    /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam
```
Alternatively, StringTie can be used to estimate transcript abundance:
```bash
stringtie \
    /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam \
    -G DL_diatom.braker4.gtf \
    -e \
    -B \
    -p 24 \
    -o DL_diatom.stringtie.gtf
```
The resulting expression estimates can be used with the functional annotation to summarize transcriptional activity across major diatom functional categories.
