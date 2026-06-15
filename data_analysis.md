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
# 12. Diatom Genome Annotation with BRAKER4 ET Mode
Gene models were generated using BRAKER4 with a diatom nuclear-enriched genome assembly and RNA-seq evidence. The genome assembly was not soft-masked before annotation; therefore, repeat masking was performed within the BRAKER4 workflow using RepeatModeler, RepeatMasker, and TRF. RNA-seq reads were aligned to the genome with STAR before BRAKER4, and the resulting coordinate-sorted BAM file was supplied to BRAKER4 as transcript evidence. The final annotation was run in ET mode, using RNA-seq evidence only.
```text
Genome FASTA
   ↓
STAR genome indexing
   ↓
STAR RNA-seq alignment
   ↓
Coordinate-sorted BAM
   ↓
BRAKER4 internal repeat masking
   ↓
BRAKER4 ET mode
   ↓
GeneMark-ET
   ↓
AUGUSTUS training and prediction
   ↓
TSEBRA refinement
   ↓
Final gene models: GTF/GFF3
   ↓
Predicted proteins and CDS
   ↓
BUSCO assessment
   ↓
Functional annotation and expression quantification
```
## 12.1 BRAKER4 Setup
BRAKER4 was cloned from GitHub, and the BRAKER container image was used for reproducible execution.
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
This path was provided in `config.ini`.
## 12.2 RNA-seq Alignment Evidence
RNA-seq reads from the diatom consortium were aligned to the nuclear-enriched genome assembly using STAR before BRAKER4. BRAKER4 was therefore supplied with a precomputed coordinate-sorted BAM file rather than raw RNA-seq FASTQ files.
The genome assembly used for STAR indexing and BRAKER4 was:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
```
A STAR genome index was generated as follows:
```bash
STAR \
    --runThreadN 24 \
    --runMode genomeGenerate \
    --genomeDir /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index \
    --genomeFastaFiles /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta \
    --genomeSAindexNbases 10
```
RNA-seq reads were aligned to the genome with STAR, and the output was written directly as a coordinate-sorted BAM file:
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
The coordinate-sorted BAM file used as RNA-seq evidence for BRAKER4 was:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam
```
The BAM file was indexed before BRAKER4 execution:
```bash
samtools index /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam
```
## 12.3 Genome Input and Repeat Masking
The genome assembly used for BRAKER4 was:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
```
The assembly contained 3,010 contigs and had a total length of approximately 82.17 Mbp.
```bash
seqkit stats /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
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
Because the genome was not pre-masked, the `genome_masked` column in `samples.csv` was left empty, and internal repeat masking was enabled in `config.ini`. BRAKER4 generated a masked genome using RepeatModeler, RepeatMasker, and TRF:
```bash
output/DL_diatom/preprocessing/genome.fa.masked
```
RepeatModeler identified 173 repeat families, and TRF identified 8,994 tandem repeat regions. External pre-masking was therefore not performed separately.
## 12.4 BRAKER4 Sample Configuration
The `samples.csv` file specified the genome assembly, RNA-seq BAM file, and BUSCO lineage. The `protein_fasta` column was left empty to force ET mode, meaning that the final annotation used RNA-seq evidence only and did not use protein evidence for GeneMark-ETP.
```bash
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4
nano samples.csv
```
```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage
DL_diatom,/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta,,,/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam,,,,,,,,stramenopiles_odb12
```
The file was checked to confirm that each row had the expected number of columns:
```bash
awk -F',' '{print NR, NF}' samples.csv
```
Expected output:
```text
1 13
2 13
```
## 12.5 BRAKER4 Configuration
The `config.ini` file specified the BRAKER container, GeneMark key, sample file, repeat masking option, and run parameters.
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
species = diatom_ET_v1
mode = et

[SLURM_ARGS]
cpus_per_task = 32
mem_of_node = 350000
max_runtime = 7200
```
The key settings were:
```ini
run_red = True
mode = et
```
`run_red = True` enabled internal repeat masking, and `mode = et` selected the RNA-seq-only BRAKER4 ET workflow.
## 12.6 Snakemake Dry Run
A Snakemake dry run was performed before launching the full BRAKER4 workflow:
```bash
snakemake \
    -s Snakefile \
    --use-singularity \
    --singularity-args "--bind /work/ebg_lab/eb/diatom_consortia/metatranscriptomics,/home/ruchita.solanki" \
    --cores 24 \
    --latency-wait 120 \
    --printshellcmds \
    --rerun-incomplete \
    -n
```
The dry run successfully built the workflow DAG and included the expected ET-mode rules:
```text
run_stringtie
bam2hints
run_genemark_et
train_augustus
run_augustus_hints
run_tsebra
busco_proteins
collect_results
```
The presence of `run_genemark_et` confirmed that BRAKER4 was configured in ET mode. The absence of `run_genemark_etp` confirmed that protein evidence was not used for GeneMark-ETP training.
## 12.7 BRAKER4 Execution
BRAKER4 was executed using the official BRAKER4 `Snakefile`:
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
BRAKER4 performed internal repeat masking, StringTie transcript reconstruction, RNA-seq hint generation, GeneMark-ET training and prediction, AUGUSTUS training and prediction, evidence integration, TSEBRA refinement, BUSCO assessment, and final result collection.
The ET-mode run produced 9,000 StringTie transcripts and 19,114 intron hints from the RNA-seq BAM file. These hints were used as transcript evidence during GeneMark-ET and AUGUSTUS prediction.
## 12.8 Final BRAKER4 Outputs
Final outputs were collected into:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET
```
The main final files were:
```text
DL_diatom.braker4.ET.gff3.gz
DL_diatom.braker4.ET.gtf.gz
DL_diatom.braker4.ET.proteins.faa.gz
DL_diatom.braker4.ET.cds.fna.gz
DL_diatom.braker4.ET.utr.gtf.gz
gene_support.tsv
software_versions.tsv
braker_report.html
braker_citations.bib
quality_control/
```
The final GFF3, protein FASTA, and CDS FASTA files can be decompressed as follows:
```bash
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET

gunzip -c DL_diatom.braker4.ET.gff3.gz > DL_diatom.braker4.ET.gff3
gunzip -c DL_diatom.braker4.ET.gtf.gz > DL_diatom.braker4.ET.gtf
gunzip -c DL_diatom.braker4.ET.proteins.faa.gz > DL_diatom.braker4.ET.proteins.faa
gunzip -c DL_diatom.braker4.ET.cds.fna.gz > DL_diatom.braker4.ET.cds.fna
```
No additional TSEBRA run is required because TSEBRA refinement was included within the BRAKER4 workflow.
## 12.9 Annotation Statistics
Basic annotation statistics were generated from the final GFF3, protein FASTA, and CDS FASTA files.
```bash
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET

grep -c $'\tgene\t' DL_diatom.braker4.ET.gff3
grep -c $'\ttranscript\t' DL_diatom.braker4.ET.gff3
grep -c $'\tCDS\t' DL_diatom.braker4.ET.gff3

grep -c "^>" DL_diatom.braker4.ET.proteins.faa
grep -c "^>" DL_diatom.braker4.ET.cds.fna

seqkit stats DL_diatom.braker4.ET.proteins.faa DL_diatom.braker4.ET.cds.fna
```
The final ET-mode annotation contained:
```text
Genes:        15,102
Transcripts:  16,947
Proteins:     16,947
CDS FASTA:    16,947
CDS features: 31,713
Exons:        31,713
Introns:      14,952
```
The final protein and CDS FASTA files had the following sequence counts:
```text
DL_diatom.braker4.ET.proteins.faa    16,947 proteins
DL_diatom.braker4.ET.cds.fna         16,947 CDS sequences
```
Predicted proteins were checked for internal stop codons:
```bash
grep -n "\*" DL_diatom.braker4.ET.proteins.faa | head
```
No internal stop codons were detected.
Very short predicted proteins can be inspected with:
```bash
seqkit fx2tab -n -l DL_diatom.braker4.ET.proteins.faa \
    | awk '$2 < 50' \
    | head
```
## 12.10 BUSCO Assessment of the Final Protein Set
The final predicted protein set was assessed with BUSCO v6 using the `stramenopiles_odb12` lineage dataset.
```bash
busco \
    -i DL_diatom.braker4.ET.proteins.faa \
    -l stramenopiles_odb12 \
    -m proteins \
    -o busco_DL_diatom_braker4_ET_proteins_odb12 \
    -c 24 \
    --download_path /work/ebg_lab/eb/diatom_consortia/databases/busco \
    --offline
```
The `stramenopiles_odb12` dataset contained 697 BUSCO markers. The final predicted protein set produced the following BUSCO result:
```text
C:84.8%[S:81.1%,D:3.7%],F:1.9%,M:13.3%,n=697
```
This corresponds to:
```text
Complete BUSCOs:              84.8%
Complete single-copy BUSCOs:  81.1%
Complete duplicated BUSCOs:    3.7%
Fragmented BUSCOs:             1.9%
Missing BUSCOs:               13.3%
Total BUSCO groups searched:   697
```
## 12.11 Annotation Acceptance
The final BRAKER4 ET annotation was accepted for downstream functional annotation because it produced a gene set of plausible size for the nuclear-enriched diatom assembly, showed no internal stop codon issues in the predicted protein FASTA, and recovered 84.8% of the `stramenopiles_odb12` BUSCO protein set with low duplication.
The final accepted annotation files are:
```text
DL_diatom.braker4.ET.gff3
DL_diatom.braker4.ET.gtf
DL_diatom.braker4.ET.proteins.faa
DL_diatom.braker4.ET.cds.fna
```
## 12.12 Functional Annotation
After accepting the BRAKER4 ET annotation, the final protein set can be used for downstream functional annotation.
```bash
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET
mkdir -p functional_annotation/eggnog
emapper.py \
    -i DL_diatom.braker4.ET.proteins.faa \
    --itype proteins \
    -m diamond \
    --cpu 24 \
    -o DL_diatom_braker4_ET \
    --output_dir functional_annotation/eggnog
```
Additional annotation tools may include:
```text
InterProScan
DIAMOND/BLASTP against UniProt or Swiss-Prot
KEGG/KO annotation
Pfam domain annotation
```
Functional annotation should focus on pathways relevant to the diatom-dominated consortium, including photosynthesis, carbon-concentrating mechanisms, silica and frustule formation, nitrogen assimilation, lipid metabolism, vitamin and cofactor metabolism, oxidative stress responses, motility, and extracellular polymeric substance production.
## 12.13 Expression Quantification Against the Final Gene Models
After accepting the final BRAKER4 annotation, RNA-seq expression can be quantified against the final GTF file.
For gene-level counts, `featureCounts` can be used:
```bash
featureCounts \
    -T 24 \
    -p \
    -t exon \
    -g gene_id \
    -a DL_diatom.braker4.ET.gtf \
    -o DL_diatom.braker4.ET.featureCounts.txt \
    /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam
```
Alternatively, StringTie can be used to estimate transcript abundance:
```bash
stringtie \
    /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam \
    -G DL_diatom.braker4.ET.gtf \
    -e \
    -B \
    -p 24 \
    -o DL_diatom.braker4.ET.stringtie.gtf
```
The resulting expression estimates can be combined with the functional annotation table to summarize transcriptional activity across major diatom functional categories.
## 12.14 Rationale for Using ET Mode Instead of ETP
BRAKER4 was initially tested in ETP mode, which combines RNA-seq evidence with protein evidence. However, the GeneMark-ETP step failed during model training. Although protein-supported alignments were generated, the GeneMark-ETP training set did not produce valid gene and transcript models for training. The failed run reported zero parsed genes and transcripts, CDS-only training entries, and phase-distribution errors, followed by division-by-zero errors in the GeneMark-ETP training scripts:
```text
genes: 0
transcripts: 0
CDS: 1724

Use of uninitialized value $ph1 in addition (+) at /opt/ETP/bin/gmes/parse_set.pl line 205.
Use of uninitialized value $ph0 in division (/) at /opt/ETP/bin/gmes/parse_set.pl line 208.
Illegal division by zero at /opt/ETP/bin/gmes/parse_set.pl line 208.
Illegal division by zero at /opt/ETP/bin/train_super.pl line 184.
ERROR: GeneMark-ETP failed, no genemark.gtf
```
Because GeneMark-ETP did not complete successfully, the annotation was rerun in ET mode using RNA-seq evidence only. This avoided the failed protein-dependent GeneMark-ETP training step while retaining transcript evidence from the coordinate-sorted STAR BAM file. The final successful workflow therefore used GeneMark-ET, AUGUSTUS, and TSEBRA, with `protein_fasta` left empty in `samples.csv` and `mode = et` specified in `config.ini`.
The ET-mode dry run confirmed the expected workflow by including `run_genemark_et` and excluding `run_genemark_etp`. The final ET-mode annotation completed successfully and was therefore used for downstream analysis.

---
# 13. Nuclear Genome Identification
After the initial BRAKER4 ET annotation, the nuclear genome was defined more explicitly from the polished whole-genome assembly. At this stage, three genome components had been recovered from the diatom-associated consortium:
```text
Whole draft genome assembly: ~84 Mbp
Chloroplast genome:          ~120 kbp
Mitogenome:                  ~38 kbp
```
The working definition used here was:
```text
nuclear-enriched genome = whole draft genome assembly - chloroplast-derived contigs - mitochondrial-derived contigs
```
No BLAST-based taxonomic filtering was performed at this stage. Instead, the nuclear genome was identified by subtractive removal of chloroplast-like and mitochondrial-like contigs from the polished whole-genome assembly.
## 13.1 Rationale
The polished whole-genome assembly contained the diatom nuclear genome as well as organelle-derived sequences. Because the chloroplast and mitochondrial genomes had already been identified, these organelle assemblies were used as internal references to identify and remove organelle-like contigs from the whole draft assembly.
The remaining contigs were retained as the nuclear-enriched diatom genome assembly. This assembly was then used as the basis for downstream genome interpretation and could be used for reannotation, comparative genomics, and functional analysis.
## 13.2 Organelle Genome References
The recovered chloroplast and mitochondrial genomes were used as references for subtractive filtering.
```text
Chloroplast genome: ~120 kbp
Mitogenome:         ~38 kbp
```
These organelle genomes were combined into a single reference FASTA file.

```bash
cat chloroplast_120kb.fasta mitogenome_38kb.fasta > organelles.fasta
```
## 13.3 Identify Organelle-like Contigs in the Whole Assembly
The polished whole-genome assembly was aligned against the combined organelle reference.
```bash
minimap2 -x asm5 \
    organelles.fasta \
    whole_genome_84M.fasta \
    > whole_vs_organelle.paf
```
Contigs were flagged as organelle-like if most of the contig aligned to the chloroplast or mitochondrial genome. A 70% contig-coverage threshold was used to identify contigs with strong organelle similarity.

```bash
awk '{
  q=$1; qlen=$2; aln=$11;
  cov[q]+=aln; len[q]=qlen
}
END {
  for (q in cov) {
    if (cov[q]/len[q] >= 0.70) print q
  }
}' whole_vs_organelle.paf > organelle_like_contigs.txt
```
The resulting file contained contig IDs classified as chloroplast-like or mitochondrial-like.
```text
organelle_like_contigs.txt
```
## 13.4 Remove Organelle-like Contigs
Organelle-like contigs were removed from the polished whole-genome assembly using `seqkit`.
```bash
seqkit grep \
    -v \
    -f organelle_like_contigs.txt \
    whole_genome_84M.fasta \
    > diatom_nuclear_genome.fasta
```
The output FASTA represented the nuclear-enriched diatom genome assembly.
```text
diatom_nuclear_genome.fasta
```
## 13.5 Assembly Statistics
Assembly statistics were calculated before and after organelle removal.
```bash
seqkit stats \
    whole_genome_84M.fasta \
    diatom_nuclear_genome.fasta
```

This comparison was used to confirm the change in assembly size after removing chloroplast-like and mitochondrial-like contigs.
## 13.6 BUSCO Assessment
BUSCO was used to assess the completeness of the nuclear-enriched genome assembly.

```bash
busco \
    -i diatom_nuclear_genome.fasta \
    -l stramenopiles_odb10 \
    -m genome \
    -o busco_diatom_nuclear_genome \
    -c 32
```
The BUSCO result was used to evaluate whether the organelle-filtered assembly retained conserved stramenopile genes expected from a diatom nuclear genome.
## 13.7 Transcriptome Context
An assembled transcriptome was also available for the diatom consortium. This transcriptome provided additional biological context for downstream annotation and interpretation of the nuclear-enriched assembly. However, transcriptome support was not used as a taxonomic filtering step during nuclear genome identification.
## 13.8 Final Nuclear Genome Definition
The final nuclear genome was defined as the polished whole-genome assembly after removal of contigs with strong similarity to the recovered chloroplast and mitochondrial genomes.
```text
Final nuclear-enriched genome:
diatom_nuclear_genome.fasta
```
This file represents the diatom nuclear-enriched genome assembly used for downstream interpretation.
## 13.9 Summary of Nuclear Genome Identification
```text
Polished whole-genome assembly
   ↓
Recovered chloroplast and mitochondrial genomes
   ↓
Combined organelle reference FASTA
   ↓
Aligned whole assembly against organelle genomes
   ↓
Identified chloroplast-like and mitochondrial-like contigs
   ↓
Removed organelle-like contigs
   ↓
Generated nuclear-enriched diatom genome FASTA
   ↓
Assessed assembly statistics and BUSCO completeness
```
In this study, the nuclear-enriched diatom genome was identified after the initial BRAKER4 annotation by subtractive filtering of the polished whole-genome assembly. The final genome was defined by removing contigs with strong similarity to the recovered chloroplast and mitochondrial genomes. No BLAST-based taxonomic filtering was performed at this stage.
## 13.8 Final Nuclear Genome Definition
The final nuclear genome was defined as the polished whole-genome assembly after removal of contigs with strong similarity to the recovered chloroplast and mitochondrial genomes.
```text
Final nuclear-enriched genome:
diatom_candidate.no_organelle.v2.fasta
```
Full path:
```bash
/work/ebg_lab/eb/diatom_consortia/nuclear_genome_filtering/diatom_candidate.no_organelle.v2.fasta
```
This file represents the organelle-filtered, nuclear-enriched diatom genome used for downstream comparative analysis.

---
# 14. Pairwise Genome Comparison with *Phaeodactylum tricornutum*
A pairwise genome comparison was performed between the candidate diatom nuclear genome and the reference genome of *Phaeodactylum tricornutum*. This analysis was used as a nucleotide-level similarity screen to identify regions of the candidate diatom genome with detectable similarity to a well-studied reference diatom genome.
This analysis should not be interpreted as a complete gene orthology analysis. BLASTN detects conserved nucleotide regions, but many homologous genes may be too diverged at the nucleotide level to be recovered. Protein-level comparison using BLASTP, DIAMOND, or OrthoFinder should be used for a more complete assessment of conserved gene content.
## 14.1 Input Files
The nuclear-enriched candidate diatom genome was used as the BLAST query:
```bash
MY_GENOME=/work/ebg_lab/eb/diatom_consortia/nuclear_genome_filtering/diatom_candidate.no_organelle.v2.fasta
```
The BRAKER4 GFF3 annotation was used only after BLASTN to identify which predicted genes overlapped *P. tricornutum*-like regions:
```bash
MY_GFF=/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/output/DL_diatom/results/braker.gff3
```
The BRAKER4 GFF3 file was generated from the broader genome assembly and therefore contained annotations on nuclear and non-nuclear contigs. Because the BLAST query was the nuclear-only genome FASTA, the GFF3 file was first filtered to retain only features present on contigs in the nuclear-enriched genome.
## 14.2 Prepare Working Directory
```bash
mkdir -p /work/ebg_lab/eb/diatom_consortia/phaeodactylum_blast
cd /work/ebg_lab/eb/diatom_consortia/phaeodactylum_blast
mkdir -p blast_out blast_db phaeodactylum
```
## 14.3 Filter the BRAKER4 GFF3 to Nuclear Contigs
A list of nuclear contigs was generated from the organelle-filtered genome FASTA:
```bash
seqkit seq -n $MY_GENOME > nuclear_contigs.txt
```
The full BRAKER4 GFF3 was then filtered to retain only rows whose contig ID matched the nuclear contig list:
```bash
awk 'BEGIN{FS=OFS="\t"}
FNR==NR {
    keep[$1]=1;
    next
}
$0 ~ /^#/ {
    if ($0 !~ /^##FASTA/) print;
    next
}
$1 in keep {
    print
}' nuclear_contigs.txt $MY_GFF \
> braker.nuclear_only.gff3
```
Feature counts were inspected:
```bash
grep -v "^#" braker.nuclear_only.gff3 \
    | cut -f3 \
    | sort \
    | uniq -c \
    | sort -nr \
    | head
```
The nuclear-only GFF3 contained 14,996 predicted gene features.
## 14.4 Convert Nuclear Gene Models to BED
The nuclear-only GFF3 was converted to BED format for overlap analysis with BLASTN alignments:
```bash
awk -F'\t' 'BEGIN{OFS="\t"} 
!/^#/ && $3=="gene" {
  id=$9;
  sub(/.*ID=/,"",id);
  sub(/;.*/,"",id);
  print $1,$4-1,$5,id,$6,$7
}' braker.nuclear_only.gff3 \
> blast_out/diatom_candidate_nuclear_genes.bed
```
The BED file contained 14,996 nuclear gene models.
## 14.5 Download the *Phaeodactylum tricornutum* Reference Genome
The *P. tricornutum* RefSeq genome assembly `GCF_000150955.2` was downloaded using NCBI Datasets:
```bash
datasets download genome accession GCF_000150955.2 \
  --include genome,gff3 \
  --filename phaeodactylum/Phaeodactylum_tricornutum_GCF_000150955.2.zip

unzip phaeodactylum/Phaeodactylum_tricornutum_GCF_000150955.2.zip -d phaeodactylum
```
The genome FASTA and GFF3 files were identified:
```bash
PT_GENOME=$(find phaeodactylum -name "*genomic.fna" | head -n 1)
PT_GFF=$(find phaeodactylum -name "*.gff" -o -name "*.gff3" | head -n 1)
echo $PT_GENOME
echo $PT_GFF
```
## 14.6 Build the *Phaeodactylum* BLAST Database
```bash
makeblastdb \
  -in $PT_GENOME \
  -dbtype nucl \
  -parse_seqids \
  -out blast_db/Phaeodactylum_tricornutum
```
## 14.7 Run Genome-vs-Genome BLASTN
The candidate diatom nuclear genome was aligned against the *P. tricornutum* genome using BLASTN:
```bash
blastn \
  -query $MY_GENOME \
  -db blast_db/Phaeodactylum_tricornutum \
  -task blastn \
  -evalue 1e-10 \
  -num_threads 32 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
  -out blast_out/diatom_candidate_vs_phaeodactylum.blastn.tsv
```
The BLASTN output was filtered using the following thresholds:
```text
Percent identity:     ≥70%
Alignment length:     ≥200 bp
E-value:              ≤1e-10
```
```bash
awk 'BEGIN{OFS="\t"} $3 >= 70 && $4 >= 200 && $11 <= 1e-10' \
  blast_out/diatom_candidate_vs_phaeodactylum.blastn.tsv \
  > blast_out/diatom_candidate_vs_phaeodactylum.filtered.tsv
```
This produced 5,246 filtered BLASTN alignments.
```bash
wc -l blast_out/diatom_candidate_vs_phaeodactylum.filtered.tsv
```
Output:
```text
5246 blast_out/diatom_candidate_vs_phaeodactylum.filtered.tsv
```
## 14.8 Convert BLASTN Hits to Candidate Diatom Genome Coordinates
Filtered BLASTN hits were converted to BED format using coordinates on the candidate diatom nuclear genome:
```bash
awk 'BEGIN{OFS="\t"} {
  qs=($7<$8?$7:$8)-1;
  qe=($7>$8?$7:$8);
  strand=($9<=$10?"+":"-");
  hit="hit_"NR;
  print $1, qs, qe, hit, $12, strand, $2, $9, $10, $3, $4, $11
}' blast_out/diatom_candidate_vs_phaeodactylum.filtered.tsv \
> blast_out/blast_hits_on_diatom_candidate.bed
```
## 14.9 Identify Diatom Genes Overlapping *Phaeodactylum*-like Regions
Filtered BLASTN alignments were intersected with nuclear BRAKER4 gene models:
```bash
bedtools intersect \
  -a blast_out/diatom_candidate_nuclear_genes.bed \
  -b blast_out/blast_hits_on_diatom_candidate.bed \
  -wa -wb \
  > blast_out/diatom_candidate_nuclear_genes_with_phaeodactylum_hits.tsv
```
The intersect produced 2,257 gene-alignment overlaps.
```bash
wc -l blast_out/diatom_candidate_nuclear_genes_with_phaeodactylum_hits.tsv
```
Output:
```text
2257 blast_out/diatom_candidate_nuclear_genes_with_phaeodactylum_hits.tsv
```
The number of unique diatom genes overlapping *P. tricornutum*-like regions was calculated:
```bash
cut -f4 blast_out/diatom_candidate_nuclear_genes_with_phaeodactylum_hits.tsv \
  | sort -u \
  | wc -l
```
Output:
```text
1487
```
The proportion of nuclear genes with detectable *P. tricornutum*-like nucleotide similarity was calculated:
```bash
TOTAL_GENES=$(wc -l < blast_out/diatom_candidate_nuclear_genes.bed)

HIT_GENES=$(cut -f4 blast_out/diatom_candidate_nuclear_genes_with_phaeodactylum_hits.tsv \
  | sort -u \
  | wc -l)

awk -v h=$HIT_GENES -v t=$TOTAL_GENES 'BEGIN{
  print "Total nuclear genes:", t
  print "Genes with Phaeodactylum-like BLAST hits:", h
  print "Percent:", (h/t)*100 "%"
}'
```
Output:
```text
Total nuclear genes: 14996
Genes with Phaeodactylum-like BLAST hits: 1487
Percent: 9.91598%
```
Thus, 1,487 of 14,996 predicted nuclear genes, corresponding to 9.92% of the nuclear gene set, overlapped filtered BLASTN alignments to *P. tricornutum*.
## 14.10 Select the Best BLASTN Hit per Diatom Gene
Because individual genes can overlap more than one BLASTN alignment, one representative best hit was selected for each diatom gene using the highest BLAST bitscore:
```bash
sort -k4,4 -k11,11nr blast_out/diatom_candidate_nuclear_genes_with_phaeodactylum_hits.tsv \
  | awk 'BEGIN{OFS="\t"} !seen[$4]++' \
  > blast_out/diatom_candidate_nuclear_genes_best_phaeodactylum_hit.tsv
```
The best-hit table contained 1,487 diatom genes.
```bash
wc -l blast_out/diatom_candidate_nuclear_genes_best_phaeodactylum_hit.tsv
```
Mean nucleotide identity and mean alignment length were calculated for the best-hit set:
```bash
awk 'BEGIN{FS=OFS="\t"}
{
  pident_sum += $16;
  aln_len_sum += $17;
  n += 1
}
END{
  print "Number of best-hit genes:", n
  print "Mean percent identity:", pident_sum/n
  print "Mean alignment length:", aln_len_sum/n
}' blast_out/diatom_candidate_nuclear_genes_best_phaeodactylum_hit.tsv
```
Output:
```text
Number of best-hit genes: 1487
Mean percent identity: 72.8159
Mean alignment length: 598.428
```
## 14.11 Create a Clean Diatom Best-hit Table
A simplified table was generated for downstream inspection:
```bash
awk 'BEGIN{FS=OFS="\t";
print "diatom_contig","diatom_gene_start","diatom_gene_end","diatom_gene_id","diatom_gene_strand","blast_hit_id","bitscore","blast_strand","pt_contig","pt_start","pt_end","pident","aln_len","evalue"
}
{
print $1,$2,$3,$4,$6,$10,$11,$12,$13,$14,$15,$16,$17,$18
}' blast_out/diatom_candidate_nuclear_genes_best_phaeodactylum_hit.tsv \
> blast_out/diatom_best_hits_to_phaeodactylum.clean.tsv
```
The resulting table links each candidate diatom gene to its best BLASTN hit, including *P. tricornutum* genomic coordinates, percent identity, alignment length, and e-value.
## 14.12 Map BLASTN Hits to Annotated *Phaeodactylum* Genes
To identify which *P. tricornutum* genes overlapped the BLASTN regions, the filtered BLASTN alignments were also converted to BED format using coordinates on the *P. tricornutum* genome:
```bash
awk 'BEGIN{OFS="\t"} {
  ss=($9<$10?$9:$10)-1;
  se=($9>$10?$9:$10);
  strand=($9<=$10?"+":"-");
  hit="hit_"NR;
  print $2, ss, se, hit, $12, strand, $1, $7, $8, $3, $4, $11
}' blast_out/diatom_candidate_vs_phaeodactylum.filtered.tsv \
> blast_out/blast_hits_on_phaeodactylum.bed
```
The *P. tricornutum* GFF3 was converted to BED format:
```bash
awk -F'\t' 'BEGIN{OFS="\t"}
function get_attr(attr,key,   n,a,i,b) {
  n=split(attr,a,";");
  for(i=1;i<=n;i++) {
    split(a[i],b,"=");
    if(b[1]==key) return b[2];
  }
  return "NA";
}
!/^#/ && $3=="gene" {
  id=get_attr($9,"ID");
  name=get_attr($9,"Name");
  gene=get_attr($9,"gene");
  locus=get_attr($9,"locus_tag");

  if(name=="NA" && gene!="NA") name=gene;
  if(name=="NA" && locus!="NA") name=locus;

  print $1,$4-1,$5,id,$6,$7,name,locus;
}' $PT_GFF \
> blast_out/phaeodactylum_genes.bed
```
The *P. tricornutum* GFF3 used contig names such as:
```text
NC_011669.1
```
whereas the BLASTN subject IDs were formatted as:
```text
ref|NC_011669.1|
```
Because `bedtools intersect` requires exact sequence ID matches, the BLASTN BED file was normalized by removing the `ref|` prefix and trailing pipe symbol:
```bash
awk 'BEGIN{OFS="\t"}
{
  seq=$1;
  if (seq ~ /\|/) {
    split(seq,a,"|");
    seq=a[2];
  }
  $1=seq;
  print
}' blast_out/blast_hits_on_phaeodactylum.bed \
> blast_out/blast_hits_on_phaeodactylum.normalized.bed
```
The normalized BLASTN intervals were intersected with the *P. tricornutum* gene BED file:
```bash
bedtools intersect \
  -a blast_out/blast_hits_on_phaeodactylum.normalized.bed \
  -b blast_out/phaeodactylum_genes.bed \
  -wa -wb \
  > blast_out/blast_hits_overlapping_phaeodactylum_genes.tsv
```
A BLAST hit to *P. tricornutum* gene mapping table was generated:
```bash
awk 'BEGIN{FS=OFS="\t";
print "blast_hit_id","pt_gene_id","pt_gene_name","pt_locus_tag"
}
{
  print $4,$16,$19,$20
}' blast_out/blast_hits_overlapping_phaeodactylum_genes.tsv \
> blast_out/blast_hit_to_phaeodactylum_gene.tsv
```
This produced 2,856 lines, corresponding to one header line and 2,855 BLAST hit to *P. tricornutum* gene overlaps.
```bash
wc -l blast_out/blast_hit_to_phaeodactylum_gene.tsv
```
Output:

```text
2856 blast_out/blast_hit_to_phaeodactylum_gene.tsv
```
## 14.13 Merge *Phaeodactylum* Gene IDs into the Diatom Best-hit Table
The *P. tricornutum* gene mapping table was merged with the clean diatom best-hit table using `blast_hit_id`:
```bash
python3 <<'PY'
import csv
from collections import defaultdict

clean_file = "blast_out/diatom_best_hits_to_phaeodactylum.clean.tsv"
map_file = "blast_out/blast_hit_to_phaeodactylum_gene.tsv"
out_file = "blast_out/diatom_best_hits_to_phaeodactylum.with_PT_genes.tsv"

hit_to_pt = defaultdict(lambda: {
    "pt_gene_id": set(),
    "pt_gene_name": set(),
    "pt_locus_tag": set()
})

with open(map_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        hit = row["blast_hit_id"]
        for key in ["pt_gene_id", "pt_gene_name", "pt_locus_tag"]:
            val = row[key]
            if val and val != "NA":
                hit_to_pt[hit][key].add(val)

with open(clean_file) as fin, open(out_file, "w") as fout:
    reader = csv.DictReader(fin, delimiter="\t")
    fieldnames = reader.fieldnames + ["pt_gene_id", "pt_gene_name", "pt_locus_tag"]
    writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    for row in reader:
        hit = row["blast_hit_id"]
        row["pt_gene_id"] = ";".join(sorted(hit_to_pt[hit]["pt_gene_id"])) or "NA"
        row["pt_gene_name"] = ";".join(sorted(hit_to_pt[hit]["pt_gene_name"])) or "NA"
        row["pt_locus_tag"] = ";".join(sorted(hit_to_pt[hit]["pt_locus_tag"])) or "NA"
        writer.writerow(row)

print("Wrote:", out_file)
PY
```
The merged output contained 1,488 lines, corresponding to one header line and 1,487 diatom genes with best BLASTN hits.

```bash
wc -l blast_out/diatom_best_hits_to_phaeodactylum.with_PT_genes.tsv
```
Output:
```text
1488 blast_out/diatom_best_hits_to_phaeodactylum.with_PT_genes.tsv
```
A final version with cleaned *P. tricornutum* gene IDs was generated by removing the `gene-` prefix from PHATRDRAFT identifiers:
```bash
awk 'BEGIN{FS=OFS="\t"}
NR==1 {print; next}
{
  gsub(/gene-/, "", $15);
  print
}' blast_out/diatom_best_hits_to_phaeodactylum.with_PT_genes.tsv \
> blast_out/diatom_best_hits_to_phaeodactylum.with_PT_genes.cleanIDs.tsv
```
## 14.14 Final Output Files
```text
blast_out/diatom_candidate_vs_phaeodactylum.blastn.tsv
```
Raw BLASTN output from the candidate diatom nuclear genome against the *P. tricornutum* genome.
```text
blast_out/diatom_candidate_vs_phaeodactylum.filtered.tsv
```
Filtered BLASTN alignments using ≥70% identity, alignment length ≥200 bp, and e-value ≤1e-10.
```text
blast_out/diatom_candidate_nuclear_genes_with_phaeodactylum_hits.tsv
```
All overlaps between nuclear diatom gene models and filtered BLASTN alignments.
```text
blast_out/diatom_candidate_nuclear_genes_best_phaeodactylum_hit.tsv
```
One best BLASTN hit per diatom gene, selected by highest BLAST bitscore.
```text
blast_out/diatom_best_hits_to_phaeodactylum.clean.tsv
```
Clean table linking each diatom gene to its best BLASTN hit and *P. tricornutum* genomic coordinates.
```text
blast_out/blast_hit_to_phaeodactylum_gene.tsv
```
Mapping file connecting BLAST hit IDs to overlapping *P. tricornutum* PHATRDRAFT gene models.
```text
blast_out/diatom_best_hits_to_phaeodactylum.with_PT_genes.cleanIDs.tsv
```
Final interpreted table linking each candidate diatom gene with its best BLASTN hit, *P. tricornutum* genomic coordinates, nucleotide identity, alignment length, e-value, and overlapping PHATRDRAFT gene model.
## 14.15 Summary
Pairwise BLASTN comparison identified 5,246 filtered nucleotide alignments between the candidate diatom nuclear genome and the *Phaeodactylum tricornutum* reference genome. These alignments overlapped 1,487 of 14,996 predicted nuclear genes, corresponding to 9.92% of the nuclear gene set. For the best hit per gene, the mean nucleotide identity was 72.8%, and the mean alignment length was 598 bp.
The best-hit regions were further mapped to annotated *P. tricornutum* PHATRDRAFT gene models. This produced a final gene-level table connecting candidate diatom genes to their best nucleotide-level *P. tricornutum* matches and the corresponding reference gene annotations.
This analysis provides a conservative nucleotide-level comparison between the candidate diatom genome and *P. tricornutum*. Because nucleotide similarity can be lost despite conservation at the protein level, this result should be treated as a genome-level similarity screen. A protein-level comparison using BLASTP, DIAMOND, or OrthoFinder is recommended for orthology and conserved gene content analysis.
