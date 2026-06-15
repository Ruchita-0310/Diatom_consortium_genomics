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
After BRAKER4 ET annotation, the genome assembly used for annotation was further filtered to define a nuclear-enriched diatom genome. This step was performed after annotation because BRAKER4 had already been run on the broader diatom genome assembly, `18_diatom.fasta`.
The input genome for this step was:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
```
The recovered chloroplast and mitochondrial genomes were used as organelle references:
```bash
/work/ebg_lab/eb/diatom_consortia/organelle/2_chloro/chloroplast_contig_1443_trimmed.fasta
/work/ebg_lab/eb/diatom_consortia/organelle/mito/diatom_candidate_mitochondrion_2contigs.fasta
```
The mitochondrial reference consisted of two contigs with a combined length of 104,526 bp. Both contigs were retained in the mitochondrial reference because the goal was to identify mitochondrial-like contigs in the genome assembly, not to require a circularized mitochondrial genome.
The working definition used here was:
```text
nuclear-enriched genome = 18_diatom.fasta - chloroplast-like contigs - mitochondrial-like contigs
```
No BLAST-based taxonomic contaminant filtering was performed during this step. Filtering was limited to removal of contigs with strong similarity to the recovered chloroplast and mitochondrial genomes.
## 13.1 Input Files
A new working directory was created for nuclear genome filtering:
```bash
mkdir -p /work/ebg_lab/eb/diatom_consortia/nuclear_genome_filtering_18_diatom
cd /work/ebg_lab/eb/diatom_consortia/nuclear_genome_filtering_18_diatom
```
Input files were defined as follows:
```bash
WHOLE=/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta

CHLORO=/work/ebg_lab/eb/diatom_consortia/organelle/2_chloro/chloroplast_contig_1443_trimmed.fasta

MITO=/work/ebg_lab/eb/diatom_consortia/organelle/mito/diatom_candidate_mitochondrion_2contigs.fasta
```
Input assembly statistics were checked with `seqkit`:
```bash
seqkit stats $WHOLE $CHLORO $MITO
```
The input genome assembly contained 3,010 contigs and had a total length of 82,172,226 bp.
```text
18_diatom.fasta                                      3,010 contigs   82,172,226 bp
chloroplast_contig_1443_trimmed.fasta                   1 contig       120,429 bp
diatom_candidate_mitochondrion_2contigs.fasta            2 contigs      104,526 bp
```
## 13.2 Combine Organelle References
The chloroplast and mitochondrial FASTA files were combined into a single organelle reference file.
```bash
cat $CHLORO $MITO > organelles_chloro_mito.fasta
```
The combined organelle reference contained three sequences with a total length of 224,955 bp.
```bash
seqkit stats organelles_chloro_mito.fasta
```
## 13.3 Align `18_diatom.fasta` Against the Organelle References
The genome assembly was aligned against the combined organelle reference using minimap2.
```bash
minimap2 -x asm5 -c \
    organelles_chloro_mito.fasta \
    $WHOLE \
    > whole_vs_organelles.paf
```
The PAF output was used to identify regions of `18_diatom.fasta` with similarity to the chloroplast or mitochondrial genomes.
## 13.4 Calculate Organelle-aligned Coverage per Contig
Query intervals from the PAF file were extracted in BED format.
```bash
awk 'BEGIN{OFS="\t"} {
    print $1, $3, $4
}' whole_vs_organelles.paf \
> whole_vs_organelles.query_intervals.bed
```
Intervals were sorted and merged to avoid double-counting overlapping alignments.
```bash
sort -k1,1 -k2,2n whole_vs_organelles.query_intervals.bed \
> whole_vs_organelles.query_intervals.sorted.bed

bedtools merge \
    -i whole_vs_organelles.query_intervals.sorted.bed \
> whole_vs_organelles.query_intervals.merged.bed
```
Contig lengths were extracted from `18_diatom.fasta`.
```bash
seqkit fx2tab -n -l $WHOLE > whole_contig_lengths.tsv
```
The total organelle-aligned length was calculated for each contig.
```bash
awk 'BEGIN{OFS="\t"} {
    aligned[$1] += ($3 - $2)
}
END {
    for (c in aligned) print c, aligned[c]
}' whole_vs_organelles.query_intervals.merged.bed \
> organelle_aligned_length_per_contig.tsv
```
The percentage of each contig covered by organelle alignments was then calculated.
```bash
awk 'BEGIN{OFS="\t"}
FNR==NR {
    len[$1]=$2;
    next
}
{
    contig=$1;
    aligned=$2;
    pct=(aligned/len[contig])*100;
    print contig, len[contig], aligned, pct
}' whole_contig_lengths.tsv organelle_aligned_length_per_contig.tsv \
> organelle_coverage_per_contig.tsv
```
The strongest organelle-like contigs were inspected.
```bash
sort -k4,4nr organelle_coverage_per_contig.tsv | head -n 50
```
## 13.5 Identify Organelle-like Contigs
Contigs were classified as organelle-like if at least 70% of the contig length aligned to the chloroplast or mitochondrial reference.
```bash
awk '$4 >= 70 {print $1}' organelle_coverage_per_contig.tsv \
> organelle_like_contigs.70pct.txt
```
This identified three organelle-like contigs:
```text
contig_1443
contig_5628
contig_1647
```
These contigs corresponded to the chloroplast-like contig and the two mitochondrial-like contigs.
## 13.6 Generate the Nuclear-enriched Genome
The three organelle-like contigs were removed from `18_diatom.fasta` using `seqkit`.
```bash
seqkit grep \
    -v \
    -f organelle_like_contigs.70pct.txt \
    $WHOLE \
    > 18_diatom_nuclear_enriched.v1.fasta
```
Assembly statistics were calculated before and after organelle removal.
```bash
seqkit stats \
    $WHOLE \
    organelles_chloro_mito.fasta \
    18_diatom_nuclear_enriched.v1.fasta
```
The final nuclear-enriched genome contained 3,007 contigs and had a total length of 81,911,772 bp.
```text
18_diatom.fasta                         3,010 contigs   82,172,226 bp
organelles_chloro_mito.fasta                3 contigs      224,955 bp
18_diatom_nuclear_enriched.v1.fasta     3,007 contigs   81,911,772 bp
```
In total, three organelle-like contigs were removed, corresponding to 260,454 bp. This represented 0.317% of the `18_diatom.fasta` assembly.
## 13.7 Relationship to BRAKER4 Annotation
BRAKER4 was not rerun after organelle filtering because only three organelle-like contigs were removed from the BRAKER4 input assembly. These contigs represented 260,454 bp, corresponding to 0.317% of the 82,172,226 bp assembly. Because this filtering removed a very small fraction of the assembly and targeted organelle-derived contigs, the existing BRAKER4 annotation was retained.
For downstream analyses requiring a nuclear-only gene set, the BRAKER4 GFF3 can be filtered to retain only gene models located on contigs present in `18_diatom_nuclear_enriched.v1.fasta`.
## 13.8 Filter BRAKER4 GFF3 to Nuclear Contigs
A list of contigs present in the nuclear-enriched genome was generated.
```bash
seqkit seq -n 18_diatom_nuclear_enriched.v1.fasta > nuclear_contigs.v1.txt
```
The full BRAKER4 GFF3 was then filtered to retain only features located on contigs present in the nuclear-enriched genome.
```bash
MY_GFF=/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/output/DL_diatom/results/braker.gff3

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
}' nuclear_contigs.v1.txt $MY_GFF \
> braker.18_diatom_nuclear_only.v1.gff3
```
The resulting file represents the BRAKER4 annotation filtered to the final nuclear-enriched genome.
```text
braker.18_diatom_nuclear_only.v1.gff3
```
The number of retained nuclear gene models can be counted as follows:
```bash
grep -c $'\tgene\t' braker.18_diatom_nuclear_only.v1.gff3
```
## 13.9 Final Nuclear Genome Files
The final nuclear-enriched genome FASTA was:
```text
18_diatom_nuclear_enriched.v1.fasta
```
Full path:
```bash
/work/ebg_lab/eb/diatom_consortia/nuclear_genome_filtering_18_diatom/18_diatom_nuclear_enriched.v1.fasta
```
The corresponding nuclear-filtered BRAKER4 annotation was:
```text
braker.18_diatom_nuclear_only.v1.gff3
```
## 13.10 Summary
```text
18_diatom.fasta
   ↓
combine recovered chloroplast and mitochondrial references
   ↓
align 18_diatom.fasta against organelle references
   ↓
calculate organelle-aligned coverage per contig
   ↓
remove contigs with ≥70% organelle-aligned coverage
   ↓
generate 18_diatom_nuclear_enriched.v1.fasta
   ↓
retain existing BRAKER4 annotation
   ↓
filter BRAKER4 GFF3 to contigs present in the nuclear-enriched genome
```
The nuclear-enriched genome was generated from `18_diatom.fasta` by removing three organelle-like contigs: `contig_1443`, `contig_5628`, and `contig_1647`. The resulting assembly, `18_diatom_nuclear_enriched.v1.fasta`, contained 3,007 contigs and had a total length of 81,911,772 bp. BRAKER4 was not rerun because organelle filtering removed only 0.317% of the assembly. Instead, the existing BRAKER4 annotation was retained and filtered to match the final nuclear-enriched genome.
