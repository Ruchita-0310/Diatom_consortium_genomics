# Diatom Consortia: Metagenomic and Metatranscriptomic Pipeline
This repository documents the workflow used to assemble, polish, bin, classify, annotate, and compare genomes and transcriptomes from a diatom-associated microbial consortium. The pipeline combines long-read metagenomic assembly, short-read polishing, metagenomic binning, contig-level taxonomic screening, organelle identification, transcriptome analysis, BRAKER4 ET gene prediction, nuclear-enriched genome generation, functional annotation, expression integration, and comparison with the reference diatom *Phaeodactylum tricornutum*. The final gene table is a clean BRAKER4 isoform-level table with one row per predicted protein isoform.
## Workflow overview
```text
Nanopore reads
   ↓
Flye metagenome assembly
   ↓
Medaka long-read polishing
   ↓
Polypolish + Pypolca short-read polishing
   ↓
Assembly assessment and read mapping
   ↓
MetaBAT2 binning
   ↓
CheckM2 + GTDB-Tk + MetaEuk contig classification
   ↓
Organelle identification using reference chloroplast and mitochondrial genomes
   ↓
BRAKER4 ET gene annotation using RNA-seq evidence
   ↓
Nuclear-enriched genome generation
   ↓
Functional annotation with Swiss-Prot, Bacillariophyta UniProtKB, InterProScan, and AntiFam
   ↓
Expression integration using best TransDecoder ORF-to-BRAKER4 mappings and Average_TPM only
   ↓
Phaeodactylum tricornutum comparison summarized as yes/no only
   ↓
Final clean BRAKER4 isoform-level gene table for pathway curation
```
## Software and environments
The workflow used Conda environments, Singularity containers, and local HPC modules depending on software availability.
| Step | Tools |
|---|---|
| Assembly and polishing | Flye, Medaka, BWA-MEM, Polypolish, Pypolca |
| Read mapping and coverage | minimap2, samtools, bedtools, seqkit |
| Assembly quality | BUSCO, QUAST/MetaQUAST |
| Binning and bin quality | MetaBAT2, CheckM2 |
| Taxonomy and abundance | GTDB-Tk, MetaEuk, CoverM |
| Organelle identification | MetaQUAST, minimap2, bedtools, seqkit, GeSeq/OGDRAW |
| Phylogenetics | Clustal Omega, TrimAl, IQ-TREE 2 |
| Transcriptomics | Nextflow, nf-core/metatdenovo, TransDecoder, Barrnap, STAR |
| Genome annotation | BRAKER4, GeneMark-ET, AUGUSTUS, TSEBRA, STAR, BUSCO/compleasm |
| Functional annotation | DIAMOND, UniProtKB/Swiss-Prot, UniProtKB Bacillariophyta, InterProScan, Pfam, PANTHER, Gene3D, CDD, SMART, SUPERFAMILY, ProSite, Python |
| Expression integration | DIAMOND, Python, pandas, TransDecoder ORFs, Average_TPM table |
| Comparative genomics | NCBI Datasets, BLASTN, bedtools, Python |
## Repository structure for scripts
Custom Python and SLURM scripts are stored in the `scripts/` directory rather than embedded directly in this markdown workflow.

```text
scripts/
├── classify_metaeuk_contigs.py
├── make_swissprot_best_hits.py
├── make_bacillariophyta_best_hits.py
├── summarize_interproscan.py
├── merge_functional_annotation_layers.py
├── run_diamond_transdecoder_vs_braker_clean.slurm
├── 07_add_BRAKER_lengths_clean.py
├── 08_make_best_ORF_to_BRAKER_mapping_clean.py
├── 09_add_ONLY_Average_TPM_clean.py
└── 10_make_FINAL_clean_BRAKER_isoform_table.py
```

Script purposes:

```text
classify_metaeuk_contigs.py
  Classifies contigs using MetaEuk ORF-level taxonomy and a priority-based contig classification scheme.

make_swissprot_best_hits.py
  Parses Swiss-Prot DIAMOND output, calculates query and subject coverage, assigns confidence classes, and writes best-hit annotation tables.

make_bacillariophyta_best_hits.py
  Parses UniProtKB Bacillariophyta DIAMOND output using the same confidence framework as the Swiss-Prot parser.

summarize_interproscan.py
  Collapses raw InterProScan TSV output to one row per predicted protein and summarizes database contributions, GO terms, and pathway annotations.

merge_functional_annotation_layers.py
  Merges Swiss-Prot, Bacillariophyta, InterProScan, and AntiFam evidence into a master BRAKER4 functional annotation table.

run_diamond_transdecoder_vs_braker_clean.slurm
  Runs DIAMOND BLASTP to align TransDecoder ORFs against BRAKER4 predicted proteins.

07_add_BRAKER_lengths_clean.py
  Parses the BRAKER4 GFF3 file to add contig ID, gene model coordinates, strand, gene length, CDS length, and protein length to each BRAKER4 isoform.

08_make_best_ORF_to_BRAKER_mapping_clean.py
  Parses the raw TransDecoder ORF versus BRAKER4 DIAMOND output, calculates coverage, and keeps one best BRAKER4 hit per TransDecoder ORF.

09_add_ONLY_Average_TPM_clean.py
  Adds `transdecoder_orf_id` and `Average_TPM` to the BRAKER4 annotation table using the best ORF-to-BRAKER mapping. It does not sum, average, or duplicate TPM values across all hits.

10_make_FINAL_clean_BRAKER_isoform_table.py
  Creates the final clean BRAKER4 isoform-level table. It keeps one row per BRAKER4 protein isoform, adds compartment from contig ID, adds only yes/no *Phaeodactylum tricornutum* status, and removes all PT detail columns.
```

Deprecated all-hit TPM scripts and GBK organelle-merging scripts were not used in the final clean table. The accepted final output remains BRAKER4 isoform-based and does not collapse gene IDs or append GenBank-derived organelle rows.

Each script can be run from the command line in the relevant working directory, as shown in the sections below.

## Analysis workflow
Click each section to expand the commands, notes, and outputs.

<details>
<summary><strong>1. Genome assembly</strong> - Flye metagenome assembly</summary>

Long-read assembly was performed with Nanopore reads basecalled with Guppy. Flye was run in metagenome mode because the sample represented a diatom-associated microbial consortium rather than an isolate genome.
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
The Flye assembly was used as the starting point for read mapping, polishing, binning, organelle screening, and genome annotation.

</details>

<details>
<summary><strong>2. Read mapping and assembly support</strong> - minimap2 and samtools</summary>

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

</details>

<details>
<summary><strong>3. Assembly polishing</strong> - Medaka, Polypolish, and Pypolca</summary>

Assembly polishing was performed using Medaka for long-read polishing, followed by Polypolish and Pypolca for short-read correction.
### 3.1 Long-read polishing with Medaka
```bash
medaka_consensus \
    -i pass_trim.fastq.gz \
    -d guppy_flye_assembly.fasta \
    -o medaka_euk_polished \
    -t 12
```
### 3.2 Short-read alignment for Polypolish
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
### 3.3 Polypolish filtering and polishing
```bash
polypolish filter \
    --in1 alignments_1.sam \
    --in2 alignments_2.sam \
    --out1 filtered_1.sam \
    --out2 filtered_2.sam

polypolish polish \
    /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/medaka_euk_polished/consensus.fasta \
    filtered_1.sam \
    filtered_2.sam \
    > sr_poly.fasta
```
### 3.4 Pypolca polishing
```bash
pypolca run \
    -a sr_poly.fasta \
    -1 /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz \
    -2 /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz \
    -t 12 \
    -o sr_pypolca_output \
    --careful
```
The final corrected assembly was used for downstream binning, organelle identification, and gene annotation.
### 3.5 BUSCO assessment of the polished assembly
```bash
busco \
    -i pypolca_corrected.fasta \
    -l busco_downloads/lineages/stramenopiles_odb10 \
    -o busco_report \
    -m genome
```

</details>

<details>
<summary><strong>4. Metagenomic binning</strong> - MetaBAT2</summary>

MetaBAT2 was used to recover genome bins from the polished assembly using Nanopore read coverage.
### 4.1 Install MetaBAT2
```bash
conda create -n metabat2_v2 -c conda-forge -c bioconda metabat2 libdeflate=1.10
conda activate metabat2_v2
```
### 4.2 Map Nanopore reads to the polished assembly
```bash
minimap2 -ax map-ont -t 16 \
    1_sr_pypolca_output/pypolca_corrected.fasta \
    pass_trim.fastq.gz \
    | samtools view -@ 16 -bS - \
    | samtools sort -@ 16 -m 10G -o aligned_reads.sorted.bam

samtools index -@ 16 aligned_reads.sorted.bam
```
### 4.3 Generate contig depth file
```bash
jgi_summarize_bam_contig_depths \
    --outputDepth depth.txt \
    --percentIdentity 85 \
    aligned_reads.sorted.bam
```
### 4.4 Run MetaBAT2
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
The resulting bins were used for quality assessment and taxonomic classification.

</details>

<details>
<summary><strong>5. Bin quality assessment</strong> - CheckM2</summary>

CheckM2 was used to estimate completeness and contamination of recovered genome bins.
```bash
checkm2 predict \
    --threads 16 \
    --input 2_metabat2_bins/ \
    --output_directory 3_checkm2_results
```

</details>

<details>
<summary><strong>6. Taxonomic classification</strong> - GTDB-Tk</summary>

GTDB-Tk was used to assign bacterial and archaeal taxonomy to recovered genome bins using the Genome Taxonomy Database.
```bash
gtdbtk classify_wf \
    --genome_dir 2_metabat2_bins/ \
    --out_dir 4_gtdbtk_output \
    --cpus 16 \
    -x fa
```
GTDB-Tk was used for bacterial and archaeal genome bins. Because the consortium also contained a dominant eukaryotic diatom, MetaEuk-based ORF taxonomy was used as an additional contig-level screen for eukaryotic, bacterial, ambiguous, and unclassified contig fractions.

</details>

<details>
<summary><strong>7. MetaEuk-based contig classification</strong> - MetaEuk and custom Python script</summary>

MetaEuk ORF-level taxonomic assignments were used to classify contigs across recovered bins. This step was added because the assembly originated from a diatom-associated microbial consortium, and bin-level bacterial taxonomy alone does not resolve eukaryotic contigs or mixed bins.
Because MetaEuk uses last common ancestor assignments, organelle-derived sequences can be assigned to bacterial lineages. To account for this, the contig classification used a priority hierarchy that grouped direct eukaryotic hits together with mitochondrial and chloroplast-derived signatures when calculating the eukaryotic score. Mitochondrial-like hits were identified using `o_Rickettsiales` or `o__Rickettsiales`, and chloroplast-like hits were identified using `p_Cyanobacteria`.
### 7.1 Input files
```text
metaeuk_output_polyp_taxonomy_tax_per_pred.tsv
contig_to_bin.txt
```
`metaeuk_output_polyp_taxonomy_tax_per_pred.tsv` contains ORF-level MetaEuk taxonomic assignments. `contig_to_bin.txt` links contig IDs to bin IDs.
The MetaEuk table was expected to contain the following columns:
```text
Contig_ID
Classification
```
### 7.2 ORF-level labels
Each predicted ORF was assigned to one of the following labels:
| Label | Rule |
|---|---|
| Eukaryota | `Classification` contains `d_Eukaryota` |
| Mitochondria-derived | `Classification` contains `o_Rickettsiales` or `o__Rickettsiales` |
| Chloroplast-derived | `Classification` contains `p_Cyanobacteria` |
| Bacteria | `Classification` contains `d_Bacteria`, excluding the organelle-derived categories above |
| Ambiguous (Cellular Org) | `Classification` is exactly `_cellular organisms` |
| Other | all other biological hits, including Archaea or viruses |
| Unclassified | no MetaEuk classification available |
### 7.3 Final contig-level priority hierarchy
For each contig, ORF-level labels were counted and assigned using a priority hierarchy.
First, the 30% eukaryotic/organelle rule was applied. A contig was classified as `Eukaryota` if more than 30% of its biological ORF assignments were eukaryotic, mitochondrial, or chloroplast-derived:
```text
(Eukaryota + Mitochondria-derived + Chloroplast-derived) / Total biological ORFs > 0.30
```
where:
```text
Total biological ORFs = Eukaryota + Mitochondria-derived + Chloroplast-derived + Bacteria + Ambiguous + Other
```
This rule has the highest priority. If the >30% threshold is met, the contig is labeled `Eukaryota` immediately, even if bacterial ORFs have a higher raw count. For example, a contig with 35% eukaryotic/organelle hits and 40% bacterial hits is still classified as `Eukaryota`.
If the contig does not meet the >30% eukaryotic/organelle threshold, the remaining rules are applied:
```text
Bacteria  = bacterial ORFs are strictly greater than the eukaryotic/organelle group and Other
Other     = Other ORFs are strictly greater than the eukaryotic/organelle group and Bacteria
Ambiguous = no group has a strict majority, including tied or unresolved cases
Unclassified = no biological MetaEuk hits are detected
```
Thus, contigs that fail the 30% rule and do not have a single strictly dominant biological category are labeled `Ambiguous (Cellular Org)`.
### 7.4 Python script
The full Python script is saved in:

```text
scripts/classify_metaeuk_contigs.py
```

Run the script from the directory containing `metaeuk_output_polyp_taxonomy_tax_per_pred.tsv` and `contig_to_bin.txt`:

```bash
python scripts/classify_metaeuk_contigs.py
```
### 7.5 Outputs
```text
contig_classification_final_priority.csv
normalized_histogram_final.png
```
The CSV file contains the bin name, contig ID, ORF-level category counts, and final contig classification. The stacked bar plot summarizes the normalized proportion of bacterial, eukaryotic, ambiguous, unclassified, and other contigs per bin.

</details>

<details>
<summary><strong>8. Genome coverage and relative abundance</strong> - CoverM</summary>

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

</details>

<details>
<summary><strong>9. 18S rRNA phylogenetic analysis</strong> - Clustal Omega, TrimAl, and IQ-TREE 2</summary>

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

</details>

<details>
<summary><strong>10. Transcriptome analysis</strong> - nf-core/metatdenovo</summary>

Transcriptome assembly and annotation were performed using the nf-core/metatdenovo workflow.
### 10.1 Java setup
```bash
module purge
module load java/openjdk-23.0.1

export JAVA_HOME=$(dirname $(dirname $(readlink -f $(which java))))
export PATH=$JAVA_HOME/bin:$PATH
```
### 10.2 nf-core/metatdenovo execution
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

</details>

<details>
<summary><strong>11. rRNA gene identification from transcriptome assemblies</strong> - Barrnap</summary>

Barrnap was used to identify eukaryotic, bacterial, and mitochondrial rRNA genes from the assembled transcriptome.
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
The resulting FASTA and GFF files were used to identify eukaryotic, bacterial, and mitochondrial rRNA transcripts.

</details>

<details>
<summary><strong>12. Organelle genome identification</strong> - MetaQUAST</summary>

Organelle contigs were identified by comparing the polished assembly against reference mitochondrial and chloroplast genomes.
```text
Mitogenome reference:   MT742552
Chloroplast reference: MT742551
```
MetaQUAST was used to compare assembly contigs against the organelle references.
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

</details>

<details>
<summary><strong>13. Diatom genome annotation with BRAKER4 ET mode</strong> - BRAKER4 ET mode</summary>

Gene models were generated with BRAKER4 using a diatom genome assembly and RNA-seq evidence. The genome was not soft-masked before annotation; therefore, repeat masking was performed internally within the BRAKER4 workflow using RepeatModeler, RepeatMasker, and TRF. The final accepted run used ET mode, meaning that gene prediction was based on RNA-seq evidence only.
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
Final gene models, proteins, CDS, and BUSCO assessment
```
### 13.1 BRAKER4 setup
```bash
git clone https://github.com/Gaius-Augustus/BRAKER4.git
cd BRAKER4

singularity pull braker3.sif docker://teambraker/braker3:latest
```
Snakemake was installed in the working environment.
```bash
conda install -c conda-forge -c bioconda snakemake
```
GeneMark requires a license key. The key was stored at:
```bash
/home/ruchita.solanki/.gm_key
```
### 13.2 STAR RNA-seq alignment evidence
RNA-seq reads were aligned to the genome using STAR. BRAKER4 was supplied with a precomputed coordinate-sorted BAM file rather than raw RNA-seq FASTQ files.
```bash
STAR \
    --runThreadN 24 \
    --runMode genomeGenerate \
    --genomeDir /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index \
    --genomeFastaFiles /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta \
    --genomeSAindexNbases 10
```
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
The BAM file was indexed before BRAKER4 execution.
```bash
samtools index /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam
```
### 13.3 Genome input and repeat masking
The genome assembly used for STAR indexing and BRAKER4 was:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
```
The assembly contained 3,010 contigs and had a total length of approximately 82.17 Mbp.
```bash
seqkit stats /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
```
The genome was checked for soft masking.
```bash
grep -v "^>" /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta \
    | grep -q '[a-z]' && echo "soft-masked" || echo "not soft-masked"
```
Output:
```text
not soft-masked
```
Because the genome was not pre-masked, the `genome_masked` column in `samples.csv` was left empty and internal repeat masking was enabled in `config.ini`.
### 13.4 BRAKER4 sample file
```bash
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4
nano samples.csv
```
```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage
DL_diatom,/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta,,,/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam,,,,,,,,stramenopiles_odb12
```
The `protein_fasta` column was left empty to force ET mode and avoid GeneMark-ETP.
```bash
awk -F',' '{print NR, NF}' samples.csv
```
Expected output:
```text
1 13
2 13
```
### 13.5 BRAKER4 configuration
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
### 13.6 Snakemake dry run
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
The dry run included the expected ET-mode rules:
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
The presence of `run_genemark_et` and absence of `run_genemark_etp` confirmed that BRAKER4 was configured in ET mode.
### 13.7 BRAKER4 execution
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
The ET-mode run produced 9,000 StringTie transcripts and 19,114 intron hints from the RNA-seq BAM file. These hints were used during GeneMark-ET and AUGUSTUS prediction.
### 13.8 Final BRAKER4 outputs
Final outputs were collected into:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET
```
Main output files:
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
Decompress final files:
```bash
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET

gunzip -c DL_diatom.braker4.ET.gff3.gz > DL_diatom.braker4.ET.gff3
gunzip -c DL_diatom.braker4.ET.gtf.gz > DL_diatom.braker4.ET.gtf
gunzip -c DL_diatom.braker4.ET.proteins.faa.gz > DL_diatom.braker4.ET.proteins.faa
gunzip -c DL_diatom.braker4.ET.cds.fna.gz > DL_diatom.braker4.ET.cds.fna
```
No additional TSEBRA run was required because TSEBRA refinement was included within BRAKER4.
### 13.9 Annotation statistics
```bash
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET

grep -c $'\tgene\t' DL_diatom.braker4.ET.gff3
grep -c $'\ttranscript\t' DL_diatom.braker4.ET.gff3
grep -c $'\tCDS\t' DL_diatom.braker4.ET.gff3

grep -c "^>" DL_diatom.braker4.ET.proteins.faa
grep -c "^>" DL_diatom.braker4.ET.cds.fna

seqkit stats DL_diatom.braker4.ET.proteins.faa DL_diatom.braker4.ET.cds.fna
```
Final ET-mode annotation statistics:
```text
Genes:        15,102
Transcripts:  16,947
Proteins:     16,947
CDS FASTA:    16,947
CDS features: 31,713
Exons:        31,713
Introns:      14,952
```
Predicted proteins were checked for internal stop codons.
```bash
grep -n "\*" DL_diatom.braker4.ET.proteins.faa | head
```
No internal stop codons were detected.
### 13.10 BUSCO assessment of the final protein set
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
The final predicted protein set produced the following BUSCO result using `stramenopiles_odb12`:
```text
C:84.8%[S:81.1%,D:3.7%],F:1.9%,M:13.3%,n=697
```
### 13.11 Annotation acceptance
The final BRAKER4 ET annotation was accepted for downstream analysis because it produced a plausible gene set for the diatom genome assembly, showed no internal stop codon issues in the predicted protein FASTA, and recovered 84.8% of the `stramenopiles_odb12` BUSCO protein set with low duplication.
Final accepted annotation files:
```text
DL_diatom.braker4.ET.gff3
DL_diatom.braker4.ET.gtf
DL_diatom.braker4.ET.proteins.faa
DL_diatom.braker4.ET.cds.fna
```
### 13.12 Rationale for ET mode instead of ETP
BRAKER4 was initially tested in ETP mode, which combines RNA-seq evidence with protein evidence. However, GeneMark-ETP failed during model training. Although protein-supported alignments were generated, the GeneMark-ETP training set did not produce valid gene and transcript models. The failed run reported zero parsed genes and transcripts, CDS-only training entries, phase-distribution errors, and division-by-zero errors in the GeneMark-ETP scripts.
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
Because GeneMark-ETP did not complete successfully, the annotation was rerun in ET mode using RNA-seq evidence only. This avoided the failed protein-dependent GeneMark-ETP training step while retaining transcript evidence from the coordinate-sorted STAR BAM file. The final successful workflow used GeneMark-ET, AUGUSTUS, and TSEBRA, with `protein_fasta` left empty in `samples.csv` and `mode = et` specified in `config.ini`.

</details>

<details>
<summary><strong>14. Functional annotation of BRAKER4-predicted proteins</strong> - DIAMOND, Swiss-Prot, Bacillariophyta, and InterProScan</summary>
   
After accepting the BRAKER4 ET annotation, the predicted protein set was used for downstream functional annotation. The final annotation strategy was designed specifically for a eukaryotic diatom genome rather than a prokaryote-centered metagenomic annotation workflow. For this reason, the workflow used three complementary annotation layers: curated Swiss-Prot homology, diatom-focused UniProtKB Bacillariophyta homology, and InterProScan domain/family annotation.
The final BRAKER4 ET protein file was used as the main input:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET/DL_diatom.braker4.ET.proteins.faa
```
The protein set contained:
```text
16,947 predicted proteins
```
A new functional annotation working directory was created:
```bash
mkdir -p /work/ebg_lab/eb/diatom_consortia/functional_annotation_swissprot
cd /work/ebg_lab/eb/diatom_consortia/functional_annotation_swissprot

mkdir -p 00_databases 01_input 02_diamond 03_best_hits 05_interproscan 06_combined_annotation logs scripts slurm
```
The BRAKER4 protein file was linked into the working directory:
```bash
ln -sfn /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET/DL_diatom.braker4.ET.proteins.faa \
    01_input/diatom_predicted_proteins.fa
```
### 14.1 Logic of the annotation strategy
The annotation workflow used multiple evidence layers because no single database provides complete and fully reliable functional annotation for a non-model diatom genome. Swiss-Prot was used as the conservative curated layer because its entries are manually reviewed, but it is expected to annotate fewer proteins because it is smaller and not diatom-rich. A UniProtKB Bacillariophyta database was added to improve diatom-specific homolog detection, while InterProScan was used to identify conserved domains, protein families, GO terms, and pathway signatures.
eggNOG, KEGG, COG, and dbCAN were not used in the final workflow. eggNOG, KEGG, and COG were avoided because the objective was a eukaryote- and diatom-focused annotation rather than broad prokaryotic orthology assignment. dbCAN was not included because specialized carbohydrate-active enzyme classification was not the central objective of this analysis.
### 14.2 Software environment for DIAMOND annotation
A conda environment was used for DIAMOND searches and parsing:
```bash
conda create -n swissprot_annot -c conda-forge -c bioconda diamond pandas seqkit wget pigz -y
conda activate swissprot_annot
```

### 14.3 Swiss-Prot database setup
Swiss-Prot was downloaded and stored in the home directory to avoid filling the project working directory:
```bash
mkdir -p $HOME/databases/swissprot/raw
mkdir -p $HOME/databases/swissprot/diamond
```
The downloaded Swiss-Prot release was:
```text
UniProtKB/Swiss-Prot Release 2026_02 of 10-Jun-2026
```

The database contained:
```text
575,503 reviewed protein sequences
```
A DIAMOND database was built:
```bash
diamond makedb \
    --in $HOME/databases/swissprot/raw/uniprot_sprot.fasta.gz \
    --db $HOME/databases/swissprot/diamond/uniprot_sprot.dmnd
```
The database directory was linked into the project:
```bash
ln -sfn $HOME/databases/swissprot 00_databases/swissprot_home
```
### 14.4 DIAMOND search against Swiss-Prot
Predicted proteins were searched against Swiss-Prot using DIAMOND BLASTP:
```bash
diamond blastp \
    --query 01_input/diatom_predicted_proteins.fa \
    --db 00_databases/swissprot_home/diamond/uniprot_sprot.dmnd \
    --out 02_diamond/DL_diatom_braker4_ET_vs_swissprot.tsv \
    --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle \
    --evalue 1e-5 \
    --max-target-seqs 5 \
    --sensitive \
    --threads 32
```
The Swiss-Prot search was used to assign conservative protein names where strong curated homologs were available. The `--sensitive` setting was used to improve detection of more distant homologs, while `--evalue 1e-5` retained candidate hits for downstream filtering. Up to five hits per query were retained so that the best hit could later be selected after calculating coverage and confidence metrics.
#### Swiss-Prot result summary
```text
Total predicted proteins: 16,947
Total DIAMOND hit lines: 36,669
Proteins with at least one Swiss-Prot hit: 8,078
Proteins with strict Swiss-Prot hit: 4,284

Percent with at least one Swiss-Prot hit: 47.67%
Percent with strict Swiss-Prot hit: 25.28%
```
Swiss-Prot confidence counts:
```text
High:                    2,037
Medium:                  2,247
Low:                     2,826
Weak domain or fragment:   968
```
### 14.5 UniProtKB Bacillariophyta database setup
A diatom-focused UniProtKB database was created using the Bacillariophyta taxonomic group. The database was stored in the home directory:
```bash
mkdir -p $HOME/databases/uniprot_bacillariophyta/raw
mkdir -p $HOME/databases/uniprot_bacillariophyta/diamond
```
The compressed FASTA file was downloaded from UniProtKB:
```bash
curl -L --retry 5 --retry-delay 10 \
    -o $HOME/databases/uniprot_bacillariophyta/raw/uniprotkb_bacillariophyta_taxid2836.fasta.gz \
    "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28taxonomy_id%3A2836%29"
```

A DIAMOND database was built:
```bash
diamond makedb \
    --in $HOME/databases/uniprot_bacillariophyta/raw/uniprotkb_bacillariophyta_taxid2836.fasta.gz \
    --db $HOME/databases/uniprot_bacillariophyta/diamond/uniprotkb_bacillariophyta_taxid2836.dmnd
```
The database was linked into the project:
```bash
ln -sfn $HOME/databases/uniprot_bacillariophyta 00_databases/uniprot_bacillariophyta_home
```
### 14.6 DIAMOND search against UniProtKB Bacillariophyta
Predicted proteins were searched against the Bacillariophyta database using DIAMOND BLASTP:
```bash
diamond blastp \
    --query 01_input/diatom_predicted_proteins.fa \
    --db 00_databases/uniprot_bacillariophyta_home/diamond/uniprotkb_bacillariophyta_taxid2836.dmnd \
    --out 02_diamond/DL_diatom_braker4_ET_vs_uniprot_bacillariophyta.tsv \
    --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle \
    --evalue 1e-5 \
    --max-target-seqs 5 \
    --sensitive \
    --threads 32
```
#### Bacillariophyta result summary
```text
Total predicted proteins: 16,947
Total DIAMOND hit lines: 65,416
Proteins with at least one Bacillariophyta hit: 13,574
Proteins with strict Bacillariophyta hit: 11,294

Percent with at least one Bacillariophyta hit: 80.10%
Percent with strict Bacillariophyta hit: 66.64%
```
Bacillariophyta confidence counts:
```text
High:                     8,557
Medium:                   2,737
Low:                      1,668
Weak domain or fragment:    612
```
### 14.7 Best-hit parsing and confidence filtering
Raw DIAMOND outputs were parsed into best-hit tables. For each predicted protein, query coverage and subject coverage were calculated, UniProt identifiers were parsed, and protein name, organism, gene name, taxon ID, and protein evidence fields were extracted from the hit title.
The best hit per query was selected using:
```text
lowest e-value
highest bitscore
highest query coverage
highest percent identity
```
Confidence classes were assigned as:
```text
High:
e-value <= 1e-20
query coverage >= 70%
percent identity >= 40%

Medium:
e-value <= 1e-10
query coverage >= 50%
percent identity >= 30%

Low:
e-value <= 1e-5
query coverage >= 30%

Weak domain or fragment:
all remaining reported hits
```

The strict filtered set was defined as:
```text
e-value <= 1e-10
query coverage >= 50%
percent identity >= 30%
```
The main parsed Swiss-Prot outputs were:
```text
03_best_hits/DL_diatom_swissprot_all_hits_with_coverage.tsv
03_best_hits/DL_diatom_swissprot_best_hits.tsv
03_best_hits/DL_diatom_swissprot_best_hits_strict.tsv
03_best_hits/DL_diatom_all_proteins_with_swissprot_annotation.tsv
03_best_hits/DL_diatom_swissprot_annotation_summary.txt
```
The main parsed Bacillariophyta outputs were:
```text
03_best_hits/DL_diatom_bacillariophyta_all_hits_with_coverage.tsv
03_best_hits/DL_diatom_bacillariophyta_best_hits.tsv
03_best_hits/DL_diatom_bacillariophyta_best_hits_strict.tsv
03_best_hits/DL_diatom_all_proteins_with_bacillariophyta_annotation.tsv
03_best_hits/DL_diatom_bacillariophyta_annotation_summary.txt
```
### 14.7.1 Swiss-Prot best-hit parsing script
The raw Swiss-Prot DIAMOND output was parsed with a custom Python script. The script calculated query and subject coverage, extracted UniProt accession identifiers, protein names, organism names, gene names, taxon IDs, and protein evidence fields, and selected the best hit per predicted protein.
The full Python script is saved in:
```text
scripts/make_swissprot_best_hits.py
```
Run the script with:
```bash
conda activate swissprot_annot
python scripts/make_swissprot_best_hits.py
```
#### Logic
The raw DIAMOND output was not used directly because some statistically significant hits only aligned to short domains or fragments. This script calculated query coverage and assigned confidence classes so that full-length or near-full-length homologs could be separated from weaker domain-only matches. The all-protein output was retained so that proteins without Swiss-Prot hits could still be included in downstream InterProScan and expression merges.

### 14.7.2 Bacillariophyta best-hit parsing script
The UniProtKB Bacillariophyta DIAMOND output was parsed using the same logic as the Swiss-Prot output. This kept the confidence framework consistent between the curated Swiss-Prot layer and the diatom-specific homolog layer.
The full Python script is saved in:
```text
scripts/make_bacillariophyta_best_hits.py
```
Run the script with:
```bash
conda activate swissprot_annot
python scripts/make_bacillariophyta_best_hits.py
```
#### Logic
The Bacillariophyta database increased diatom-specific homolog detection, but many entries are unreviewed. Therefore, hits were filtered and ranked using the same e-value, coverage, and identity criteria used for Swiss-Prot. This allowed the Bacillariophyta output to be used as a homolog-support layer rather than as an unfiltered functional naming source.
### 14.8 InterProScan setup and annotation
InterProScan was installed in the home tools directory and linked into the project:
```bash
ln -sfn $HOME/tools/interproscan/current 00_databases/interproscan_home
```

InterProScan was run in a Java 11 conda environment:
```bash
conda create -n interproscan_env -c conda-forge openjdk=11 -y
conda activate interproscan_env
```

The full BRAKER4 ET protein set was submitted to InterProScan using SLURM:
```bash
sbatch slurm/interproscan_diatom_full.sh
```
The SLURM script used was:
```bash
#!/bin/bash
set -euo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate interproscan_env

cd /work/ebg_lab/eb/diatom_consortia/functional_annotation_swissprot

mkdir -p 05_interproscan
mkdir -p logs
mkdir -p temp

00_databases/interproscan_home/interproscan.sh \
    -i 01_input/diatom_predicted_proteins.fa \
    -f TSV \
    -dp \
    -goterms \
    -pa \
    -exclappl MobiDBLite \
    -cpu 32 \
    -o 05_interproscan/DL_diatom_braker4_ET_interproscan.tsv
```
MobiDBLite was excluded after the initial run failed because of a Python compatibility error in the bundled MobiDBLite script. Because MobiDBLite predicts intrinsically disordered regions and was not central to the functional annotation goals, it was excluded while retaining the main protein family, domain, and GO annotation resources.
The completed InterProScan run included:
```text
AntiFam
CDD
Coils
FunFam
Gene3D
Hamap
NCBIfam
PANTHER
Pfam
PIRSF
PIRSR
PRINTS
ProSitePatterns
ProSiteProfiles
SFLD
SMART
SUPERFAMILY
```
The `-goterms` option was used to report GO terms, and `-pa` was used to include pathway annotations where available. The `-dp` option disabled the pre-calculated match lookup service, making the job independent of internet access from the compute node.
#### InterProScan result summary
InterProScan completed successfully and produced:
```text
Raw InterProScan rows:              102,153
Proteins with InterProScan hits:     13,106 / 16,947
Percent with InterProScan hits:      77.34%
```

The raw InterProScan file was:
```text
05_interproscan/DL_diatom_braker4_ET_interproscan.tsv
```
The raw file contains multiple rows per protein because a single protein can contain multiple domains, repeats, sites, signatures, GO terms, or pathway annotations. Therefore, the raw output was summarized to one row per predicted protein before integration with Swiss-Prot and Bacillariophyta annotations.
### 14.9 InterProScan summary by protein
The raw InterProScan TSV was collapsed into one row per protein using a custom Python script:
```text
scripts/summarize_interproscan.py
```
Run the script with:
```bash
conda activate swissprot_annot
python scripts/summarize_interproscan.py
```
The output files were:
```text
06_combined_annotation/DL_diatom_interproscan_summary_by_protein.tsv
06_combined_annotation/DL_diatom_all_proteins_with_interproscan_summary.tsv
06_combined_annotation/DL_diatom_interproscan_analysis_counts.tsv
06_combined_annotation/DL_diatom_interproscan_summary_stats.txt
```
The summary table includes:
```text
protein_id
protein_length
number of InterProScan rows
analyses
signature accessions
signature descriptions
InterPro accessions
InterPro descriptions
GO terms
pathway annotations
```
The InterProScan result was also summarized by analysis source:
```text
analysis         raw_rows
Pfam             17,891
Gene3D           17,698
SUPERFAMILY      13,811
PANTHER           9,840
SMART             7,367
PRINTS            7,354
ProSiteProfiles   7,018
Coils             5,536
CDD               5,075
FunFam            3,292
NCBIfam           2,853
ProSitePatterns   2,748
Hamap               901
PIRSF               537
SFLD                230
AntiFam               2
```
Pfam, Gene3D, SUPERFAMILY, PANTHER, SMART, PRINTS, and ProSiteProfiles contributed the largest number of raw annotation rows. These counts represent raw annotation rows rather than unique proteins, because individual proteins can contain multiple domains or match multiple databases.
#### AntiFam screening
InterProScan reported two AntiFam matches:
```text
protein_id   AntiFam accession   AntiFam description
g10893.t1    ANF00012            tRNA
g11404.t1    ANF00005            Antisense to 23S rRNA
```

AntiFam matches were treated as warning flags rather than functional annotations. These two proteins were marked as potentially spurious predictions or RNA-associated ORFs and should not be used for biological interpretation of protein function.

The AntiFam hits were inspected with:
```bash
awk -F '\t' '$4=="AntiFam"' 05_interproscan/DL_diatom_braker4_ET_interproscan.tsv | column -t
```
### 14.10 Master functional annotation table
After Swiss-Prot, Bacillariophyta, and InterProScan outputs were complete, all evidence layers were merged into a single master annotation table. The merge also incorporated the two AntiFam warning flags.
The AntiFam flag file was created with explicit tab separators:
```bash
printf "protein_id\tantifam_accession\tantifam_description\tinterpretation\n" \
> 06_combined_annotation/DL_diatom_antifam_flagged_proteins.tsv

printf "g10893.t1\tANF00012\ttRNA\tpotential_spurious_or_RNA_associated_ORF\n" \
>> 06_combined_annotation/DL_diatom_antifam_flagged_proteins.tsv

printf "g11404.t1\tANF00005\tAntisense to 23S rRNA\tpotential_spurious_or_RNA_associated_ORF\n" \
>> 06_combined_annotation/DL_diatom_antifam_flagged_proteins.tsv
```
The final annotation layers were merged using:
```text
scripts/merge_functional_annotation_layers.py
```
Run the script with:
```bash
conda activate swissprot_annot
python scripts/merge_functional_annotation_layers.py
```
The script merged:
```text
03_best_hits/DL_diatom_all_proteins_with_swissprot_annotation.tsv
03_best_hits/DL_diatom_all_proteins_with_bacillariophyta_annotation.tsv
06_combined_annotation/DL_diatom_all_proteins_with_interproscan_summary.tsv
06_combined_annotation/DL_diatom_antifam_flagged_proteins.tsv
```
The merged master table was checked with:
```bash
wc -l 06_combined_annotation/DL_diatom_master_functional_annotation.tsv
```
The final master table contained:
```text
16,948 lines = 16,947 predicted proteins + 1 header line
```
AntiFam-flagged proteins were checked in the merged table:
```bash
grep -E 'g10893.t1|g11404.t1' \
  06_combined_annotation/DL_diatom_master_functional_annotation.tsv \
  | column -t -s $'\t'
```
Both AntiFam proteins were retained in the master table but flagged with `has_antifam_flag = yes`. This preserves the annotation evidence while allowing these entries to be excluded from biological interpretation.
### 14.11 AntiFam-filtered interpretation tables
Two downstream interpretation files were generated from the master annotation table. One file contains only AntiFam-flagged proteins, and the other removes AntiFam-flagged proteins for cleaner downstream interpretation.
```bash
awk -F'\t' 'NR==1 || $9=="yes"' \
  06_combined_annotation/DL_diatom_master_functional_annotation.tsv \
  > 06_combined_annotation/DL_diatom_master_functional_annotation_AntiFam_flagged_only.tsv
```
This file contains the two proteins flagged by AntiFam plus the table header.
```bash
awk -F'\t' 'NR==1 || $9!="yes"' \
  06_combined_annotation/DL_diatom_master_functional_annotation.tsv \
  > 06_combined_annotation/DL_diatom_master_functional_annotation_no_AntiFam.tsv
```
This file excludes AntiFam-flagged proteins and can be used for downstream functional summaries where potential RNA-associated ORFs should not be interpreted as protein-coding genes.
### 14.12 Final functional annotation outputs
The final functional annotation outputs were:
```text
06_combined_annotation/DL_diatom_interproscan_summary_by_protein.tsv
06_combined_annotation/DL_diatom_all_proteins_with_interproscan_summary.tsv
06_combined_annotation/DL_diatom_interproscan_analysis_counts.tsv
06_combined_annotation/DL_diatom_interproscan_summary_stats.txt
06_combined_annotation/DL_diatom_antifam_flagged_proteins.tsv
06_combined_annotation/DL_diatom_master_functional_annotation.tsv
06_combined_annotation/DL_diatom_master_functional_annotation_for_manual_categories.tsv
06_combined_annotation/DL_diatom_master_functional_annotation_AntiFam_flagged_only.tsv
06_combined_annotation/DL_diatom_master_functional_annotation_no_AntiFam.tsv
06_combined_annotation/DL_diatom_master_functional_annotation_summary.txt
```

The master annotation table contains one row per predicted protein and includes:
```text
protein_id
recommended annotation
recommended annotation source
recommended annotation confidence
Swiss-Prot annotation fields
Bacillariophyta annotation fields
InterProScan domains and signatures
GO terms
pathway annotations
AntiFam flag
```
#### Logic
The merged table combines conservative curated protein names, diatom-specific homolog support, domain-based functional evidence, GO terms, pathway annotations, and AntiFam warning flags. This makes the annotation more interpretable than any single database alone while retaining proteins without homology hits for later expression integration and manual category assignment.
### 14.13 Expression integration status

Expression integration was completed after the master functional annotation table was generated. Transcriptome ORFs predicted by TransDecoder were aligned to the BRAKER4 ET protein set using DIAMOND BLASTP. One best BRAKER4 hit was retained per TransDecoder ORF to prevent the same transcript-level `Average_TPM` value from being duplicated across multiple BRAKER4 proteins.

The final expression merge used only the `Average_TPM` column from the trusted transcriptome table. No all-hit TPM sums, means, hit counts, or mapped-ORF lists were carried into the final table. If multiple TransDecoder ORFs mapped best to the same BRAKER4 protein, the ORF with the highest valid `Average_TPM` was used as the representative expression-supported ORF. The corresponding `transdecoder_orf_id` was retained next to the BRAKER4 `gene_id` for traceability.

The detailed expression integration and final gene-table construction workflow is documented in Section 17.

### 14.14 Current status of functional annotation
Completed:
```text
Swiss-Prot database download and DIAMOND database construction
Swiss-Prot DIAMOND search
Swiss-Prot best-hit parsing
UniProtKB Bacillariophyta database download and DIAMOND database construction
Bacillariophyta DIAMOND search
Bacillariophyta best-hit parsing
InterProScan installation and Java 11 setup
InterProScan full run completed
InterProScan summary by protein completed
InterProScan database contribution summary completed
AntiFam screen completed
Master functional annotation table generated
AntiFam-flagged interpretation table generated
AntiFam-filtered interpretation table generated
Manual functional category scaffold generated
Expression integration with TransDecoder ORFs and Average_TPM completed
Final clean BRAKER4 isoform-level gene table generated
```

</details>

<details>
<summary><strong>15. Nuclear-enriched genome generation</strong> - organelle contig removal and BRAKER4 annotation filtering</summary>

### 15.1 Input files
```bash
mkdir -p /work/ebg_lab/eb/diatom_consortia/nuclear_genome_filtering_18_diatom
cd /work/ebg_lab/eb/diatom_consortia/nuclear_genome_filtering_18_diatom

WHOLE=/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
CHLORO=/work/ebg_lab/eb/diatom_consortia/organelle/2_chloro/chloroplast_contig_1443_trimmed.fasta
MITO=/work/ebg_lab/eb/diatom_consortia/organelle/mito/diatom_candidate_mitochondrion_2contigs.fasta
```
This command creates a dedicated working directory for nuclear genome filtering and moves into it. The variables define the whole diatom assembly, the chloroplast reference, and the mitochondrial reference used for organelle contig detection.
Input assembly statistics:
```bash
seqkit stats $WHOLE $CHLORO $MITO
```
This command summarizes the number of sequences and total length of the input assembly and organelle references. It provides a starting point for checking how much sequence is removed during nuclear-enrichment filtering.
```text
18_diatom.fasta                              3,010 contigs   82,172,226 bp
chloroplast_contig_1443_trimmed.fasta           1 contig       120,429 bp
diatom_candidate_mitochondrion_2contigs.fasta    2 contigs      104,526 bp
```
### 15.2 Combine organelle references
```bash
cat $CHLORO $MITO > organelles_chloro_mito.fasta
seqkit stats organelles_chloro_mito.fasta
```
This command combines the chloroplast and mitochondrial FASTA files into a single organelle reference file. The `seqkit stats` command confirms the number of organelle sequences and total reference length.
The combined organelle reference contained three sequences with a total length of 224,955 bp.
### 15.3 Align the genome against organelle references
```bash
minimap2 -x asm5 -c \
    organelles_chloro_mito.fasta \
    $WHOLE \
    > whole_vs_organelles.paf
```
This command aligns the whole diatom assembly against the combined chloroplast and mitochondrial references. The `asm5` preset is appropriate for closely related assembled sequences, and the PAF output records which assembly contigs align to organelle references.
The PAF output was used to identify regions of `18_diatom.fasta` with similarity to the chloroplast or mitochondrial genomes.
### 15.4 Calculate organelle-aligned coverage per contig
```bash
awk 'BEGIN{OFS="\t"} {print $1, $3, $4}' whole_vs_organelles.paf \
    > whole_vs_organelles.query_intervals.bed

sort -k1,1 -k2,2n whole_vs_organelles.query_intervals.bed \
    > whole_vs_organelles.query_intervals.sorted.bed

bedtools merge \
    -i whole_vs_organelles.query_intervals.sorted.bed \
    > whole_vs_organelles.query_intervals.merged.bed

seqkit fx2tab -n -l $WHOLE > whole_contig_lengths.tsv

awk 'BEGIN{OFS="\t"} {aligned[$1] += ($3 - $2)} END {for (c in aligned) print c, aligned[c]}' \
    whole_vs_organelles.query_intervals.merged.bed \
    > organelle_aligned_length_per_contig.tsv

awk 'BEGIN{OFS="\t"}
FNR==NR {len[$1]=$2; next}
{
    contig=$1;
    aligned=$2;
    pct=(aligned/len[contig])*100;
    print contig, len[contig], aligned, pct
}' whole_contig_lengths.tsv organelle_aligned_length_per_contig.tsv \
    > organelle_coverage_per_contig.tsv

sort -k4,4nr organelle_coverage_per_contig.tsv | head -n 50
```
The first `awk` command extracts query contig coordinates from the PAF file and converts them into BED-like intervals. The `sort` and `bedtools merge` commands merge overlapping aligned intervals so that aligned regions are not counted multiple times.
`seqkit fx2tab` records the total length of each contig in the whole assembly. The following `awk` commands calculate the total organelle-aligned length per contig and convert this to percent coverage, while the final `sort` command displays the strongest organelle-like candidates.
### 15.5 Remove organelle-like contigs
Contigs were classified as organelle-like if at least 70% of the contig length aligned to the chloroplast or mitochondrial reference.
```bash
awk '$4 >= 70 {print $1}' organelle_coverage_per_contig.tsv \
    > organelle_like_contigs.70pct.txt
```
This command selects contigs for removal if at least 70% of their length aligns to the organelle reference. The resulting text file contains only the names of contigs classified as chloroplast-like or mitochondrial-like.
This identified three organelle-like contigs:
```text
contig_1443
contig_5628
contig_1647
```
These corresponded to the chloroplast-like contig and two mitochondrial-like contigs.
```bash
seqkit grep \
    -v \
    -f organelle_like_contigs.70pct.txt \
    $WHOLE \
    > 18_diatom_nuclear_enriched.v1.fasta
```
This command removes the organelle-like contigs from the whole diatom assembly. The `-v` option keeps all sequences except those listed in `organelle_like_contigs.70pct.txt`, producing the nuclear-enriched genome FASTA.
Assembly statistics were calculated before and after organelle removal.
```bash
seqkit stats \
    $WHOLE \
    organelles_chloro_mito.fasta \
    18_diatom_nuclear_enriched.v1.fasta
```
This command compares the original assembly, the combined organelle reference, and the final nuclear-enriched genome. It confirms the number of contigs and total assembly length after organelle-like contigs were removed.
Final nuclear-enriched genome statistics:
```text
18_diatom.fasta                         3,010 contigs   82,172,226 bp
organelles_chloro_mito.fasta                3 contigs      224,955 bp
18_diatom_nuclear_enriched.v1.fasta     3,007 contigs   81,911,772 bp
```
Three organelle-like contigs were removed, corresponding to 260,454 bp or 0.317% of the `18_diatom.fasta` assembly.
### 15.6 Filter BRAKER4 annotation to nuclear contigs
BRAKER4 was not rerun after organelle filtering because only three organelle-like contigs were removed from the BRAKER4 input assembly. Instead, the existing BRAKER4 annotation was filtered to retain only features located on contigs present in the nuclear-enriched genome.
```bash
seqkit seq -n 18_diatom_nuclear_enriched.v1.fasta > nuclear_contigs.v1.txt
MY_GFF=/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET/DL_diatom.braker4.ET.gff3
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
The `seqkit seq -n` command extracts the names of contigs retained in the nuclear-enriched genome, and `MY_GFF` defines the original BRAKER4 annotation file. The `awk` command then keeps only GFF3 features located on retained nuclear contigs while excluding the embedded FASTA section.
Count retained nuclear gene models:
```bash
grep -c $'\tgene\t' braker.18_diatom_nuclear_only.v1.gff3
```
This command counts the number of gene features retained in the nuclear-filtered BRAKER4 annotation. It verifies the final number of predicted gene models associated with the nuclear-enriched genome.
Final nuclear genome files:
```text
18_diatom_nuclear_enriched.v1.fasta
braker.18_diatom_nuclear_only.v1.gff3
```
</details>

<details>
<summary><strong>16. Pairwise genome comparison with <em>Phaeodactylum tricornutum</em></strong> - BLASTN and bedtools</summary>
   
A pairwise genome comparison was performed between the updated nuclear-enriched diatom genome and the reference genome of *Phaeodactylum tricornutum*. This analysis was used as a nucleotide-level similarity screen and was not intended to define complete gene orthology. Protein-level comparison using BLASTP, DIAMOND, or OrthoFinder is recommended for a more complete assessment of conserved gene content.

### 16.1 Input files and working directory
```bash
MY_GENOME=/work/ebg_lab/eb/diatom_consortia/nuclear_genome_filtering_18_diatom/18_diatom_nuclear_enriched.v1.fasta
MY_GFF=/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET/DL_diatom.braker4.ET.gff3

mkdir -p /work/ebg_lab/eb/diatom_consortia/phaeodactylum_blast_18_diatom_v1
cd /work/ebg_lab/eb/diatom_consortia/phaeodactylum_blast_18_diatom_v1
mkdir -p blast_out blast_db phaeodactylum
```
These variables define the nuclear-enriched diatom genome and the original BRAKER4 annotation file used in the comparison. The working directory and subfolders organize the BLAST outputs, BLAST database files, and downloaded *Phaeodactylum* reference files.
### 16.2 Filter BRAKER4 GFF3 to nuclear contigs
```bash
seqkit seq -n $MY_GENOME > nuclear_contigs.v1.txt

awk 'BEGIN{FS=OFS="\t"}
FNR==NR {keep[$1]=1; next}
$0 ~ /^#/ {
    if ($0 !~ /^##FASTA/) print;
    next
}
$1 in keep {print}' nuclear_contigs.v1.txt $MY_GFF \
    > braker.18_diatom_nuclear_only.v1.gff3
```
This step extracts the contig names retained in the nuclear-enriched genome and uses them to filter the BRAKER4 GFF3 file. Only gene models located on retained nuclear contigs are kept, while features on removed organelle-like contigs are excluded.
The number of genes before and after filtering was checked.
```bash
grep -c $'\tgene\t' $MY_GFF
grep -c $'\tgene\t' braker.18_diatom_nuclear_only.v1.gff3
```
These commands count gene features before and after nuclear-contig filtering. The comparison confirms how many BRAKER4 gene models were removed because they were located on organelle-like contigs.
```text
Full BRAKER4 genes:              15102
Nuclear-filtered BRAKER4 genes:  15048
```
Thus, 54 gene models located on removed organelle-like contigs were excluded from the nuclear gene set.
### 16.3 Convert nuclear gene models to BED
```bash
MY_NUCLEAR_GFF=braker.18_diatom_nuclear_only.v1.gff3
awk -F'\t' 'BEGIN{OFS="\t"}
!/^#/ && $3=="gene" {
  id=$9;
  sub(/.*ID=/,"",id);
  sub(/;.*/,"",id);
  print $1,$4-1,$5,id,$6,$7
}' $MY_NUCLEAR_GFF \
    > blast_out/18_diatom_nuclear_v1_genes.bed
```
This command converts nuclear-filtered BRAKER4 gene coordinates from GFF3 format to BED format. The BED file is required for `bedtools intersect`, which is used later to overlap predicted genes with BLASTN hit regions.
The resulting BED file contained 15,048 nuclear gene models.
```bash
wc -l blast_out/18_diatom_nuclear_v1_genes.bed
```
This command counts the number of nuclear gene models in the BED file. It verifies that the BED conversion retained the expected 15,048 nuclear-filtered genes.
### 16.4 Download the *Phaeodactylum tricornutum* reference genome
```bash
datasets download genome accession GCF_000150955.2 \
    --include genome,gff3 \
    --filename phaeodactylum/Phaeodactylum_tricornutum_GCF_000150955.2.zip

unzip phaeodactylum/Phaeodactylum_tricornutum_GCF_000150955.2.zip -d phaeodactylum

PT_GENOME=$(find phaeodactylum -name "*genomic.fna" | head -n 1)
PT_GFF=$(find phaeodactylum -name "*.gff" -o -name "*.gff3" | head -n 1)

echo $PT_GENOME
echo $PT_GFF
```
The NCBI Datasets command downloads the *Phaeodactylum tricornutum* reference genome and GFF3 annotation. The `find` commands automatically locate the downloaded genome FASTA and annotation files for use in BLAST and gene-overlap analyses.
### 16.5 Build the BLAST database
```bash
makeblastdb \
    -in $PT_GENOME \
    -dbtype nucl \
    -parse_seqids \
    -out blast_db/Phaeodactylum_tricornutum

PT_DB=blast_db/Phaeodactylum_tricornutum
```
This command builds a nucleotide BLAST database from the *P. tricornutum* reference genome. The database is then assigned to `PT_DB` so it can be used as the target in the genome-vs-genome BLASTN search.
### 16.6 Run genome-vs-genome BLASTN
```bash
blastn \
    -query $MY_GENOME \
    -db $PT_DB \
    -task blastn \
    -evalue 1e-10 \
    -num_threads 32 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
    -out blast_out/18_diatom_nuclear_v1_vs_phaeodactylum.blastn.tsv
```
This command compares the nuclear-enriched diatom genome against the *P. tricornutum* reference genome at the nucleotide level. The output reports alignment coordinates, percent identity, e-values, bitscores, sequence lengths, and query coverage for downstream filtering.
The raw BLASTN output contained 34,991 alignments.
```bash
wc -l blast_out/18_diatom_nuclear_v1_vs_phaeodactylum.blastn.tsv
```
This command counts the number of raw BLASTN alignments. It provides the starting number of genome-level nucleotide matches before applying quality filters.
### 16.7 Filter BLASTN alignments
Alignments were filtered using the following thresholds:
```text
Percent identity:  ≥70%
Alignment length:  ≥200 bp
E-value:           ≤1e-10
```
```bash
awk 'BEGIN{OFS="\t"} $3 >= 70 && $4 >= 200 && $11 <= 1e-10' \
    blast_out/18_diatom_nuclear_v1_vs_phaeodactylum.blastn.tsv \
    > blast_out/18_diatom_nuclear_v1_vs_phaeodactylum.filtered.tsv
```
This command keeps only stronger BLASTN alignments with at least 70% identity, at least 200 bp alignment length, and e-value ≤ 1e-10. The filtered file reduces weak or very short matches before gene-overlap analysis.
This produced 4,895 filtered BLASTN alignments.
```bash
wc -l blast_out/18_diatom_nuclear_v1_vs_phaeodactylum.filtered.tsv
```
This command counts the number of BLASTN alignments that passed the filtering thresholds. It confirms that 4,895 nucleotide alignments were retained for downstream interpretation.
### 16.8 Identify diatom genes overlapping *Phaeodactylum*-like regions
Filtered BLASTN hits were converted to BED format using coordinates on the updated nuclear-enriched diatom genome.
```bash
awk 'BEGIN{OFS="\t"} {
  qs=($7<$8?$7:$8)-1;
  qe=($7>$8?$7:$8);
  strand=($9<=$10?"+":"-");
  hit="hit_"NR;
  print $1, qs, qe, hit, $12, strand, $2, $9, $10, $3, $4, $11
}' blast_out/18_diatom_nuclear_v1_vs_phaeodactylum.filtered.tsv \
    > blast_out/blast_hits_on_18_diatom_nuclear_v1.bed
```
This command converts filtered BLASTN hits into BED-like coordinates on the nuclear-enriched diatom genome. The coordinate conversion allows BLASTN regions to be intersected with predicted diatom gene models.
Filtered BLASTN alignments were intersected with nuclear-filtered BRAKER4 gene models.
```bash
bedtools intersect \
    -a blast_out/18_diatom_nuclear_v1_genes.bed \
    -b blast_out/blast_hits_on_18_diatom_nuclear_v1.bed \
    -wa -wb \
    > blast_out/18_diatom_nuclear_v1_genes_with_phaeodactylum_hits.tsv
```
This command identifies nuclear diatom genes that overlap filtered *Phaeodactylum*-like BLASTN regions. The output contains both the diatom gene coordinates and the overlapping BLAST hit information.
This produced 2,269 gene-alignment overlaps. The number of unique diatom genes with *P. tricornutum*-like BLASTN hits was 1,492.
```bash
cut -f4 blast_out/18_diatom_nuclear_v1_genes_with_phaeodactylum_hits.tsv \
    | sort -u \
    | wc -l
```
This command counts the number of unique diatom gene IDs that overlap filtered BLASTN hits. It avoids counting the same gene multiple times if it overlaps more than one BLAST alignment.
The percentage of nuclear genes with detectable nucleotide similarity to *P. tricornutum* was calculated.
```bash
TOTAL_GENES=$(wc -l < blast_out/18_diatom_nuclear_v1_genes.bed)

HIT_GENES=$(cut -f4 blast_out/18_diatom_nuclear_v1_genes_with_phaeodactylum_hits.tsv \
    | sort -u \
    | wc -l)

awk -v h=$HIT_GENES -v t=$TOTAL_GENES 'BEGIN{
  print "Total nuclear genes:", t
  print "Genes with Phaeodactylum-like BLAST hits:", h
  print "Percent:", (h/t)*100 "%"
}'
```
These commands calculate the proportion of nuclear-filtered diatom genes with detectable nucleotide similarity to *P. tricornutum*. The result provides a gene-level summary of how much of the candidate nuclear gene set has BLASTN support against the reference diatom genome.
```text
Total nuclear genes: 15048
Genes with Phaeodactylum-like BLAST hits: 1492
Percent: 9.91494%
```
### 16.9 Select the best BLASTN hit per diatom gene
```bash
sort -k4,4 -k11,11nr blast_out/18_diatom_nuclear_v1_genes_with_phaeodactylum_hits.tsv \
    | awk 'BEGIN{OFS="\t"} !seen[$4]++' \
    > blast_out/18_diatom_nuclear_v1_genes_best_phaeodactylum_hit.tsv
```
This command sorts gene-overlap records by diatom gene ID and bitscore, then keeps the top-scoring BLASTN hit per gene. The resulting file provides one representative best nucleotide match for each diatom gene with a *Phaeodactylum*-like hit.
Mean nucleotide identity and mean alignment length were calculated for the best-hit set.
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
}' blast_out/18_diatom_nuclear_v1_genes_best_phaeodactylum_hit.tsv
```
This command calculates the mean percent identity and mean alignment length across the best-hit gene set. These summary statistics describe the strength and span of the best nucleotide-level matches to *P. tricornutum*.
```text
Number of best-hit genes: 1492
Mean percent identity: 72.9402
Mean alignment length: 596.025
```
### 16.10 Create final interpreted BLAST table
A simplified best-hit table was generated.
```bash
awk 'BEGIN{FS=OFS="\t";
print "diatom_contig","diatom_gene_start","diatom_gene_end","diatom_gene_id","diatom_gene_strand","blast_hit_id","bitscore","blast_strand","pt_contig","pt_start","pt_end","pident","aln_len","evalue"
}
{
print $1,$2,$3,$4,$6,$10,$11,$12,$13,$14,$15,$16,$17,$18
}' blast_out/18_diatom_nuclear_v1_genes_best_phaeodactylum_hit.tsv \
    > blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.clean.tsv
```
This command creates a cleaner best-hit table with the main diatom gene coordinates, BLAST hit coordinates, percent identity, alignment length, and e-value. The simplified table is easier to inspect and use for downstream interpretation than the full `bedtools intersect` output.
The *P. tricornutum* GFF3 was then converted to BED and intersected with BLASTN regions to link diatom best hits to PHATRDRAFT gene models.
```bash
awk 'BEGIN{OFS="\t"} {
  ss=($9<$10?$9:$10)-1;
  se=($9>$10?$9:$10);
  strand=($9<=$10?"+":"-");
  hit="hit_"NR;
  print $2, ss, se, hit, $12, strand, $1, $7, $8, $3, $4, $11
}' blast_out/18_diatom_nuclear_v1_vs_phaeodactylum.filtered.tsv \
    > blast_out/blast_hits_on_phaeodactylum.bed
```
This command converts filtered BLASTN alignments into BED-like coordinates on the *P. tricornutum* reference genome. These coordinates are needed to identify which annotated *Phaeodactylum* genes overlap the BLASTN hit regions.
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
This command converts *P. tricornutum* gene annotations from GFF3 format to BED format. It extracts gene IDs, gene names, and locus tags so BLAST hits can be linked to annotated PHATRDRAFT gene models.
Sequence IDs were normalized before `bedtools intersect`.
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
This command normalizes *Phaeodactylum* sequence IDs so that the BLAST hit coordinates and GFF-derived gene coordinates use matching chromosome or contig names. This prevents missed overlaps caused only by differences in sequence ID formatting.
```bash
bedtools intersect \
    -a blast_out/blast_hits_on_phaeodactylum.normalized.bed \
    -b blast_out/phaeodactylum_genes.bed \
    -wa -wb \
    > blast_out/blast_hits_overlapping_phaeodactylum_genes.tsv
```
This command intersects BLASTN hit regions on the *Phaeodactylum* genome with annotated *Phaeodactylum* genes. The output links each BLAST hit to any overlapping PHATRDRAFT gene model.
A BLAST hit to *P. tricornutum* gene mapping table was generated.
```bash
awk 'BEGIN{FS=OFS="\t";
print "blast_hit_id","pt_gene_id","pt_gene_name","pt_locus_tag"
}
{
  print $4,$16,$19,$20
}' blast_out/blast_hits_overlapping_phaeodactylum_genes.tsv \
    > blast_out/blast_hit_to_phaeodactylum_gene.tsv
```
This command extracts a compact mapping table between each BLAST hit ID and the overlapping *Phaeodactylum* gene identifiers. The table is later merged with the diatom best-hit table to add PHATRDRAFT gene annotations.
The mapping table was merged with the clean best-hit table using Python.
```bash
python scripts/merge_phaeodactylum_blast_hits.py
```
This Python script merges the clean diatom best-hit table with the *Phaeodactylum* gene mapping table. The output adds PHATRDRAFT gene IDs, gene names, and locus tags to the best BLASTN hit table.
A final cleaned version was generated by removing the `gene-` prefix from PHATRDRAFT identifiers and correcting the header.
```bash
awk 'BEGIN{FS=OFS="\t"}
NR==1 {print; next}
{
  gsub(/gene-/, "", $15);
  print
}' blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.with_PT_genes.tsv \
    > blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.with_PT_genes.cleanIDs.tsv

{
  echo -e "diatom_contig\tdiatom_gene_start\tdiatom_gene_end\tdiatom_gene_id\tdiatom_gene_strand\tblast_hit_id\tbitscore\tblast_strand\tpt_contig\tpt_start\tpt_end\tpident\taln_len\tevalue\tpt_gene_id\tpt_gene_name\tpt_locus_tag"
  tail -n +2 blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.with_PT_genes.cleanIDs.tsv
} > blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.with_PT_genes.cleanIDs.fixed.tsv
```
This command removes the `gene-` prefix from PHATRDRAFT identifiers and rewrites the table header. The final output is a clean, interpreted BLASTN best-hit table linking diatom genes to overlapping *P. tricornutum* gene models.
### 16.11 Final comparative-genomics outputs
```text
blast_out/18_diatom_nuclear_v1_vs_phaeodactylum.blastn.tsv
blast_out/18_diatom_nuclear_v1_vs_phaeodactylum.filtered.tsv
blast_out/18_diatom_nuclear_v1_genes_with_phaeodactylum_hits.tsv
blast_out/18_diatom_nuclear_v1_genes_best_phaeodactylum_hit.tsv
blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.clean.tsv
blast_out/blast_hit_to_phaeodactylum_gene.tsv
blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.with_PT_genes.cleanIDs.fixed.tsv
```
### 16.12 Summary
Pairwise BLASTN comparison between the updated nuclear-enriched diatom genome and the *Phaeodactylum tricornutum* reference genome identified 4,895 filtered nucleotide alignments using thresholds of ≥70% identity, alignment length ≥200 bp, and e-value ≤1e-10. These alignments overlapped 1,492 of 15,048 predicted nuclear genes, corresponding to 9.91% of the nuclear gene set. For the best hit per gene, the mean nucleotide identity was 72.94%, and the mean alignment length was 596 bp.
The best-hit regions were further mapped to annotated *P. tricornutum* PHATRDRAFT gene models. Of the 1,492 best-hit diatom genes, 1,446 overlapped annotated *P. tricornutum* gene models, representing 1,314 unique PHATRDRAFT genes.
This analysis provides a conservative nucleotide-level comparison between the candidate diatom genome and *P. tricornutum*. Because nucleotide similarity can be lost despite conservation at the protein level, this result should be treated as a genome-level similarity screen rather than a full orthology analysis.
</details> 

<details>
<summary><strong>17. Clean BRAKER4 isoform-level gene table construction</strong> - TransDecoder, DIAMOND, Average_TPM, and Phaeodactylum yes/no</summary>

This section describes the final clean table construction used after rebuilding the annotation workflow from raw input files. The final output is a BRAKER4 isoform-level table with one row per predicted protein isoform. It does not collapse BRAKER4 IDs, does not append GenBank-derived organelle rows, does not use all-hit TPM summaries, and does not retain detailed *Phaeodactylum tricornutum* hit columns.
The final table keeps:
```text
gene_id
transdecoder_orf_id
compartment
gene_length_bp
recommended_annotation
recommended_annotation_source
recommended_annotation_confidence
Average_TPM
in_Phaeodactylum_tricornutum
Swiss-Prot evidence
Bacillariophyta evidence
InterProScan evidence
GO terms
pathway annotations
AntiFam flags
BRAKER4 coordinates and length fields
```
### 17.1 Clean rebuild directory
A clean working directory was created to avoid carrying forward columns from older all-hit and GenBank-merging workflows.
```bash
mkdir -p /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/transdecoder_to_braker_ID_bridge/CLEAN_REBUILD_FROM_RAW
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/transdecoder_to_braker_ID_bridge/CLEAN_REBUILD_FROM_RAW

mkdir -p 01_input 02_diamond 03_best_hits 05_interproscan 06_combined_annotation 07_expression 08_phaeodactylum 09_final scripts logs slurm
```
### 17.2 Input files
The rebuild used the accepted BRAKER4 ET proteins and GFF3 file:
```text
01_input/diatom_predicted_proteins.fa
01_input/DL_diatom.braker4.ET.proteins.faa
01_input/DL_diatom.braker4.ET.gff3
```
The functional annotation layers were rebuilt from raw or core annotation files:
```text
02_diamond/DL_diatom_braker4_ET_vs_swissprot.tsv
02_diamond/DL_diatom_braker4_ET_vs_uniprot_bacillariophyta.tsv
05_interproscan/DL_diatom_braker4_ET_interproscan.tsv
06_combined_annotation/DL_diatom_antifam_flagged_proteins.tsv
```
The expression integration used cleaned TransDecoder peptides and the trusted transcriptome expression table:
```text
01_input/transcriptome_orfs.transdecoder.clean.pep
07_expression/master_with_custom_broad_categories.csv
```
The *Phaeodactylum tricornutum* comparison used only the raw gene-overlap file for yes/no lookup:
```text
08_phaeodactylum/18_diatom_nuclear_v1_genes_with_phaeodactylum_hits.tsv
```
The file `blast_hit_to_phaeodactylum_gene.tsv` was not required for the final clean table because no PHATRDRAFT IDs, gene names, bitscores, percent identities, or e-values were retained.
### 17.3 Rebuild functional annotation layers
Swiss-Prot, Bacillariophyta, and InterProScan outputs were parsed using the same scripts described in Section 14.
```bash
python scripts/make_swissprot_best_hits.py 2>&1 | tee swissprot_rebuild.log
python scripts/make_bacillariophyta_best_hits.py 2>&1 | tee bacillariophyta_rebuild.log
python scripts/summarize_interproscan.py 2>&1 | tee interproscan_rebuild.log
```
Each all-protein annotation file was checked to confirm one row per BRAKER4 protein isoform plus one header line:
```bash
wc -l 03_best_hits/DL_diatom_all_proteins_with_swissprot_annotation.tsv
wc -l 03_best_hits/DL_diatom_all_proteins_with_bacillariophyta_annotation.tsv
wc -l 06_combined_annotation/DL_diatom_all_proteins_with_interproscan_summary.tsv
```
Expected output:
```text
16948
16948
16948
```
The two AntiFam flags were retained as warning flags:
```text
g10893.t1    ANF00012    tRNA
g11404.t1    ANF00005    Antisense to 23S rRNA
```
The final functional annotation merge was run with:
```bash
python scripts/merge_functional_annotation_layers.py 2>&1 | tee merge_functional_annotation_rebuild.log
```
Output:
```text
06_combined_annotation/DL_diatom_master_functional_annotation.tsv
06_combined_annotation/DL_diatom_master_functional_annotation_for_manual_categories.tsv
06_combined_annotation/DL_diatom_master_functional_annotation_summary.txt
```
The master annotation table contained:
```text
16,947 BRAKER4 predicted proteins + 1 header line = 16,948 lines
```
### 17.4 Add BRAKER4 gene, CDS, and protein lengths
BRAKER4 coordinates and length fields were added from the accepted GFF3 file.
```bash
python scripts/07_add_BRAKER_lengths_clean.py 2>&1 | tee add_lengths_rebuild.log
```
Output:
```text
07_expression/DL_diatom_master_functional_annotation_with_lengths.tsv
```
Observed counts:

```text
Input annotation rows:      16,947
GFF3 transcript rows parsed: 16,947
Rows with gene length:       16,947
Rows missing gene length:         0
```
### 17.5 Align TransDecoder ORFs to BRAKER4 proteins
Cleaned TransDecoder peptides were aligned against the BRAKER4 ET protein set with DIAMOND BLASTP.
```bash
sbatch scripts/run_diamond_transdecoder_vs_braker_clean.slurm
```
The DIAMOND command retained up to 10 hits per TransDecoder ORF:
```bash
diamond blastp     -q 01_input/transcriptome_orfs.transdecoder.clean.pep     -d 07_expression/braker4_ET.proteins     -o 07_expression/transcriptome_ORFs_vs_BRAKER4_ET_proteins.tsv     --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore     --max-target-seqs 10     --evalue 1e-5     --threads 16
```
Raw DIAMOND output:
```text
118,472 hit rows
39,935 unique TransDecoder ORFs with at least one BRAKER hit
14,473 unique BRAKER proteins hit by at least one reported ORF hit
```
### 17.6 Create one best ORF-to-BRAKER mapping per TransDecoder ORF
The raw DIAMOND output was parsed to calculate ORF coverage and BRAKER protein coverage. One best BRAKER4 hit was retained per TransDecoder ORF.
```bash
python scripts/08_make_best_ORF_to_BRAKER_mapping_clean.py 2>&1 | tee make_best_mapping_rebuild.log
```
Output files:
```text
07_expression/transdecoder_ORFs_vs_BRAKER4_all_hits_with_coverage.tsv
07_expression/transdecoder_ORFs_to_BRAKER4_best_hit_per_ORF.tsv
```
Observed counts:
```text
Raw DIAMOND hit rows: 118,472
Unique TransDecoder ORFs with at least one BRAKER hit: 39,935
Unique BRAKER proteins hit by best ORF mappings: 12,277
```
The best-hit file contained:
```text
39,935 ORF mappings + 1 header line = 39,936 lines
```
### 17.7 Add TransDecoder ORF ID and Average_TPM only
The trusted transcriptome table was copied into the clean rebuild folder:
```bash
cp ../master_with_custom_broad_categories.csv 07_expression/
```
Only two columns were used from this file:
```text
orf
Average_TPM
```
Expression values were merged using the best ORF-to-BRAKER mapping. If multiple TransDecoder ORFs mapped best to the same BRAKER4 protein, the ORF with the highest valid `Average_TPM` was retained. ORFs without a valid numeric `Average_TPM` were excluded from the expression-supported ORF column so that `transdecoder_orf_id` and `Average_TPM` remain consistent.
```bash
python scripts/09_add_ONLY_Average_TPM_clean.py 2>&1 | tee add_average_tpm_rebuild.log
```
Output:
```text
07_expression/DL_diatom_master_functional_annotation_lengths_Average_TPM.tsv
```
Final expression integration counts:
```text
Annotation rows: 16,947
Best ORF-to-BRAKER rows: 39,935
BRAKER proteins with Average_TPM: 12,276
BRAKER proteins without Average_TPM: 4,671
BRAKER proteins with TransDecoder ORF ID: 12,276
BRAKER proteins without TransDecoder ORF ID: 4,671
```
The output retained all 16,947 BRAKER4 isoforms and added:
```text
transdecoder_orf_id
Average_TPM
```
No all-hit TPM sums, means, hit counts, or mapped ORF lists were retained.
### 17.8 Add Phaeodactylum yes/no only
The raw *Phaeodactylum tricornutum* comparison file contains detailed BLASTN overlap fields, but the final clean table used it only as a yes/no lookup.
```text
08_phaeodactylum/18_diatom_nuclear_v1_genes_with_phaeodactylum_hits.tsv
```
The file contains gene roots such as:
```text
g8700
g8702
g14082
```
The final BRAKER4 table keeps isoform IDs such as:
```text
g8700.t1
g8702.t1
```
Therefore, the script used the root only for lookup while preserving the full BRAKER4 isoform ID in the final table. The only PT-derived output column is:
```text
in_Phaeodactylum_tricornutum
```
Detailed PT columns such as `blast_hit_id`, `pt_gene_id`, `pt_gene_name`, `bitscore`, `pident`, and `evalue` were not retained.
### 17.9 Create final clean BRAKER4 isoform-level table
```bash
python scripts/10_make_FINAL_clean_BRAKER_isoform_table.py 2>&1 | tee make_final_clean_table.log
```
Output files:
```text
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table.tsv
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table_sorted_by_Average_TPM.tsv
```
Both final files contained:
```text
16,947 BRAKER4 isoform rows + 1 header line = 16,948 lines
```
Final compartment counts:
```text
nuclear        16,873
plastid_like       72
mito_like           2
```
The compartment labels were assigned from BRAKER4 contig IDs only. No GenBank organelle rows were appended.
Final *Phaeodactylum tricornutum* yes/no counts:
```text
no     15,225
yes     1,722
```
The script also confirmed:
```text
PT detail columns found: []
```
### 17.10 Final table structure
The final clean table begins with:
```text
gene_id
transdecoder_orf_id
compartment
gene_length_bp
recommended_annotation
recommended_annotation_source
recommended_annotation_confidence
Average_TPM
in_Phaeodactylum_tricornutum
```
`gene_id` is the original BRAKER4 protein isoform ID and is not collapsed to the gene root. `transdecoder_orf_id` is the TransDecoder ORF that contributed the expression value for that BRAKER4 isoform. This column is included for traceability between the transcriptome assembly and the BRAKER4 gene model.

The `transdecoder_orf_id` column is populated only when the mapped ORF has a valid numeric `Average_TPM`. If a TransDecoder ORF maps to a BRAKER4 protein but does not have a usable `Average_TPM`, both `transdecoder_orf_id` and `Average_TPM` are left blank in the final output. This prevents the table from showing an ORF ID that did not contribute an expression value.

If multiple TransDecoder ORFs map best to the same BRAKER4 isoform, the ORF with the highest valid `Average_TPM` is retained. TPM values are not summed, averaged, or duplicated across all DIAMOND hits.
The final table therefore reports:
```text
one BRAKER4 isoform row
one representative TransDecoder ORF ID when expression support is available
one Average_TPM value
one Phaeodactylum tricornutum yes/no value
```
The remaining columns retain the functional annotation evidence and BRAKER4 coordinate fields:
```text
has_any_annotation
has_swissprot_hit
has_bacillariophyta_hit
has_interproscan_hit
has_antifam_flag
antifam_accession
antifam_description
antifam_flag
Swiss-Prot fields
Bacillariophyta fields
InterProScan fields
GO terms
pathway annotations
contig_id
gene_model_start
gene_model_end
strand
cds_length_bp
protein_length_aa
```
### 17.11 Final output files
The main final table is:
```text
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table.tsv
```
The expression-sorted version is:
```text
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table_sorted_by_Average_TPM.tsv
```
The sorted file is useful for manual inspection and pathway curation because highly expressed genes appear first. The unsorted file should be retained as the master output.
#### Logic
This clean workflow links BRAKER4-predicted diatom proteins to functional annotation, transcriptome expression, and comparative genomic evidence without inflating expression values or adding unnecessary PT detail columns. `Average_TPM` is transferred through a best TransDecoder ORF-to-BRAKER4 mapping, and the TransDecoder ORF ID used for expression is retained next to the BRAKER4 gene ID. *Phaeodactylum tricornutum* similarity is summarized only as yes/no, and the final table remains one row per BRAKER4 protein isoform.

</details>
