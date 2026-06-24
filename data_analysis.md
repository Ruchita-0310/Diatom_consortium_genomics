# Diatom Consortia: Metagenomic and Metatranscriptomic Pipeline
This repository documents the workflow used to assemble, polish, bin, classify, annotate, and compare genomes and transcriptomes from a diatom-associated microbial consortium. The pipeline combines long-read metagenomic assembly, short-read polishing, metagenomic binning, contig-level taxonomic screening, organelle identification, transcriptome analysis, BRAKER4 ET gene prediction, nuclear-enriched genome generation, functional annotation, expression integration, comparison with the reference diatom *Phaeodactylum tricornutum*, and Hi-C read mapping/contact-network analysis.
The final gene table is a clean BRAKER4 isoform-level table with one row per predicted protein isoform.

---

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
   ↓
Hi-C read mapping to polished whole assembly
   ↓
Contig-level Hi-C representation summary
   ↓
Hi-C proximal-ligation contact network
```
---

## Software and environments
The workflow used Conda environments, Singularity containers, and local HPC modules depending on software availability.

| Step                             | Tools                                                                                                                                   |
| -------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| Assembly and polishing           | Flye, Medaka, BWA-MEM, Polypolish, Pypolca                                                                                              |
| Read mapping and coverage        | minimap2, samtools, bedtools, seqkit                                                                                                    |
| Assembly quality                 | BUSCO, QUAST/MetaQUAST                                                                                                                  |
| Binning and bin quality          | MetaBAT2, CheckM2                                                                                                                       |
| Taxonomy and abundance           | GTDB-Tk, MetaEuk, CoverM                                                                                                                |
| Organelle identification         | MetaQUAST, minimap2, bedtools, seqkit, GeSeq/OGDRAW                                                                                     |
| Phylogenetics                    | Clustal Omega, TrimAl, IQ-TREE 2                                                                                                        |
| Transcriptomics                  | Nextflow, nf-core/metatdenovo, TransDecoder, Barrnap, STAR                                                                              |
| Genome annotation                | BRAKER4, GeneMark-ET, AUGUSTUS, TSEBRA, STAR, BUSCO/compleasm                                                                           |
| Functional annotation            | DIAMOND, UniProtKB/Swiss-Prot, UniProtKB Bacillariophyta, InterProScan, Pfam, PANTHER, Gene3D, CDD, SMART, SUPERFAMILY, ProSite, Python |
| Expression integration           | DIAMOND, Python, pandas, TransDecoder ORFs, Average_TPM table                                                                           |
| Comparative genomics             | NCBI Datasets, BLASTN, bedtools, Python                                                                                                 |
| Hi-C mapping and contact network | FastQC, MultiQC, BWA-MEM, samtools, seqkit, YaHS, awk, Python                                                                           |

---

## Repository structure for scripts
Custom Python scripts are stored in the scripts/ directory rather than embedded directly in this markdown workflow. Scripts are numbered in the order they are used in the analysis.
```text
scripts/
├── 01_classify_metaeuk_contigs.py
├── 02_make_swissprot_best_hits.py
├── 03_make_bacillariophyta_best_hits.py
├── 04_summarize_interproscan.py
├── 05_merge_functional_annotation_layers.py
├── 06_merge_phaeodactylum_blast_hits.py
├── 07_add_BRAKER_lengths_clean.py
├── 08_make_best_ORF_to_BRAKER_mapping_clean.py
├── 09_add_ONLY_Average_TPM_clean.py
├── 10_make_FINAL_clean_BRAKER_isoform_table.py
└── 11_make_hic_network_files.py
```
Script purposes:
```text
01_classify_metaeuk_contigs.py
  Classifies contigs using MetaEuk ORF-level taxonomy and assigns each contig to a final category.

02_make_swissprot_best_hits.py
  Parses Swiss-Prot DIAMOND output, calculates coverage, assigns confidence classes, and writes best-hit annotation tables.

03_make_bacillariophyta_best_hits.py
  Parses UniProtKB Bacillariophyta DIAMOND output using the same confidence framework as the Swiss-Prot parser.

04_summarize_interproscan.py
  Collapses raw InterProScan TSV output to one row per predicted protein and summarizes domains, GO terms, pathways, and database sources.

05_merge_functional_annotation_layers.py
  Merges Swiss-Prot, Bacillariophyta, InterProScan, and AntiFam evidence into one BRAKER4 functional annotation table.

06_merge_phaeodactylum_blast_hits.py
  Adds overlapping Phaeodactylum gene information to the cleaned BLASTN best-hit table.

07_add_BRAKER_lengths_clean.py
  Adds BRAKER4 coordinates, contig IDs, strand, gene length, CDS length, and protein length to each isoform.

08_make_best_ORF_to_BRAKER_mapping_clean.py
  Parses TransDecoder ORF versus BRAKER4 DIAMOND output and keeps one best BRAKER4 hit per TransDecoder ORF.

09_add_ONLY_Average_TPM_clean.py
  Adds the matching TransDecoder ORF ID and Average_TPM value to BRAKER4 isoforms using the best ORF-to-BRAKER mapping.

10_make_FINAL_clean_BRAKER_isoform_table.py
  Creates the final clean BRAKER4 isoform-level table with one row per protein isoform.

11_make_hic_network_files.py
  Converts the Hi-C contig-contact table into GEXF and GraphML network files.
```
---

# Analysis workflow
Click each section to expand the commands, notes, and outputs.

---

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
This command assembled Nanopore reads using Flye in metagenome mode. The output assembly was used as the starting point for polishing, read mapping, binning, organelle screening, and genome annotation.

</details>

---

<details>
<summary><strong>2. Read mapping and assembly support</strong> - minimap2 and samtools</summary>

Short reads were mapped to the Nanopore assembly to assess read support, mapping rate, and contig-level coverage.
```bash
minimap2 -ax sr \
    /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/guppy_flye_assembly.fasta \
    Diatoms_merged.fastq.gz \
    > sr_alignment.sam
```
This command maps short reads to the Nanopore assembly using the short-read preset in minimap2.
```bash
samtools view -S -b sr_alignment.sam > alignment.bam
samtools sort alignment.bam -o alignment_sorted.bam
samtools index alignment_sorted.bam
```
These commands convert the SAM file to BAM, sort the alignments by coordinate, and index the sorted BAM file.
```bash
samtools flagstat alignment_sorted.bam > mapping_stats.txt
```
This command summarizes the number and proportion of reads that mapped to the assembly.
```bash
samtools idxstats alignment_sorted.bam \
    | sort -k3,3rn \
    > sr_all_nanopore_hits.tsv
```
This command reports mapped-read counts per contig and sorts contigs by the number of mapped reads.
```bash
samtools depth alignment_sorted.bam > sr_depth.txt
```
This command calculates per-base read depth across the assembly.
```bash
awk '{sum[$1]+=$3; count[$1]++} END {for (c in sum) print c, sum[c]/count[c]}' sr_depth.txt \
    | sort -k2,2nr \
    > sr_mean_depth.tsv
```
This command calculates mean read depth per contig and sorts contigs from highest to lowest coverage.
The resulting files were used to evaluate mapping rate, contig-level coverage, and short-read support across the assembly.

</details>

---

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
This command uses Nanopore reads to correct consensus errors in the Flye assembly.
### 3.2 Short-read alignment for Polypolish
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
These commands align each short-read file separately to the Medaka-polished assembly. The resulting SAM files were used as input for Polypolish.
### 3.3 Polypolish filtering and polishing
```bash
polypolish filter \
    --in1 alignments_1.sam \
    --in2 alignments_2.sam \
    --out1 filtered_1.sam \
    --out2 filtered_2.sam
```
This command filters paired-end short-read alignments into the format expected by Polypolish.
```bash
polypolish polish \
    /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/medaka_euk_polished/consensus.fasta \
    filtered_1.sam \
    filtered_2.sam \
    > sr_poly.fasta
```
This command uses short-read alignments to polish the Medaka-corrected assembly.
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
This command performs an additional short-read polishing step using Pypolca. The final corrected assembly was used for downstream binning, organelle identification, and gene annotation.
### 3.5 BUSCO assessment of the polished assembly
```bash
busco \
    -i pypolca_corrected.fasta \
    -l busco_downloads/lineages/stramenopiles_odb10 \
    -o busco_report \
    -m genome
```
This command evaluates assembly completeness using the stramenopile BUSCO marker set.

</details>

---

<details>
<summary><strong>4. Metagenomic binning</strong> - MetaBAT2</summary>

MetaBAT2 was used to recover genome bins from the polished assembly using Nanopore read coverage.
### 4.1 Map Nanopore reads to the polished assembly
```bash
minimap2 -ax map-ont -t 16 \
    1_sr_pypolca_output/pypolca_corrected.fasta \
    pass_trim.fastq.gz \
    | samtools view -@ 16 -bS - \
    | samtools sort -@ 16 -m 10G -o aligned_reads.sorted.bam
```
This command maps Nanopore reads back to the polished assembly and creates a coordinate-sorted BAM file for coverage estimation.
```bash
samtools index -@ 16 aligned_reads.sorted.bam
```
This command indexes the sorted BAM file so it can be used by downstream coverage tools.
### 4.2 Generate contig depth file
```bash
jgi_summarize_bam_contig_depths \
    --outputDepth depth.txt \
    --percentIdentity 85 \
    aligned_reads.sorted.bam
```
This command calculates contig-level depth from the Nanopore read mapping file. The depth file was used by MetaBAT2 for binning.
### 4.3 Run MetaBAT2
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
This command bins the polished metagenomic assembly using contig sequence composition and read-depth information.
The resulting bins were used for quality assessment and taxonomic classification.

</details>

---

<details>
<summary><strong>5. Bin quality assessment</strong> - CheckM2</summary>

CheckM2 was used to estimate completeness and contamination of recovered genome bins.
```bash
checkm2 predict \
    --threads 16 \
    --input 2_metabat2_bins/ \
    --output_directory 3_checkm2_results
```
This command estimates genome-bin completeness and contamination. The output was used to assess bin quality before downstream taxonomy and interpretation.

</details>

---

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
This command classifies bacterial and archaeal genome bins using the GTDB-Tk workflow.
GTDB-Tk was used for bacterial and archaeal genome bins. Because the consortium also contained a dominant eukaryotic diatom, MetaEuk-based ORF taxonomy was used as an additional contig-level screen for eukaryotic, bacterial, ambiguous, and unclassified contig fractions.

</details>

---

<details>
<summary><strong>7. MetaEuk-based contig classification</strong> - MetaEuk and custom Python script</summary>

MetaEuk ORF-level taxonomic assignments were used to classify contigs across recovered bins. This step was added because bin-level bacterial taxonomy alone does not resolve eukaryotic contigs or mixed bins in a diatom-associated consortium.

Because MetaEuk uses last common ancestor assignments, organelle-derived sequences can be assigned to bacterial lineages. To account for this, the contig classification grouped direct eukaryotic hits together with mitochondrial and chloroplast-derived signatures when calculating the eukaryotic score.

Mitochondrial-like hits were identified using:
```text
o_Rickettsiales
o__Rickettsiales
```
Chloroplast-like hits were identified using:
```text
p_Cyanobacteria
```
### 7.1 Input files
```text
metaeuk_output_polyp_taxonomy_tax_per_pred.tsv
contig_to_bin.txt
```
`metaeuk_output_polyp_taxonomy_tax_per_pred.tsv` contains ORF-level MetaEuk taxonomic assignments. `contig_to_bin.txt` links contig IDs to bin IDs.
The MetaEuk table was expected to contain:
```text
Contig_ID
Classification
```
### 7.2 ORF-level labels

Each predicted ORF was assigned to one of the following labels:

| Label                | Rule                                                                           |
| -------------------- | ------------------------------------------------------------------------------ |
| Eukaryota            | `Classification` contains `d_Eukaryota`                                        |
| Mitochondria-derived | `Classification` contains `o_Rickettsiales` or `o__Rickettsiales`              |
| Chloroplast-derived  | `Classification` contains `p_Cyanobacteria`                                    |
| Bacteria             | `Classification` contains `d_Bacteria`, excluding organelle-derived categories |
| Ambiguous            | `Classification` is exactly `_cellular organisms`                              |
| Other                | all other biological hits, including Archaea or viruses                        |
| Unclassified         | no MetaEuk classification available                                            |

### 7.3 Final contig-level classification rule
A contig was classified as `Eukaryota` if more than 30% of its biological ORF assignments were eukaryotic, mitochondrial, or chloroplast-derived:
```text
(Eukaryota + Mitochondria-derived + Chloroplast-derived) / Total biological ORFs > 0.30
```
where:
```text
Total biological ORFs = Eukaryota + Mitochondria-derived + Chloroplast-derived + Bacteria + Ambiguous + Other
```

If the contig did not meet the >30% eukaryotic/organelle threshold, the remaining labels were assigned based on the dominant biological category. Contigs without a clear dominant category were labeled ambiguous.
### 7.4 Python script
The full Python script is saved in:
```text
scripts/01_classify_metaeuk_contigs.py
```
Run the script from the directory containing the MetaEuk taxonomy table and `contig_to_bin.txt`. 

This script classifies each contig using ORF-level MetaEuk assignments and writes the final contig-level classification table.
### 7.5 Output
```text
contig_classification_final_priority.csv
```
The output table contains the bin name, contig ID, ORF-level category counts, and final contig classification.

</details>

---

<details>
<summary><strong>8. Genome coverage and relative abundance</strong> - CoverM</summary>

CoverM was used to calculate coverage and relative abundance of bacterial genome bins using paired-end short reads.
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
This command maps short reads to bacterial genome bins and calculates mean coverage, relative abundance, and covered fraction for each bin.
The output table was:
```text
bac_output_coverm.tsv
```

</details>

---

<details>
<summary><strong>9. 18S rRNA phylogenetic analysis</strong> - Clustal Omega, TrimAl, and IQ-TREE 2</summary>
A phylogenetic tree was generated from 18S rRNA sequences using Clustal Omega, TrimAl, and IQ-TREE 2.
```bash
cat *.fasta > 18S_new.fasta
```
This command combines individual 18S FASTA files into one input file for alignment.

```bash
clustalo \
    -i 18S_new.fasta \
    -o 18S_aligned.fasta
```
This command aligns the combined 18S rRNA sequences.
```bash
trimal \
    -in 18S_aligned.fasta \
    -out 18S_trimmed.fasta \
    -automated1
```
This command trims poorly aligned regions from the 18S alignment.
```bash
/home/ruchita.solanki/iqtree-2.2.2.7-Linux/bin/iqtree2 \
    -s 18S_trimmed.fasta \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -nt AUTO
```
This command infers a phylogenetic tree, selects the best-fit model using ModelFinder, and estimates branch support using ultrafast bootstrap and SH-aLRT values.

</details>

---

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
This command runs the nf-core/metatdenovo transcriptome workflow using SPAdes for transcript assembly and TransDecoder for ORF prediction.

The workflow generated transcript assemblies, predicted ORFs, and transcript-level annotation files that were later used for expression integration.

</details>

---

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
```
This command identifies eukaryotic rRNA transcripts from the transcriptome assembly.

```bash
barrnap \
    --kingdom bac \
    spades.transcripts.fa \
    --outseq bac_transcript_rRNA.fna \
    > diatom_bac_rRNA.gff
```
This command identifies bacterial rRNA transcripts from the transcriptome assembly.

```bash
barrnap \
    --kingdom mito \
    spades.transcripts.fa \
    --outseq mito_transcript_rRNA.fna \
    > diatom_mito_rRNA.gff
```
This command identifies mitochondrial rRNA transcripts from the transcriptome assembly.

The resulting FASTA and GFF files were used to inspect rRNA transcript origin in the consortium transcriptome.

</details>

---

<details>
<summary><strong>12. Organelle genome identification</strong> - MetaQUAST</summary>

Organelle contigs were identified by comparing the polished assembly against reference mitochondrial and chloroplast genomes.
```text
Mitogenome reference:   MT742552
Chloroplast reference:  MT742551
```
### 12.1 MetaQUAST comparison of a candidate bin
```bash
metaquast.py \
    8_diatom.fasta \
    -R /work/ebg_lab/eb/diatom_consortia/organelle/ref/ \
    -o ./8_metaquast_output
```
This command compares the candidate diatom bin against the chloroplast and mitochondrial reference genomes.
### 12.2 MetaQUAST comparison of the polished whole assembly
```bash
metaquast.py \
    /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_output/pypolca_corrected.fasta \
    -R /work/ebg_lab/eb/diatom_consortia/organelle/ref/ \
    -o ./whole_metaquast_output
```
This command compares the polished whole assembly against the organelle references to identify candidate chloroplast and mitochondrial contigs outside the binned assembly.

The MetaQUAST output was used to identify candidate chloroplast and mitochondrial contigs for downstream organelle genome refinement and annotation.

</details>

---

<details>
<summary><strong>13. Diatom genome annotation with BRAKER4 ET mode</strong> - BRAKER4, STAR, GeneMark-ET, AUGUSTUS, and TSEBRA</summary>

Gene models were generated with BRAKER4 using a diatom genome assembly and RNA-seq evidence. The genome was not soft-masked before annotation; therefore, repeat masking was performed internally within the BRAKER4 workflow using RepeatModeler, RepeatMasker, and TRF.
The final accepted run used ET mode, meaning that gene prediction was based on RNA-seq evidence only.
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
These commands downloaded the BRAKER4 workflow and pulled the BRAKER container used for the annotation run.
GeneMark requires a license key. The key was stored at:
```bash
/home/ruchita.solanki/.gm_key
```
### 13.2 STAR genome indexing
```bash
STAR \
    --runThreadN 24 \
    --runMode genomeGenerate \
    --genomeDir /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index \
    --genomeFastaFiles /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta \
    --genomeSAindexNbases 10
```
This command builds a STAR genome index for the diatom assembly so that RNA-seq reads can be aligned to the genome.
### 13.3 STAR RNA-seq alignment
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
This command aligns RNA-seq reads to the indexed diatom genome and writes a coordinate-sorted BAM file for BRAKER4.
```bash
samtools index /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam
```
This command indexes the STAR BAM file so BRAKER4 can read the RNA-seq alignments.
### 13.4 Genome input and repeat-masking check
The genome assembly used for STAR indexing and BRAKER4 was:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
```
The assembly contained 3,010 contigs and had a total length of approximately 82.17 Mbp.
```bash
seqkit stats /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
```
This command reports the number of contigs and total genome length before annotation.
```bash
grep -v "^>" /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta \
    | grep -q '[a-z]' && echo "soft-masked" || echo "not soft-masked"
```
This command checks whether the genome FASTA contains lowercase bases, which would indicate soft masking.
Output:
```text
not soft-masked
```
Because the genome was not pre-masked, the `genome_masked` column in `samples.csv` was left empty and internal repeat masking was enabled in `config.ini`.
### 13.5 BRAKER4 sample file
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
This command checks that the header and sample row contain the same number of comma-separated fields.
Expected output:
```text
1 13
2 13
```
### 13.6 BRAKER4 configuration
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
### 13.7 Snakemake dry run
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
This command performs a dry run of the BRAKER4 workflow without executing jobs. It was used to confirm that the expected ET-mode rules would run.
The dry run included:
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
### 13.8 BRAKER4 execution
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
This command runs the full BRAKER4 ET workflow. BRAKER4 performed internal repeat masking, RNA-seq hint generation, GeneMark-ET prediction, AUGUSTUS training and prediction, TSEBRA refinement, BUSCO assessment, and final result collection.

The ET-mode run produced:
```text
StringTie transcripts: 9,000
RNA-seq intron hints: 19,114
```
### 13.9 Final BRAKER4 outputs
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
These commands create uncompressed GFF3, GTF, protein FASTA, and CDS FASTA files for downstream annotation and table construction.

No additional TSEBRA run was required because TSEBRA refinement was included within BRAKER4.
### 13.10 Annotation statistics
```bash
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET

grep -c $'\tgene\t' DL_diatom.braker4.ET.gff3
grep -c $'\ttranscript\t' DL_diatom.braker4.ET.gff3
grep -c $'\tCDS\t' DL_diatom.braker4.ET.gff3

grep -c "^>" DL_diatom.braker4.ET.proteins.faa
grep -c "^>" DL_diatom.braker4.ET.cds.fna

seqkit stats DL_diatom.braker4.ET.proteins.faa DL_diatom.braker4.ET.cds.fna
```
These commands count genes, transcripts, CDS features, protein sequences, and CDS sequences in the final BRAKER4 output.
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
Predicted proteins were checked for internal stop codons:
```bash
grep -n "\*" DL_diatom.braker4.ET.proteins.faa | head
```
This command searches the predicted protein FASTA for internal stop codons.
No internal stop codons were detected.
### 13.11 BUSCO assessment of the final protein set
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
This command evaluates completeness of the predicted protein set using the `stramenopiles_odb12` BUSCO marker set.
The final predicted protein set produced:
```text
C:84.8%[S:81.1%,D:3.7%],F:1.9%,M:13.3%,n=697
```
### 13.12 Annotation acceptance
The final BRAKER4 ET annotation was accepted for downstream analysis because it produced a plausible gene set for the diatom genome assembly, showed no internal stop codon issues in the predicted protein FASTA, and recovered 84.8% of the `stramenopiles_odb12` BUSCO protein set with low duplication.
Final accepted annotation files:
```text
DL_diatom.braker4.ET.gff3
DL_diatom.braker4.ET.gtf
DL_diatom.braker4.ET.proteins.faa
DL_diatom.braker4.ET.cds.fna
```
### 13.13 Rationale for ET mode instead of ETP
BRAKER4 was initially tested in ETP mode, which combines RNA-seq evidence with protein evidence. However, GeneMark-ETP failed during model training. Although protein-supported alignments were generated, the GeneMark-ETP training set did not produce valid gene and transcript models.
The failed run reported:
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
Because GeneMark-ETP did not complete successfully, the annotation was rerun in ET mode using RNA-seq evidence only. This avoided the failed protein-dependent GeneMark-ETP training step while retaining transcript evidence from the coordinate-sorted STAR BAM file.

The final successful workflow used GeneMark-ET, AUGUSTUS, and TSEBRA, with `protein_fasta` left empty in `samples.csv` and `mode = et` specified in `config.ini`.

</details>

---

<details>
<summary><strong>14. Functional annotation of BRAKER4-predicted proteins</strong> - DIAMOND, Swiss-Prot, Bacillariophyta, InterProScan, and AntiFam</summary>

After accepting the BRAKER4 ET annotation, the predicted protein set was used for downstream functional annotation. The final annotation strategy used three complementary evidence layers: curated Swiss-Prot homology, diatom-focused UniProtKB Bacillariophyta homology, and InterProScan domain/family annotation.
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
This command makes the accepted BRAKER4 protein set available in the functional annotation directory without copying the original file.
### 14.1 Annotation strategy
Swiss-Prot was used as the conservative curated layer because entries are manually reviewed, although fewer proteins are expected to receive hits. A UniProtKB Bacillariophyta database was added to improve detection of diatom-specific homologs. InterProScan was used to identify conserved domains, protein families, GO terms, and pathway signatures.
eggNOG, KEGG, COG, and dbCAN were not used in the final workflow. The final strategy prioritized eukaryote- and diatom-focused annotation rather than broad prokaryotic orthology or specialized carbohydrate-active enzyme classification.
### 14.2 Swiss-Prot database setup
Swiss-Prot was stored in the home directory to avoid filling the project working directory.
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
```bash
diamond makedb \
    --in $HOME/databases/swissprot/raw/uniprot_sprot.fasta.gz \
    --db $HOME/databases/swissprot/diamond/uniprot_sprot.dmnd
```
This command builds a DIAMOND database from the Swiss-Prot FASTA file.
```bash
ln -sfn $HOME/databases/swissprot 00_databases/swissprot_home
```
This command links the home-directory Swiss-Prot database into the project directory.
### 14.3 DIAMOND search against Swiss-Prot
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
This command searches BRAKER4-predicted proteins against Swiss-Prot using DIAMOND BLASTP. Up to five candidate hits per protein were retained for downstream best-hit selection and confidence filtering.

Swiss-Prot result summary:

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
### 14.4 UniProtKB Bacillariophyta database setup
A diatom-focused UniProtKB database was created using the Bacillariophyta taxonomic group.
```bash
mkdir -p $HOME/databases/uniprot_bacillariophyta/raw
mkdir -p $HOME/databases/uniprot_bacillariophyta/diamond
```
```bash
curl -L --retry 5 --retry-delay 10 \
    -o $HOME/databases/uniprot_bacillariophyta/raw/uniprotkb_bacillariophyta_taxid2836.fasta.gz \
    "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28taxonomy_id%3A2836%29"
```
This command downloads UniProtKB protein sequences assigned to Bacillariophyta.
```bash
diamond makedb \
    --in $HOME/databases/uniprot_bacillariophyta/raw/uniprotkb_bacillariophyta_taxid2836.fasta.gz \
    --db $HOME/databases/uniprot_bacillariophyta/diamond/uniprotkb_bacillariophyta_taxid2836.dmnd
```
This command builds a DIAMOND database from the Bacillariophyta UniProtKB FASTA file.
```bash
ln -sfn $HOME/databases/uniprot_bacillariophyta 00_databases/uniprot_bacillariophyta_home
```
This command links the home-directory Bacillariophyta database into the project directory.
### 14.5 DIAMOND search against UniProtKB Bacillariophyta
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
This command searches BRAKER4-predicted proteins against the diatom-focused Bacillariophyta UniProtKB database.

Bacillariophyta result summary:
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
### 14.6 Best-hit parsing and confidence filtering
Raw DIAMOND outputs were parsed into best-hit tables. For each predicted protein, query coverage and subject coverage were calculated, UniProt identifiers were parsed, and protein name, organism, gene name, taxon ID, and protein evidence fields were extracted.
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
### 14.7 Swiss-Prot best-hit parsing script
The full Python script is saved in:
```text
scripts/02_make_swissprot_best_hits.py
```
Run the script with:
```bash
conda activate swissprot_annot
python scripts/02_make_swissprot_best_hits.py
```
This script parses the Swiss-Prot DIAMOND output, calculates coverage, assigns confidence classes, selects the best hit per protein, and writes both best-hit and all-protein annotation tables.

Main parsed Swiss-Prot outputs:
```text
03_best_hits/DL_diatom_swissprot_all_hits_with_coverage.tsv
03_best_hits/DL_diatom_swissprot_best_hits.tsv
03_best_hits/DL_diatom_swissprot_best_hits_strict.tsv
03_best_hits/DL_diatom_all_proteins_with_swissprot_annotation.tsv
03_best_hits/DL_diatom_swissprot_annotation_summary.txt
```
### 14.8 Bacillariophyta best-hit parsing script
The full Python script is saved in:
```text
scripts/03_make_bacillariophyta_best_hits.py
```
Run the script with:
```bash
conda activate swissprot_annot
python scripts/03_make_bacillariophyta_best_hits.py
```
This script parses the Bacillariophyta DIAMOND output using the same coverage, ranking, and confidence framework used for Swiss-Prot.
Main parsed Bacillariophyta outputs:
```text
03_best_hits/DL_diatom_bacillariophyta_all_hits_with_coverage.tsv
03_best_hits/DL_diatom_bacillariophyta_best_hits.tsv
03_best_hits/DL_diatom_bacillariophyta_best_hits_strict.tsv
03_best_hits/DL_diatom_all_proteins_with_bacillariophyta_annotation.tsv
03_best_hits/DL_diatom_bacillariophyta_annotation_summary.txt
```
### 14.9 InterProScan annotation
InterProScan was installed in the home tools directory and linked into the project:
```bash
ln -sfn $HOME/tools/interproscan/current 00_databases/interproscan_home
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
This script runs InterProScan on the full BRAKER4 protein set and writes a TSV output file. GO terms and pathway annotations were requested where available.

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
InterProScan result summary:
```text
Raw InterProScan rows:              102,153
Proteins with InterProScan hits:     13,106 / 16,947
Percent with InterProScan hits:      77.34%
```
The raw InterProScan file was:
```text
05_interproscan/DL_diatom_braker4_ET_interproscan.tsv
```
### 14.10 InterProScan summary by protein
The raw InterProScan TSV was collapsed into one row per protein using a custom Python script:
```text
scripts/04_summarize_interproscan.py
```
Run the script with:
```bash
conda activate swissprot_annot
python scripts/04_summarize_interproscan.py
```
This script summarizes raw InterProScan matches by protein so that domain, family, GO, and pathway evidence can be merged with the homology-based annotation tables.
Output files:
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
InterProScan database contribution summary:
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
These counts represent raw annotation rows rather than unique proteins, because individual proteins can contain multiple domains or match multiple databases.
### 14.11 AntiFam screening
InterProScan reported two AntiFam matches:
```text
protein_id   AntiFam accession   AntiFam description
g10893.t1    ANF00012            tRNA
g11404.t1    ANF00005            Antisense to 23S rRNA
```
AntiFam matches were treated as warning flags rather than functional annotations. These two proteins were marked as potentially spurious predictions or RNA-associated ORFs.
The AntiFam hits were inspected with:
```bash
awk -F '\t' '$4=="AntiFam"' 05_interproscan/DL_diatom_braker4_ET_interproscan.tsv | column -t
```
This command extracts AntiFam rows from the InterProScan output for manual inspection.
### 14.12 Master functional annotation table
The AntiFam flag file was created with explicit tab separators:
```bash
printf "protein_id\tantifam_accession\tantifam_description\tinterpretation\n" \
> 06_combined_annotation/DL_diatom_antifam_flagged_proteins.tsv

printf "g10893.t1\tANF00012\ttRNA\tpotential_spurious_or_RNA_associated_ORF\n" \
>> 06_combined_annotation/DL_diatom_antifam_flagged_proteins.tsv

printf "g11404.t1\tANF00005\tAntisense to 23S rRNA\tpotential_spurious_or_RNA_associated_ORF\n" \
>> 06_combined_annotation/DL_diatom_antifam_flagged_proteins.tsv
```
These commands create a small AntiFam flag table so the flagged proteins can be retained but marked in downstream annotation files.
The final annotation layers were merged using:
```text
scripts/05_merge_functional_annotation_layers.py
```
Run the script with:
```bash
conda activate swissprot_annot
python scripts/05_merge_functional_annotation_layers.py
```
This script merges Swiss-Prot, Bacillariophyta, InterProScan, and AntiFam evidence into one master annotation table.
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
This command confirms that the master annotation table contains one row per predicted protein plus one header line.
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
This command verifies that the two AntiFam-flagged proteins were retained and correctly marked in the master table.
### 14.13 AntiFam-filtered interpretation tables
```bash
awk -F'\t' 'NR==1 || $9=="yes"' \
  06_combined_annotation/DL_diatom_master_functional_annotation.tsv \
  > 06_combined_annotation/DL_diatom_master_functional_annotation_AntiFam_flagged_only.tsv
```
This command creates a table containing only AntiFam-flagged proteins plus the header.
```bash
awk -F'\t' 'NR==1 || $9!="yes"' \
  06_combined_annotation/DL_diatom_master_functional_annotation.tsv \
  > 06_combined_annotation/DL_diatom_master_functional_annotation_no_AntiFam.tsv
```
This command creates an AntiFam-filtered version of the master annotation table for downstream interpretation.
### 14.14 Final functional annotation outputs
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
### 14.15 Expression integration status
Expression integration was completed after the master functional annotation table was generated. Transcriptome ORFs predicted by TransDecoder were aligned to the BRAKER4 ET protein set using DIAMOND BLASTP.

One best BRAKER4 hit was retained per TransDecoder ORF to prevent the same transcript-level `Average_TPM` value from being duplicated across multiple BRAKER4 proteins. The final expression merge used only the `Average_TPM` column from the trusted transcriptome table. No all-hit TPM sums, means, hit counts, or mapped-ORF lists were carried into the final table.
The detailed expression integration and final gene-table construction workflow is documented in Section 17.
### 14.16 Functional annotation status
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

---

<details>
<summary><strong>15. Nuclear-enriched genome generation</strong> - organelle contig removal and BRAKER4 annotation filtering</summary>

The nuclear-enriched genome was generated by removing contigs with strong chloroplast or mitochondrial similarity from the whole diatom assembly.
### 15.1 Input files
```bash
mkdir -p /work/ebg_lab/eb/diatom_consortia/nuclear_genome_filtering_18_diatom
cd /work/ebg_lab/eb/diatom_consortia/nuclear_genome_filtering_18_diatom

WHOLE=/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
CHLORO=/work/ebg_lab/eb/diatom_consortia/organelle/2_chloro/chloroplast_contig_1443_trimmed.fasta
MITO=/work/ebg_lab/eb/diatom_consortia/organelle/mito/diatom_candidate_mitochondrion_2contigs.fasta
```
These variables define the whole diatom assembly, the chloroplast sequence, and the mitochondrial sequence used for organelle-contig detection.
```bash
seqkit stats $WHOLE $CHLORO $MITO
```
This command summarizes the number of sequences and total length of the input assembly and organelle references.
Input assembly statistics:
```text
18_diatom.fasta                               3,010 contigs   82,172,226 bp
chloroplast_contig_1443_trimmed.fasta            1 contig        120,429 bp
diatom_candidate_mitochondrion_2contigs.fasta    2 contigs       104,526 bp
```
### 15.2 Combine organelle references
```bash
cat $CHLORO $MITO > organelles_chloro_mito.fasta
seqkit stats organelles_chloro_mito.fasta
```
These commands combine the chloroplast and mitochondrial FASTA files into one organelle reference file and check its total length.
The combined organelle reference contained three sequences with a total length of 224,955 bp.
### 15.3 Align the genome against organelle references
```bash
minimap2 -x asm5 -c \
    organelles_chloro_mito.fasta \
    $WHOLE \
    > whole_vs_organelles.paf
```
This command aligns the whole diatom assembly against the combined chloroplast and mitochondrial reference file. The PAF output records assembly regions with organelle similarity.
### 15.4 Calculate organelle-aligned coverage per contig
```bash
awk 'BEGIN{OFS="\t"} {print $1, $3, $4}' whole_vs_organelles.paf \
    > whole_vs_organelles.query_intervals.bed
```
This command extracts aligned query intervals from the PAF file.
```bash
sort -k1,1 -k2,2n whole_vs_organelles.query_intervals.bed \
    > whole_vs_organelles.query_intervals.sorted.bed

bedtools merge \
    -i whole_vs_organelles.query_intervals.sorted.bed \
    > whole_vs_organelles.query_intervals.merged.bed
```
These commands sort and merge overlapping aligned intervals so aligned bases are not counted multiple times.
```bash
seqkit fx2tab -n -l $WHOLE > whole_contig_lengths.tsv
```
This command writes contig names and contig lengths for the whole assembly.

```bash
awk 'BEGIN{OFS="\t"} {aligned[$1] += ($3 - $2)} END {for (c in aligned) print c, aligned[c]}' \
    whole_vs_organelles.query_intervals.merged.bed \
    > organelle_aligned_length_per_contig.tsv
```
This command calculates the total organelle-aligned length for each assembly contig.
```bash
awk 'BEGIN{OFS="\t"}
FNR==NR {len[$1]=$2; next}
{
    contig=$1;
    aligned=$2;
    pct=(aligned/len[contig])*100;
    print contig, len[contig], aligned, pct
}' whole_contig_lengths.tsv organelle_aligned_length_per_contig.tsv \
    > organelle_coverage_per_contig.tsv
```
This command calculates the percentage of each contig covered by organelle-like alignments.
```bash
sort -k4,4nr organelle_coverage_per_contig.tsv | head -n 50
```
This command lists the strongest organelle-like candidate contigs.
### 15.5 Remove organelle-like contigs
Contigs were classified as organelle-like if at least 70% of the contig length aligned to the chloroplast or mitochondrial reference.
```bash
awk '$4 >= 70 {print $1}' organelle_coverage_per_contig.tsv \
    > organelle_like_contigs.70pct.txt
```
This command selects contigs for removal if at least 70% of their length aligns to the organelle reference.

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
This command removes the organelle-like contigs from the whole diatom assembly and writes the nuclear-enriched genome FASTA.
```bash
seqkit stats \
    $WHOLE \
    organelles_chloro_mito.fasta \
    18_diatom_nuclear_enriched.v1.fasta
```
This command compares the original assembly, the combined organelle reference, and the final nuclear-enriched genome.

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
```
This command extracts the names of contigs retained in the nuclear-enriched genome.
```bash
MY_GFF=/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET/DL_diatom.braker4.ET.gff3
```
This variable defines the original BRAKER4 GFF3 annotation file.
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
}' nuclear_contigs.v1.txt $MY_GFF \
    > braker.18_diatom_nuclear_only.v1.gff3
```
This command keeps only BRAKER4 features located on retained nuclear contigs and removes the embedded FASTA section from the GFF3.
```bash
grep -c $'\tgene\t' braker.18_diatom_nuclear_only.v1.gff3
```
This command counts the number of gene features retained in the nuclear-filtered BRAKER4 annotation.
Final nuclear genome files:
```text
18_diatom_nuclear_enriched.v1.fasta
braker.18_diatom_nuclear_only.v1.gff3
```
</details>

---

<details>
<summary><strong>16. Pairwise genome comparison with <em>Phaeodactylum tricornutum</em></strong> - BLASTN and bedtools</summary>

A pairwise genome comparison was performed between the updated nuclear-enriched diatom genome and the reference genome of *Phaeodactylum tricornutum*. This analysis was used as a nucleotide-level similarity screen and was not intended to define complete gene orthology.
### 16.1 Input files and working directory
```bash
MY_GENOME=/work/ebg_lab/eb/diatom_consortia/nuclear_genome_filtering_18_diatom/18_diatom_nuclear_enriched.v1.fasta
MY_GFF=/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET/DL_diatom.braker4.ET.gff3

mkdir -p /work/ebg_lab/eb/diatom_consortia/phaeodactylum_blast_18_diatom_v1
cd /work/ebg_lab/eb/diatom_consortia/phaeodactylum_blast_18_diatom_v1
mkdir -p blast_out blast_db phaeodactylum
```
These variables define the nuclear-enriched diatom genome and the original BRAKER4 annotation file used for comparison.
### 16.2 Filter BRAKER4 GFF3 to nuclear contigs
```bash
seqkit seq -n $MY_GENOME > nuclear_contigs.v1.txt
```
This command extracts the contig names retained in the nuclear-enriched genome.
```bash
awk 'BEGIN{FS=OFS="\t"}
FNR==NR {keep[$1]=1; next}
$0 ~ /^#/ {
    if ($0 !~ /^##FASTA/) print;
    next
}
$1 in keep {print}' nuclear_contigs.v1.txt $MY_GFF \
    > braker.18_diatom_nuclear_only.v1.gff3
```
This command filters the BRAKER4 GFF3 file to keep only features on nuclear-enriched contigs.
```bash
grep -c $'\tgene\t' $MY_GFF
grep -c $'\tgene\t' braker.18_diatom_nuclear_only.v1.gff3
```
These commands count gene features before and after nuclear-contig filtering.
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
This command converts nuclear-filtered BRAKER4 gene coordinates from GFF3 to BED format for downstream overlap analysis.
```bash
wc -l blast_out/18_diatom_nuclear_v1_genes.bed
```
This command counts the number of nuclear gene models in the BED file.
The resulting BED file contained 15,048 nuclear gene models.
### 16.4 Download the *Phaeodactylum tricornutum* reference genome
```bash
datasets download genome accession GCF_000150955.2 \
    --include genome,gff3 \
    --filename phaeodactylum/Phaeodactylum_tricornutum_GCF_000150955.2.zip
```
This command downloads the *P. tricornutum* reference genome and annotation from NCBI.

```bash
unzip phaeodactylum/Phaeodactylum_tricornutum_GCF_000150955.2.zip -d phaeodactylum
```
This command extracts the downloaded reference genome package.
```bash
PT_GENOME=$(find phaeodactylum -name "*genomic.fna" | head -n 1)
PT_GFF=$(find phaeodactylum -name "*.gff" -o -name "*.gff3" | head -n 1)

echo $PT_GENOME
echo $PT_GFF
```
These commands locate the downloaded *Phaeodactylum* genome FASTA and GFF annotation files.
### 16.5 Build the BLAST database
```bash
makeblastdb \
    -in $PT_GENOME \
    -dbtype nucl \
    -parse_seqids \
    -out blast_db/Phaeodactylum_tricornutum

PT_DB=blast_db/Phaeodactylum_tricornutum
```
This command builds a nucleotide BLAST database from the *P. tricornutum* reference genome.
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
This command compares the nuclear-enriched diatom genome against the *P. tricornutum* reference genome at the nucleotide level.
```bash
wc -l blast_out/18_diatom_nuclear_v1_vs_phaeodactylum.blastn.tsv
```
This command counts the number of raw BLASTN alignments.
The raw BLASTN output contained:
```text
34,991 alignments
```
### 16.7 Filter BLASTN alignments
Alignments were filtered using:
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
This command keeps stronger BLASTN alignments and removes weak or very short matches before gene-overlap analysis.
```bash
wc -l blast_out/18_diatom_nuclear_v1_vs_phaeodactylum.filtered.tsv
```This command counts the number of BLASTN alignments that passed the filtering thresholds.

This produced:
```text
4,895 filtered BLASTN alignments
```
### 16.8 Identify diatom genes overlapping *Phaeodactylum*-like regions
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
This command converts filtered BLASTN hits into BED-like coordinates on the nuclear-enriched diatom genome.
```bash
bedtools intersect \
    -a blast_out/18_diatom_nuclear_v1_genes.bed \
    -b blast_out/blast_hits_on_18_diatom_nuclear_v1.bed \
    -wa -wb \
    > blast_out/18_diatom_nuclear_v1_genes_with_phaeodactylum_hits.tsv
```
This command identifies nuclear diatom genes that overlap filtered *Phaeodactylum*-like BLASTN regions.

This produced:
```text
2,269 gene-alignment overlaps
1,492 unique diatom genes with P. tricornutum-like BLASTN hits
```
```bash
cut -f4 blast_out/18_diatom_nuclear_v1_genes_with_phaeodactylum_hits.tsv \
    | sort -u \
    | wc -l
```
This command counts the number of unique diatom gene IDs overlapping filtered BLASTN hits.
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
These commands calculate the proportion of nuclear-filtered diatom genes with detectable nucleotide similarity to *P. tricornutum*.

Output:
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
This command sorts gene-overlap records by diatom gene ID and bitscore, then keeps the top-scoring BLASTN hit per gene.
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
This command calculates mean percent identity and mean alignment length across the best-hit gene set.
Output:
```text
Number of best-hit genes: 1492
Mean percent identity: 72.9402
Mean alignment length: 596.025
```
### 16.10 Create final interpreted BLAST table
```bash
awk 'BEGIN{FS=OFS="\t";
print "diatom_contig","diatom_gene_start","diatom_gene_end","diatom_gene_id","diatom_gene_strand","blast_hit_id","bitscore","blast_strand","pt_contig","pt_start","pt_end","pident","aln_len","evalue"
}
{
print $1,$2,$3,$4,$6,$10,$11,$12,$13,$14,$15,$16,$17,$18
}' blast_out/18_diatom_nuclear_v1_genes_best_phaeodactylum_hit.tsv \
    > blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.clean.tsv
```
This command creates a cleaner best-hit table containing the main diatom gene coordinates and BLASTN hit metrics.
### 16.11 Link BLAST regions to annotated *Phaeodactylum* genes
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
This command converts filtered BLASTN alignments into BED-like coordinates on the *P. tricornutum* reference genome.
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
This command converts *P. tricornutum* gene annotations from GFF3 to BED format and extracts gene IDs, names, and locus tags.
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
This command normalizes *Phaeodactylum* sequence IDs so BLAST hit coordinates and GFF-derived gene coordinates use matching names.
```bash
bedtools intersect \
    -a blast_out/blast_hits_on_phaeodactylum.normalized.bed \
    -b blast_out/phaeodactylum_genes.bed \
    -wa -wb \
    > blast_out/blast_hits_overlapping_phaeodactylum_genes.tsv
```
This command identifies *Phaeodactylum* genes overlapping the filtered BLASTN hit regions.

```bash
awk 'BEGIN{FS=OFS="\t";
print "blast_hit_id","pt_gene_id","pt_gene_name","pt_locus_tag"
}
{
  print $4,$16,$19,$20
}' blast_out/blast_hits_overlapping_phaeodactylum_genes.tsv \
    > blast_out/blast_hit_to_phaeodactylum_gene.tsv
```

This command creates a compact mapping table between each BLAST hit ID and overlapping *Phaeodactylum* gene identifiers.

The mapping table was merged with the clean best-hit table using Python:

```bash
python scripts/06_merge_phaeodactylum_blast_hits.py
```

This script adds PHATRDRAFT gene IDs, gene names, and locus tags to the clean diatom best-hit table.

```bash
awk 'BEGIN{FS=OFS="\t"}
NR==1 {print; next}
{
  gsub(/gene-/, "", $15);
  print
}' blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.with_PT_genes.tsv \
    > blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.with_PT_genes.cleanIDs.tsv
```

This command removes the `gene-` prefix from PHATRDRAFT identifiers.

```bash
{
  echo -e "diatom_contig\tdiatom_gene_start\tdiatom_gene_end\tdiatom_gene_id\tdiatom_gene_strand\tblast_hit_id\tbitscore\tblast_strand\tpt_contig\tpt_start\tpt_end\tpident\taln_len\tevalue\tpt_gene_id\tpt_gene_name\tpt_locus_tag"
  tail -n +2 blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.with_PT_genes.cleanIDs.tsv
} > blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.with_PT_genes.cleanIDs.fixed.tsv
```

This command rewrites the table header and creates the final cleaned BLASTN best-hit table.

### 16.12 Final comparative-genomics outputs

```text
blast_out/18_diatom_nuclear_v1_vs_phaeodactylum.blastn.tsv
blast_out/18_diatom_nuclear_v1_vs_phaeodactylum.filtered.tsv
blast_out/18_diatom_nuclear_v1_genes_with_phaeodactylum_hits.tsv
blast_out/18_diatom_nuclear_v1_genes_best_phaeodactylum_hit.tsv
blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.clean.tsv
blast_out/blast_hit_to_phaeodactylum_gene.tsv
blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.with_PT_genes.cleanIDs.fixed.tsv
```
### 16.13 Summary
Pairwise BLASTN comparison between the updated nuclear-enriched diatom genome and the *Phaeodactylum tricornutum* reference genome identified 4,895 filtered nucleotide alignments using thresholds of ≥70% identity, alignment length ≥200 bp, and e-value ≤1e-10.

These alignments overlapped 1,492 of 15,048 predicted nuclear genes, corresponding to 9.91% of the nuclear gene set. For the best hit per gene, the mean nucleotide identity was 72.94%, and the mean alignment length was 596 bp.

The best-hit regions were further mapped to annotated *P. tricornutum* PHATRDRAFT gene models. Of the 1,492 best-hit diatom genes, 1,446 overlapped annotated *P. tricornutum* gene models, representing 1,314 unique PHATRDRAFT genes.

This analysis provides a conservative nucleotide-level comparison between the candidate diatom genome and *P. tricornutum*. Because nucleotide similarity can be lost despite conservation at the protein level, this result should be treated as a genome-level similarity screen rather than a full orthology analysis.

</details>

---

<details>
<summary><strong>17. Clean BRAKER4 isoform-level gene table construction</strong> - TransDecoder, DIAMOND, Average_TPM, and Phaeodactylum yes/no</summary>

This section describes the final clean table construction used after rebuilding the annotation workflow from raw input files. The final output is a BRAKER4 isoform-level table with one row per predicted protein isoform.

The final table does not collapse BRAKER4 IDs, does not append GenBank-derived organelle rows, does not use all-hit TPM summaries, and does not retain detailed *Phaeodactylum tricornutum* hit columns.
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

```bash
mkdir -p /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/transdecoder_to_braker_ID_bridge/CLEAN_REBUILD_FROM_RAW
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/transdecoder_to_braker_ID_bridge/CLEAN_REBUILD_FROM_RAW

mkdir -p 01_input 02_diamond 03_best_hits 05_interproscan 06_combined_annotation 07_expression 08_phaeodactylum 09_final scripts logs slurm
```
This clean rebuild directory was used to avoid carrying forward columns from older all-hit and GenBank-merging workflows.

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
```bash
python scripts/02_make_swissprot_best_hits.py 2>&1 | tee swissprot_rebuild.log
```
This script rebuilds the Swiss-Prot best-hit and all-protein annotation tables from the raw DIAMOND output.
```bash
python scripts/03_make_bacillariophyta_best_hits.py 2>&1 | tee bacillariophyta_rebuild.log
```
This script rebuilds the Bacillariophyta best-hit and all-protein annotation tables from the raw DIAMOND output.
```bash
python scripts/04_summarize_interproscan.py 2>&1 | tee interproscan_rebuild.log
```
This script summarizes the raw InterProScan TSV output to one row per predicted protein.

Each all-protein annotation file was checked to confirm one row per BRAKER4 protein isoform plus one header line:
```bash
wc -l 03_best_hits/DL_diatom_all_proteins_with_swissprot_annotation.tsv
wc -l 03_best_hits/DL_diatom_all_proteins_with_bacillariophyta_annotation.tsv
wc -l 06_combined_annotation/DL_diatom_all_proteins_with_interproscan_summary.tsv
```
These commands verify that the annotation tables contain the expected 16,947 protein rows plus one header.
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
```bash
python scripts/05_merge_functional_annotation_layers.py 2>&1 | tee merge_functional_annotation_rebuild.log
```
This script merges Swiss-Prot, Bacillariophyta, InterProScan, and AntiFam evidence into the rebuilt master functional annotation table.
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
```bash
python scripts/07_add_BRAKER_lengths_clean.py 2>&1 | tee add_lengths_rebuild.log
```
This script adds contig ID, gene coordinates, strand, gene length, CDS length, and protein length from the BRAKER4 GFF3 file.
Output:
```text
07_expression/DL_diatom_master_functional_annotation_with_lengths.tsv
```
Observed counts:
```text
Input annotation rows:       16,947
GFF3 transcript rows parsed: 16,947
Rows with gene length:       16,947
Rows missing gene length:         0
```
### 17.5 Align TransDecoder ORFs to BRAKER4 proteins
```bash
set -euo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate swissprot_annot

cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/transdecoder_to_braker_ID_bridge/CLEAN_REBUILD_FROM_RAW

mkdir -p 07_expression logs
```
```bash
diamond makedb \
    --in 01_input/diatom_predicted_proteins.fa \
    -d 07_expression/braker4_ET.proteins
```
This command builds a DIAMOND database from the accepted BRAKER4 ET protein set.
```bash
diamond blastp \
    -q 01_input/transcriptome_orfs.transdecoder.clean.pep \
    -d 07_expression/braker4_ET.proteins \
    -o 07_expression/transcriptome_ORFs_vs_BRAKER4_ET_proteins.tsv \
    --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore \
    --max-target-seqs 10 \
    --evalue 1e-5 \
    --threads 16
```
This command aligns TransDecoder-predicted ORFs against the BRAKER4 protein set and retains up to 10 candidate hits per ORF.
Raw DIAMOND output:
```text
118,472 hit rows
39,935 unique TransDecoder ORFs with at least one BRAKER hit
14,473 unique BRAKER proteins hit by at least one reported ORF hit
```
### 17.6 Create one best ORF-to-BRAKER mapping per TransDecoder ORF
```bash
python scripts/08_make_best_ORF_to_BRAKER_mapping_clean.py 2>&1 | tee make_best_mapping_rebuild.log
```
This script calculates ORF and BRAKER protein coverage from the raw DIAMOND output and keeps one best BRAKER4 hit per TransDecoder ORF.
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
This command copies the trusted transcriptome expression table into the clean rebuild directory.

Only two columns were used from this file:
```text
orf
Average_TPM
```
```bash
python scripts/09_add_ONLY_Average_TPM_clean.py 2>&1 | tee add_average_tpm_rebuild.log
```
This script adds `transdecoder_orf_id` and `Average_TPM` to BRAKER4 isoforms using the best ORF-to-BRAKER mapping.

If multiple TransDecoder ORFs mapped best to the same BRAKER4 protein, the ORF with the highest valid `Average_TPM` was retained. ORFs without valid numeric `Average_TPM` were excluded from the expression-supported ORF column.
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
### 17.8 Add *Phaeodactylum tricornutum* yes/no only
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
Therefore, the script used the root only for lookup while preserving the full BRAKER4 isoform ID in the final table. The only *Phaeodactylum*-derived output column is:
```text
in_Phaeodactylum_tricornutum
```
Detailed columns such as `blast_hit_id`, `pt_gene_id`, `pt_gene_name`, `bitscore`, `pident`, and `evalue` were not retained.
### 17.9 Create final clean BRAKER4 isoform-level table
```bash
python scripts/10_make_FINAL_clean_BRAKER_isoform_table.py 2>&1 | tee make_final_clean_table.log
```
This script creates the final clean BRAKER4 isoform-level table by combining functional annotation, BRAKER4 lengths, expression, compartment labels, and *Phaeodactylum* yes/no status.
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
`gene_id` is the original BRAKER4 protein isoform ID and is not collapsed to the gene root. `transdecoder_orf_id` is the TransDecoder ORF that contributed the expression value for that BRAKER4 isoform.

The `transdecoder_orf_id` column is populated only when the mapped ORF has a valid numeric `Average_TPM`. If a TransDecoder ORF maps to a BRAKER4 protein but does not have a usable `Average_TPM`, both `transdecoder_orf_id` and `Average_TPM` are left blank in the final output.

If multiple TransDecoder ORFs map best to the same BRAKER4 isoform, the ORF with the highest valid `Average_TPM` is retained. TPM values are not summed, averaged, or duplicated across all DIAMOND hits.
The final table therefore reports:
```text
one BRAKER4 isoform row
one representative TransDecoder ORF ID when expression support is available
one Average_TPM value
one Phaeodactylum tricornutum yes/no value
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

</details>

---

<details>
<summary><strong>18. Hi-C read mapping and contig-level proximity-ligation network</strong> - BWA-MEM, samtools, awk, YaHS, and Python</summary>

Hi-C paired-end reads were incorporated after the main assembly, annotation, expression, and comparative-genomics workflow. The goal was to assess how broadly the polished whole assembly was represented in the proximity-ligation dataset and to identify contigs connected by Hi-C read pairs.

The Hi-C analysis was performed on the polished whole assembly rather than the nuclear-enriched subset because the initial goal was to evaluate contig-level representation and proximity-ligation links across the whole assembly.
### 18.1 Input files and working directories
```bash
cd /work/ebg_lab/eb/diatom_consortia

mkdir -p hi-c_diatoms/01_qc
mkdir -p hi-c_diatoms/02_map_to_whole_assembly
mkdir -p hi-c_diatoms/03_yahs_scaffolding
mkdir -p hi-c_diatoms/04_contact_maps
```
Input Hi-C reads:
```text
/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/1574499_S6_L001_R1_001.fastq.gz
/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/1574499_S6_L001_R2_001.fastq.gz
```

Input polished whole assembly:

```text
/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_output/pypolca_corrected.fasta
```
### 18.2 Map Hi-C reads to the polished whole assembly
The original assembly directory was not writable by the Hi-C job, so the assembly was linked into the Hi-C working directory and indexed there.
```bash
cd /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/02_map_to_whole_assembly

ASM_ORIG=/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_output/pypolca_corrected.fasta
ASM_LOCAL=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/02_map_to_whole_assembly/pypolca_corrected.hic_input.fasta

ln -sfn $ASM_ORIG $ASM_LOCAL
```
This command links the polished whole assembly into the Hi-C mapping directory.
```bash
bwa index $ASM_LOCAL
samtools faidx $ASM_LOCAL
cut -f1,2 ${ASM_LOCAL}.fai > pypolca_corrected.chrom.sizes
```
These commands index the assembly for BWA and samtools, then save contig lengths for downstream contact-table annotation.
```bash
R1=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/1574499_S6_L001_R1_001.fastq.gz
R2=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/1574499_S6_L001_R2_001.fastq.gz

bwa mem -SP5M -t 24 $ASM_LOCAL $R1 $R2 \
    | samtools view -@ 8 -bS - \
    > hic_to_whole_assembly.raw.bam
```
This command maps paired Hi-C reads to the polished whole assembly using BWA-MEM with Hi-C-compatible settings.
```bash
samtools sort -@ 24 -n \
    -o hic_to_whole_assembly.name_sorted.bam \
    hic_to_whole_assembly.raw.bam
```
This command creates a name-sorted BAM file for paired-read contact extraction.
```bash
samtools sort -@ 24 \
    -o hic_to_whole_assembly.coord_sorted.bam \
    hic_to_whole_assembly.raw.bam

samtools index hic_to_whole_assembly.coord_sorted.bam
```
These commands create and index a coordinate-sorted BAM file for mapping summaries.
```bash
samtools flagstat hic_to_whole_assembly.coord_sorted.bam \
    > hic_to_whole_assembly.flagstat.txt
```
This command summarizes read-level mapping statistics.
```bash
samtools idxstats hic_to_whole_assembly.coord_sorted.bam \
    > hic_to_whole_assembly.idxstats.txt
```
This command reports mapped and unmapped Hi-C read counts per contig.
Read-level mapping summary:
```text
Primary reads:              890,810
Primary mapped reads:       622,967
Primary mapping rate:       69.93%
All mapped alignments:      803,799 / 1,071,642 = 75.01%
Read pairs:                 445,405
Singletons:                 54,019 reads = 6.06%
```
### 18.3 Summarize contig-level Hi-C representation
```bash
awk 'BEGIN {
    OFS="\t";
    print "contig","length_bp","mapped_HiC_reads","unmapped_HiC_reads","HiC_mapped"
}
$1!="*" {
    status = ($3 > 0 ? "yes" : "no");
    print $1,$2,$3,$4,status
}' hic_to_whole_assembly.idxstats.txt \
> hic_contig_mapping_presence.tsv
```
This command classifies each assembly contig as represented in the Hi-C dataset if at least one Hi-C read mapped to it.
```bash
{
    echo "Hi-C contig mapping summary"
    echo "Date: $(date)"
    echo
    echo "Input assembly:"
    echo "/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_output/pypolca_corrected.fasta"
    echo
    echo "Hi-C reads:"
    echo "/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/1574499_S6_L001_R1_001.fastq.gz"
    echo "/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/1574499_S6_L001_R2_001.fastq.gz"
    echo
    awk '
    $1!="*" {
        total++;
        if ($3 > 0) mapped++;
    }
    END {
        unmapped = total - mapped;
        printf "Total contigs: %d\n", total;
        printf "Contigs with >=1 Hi-C read mapped: %d\n", mapped;
        printf "Contigs with 0 Hi-C reads mapped: %d\n", unmapped;
        printf "Percent contigs with Hi-C reads mapped: %.2f%%\n", (mapped/total)*100;
    }
    ' hic_to_whole_assembly.idxstats.txt
    echo
    echo "Read-level mapping summary:"
    cat hic_to_whole_assembly.flagstat.txt
} > hic_contig_mapping_summary.txt
```
This command creates one readable summary file combining contig-level Hi-C representation and read-level mapping statistics.

Contig-level Hi-C representation:
```text
Total contigs: 4,925
Contigs with >=1 Hi-C read mapped: 4,010
Contigs with 0 Hi-C reads mapped: 915
Percent contigs with Hi-C reads mapped: 81.42%
```
### 18.4 Exploratory whole-assembly Hi-C scaffolding with YaHS
Whole-assembly Hi-C scaffolding was tested with YaHS as an exploratory step. Because the assembly represents a consortium, this result was treated cautiously and was not used as the final Hi-C integration output.
```bash
cd /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/03_yahs_scaffolding

ASM=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/02_map_to_whole_assembly/pypolca_corrected.hic_input.fasta
BAM=/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/02_map_to_whole_assembly/hic_to_whole_assembly.name_sorted.bam

yahs -o DL_diatom_whole_hic_yahs $ASM $BAM
```
This command tests whether Hi-C read pairs could scaffold the polished whole assembly.
```bash
seqkit stats \
    $ASM \
    DL_diatom_whole_hic_yahs_scaffolds_final.fa \
    > DL_diatom_whole_hic_yahs.seqkit_stats.txt
```
This command compares the polished input assembly with the YaHS scaffolded output.
YaHS scaffold summary:
```text
Input assembly:
4,925 contigs
189,915,395 bp
maximum contig length: 5,424,378 bp

YaHS output:
5,032 scaffolds
189,915,395 bp
maximum scaffold length: 5,424,378 bp
```
The YaHS run did not increase maximum scaffold length and increased the number of sequences. Therefore, the whole-assembly YaHS output was treated as exploratory rather than as a final scaffolded assembly.
### 18.5 Extract all-primary inter-contig Hi-C contacts
The final contig-contact analysis used primary mapped Hi-C read pairs without applying a MAPQ cutoff. Unmapped reads, mate-unmapped reads, secondary alignments, and supplementary alignments were excluded. Each read pair was counted once if the two mates mapped to different contigs.
```bash
cd /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/02_map_to_whole_assembly
BAM=hic_to_whole_assembly.name_sorted.bam
samtools view -@ 8 -f 1 -F 2316 $BAM \
    | awk 'BEGIN{OFS="\t"}
    {
        q=$1
        r=$3

        if (prev != "" && q != prev) {
            if (n == 2 && contig[1] != contig[2]) {
                a=contig[1]
                b=contig[2]

                if (a > b) {
                    tmp=a
                    a=b
                    b=tmp
                }

                print a,b
            }

            delete contig
            n=0
        }

        prev=q
        n++

        if (n <= 2) {
            contig[n]=r
        }
    }
    END {
        if (n == 2 && contig[1] != contig[2]) {
            a=contig[1]
            b=contig[2]

            if (a > b) {
                tmp=a
                a=b
                b=tmp
            }

            print a,b
        }
    }' \
    | sort -S 20G \
    | uniq -c \
    | awk 'BEGIN{OFS="\t"; print "contig_A","contig_B","HiC_contact_pairs"}
           {print $2,$3,$1}' \
    | sort -k3,3nr \
    > hic_intercontig_contacts_all_primary_pairs.tsv
```
This command extracts inter-contig Hi-C contacts from primary mapped read pairs and counts the number of read pairs supporting each contig-to-contig link.
```bash
awk 'BEGIN{OFS="\t"}
NR==FNR {
    len[$1]=$2
    next
}
FNR==1 {
    print $0,"len_A","len_B"
    next
}
{
    print $0,len[$1],len[$2]
}' pypolca_corrected.chrom.sizes \
   hic_intercontig_contacts_all_primary_pairs.tsv \
   > hic_intercontig_contacts_all_primary_pairs.with_lengths.tsv
```
This command adds contig lengths to the full inter-contig Hi-C contact table.
The main contact table reports:
```text
contig_A
contig_B
Hi-C contact pairs
length of contig_A
length of contig_B
```
Each row represents one pair of contigs connected by Hi-C proximity-ligation evidence.
### 18.6 Summarize connected contigs

```bash
awk 'NR>1 {print $1; print $2}' \
    hic_intercontig_contacts_all_primary_pairs.tsv \
    | sort -u \
    > hic_connected_contigs_all_primary_pairs.txt
```

This command lists every contig involved in at least one inter-contig Hi-C contact.

```bash
{
    echo -e "contig\tnumber_of_connected_contigs\ttotal_intercontig_HiC_pairs\tstrongest_single_contact_pairs"

    awk 'BEGIN{OFS="\t"}
    NR==1 {next}
    {
        a=$1
        b=$2
        pairs=$3

        if (!(a SUBSEP b in seen)) {
            seen[a SUBSEP b]=1
            degree[a]++
            degree[b]++
        }

        total_pairs[a]+=pairs
        total_pairs[b]+=pairs

        if (pairs > max_pair[a]) max_pair[a]=pairs
        if (pairs > max_pair[b]) max_pair[b]=pairs
    }
    END {
        for (c in total_pairs) {
            print c,degree[c],total_pairs[c],max_pair[c]
        }
    }' hic_intercontig_contacts_all_primary_pairs.tsv \
    | sort -k3,3nr
} > hic_contig_connectivity_summary_all_primary_pairs.tsv
```
This command summarizes each connected contig by number of partner contigs, total inter-contig Hi-C pairs, and strongest single contig-to-contig contact.
Final all-primary contact-network summary:
```text
Connected contigs: 3,770
Inter-contig Hi-C links: 75,703
```

### 18.7 Convert the contact table to network files
The full contig-contact table was converted into GEXF and GraphML network files using a small helper Python script.
The script is saved as:
```text
scripts/11_make_hic_network_files.py
```
Run the script with:
```bash
cd /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/02_map_to_whole_assembly
python scripts/11_make_hic_network_files.py
```
This script reads the Hi-C inter-contig contact table and exports the contact network in GEXF and GraphML formats.
The script generated:
```text
hic_contig_network_all_primary_pairs.gexf
hic_contig_network_all_primary_pairs.graphml
```
The final network contained:
```text
Nodes: 3,770 contigs
Edges: 75,703 inter-contig Hi-C links
```
In this network:
```text
Node = assembly contig
Edge = Hi-C proximity-ligation contact between two contigs
Edge weight = number of Hi-C read pairs supporting the contig-to-contig connection
```
### 18.8 Organize final Hi-C outputs
After generating the final mapping and contact-network outputs, files were organized into final mapping and contact-map folders.
```bash
cd /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms
mkdir -p 02_map_to_whole_assembly/final_mapping
mkdir -p 04_contact_maps/tables
mkdir -p 04_contact_maps/network_files
mkdir -p 04_contact_maps/scripts
```
Final mapping files:
```text
02_map_to_whole_assembly/final_mapping/
├── hic_to_whole_assembly.flagstat.txt
├── hic_to_whole_assembly.idxstats.txt
├── hic_contig_mapping_summary.txt
├── hic_contig_mapping_presence.tsv
├── hic_to_whole_assembly.name_sorted.bam
├── hic_to_whole_assembly.coord_sorted.bam
├── hic_to_whole_assembly.coord_sorted.bam.bai
└── pypolca_corrected.chrom.sizes
```
Final contact-network files:
```text
04_contact_maps/
├── tables/
│   ├── hic_intercontig_contacts_all_primary_pairs.tsv
│   ├── hic_intercontig_contacts_all_primary_pairs.with_lengths.tsv
│   ├── hic_connected_contigs_all_primary_pairs.txt
│   └── hic_contig_connectivity_summary_all_primary_pairs.tsv
├── network_files/
│   ├── hic_contig_network_all_primary_pairs.gexf
│   └── hic_contig_network_all_primary_pairs.graphml
└── scripts/
    └── 11_make_hic_network_files.py
```
### 18.9 Final Hi-C analysis summary
Hi-C reads mapped to 4,010 of 4,925 contigs in the polished whole assembly, corresponding to 81.42% of assembly contigs. At the read level, 622,967 of 890,810 primary reads mapped to the assembly, corresponding to a primary mapping rate of 69.93%.

Inter-contig proximity-ligation contacts were extracted from primary mapped Hi-C read pairs without applying a MAPQ cutoff. The final all-primary contig-contact network contained 3,770 contig nodes and 75,703 inter-contig Hi-C links.

</details>
