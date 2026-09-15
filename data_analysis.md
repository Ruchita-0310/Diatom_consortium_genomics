# Diatom Consortia: Metagenomic and Metatranscriptomic Pipeline
This repository documents the workflow used to assemble, polish, bin, classify, annotate, and compare genomes and transcriptomes from a diatom associated microbial consortium. The current workflow combines long read metagenomic assembly, short read polishing, metagenomic binning, contig level taxonomic screening, organelle identification, marker based phylogenetic analyses, transcriptome analysis, BRAKER4 ET gene prediction, nuclear enriched genome generation, functional annotation, expression integration, repeat aware gene curation, four genome nucleotide comparison, secondary protein orthology analysis, ORF level comparative metatranscriptomics, and Hi C contact network analysis.

The current manuscript level comparative gene analysis uses a repeat and plastid quality controlled set of 14,941 Deer Lake nuclear genes and compares them by `dc-megablast` with *Nitzschia inconspicua*, *Seminavis robusta*, *Phaeodactylum tricornutum*, and *Thalassiosira pseudonana*. OrthoFinder is retained as a secondary protein level analysis. Historical PT and TP pairwise BLASTN sections remain below for provenance, but the final four comparator BLASTN workflow is documented in Section 21.

A separate ORF level comparative metatranscriptomic analysis was subsequently performed from the Deer Lake **nf-core/metatdenovo** output. This analysis uses expressed TransDecoder ORFs, reciprocal best BLASTP hits against the same four reference diatoms, a within Deer Lake expression percentile, and integrated KOfam, EggNOG, and direct Pfam/HMMER annotation. The new workflow is documented in Section 23 and is complementary to, rather than a replacement for, the 14,941 nuclear gene BLASTN analysis.

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
Organelle identification
   ↓
18S rRNA, plastid 16S rRNA, and rbcL phylogenies
   ↓
BRAKER4 ET annotation using RNA-seq evidence
   ↓
Nuclear enriched genome generation
   ↓
Functional annotation + Average_TPM integration
   ↓
Representative nuclear gene curation
   ↓
RepeatModeler + RepeatMasker
   ↓
Removal of 38 high confidence TE associated models
   ↓
Plastid derived contig QC
   ↓
Final 14,941 gene nucleotide query
   ↓
Four genome dc-megablast: NI, SR, PT, TP
   ↓
Master BLASTN + Average_TPM table
   ↓
BLASTN unmatched and expression subsets
   ↓
Secondary OrthoFinder protein comparison
   ↓
Expressed Deer Lake TransDecoder ORF BLASTP RBH comparison
   ↓
Deer Lake expression percentile + KOfam/EggNOG/Pfam integration
   ↓
Manual curation of six functional systems
   ↓
Hi C mapping to polished whole assembly
   ↓
High confidence contig contact processing
   ↓
Whole assembly network visualization
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
| Organelle identification         | MetaQUAST, minimap2, bedtools, seqkit                                                                                                   |
| Phylogenetics                    | Barrnap, BLAST+, bedtools, seqkit, Clustal Omega, TrimAl, IQ-TREE 2, Python                                                            |
| Transcriptomics                  | Nextflow, nf-core/metatdenovo, TransDecoder, Barrnap, STAR                                                                              |
| Comparative metatranscriptomics    | BLAST+, KOfam, EggNOG mapper, HMMER/Pfam, Python, pandas, Matplotlib                                                                    |
| Genome annotation                | BRAKER4, GeneMark-ET, AUGUSTUS, TSEBRA, STAR, BUSCO/compleasm                                                                           |
| Functional annotation            | DIAMOND, UniProtKB/Swiss-Prot, UniProtKB Bacillariophyta, InterProScan, Pfam, PANTHER, Gene3D, CDD, SMART, SUPERFAMILY, ProSite, Python |
| Expression integration           | DIAMOND, Python, pandas, TransDecoder ORFs, Average_TPM table                                                                           |
| Comparative genomics             | NCBI RefSeq/FTP, BLASTN, bedtools, seqkit, RepeatModeler, RepeatMasker, OrthoFinder, DIAMOND, FAMSA, FastTree, Python                 |
| Hi-C mapping and contact network | FastQC, MultiQC, BWA-MEM, samtools, seqkit, YaHS, awk, Python                                                                           |

---

## Repository structure for Python scripts
Custom Python scripts are stored in `scripts/` and are numbered according to their role in the analysis. Runnable shell and SLURM workflows are documented directly within the relevant workflow sections rather than listed as separate repository scripts.

```text
scripts/
├── 01_classify_metaeuk_contigs.py
├── 02_make_swissprot_best_hits.py
├── 03_make_bacillariophyta_best_hits.py
├── 04_summarize_interproscan.py
├── 05_merge_functional_annotation_layers.py
├── 06_merge_phaeodactylum_blast_hits.py
├── 07_merge_thalassiosira_blast_hits.py
├── 08_add_BRAKER_lengths_clean.py
├── 09_make_best_ORF_to_BRAKER_mapping_clean.py
├── 10_add_ONLY_Average_TPM_clean.py
├── 11_make_FINAL_clean_BRAKER_isoform_table.py
├── 12_make_boss_review_gene_table_PTredo.py
├── 13_add_thalassiosira_yes_no.py
├── 14_make_hic_network_files.py
├── 15_make_hic_primary_mapq30_pid95_tables.py
├── 16_make_hic_pair_type_tables.py
├── 17_make_hic_simple_mixed_read_table.py
├── 18_make_DL_nuclear_representative_proteome.py
├── 19_calculate_DL_CDS_repeat_overlap.py
├── 20_classify_DL_TE_candidates.py
├── 21_prepare_reference_proteome.py
├── 22_make_final_orthofinder_gene_table.py
├── 23_prepare_DL_BLASTN_query.py
├── 24_build_BLASTN_master_table.py
├── 25_make_BLASTN_subsets.py
├── 26_plot_BLASTN_TPM_panelA.py
├── 27_make_HiC_undirected_min2_edges.py
├── 28_plot_HiC_whole_assembly_network.py
├── 29_prepare_DL_expressed_ORFs.py
├── 30_extract_unique_top_RBH.py
├── 31_build_5species_RBH_expression_master.py
├── 32_add_DL_expression_percentile.py
├── 33_integrate_transcriptome_annotations.py
├── 34_extract_six_system_candidates.py
├── 35_make_curated_6systems_expression_conservation.py
├── 36_make_final_24genes_expression_conservation.py
├── 37_make_final_30genes_expression_conservation.py
└── 38_plot_DL_expression_RBH_matrix.py
```

The scripts are numbered sequentially from `01` to `38`. Custom Python logic is stored in `scripts/`. The comparative metatranscriptomics BLAST jobs added in Section 23 are stored as plain text SLURM files in `slurm/`.

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
  Collapses Phaeodactylum and Deer Lake gene overlaps to one row per raw BLASTN hit.

07_merge_thalassiosira_blast_hits.py
  Collapses Thalassiosira and Deer Lake gene overlaps to one row per raw BLASTN hit.

08_add_BRAKER_lengths_clean.py
  Adds BRAKER4 coordinates, contig IDs, strand, gene length, CDS length, and protein length to each isoform.

09_make_best_ORF_to_BRAKER_mapping_clean.py
  Parses TransDecoder ORF versus BRAKER4 DIAMOND output and keeps one best BRAKER4 hit per TransDecoder ORF.

10_add_ONLY_Average_TPM_clean.py
  Adds the selected TransDecoder ORF ID and Average_TPM value to BRAKER4 isoforms.

11_make_FINAL_clean_BRAKER_isoform_table.py
  Creates the final clean BRAKER4 isoform-level table with one row per predicted protein isoform.

12_make_boss_review_gene_table_PTredo.py
  Creates the simplified review table containing the core annotation, expression, compartment, and Phaeodactylum comparison fields.

13_add_thalassiosira_yes_no.py
  Preserves the Phaeodactylum comparison field and adds present_in_Thalassiosira_pseudonana to the final isoform and review tables.

14_make_hic_network_files.py
  Converts the Hi-C contig-contact table into GEXF and GraphML network files.

15_make_hic_primary_mapq30_pid95_tables.py
  Parses separate Hi-C R1 and R2 BAM files, removes non-primary alignments, and keeps MAPQ >= 30 and percent identity >= 95 alignments.

16_make_hic_pair_type_tables.py
  Joins read 1 and read 2 by read ID, assigns contig types using the diatom draft genome, and creates high-confidence Hi-C pair-type tables.

17_make_hic_simple_mixed_read_table.py
  Creates the final simplified read-level table for mixed diatom-bacterial Hi-C pairs.

18_make_DL_nuclear_representative_proteome.py
  Selects one representative nuclear BRAKER4 protein per Deer Lake gene, excludes specified AntiFam gene roots, and removes only unsupported proteins shorter than 50 aa.

19_calculate_DL_CDS_repeat_overlap.py
  Intersects final representative CDS coordinates with RepeatMasker annotations and summarizes repeat overlap without double counting overlapping repeat hits.

20_classify_DL_TE_candidates.py
  Applies a conservative annotation-supported review scheme to repeat-overlapping Deer Lake protein models and writes the Tier 1 TE removal list.

21_prepare_reference_proteome.py
  Creates one representative protein per NCBI locus tag and optionally excludes organelle contigs; used for PT, TP, SR, and NI preparation.

22_make_final_orthofinder_gene_table.py
  Builds the historical OrthoFinder gene table with functional annotation, Average_TPM, and four orthogroup-sharing fields.

23_prepare_DL_BLASTN_query.py
  Extracts the final nucleotide gene query using gene roots represented in the repeat and plastid quality controlled Deer Lake protein set.

24_build_BLASTN_master_table.py
  Collapses four dc-megablast result files to one best alignment per Deer Lake gene per comparator and merges annotation, compartment, and Average_TPM metadata.

25_make_BLASTN_subsets.py
  Writes no-comparator-hit, expressed no-hit, and expressed annotated no-hit subsets without applying a gene-length cutoff.

26_plot_BLASTN_TPM_panelA.py
  Generates the color blind aware expression violin plot used for Panel A.

27_make_HiC_undirected_min2_edges.py
  Removes self contacts, retains oriented rows with at least two Hi C read pairs, and collapses reciprocal contig pairs to unique undirected edges.

28_plot_HiC_whole_assembly_network.py
  Plots all assembly contigs, the high confidence connected network, isolated contigs, and highlighted organelle-associated contigs.

29_prepare_DL_expressed_ORFs.py
  Matches nf-core/metatdenovo TransDecoder ORFs to the ORF level TPM table, calculates Average_TPM, and writes the expressed Deer Lake peptide set.

30_extract_unique_top_RBH.py
  Identifies unique top BLASTP hits in both directions and writes reciprocal best hit pairs; top bitscore ties are treated as ambiguous.

31_build_5species_RBH_expression_master.py
  Merges Deer Lake Average_TPM, four RBH tables, reference protein metadata, RBH pattern, and RBH count into one ORF level master table.

32_add_DL_expression_percentile.py
  Calculates within Deer Lake Average_TPM percentiles across the 87,621 quantified ORFs.

33_integrate_transcriptome_annotations.py
  Joins KOfam, EggNOG, and direct Pfam/HMMER evidence to the expression and RBH master table.

34_extract_six_system_candidates.py
  Performs the intentionally broad first pass screen for six functional systems prior to manual review.

35_make_curated_6systems_expression_conservation.py
  Writes the manually reviewed 36 gene table with six representative genes per functional system.

36_make_final_24genes_expression_conservation.py
  Writes the first compact 24 gene figure table and marks confirmed plastid encoded genes as not assessed for nuclear proteome RBH.

37_make_final_30genes_expression_conservation.py
  Prepared next step that expands the compact table to five genes per functional system; not yet executed at the documented stopping point.

38_plot_DL_expression_RBH_matrix.py
  Prepared plotting script for the compact expression percentile plus NI/SR/PT/TP RBH matrix; intended for the 30 gene table after final organelle checks.
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
<summary><strong>9. Phylogenetic analyses</strong> - 18S rRNA, plastid 16S rRNA, and rbcL</summary>

Three marker-based phylogenetic analyses were used to evaluate the placement of the Deer Lake diatom. The 18S rRNA sequence was recovered from the metatranscriptomic assembly, whereas plastid 16S rRNA and rbcL were obtained from the reconstructed chloroplast contig. All three trees were aligned with Clustal Omega, trimmed with TrimAl, and inferred with IQ-TREE 2 using ModelFinder, 1,000 ultrafast bootstrap replicates, and 1,000 SH-aLRT replicates.

### 9.1 18S rRNA phylogeny

The 18S rRNA phylogeny was generated from the Deer Lake 18S sequence and selected diatom reference sequences.

```bash
cat *.fasta > 18S_new.fasta
```

This command combines individual 18S FASTA files into one input file for alignment.

```bash
clustalo \
    -i 18S_new.fasta \
    -o 18S_aligned.fasta
```

The alignment was trimmed with TrimAl:

```bash
trimal \
    -in 18S_aligned.fasta \
    -out 18S_trimmed.fasta \
    -automated1
```

The maximum-likelihood tree was inferred with IQ-TREE 2:

```bash
/home/ruchita.solanki/iqtree-2.2.2.7-Linux/bin/iqtree2 \
    -s 18S_trimmed.fasta \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -nt AUTO
```

ModelFinder selected the best-fit nucleotide substitution model, and branch support was estimated using ultrafast bootstrap and SH-aLRT values.

---

### 9.2 Plastid 16S rRNA phylogeny

The plastid 16S rRNA sequence was identified from the reconstructed Deer Lake chloroplast contig:

```text
chloroplast_contig_1443_trimmed.fasta
```

The chloroplast sequence was 120,429 bp. Barrnap identified two complete 16S rRNA copies on the plastid contig:

```text
copy 1: contig_1443  46638-48119    strand -
copy 2: contig_1443 114897-116378   strand +
```

Both copies were 1,482 bp and were identical. The two loci were extracted with bedtools:

```bash
printf "contig_1443\t46637\t48119\tchloroplast_16S_copy1\t.\t-\ncontig_1443\t114896\t116378\tchloroplast_16S_copy2\t.\t+\n" \
    > chloroplast_16S.bed

bedtools getfasta \
    -fi chloroplast_contig_1443_trimmed.fasta \
    -bed chloroplast_16S.bed \
    -s \
    -name \
    > DeerLake_chloroplast_16S_copies.fasta
```

Because the two plastid copies were identical, one copy was retained for phylogenetic analysis:

```bash
seqkit grep \
    -r \
    -p 'chloroplast_16S_copy2' \
    DeerLake_chloroplast_16S_copies.fasta \
    | sed 's/^>.*/>Deer_Lake_diatom_16S/' \
    > DeerLake_diatom_16S.fasta
```

Reference plastid 16S sequences were selected from NCBI searches. Records labelled as uncultured were excluded, and duplicate inverted-repeat copies from the same plastid accession were collapsed to one representative sequence. The curated reference dataset contained 27 sequences. Addition of the Deer Lake plastid 16S sequence produced a 28-sequence dataset.

```bash
cat DeerLake_diatom_16S.fasta \
    16S_references_large.fasta \
    > all_16S_large.fasta
```

The sequences were aligned with Clustal Omega:

```bash
clustalo \
    -i all_16S_large.fasta \
    -o all_16S_large_aligned.fasta \
    --threads=32 \
    --force
```

The alignment was trimmed with TrimAl:

```bash
trimal \
    -in all_16S_large_aligned.fasta \
    -out all_16S_large_aligned_trimmed.fasta \
    -automated1
```

The final trimmed plastid 16S alignment contained:

```text
Sequences:          28
Alignment length:   1,475 bp
```

The maximum-likelihood tree was inferred with IQ-TREE 2:

```bash
/home/ruchita.solanki/iqtree-2.2.2.7-Linux/bin/iqtree2 \
    -s all_16S_large_aligned_trimmed.fasta \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -nt 32 \
    -pre DeerLake_16S_large \
    -redo
```

The main output tree is:

```text
DeerLake_16S_large.treefile
```

---

### 9.3 Plastid rbcL phylogeny

The rbcL locus from `contig_1443` was used to identify related diatom rbcL sequences in NCBI. Matching rbcL sequences were downloaded as a FASTA file for phylogenetic analysis.

The initial NCBI download contained:

```text
100 rbcL reference sequences
```

All sequences assigned to *Durinskia* were removed before tree construction. Exact nucleotide duplicates were then screened so that only literal duplicate sequences would be removed; different strains or accessions from the same species were retained when their rbcL sequences differed.

The curation step produced:

```text
Starting references:              100
After removing Durinskia:          87
After exact-sequence filtering:    87
```

The following Python workflow was used:

```python
INPUT = "seqdump.txt"
OUTPUT = "rbcL_no_Durinskia_deduplicated.fasta"

def read_fasta(path):
    records = []
    header = None
    seq = []

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))

                header = line[1:]
                seq = []

            else:
                seq.append(line)

    if header is not None:
        records.append((header, "".join(seq)))

    return records


records = read_fasta(INPUT)

no_durinskia = [
    (header, seq)
    for header, seq in records
    if "durinskia" not in header.lower()
]

seen_sequences = set()
clean = []

for header, seq in no_durinskia:
    seq_upper = seq.upper()

    if seq_upper in seen_sequences:
        continue

    seen_sequences.add(seq_upper)
    clean.append((header, seq))

with open(OUTPUT, "w") as out:
    for header, seq in clean:
        out.write(">" + header + "\n")

        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")

print("Starting sequences:", len(records))
print("After removing Durinskia:", len(no_durinskia))
print("After exact-sequence deduplication:", len(clean))
```

The curated reference dataset had the following length distribution:

```text
Reference sequences:  87
Minimum length:        1,386 bp
Average length:        1,446.9 bp
Maximum length:        1,500 bp
```

The Deer Lake rbcL sequence was added to the 87-reference dataset:

```bash
cat DeerLake_diatom_rbcL.fasta \
    rbcL_no_Durinskia_deduplicated.fasta \
    > all_rbcL_for_tree.fasta
```

This produced a dataset of 88 sequences.

#### rbcL sequence orientation

Some rbcL sequences downloaded from complete plastid genomes were stored in the opposite orientation. Before alignment, sequence orientation was checked by comparing the Deer Lake rbcL sequence against the complete rbcL dataset with local BLASTN.

```bash
makeblastdb \
    -in all_rbcL_for_tree.fasta \
    -dbtype nucl \
    -out rbcL_orientation_db
```

```bash
blastn \
    -query DeerLake_diatom_rbcL.fasta \
    -db rbcL_orientation_db \
    -max_target_seqs 200 \
    -max_hsps 1 \
    -outfmt "6 sseqid sstart send pident length" \
    > rbcL_orientation.tsv
```

For each subject sequence:

```text
sstart < send   = same orientation as Deer Lake
sstart > send   = reverse orientation
```

Sequences with `sstart > send` were reverse complemented before alignment. The resulting orientation-corrected file was:

```text
all_rbcL_oriented.fasta
```

#### rbcL alignment and trimming

The 88 orientation-corrected sequences were aligned with Clustal Omega:

```bash
clustalo \
    -i all_rbcL_oriented.fasta \
    -o all_rbcL_oriented_aligned.fasta \
    --threads=32 \
    --force
```

The aligned dataset contained:

```text
Sequences:          88
Alignment length:   1,509 bp
```

The alignment was trimmed with TrimAl:

```bash
trimal \
    -in all_rbcL_oriented_aligned.fasta \
    -out all_rbcL_oriented_aligned_trimmed.fasta \
    -automated1
```

The final trimmed rbcL alignment contained:

```text
Sequences:          88
Alignment length:   1,473 bp
```

Coordinates appended to some NCBI FASTA identifiers were removed before tree inference:

```bash
sed -E 's/^>([^: ]+):[0-9]+-[0-9]+ />\1 /' \
    all_rbcL_oriented_aligned_trimmed.fasta \
    > all_rbcL_FINAL_alignment.fasta
```

The maximum-likelihood rbcL tree was inferred with IQ-TREE 2:

```bash
/home/ruchita.solanki/iqtree-2.2.2.7-Linux/bin/iqtree2 \
    -s all_rbcL_FINAL_alignment.fasta \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -nt 32 \
    -pre DeerLake_rbcL \
    -redo
```

The main output tree is:

```text
DeerLake_rbcL.treefile
```

Across all three marker analyses, phylogenetic placement was evaluated from the inferred tree topology and branch-support values rather than assigning species identity from BLAST similarity alone.

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
### 12.1 MetaQUAST comparison of a candidate draft genome
```bash
metaquast.py \
    18_diatom.fasta \
    -R /work/ebg_lab/eb/diatom_consortia/organelle/ref/ \
    -o ./18_metaquast_output
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
### 13.5 BUSCO assessment of the BRAKER4 genome input
Genome completeness was assessed before gene prediction using BUSCO v6.0.0 and the stramenopiles_odb12 lineage dataset.
```bash
busco \
    -i /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta \
    -l stramenopiles_odb12 \
    -m genome \
    -o busco_18_diatom_stramenopiles_odb12 \
    -c 24
```
This command evaluates the completeness of the diatom genome assembly using conserved single-copy orthologues from stramenopiles. BUSCO was run in eukaryotic genome mode using MetaEuk as the gene predictor.

BUSCO configuration:
```text
BUSCO version:          6.0.0
Lineage dataset:        stramenopiles_odb12
Reference genomes:      55
BUSCO groups:           697
Analysis mode:          euk_genome_met
Gene predictor:         MetaEuk 7.bba0d80
HMMER version:          3.4
```
BUSCO completeness:
```text
C:86.1%[S:82.4%,D:3.7%],F:2.2%,M:11.8%,n:697

Complete BUSCOs:                    600
Complete and single-copy BUSCOs:   574
Complete and duplicated BUSCOs:     26
Fragmented BUSCOs:                  15
Missing BUSCOs:                     82
Total BUSCO groups searched:       697
```
Assembly statistics reported by BUSCO:
```text
Scaffolds:       3,010
Contigs:         3,010
Total length:    82,172,226 bp
Sequence gaps:   0.000%
Scaffold N50:    48 kbp
Contig N50:      48 kbp
```
The assembly recovered 86.1% of the conserved stramenopiles_odb12 BUSCO set, with 82.4% present as single-copy genes and 3.7% duplicated.
### 13.6 BRAKER4 sample file
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
### 13.7 BRAKER4 configuration
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
### 13.8 Snakemake dry run
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
### 13.9 BRAKER4 execution
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
### 13.10 Final BRAKER4 outputs
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
### 13.11 Annotation statistics
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
### 13.12 BUSCO assessment of the final protein set
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
### 13.13 Annotation acceptance
The final BRAKER4 ET annotation was accepted for downstream analysis because it produced a plausible gene set for the diatom genome assembly, showed no internal stop codon issues in the predicted protein FASTA, and recovered 84.8% of the `stramenopiles_odb12` BUSCO protein set with low duplication.
Final accepted annotation files:
```text
DL_diatom.braker4.ET.gff3
DL_diatom.braker4.ET.gtf
DL_diatom.braker4.ET.proteins.faa
DL_diatom.braker4.ET.cds.fna
```
### 13.14 Rationale for ET mode instead of ETP
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
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=120:00:00
#SBATCH --mem=100G
#SBATCH --partition=cpu2025
####### Run your script #########################
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
<summary><strong>16. Historical pairwise genome comparison with <em>Phaeodactylum tricornutum</em></strong> - BLASTN, GFF3, and bedtools; superseded by Section 21</summary>

A pairwise genome comparison was performed between the BRAKER4-annotated diatom genome and the reference genome of *Phaeodactylum tricornutum*. This analysis was used as a nucleotide-level similarity screen and was not treated as a full orthology analysis. Raw BLASTN hits were retained without filtering and then linked to overlapping gene models in both genomes.
### 16.1 Working directory and input files
```bash
cd /work/ebg_lab/eb/diatom_consortia

mkdir -p phaeodactylum_to_diatom_blastn_redo/{00_inputs,01_db,02_blast,03_filtered,04_summary,logs,scripts}

cd phaeodactylum_to_diatom_blastn_redo
```
The diatom genome used for this comparison was:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
```
The accepted BRAKER4 ET annotation file was:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET/DL_diatom.braker4.ET.gff3
```
The diatom genome was linked into the input directory:
```bash
ln -sf /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta \
    00_inputs/diatom_genome.fasta
```
### 16.2 Download the *Phaeodactylum tricornutum* reference genome and GFF3
The *P. tricornutum* reference genome and annotation were downloaded from the NCBI RefSeq assembly `GCF_000150955.2_ASM15095v2`.
```bash
wget -O 00_inputs/Phaeodactylum_tricornutum_ASM15095v2_GCF_000150955.2_genomic.fna.gz \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/955/GCF_000150955.2_ASM15095v2/GCF_000150955.2_ASM15095v2_genomic.fna.gz
```
This command downloads the *P. tricornutum* genome FASTA file.

```bash
gunzip -c 00_inputs/Phaeodactylum_tricornutum_ASM15095v2_GCF_000150955.2_genomic.fna.gz \
    > 00_inputs/phaeodactylum_genome.fna
```
This command creates an uncompressed FASTA file while keeping the original compressed file.
```bash
wget -O 00_inputs/phaeodactylum_ASM15095v2.gff3.gz \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/955/GCF_000150955.2_ASM15095v2/GCF_000150955.2_ASM15095v2_genomic.gff.gz
```
This command downloads the matching *P. tricornutum* GFF3 annotation file.
```bash
gunzip -c 00_inputs/phaeodactylum_ASM15095v2.gff3.gz \
    > 00_inputs/phaeodactylum_ASM15095v2.gff3
```
This command creates an uncompressed GFF3 annotation file.
### 16.3 Check genome FASTA files
```bash
seqkit stats 00_inputs/diatom_genome.fasta 00_inputs/phaeodactylum_genome.fna
```
This command checks the number of sequences and total genome length for both FASTA files.
Output:
```text
file                                format  type  num_seqs     sum_len  min_len   avg_len    max_len
00_inputs/diatom_genome.fasta       FASTA   DNA      3,010  82,172,226      498  27,299.7    278,139
00_inputs/phaeodactylum_genome.fna  FASTA   DNA         88  27,450,724      450   311,940  2,535,400
```
### 16.4 Build the diatom BLAST database
```bash
makeblastdb \
    -in 00_inputs/diatom_genome.fasta \
    -dbtype nucl \
    -parse_seqids \
    -out 01_db/diatom_genome_blastdb \
    -title "DL_diatom_genome"
```
This command builds a nucleotide BLAST database from the diatom genome.
Output:
```text
Adding sequences from FASTA; added 3010 sequences
```
### 16.5 Run *Phaeodactylum* versus diatom BLASTN

```bash
blastn \
    -task dc-megablast \
    -query 00_inputs/phaeodactylum_genome.fna \
    -db 01_db/diatom_genome_blastdb \
    -out 02_blast/phaeodactylum_vs_diatom_dcmegablast.tsv \
    -evalue 1e-10 \
    -perc_identity 60 \
    -num_threads 16 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"
```
This command uses the *P. tricornutum* genome as the query and the diatom genome as the BLAST database. The raw BLASTN output was retained without additional filtering.
```bash
wc -l 02_blast/phaeodactylum_vs_diatom_dcmegablast.tsv
```
Output:
```text
28,934 raw BLASTN hits
```
A headered version of the BLASTN output was also created:
```bash
printf "pt_contig\tdiatom_contig\tpident\taln_len\tmismatch\tgapopen\tpt_start\tpt_end\tdiatom_start\tdiatom_end\tevalue\tbitscore\tpt_len\tdiatom_len\tqcovs\n" \
> 02_blast/phaeodactylum_vs_diatom_dcmegablast.header.tsv

cat 02_blast/phaeodactylum_vs_diatom_dcmegablast.tsv \
>> 02_blast/phaeodactylum_vs_diatom_dcmegablast.header.tsv
```
This command adds readable column names to the raw BLASTN output.
### 16.6 Convert *Phaeodactylum* genes to BED
```bash
awk -F'\t' '
BEGIN{OFS="\t"}
$3=="gene"{
  id=name=gene=locus="NA"
  n=split($9,a,";")
  for(i=1;i<=n;i++){
    split(a[i],b,"=")
    if(b[1]=="ID") id=b[2]
    else if(b[1]=="Name") name=b[2]
    else if(b[1]=="gene") gene=b[2]
    else if(b[1]=="locus_tag") locus=b[2]
  }
  gsub(/^gene-/,"",id)
  print $1, $4-1, $5, id, name, gene, locus, $7
}' 00_inputs/phaeodactylum_ASM15095v2.gff3 \
> 00_inputs/phaeodactylum_genes.bed
```
This command converts *P. tricornutum* gene coordinates from GFF3 to BED format and extracts gene IDs, gene names, gene symbols, locus tags, and strand information.
```bash
wc -l 00_inputs/phaeodactylum_genes.bed
```
Output:
```text
10,392 Phaeodactylum genes
```
### 16.7 Convert diatom BRAKER4 ET genes to BED
```bash
awk -F'\t' '
BEGIN{OFS="\t"}
$3=="gene"{
  id=gene_id="NA"
  n=split($9,a,";")
  for(i=1;i<=n;i++){
    split(a[i],b,"=")
    if(b[1]=="ID") id=b[2]
    else if(b[1]=="gene_id") gene_id=b[2]
  }
  if(gene_id=="NA") gene_id=id
  print $1, $4-1, $5, gene_id, id, $7
}' /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET/DL_diatom.braker4.ET.gff3 \
> 00_inputs/diatom_BRAKER_ET_genes.bed
```
This command converts the accepted BRAKER4 ET diatom gene models from GFF3 to BED format.

```bash
wc -l 00_inputs/diatom_BRAKER_ET_genes.bed
```
Output:

```text
15,102 BRAKER4 ET genes
```
### 16.8 Convert raw BLASTN hits to BED intervals
```bash
awk -F'\t' '
BEGIN{OFS="\t"}
{
  hit=sprintf("hit_%06d", NR)

  qstart=$7; qend=$8
  if(qstart <= qend){qbstart=qstart-1; qbend=qend; qstrand="+"}
  else{qbstart=qend-1; qbend=qstart; qstrand="-"}

  sstart=$9; send=$10
  if(sstart <= send){sbstart=sstart-1; sbend=send; sstrand="+"}
  else{sbstart=send-1; sbend=sstart; sstrand="-"}

  print $1, qbstart, qbend, hit, $3, $4, $11, $12, $13, $2, $9, $10, qstrand > "02_blast/phaeodactylum_blast_intervals.bed"
  print $2, sbstart, sbend, hit, $3, $4, $11, $12, $14, $1, $7, $8, sstrand > "02_blast/diatom_blast_intervals.bed"
}' 02_blast/phaeodactylum_vs_diatom_dcmegablast.tsv
```
This command converts each raw BLASTN hit into BED-like intervals on both the *Phaeodactylum* query genome and the diatom subject genome. Each BLASTN alignment was assigned a stable hit ID.
```bash
wc -l 02_blast/phaeodactylum_blast_intervals.bed
wc -l 02_blast/diatom_blast_intervals.bed
```
Output:
```text
28,934 Phaeodactylum-side BLAST intervals
28,934 diatom-side BLAST intervals
```
### 16.9 Link BLASTN hits to *Phaeodactylum* genes
```bash
bedtools intersect \
    -a 02_blast/phaeodactylum_blast_intervals.bed \
    -b 00_inputs/phaeodactylum_genes.bed \
    -wa -wb -loj \
> 03_filtered/phaeodactylum_hits_with_PT_genes.tsv
```
This command identifies *P. tricornutum* genes that overlap each BLASTN interval. The `-loj` option keeps all BLASTN hits, including hits that do not overlap an annotated *Phaeodactylum* gene.
```bash
wc -l 03_filtered/phaeodactylum_hits_with_PT_genes.tsv
```
Output:

```text
28,976 rows
```
The row count is slightly higher than the raw BLASTN hit count because some BLASTN intervals overlap more than one *Phaeodactylum* gene.
### 16.10 Link BLASTN hits to diatom BRAKER4 ET genes
```bash
bedtools intersect \
    -a 02_blast/diatom_blast_intervals.bed \
    -b 00_inputs/diatom_BRAKER_ET_genes.bed \
    -wa -wb -loj \
> 03_filtered/diatom_hits_with_BRAKER_ET_genes.tsv
```
This command identifies diatom BRAKER4 ET genes that overlap each BLASTN interval. The `-loj` option keeps all BLASTN hits, including hits that do not overlap an annotated diatom gene.
```bash
wc -l 03_filtered/diatom_hits_with_BRAKER_ET_genes.tsv
```
Output:
```text
29,005 rows
```
The row count is slightly higher than the raw BLASTN hit count because some BLASTN intervals overlap more than one diatom gene.
### 16.11 Create final BLASTN gene-linked table
The raw BLASTN output, the *Phaeodactylum* gene-overlap table, and the diatom BRAKER4 ET gene-overlap table were merged using a custom Python script.
```bash
python scripts/06_merge_phaeodactylum_blast_hits.py
```
This script keeps one row per raw BLASTN hit and adds overlapping gene information from both genomes. When multiple genes overlap the same BLASTN interval, gene IDs are collapsed into semicolon-separated fields.
Final output:
```text
04_summary/phaeodactylum_vs_diatom_BLASTN_with_PT_and_BRAKER_ET_genes.tsv
```
The final table contained:
```text
28,934 BLASTN hits
24 columns
```
Main columns:
```text
blast_hit_id
pt_contig
diatom_contig
pident
aln_len
mismatch
gapopen
pt_start
pt_end
diatom_start
diatom_end
evalue
bitscore
pt_len
diatom_len
qcovs
pt_gene_id
pt_gene_name
pt_gene_symbol
pt_locus_tag
pt_gene_strand
diatom_gene_id
diatom_gene_attr_id
diatom_gene_strand
```
### 16.12 Summary of final BLASTN gene-linked table
The final merged table was summarized to count how many raw BLASTN hits overlapped annotated genes in *Phaeodactylum*, BRAKER4 ET genes in the diatom genome, or genes on both sides.
Summary:
```text
total_blast_hits                         28,934
hits_with_PT_gene                        8,529
hits_with_diatom_BRAKER_ET_gene          6,693
hits_with_both_PT_and_diatom_gene        6,173
unique_PT_genes_hit                      3,202
unique_diatom_BRAKER_ET_genes_hit        3,386
```
### 16.13 Final comparative-genomics outputs
```text
00_inputs/phaeodactylum_genome.fna
00_inputs/phaeodactylum_ASM15095v2.gff3
00_inputs/phaeodactylum_genes.bed
00_inputs/diatom_BRAKER_ET_genes.bed

01_db/diatom_genome_blastdb.*

02_blast/phaeodactylum_vs_diatom_dcmegablast.tsv
02_blast/phaeodactylum_vs_diatom_dcmegablast.header.tsv
02_blast/phaeodactylum_blast_intervals.bed
02_blast/diatom_blast_intervals.bed

03_filtered/phaeodactylum_hits_with_PT_genes.tsv
03_filtered/diatom_hits_with_BRAKER_ET_genes.tsv

04_summary/phaeodactylum_vs_diatom_BLASTN_with_PT_and_BRAKER_ET_genes.tsv
```
### 16.14 Interpretation

The raw pairwise BLASTN comparison identified 28,934 nucleotide alignments between the *P. tricornutum* reference genome and the diatom genome. These raw alignments were retained without filtering to preserve the full nucleotide-level comparison.

Of these raw BLASTN hits, 8,529 overlapped annotated *P. tricornutum* genes, 6,693 overlapped BRAKER4 ET gene models in the diatom genome, and 6,173 overlapped genes on both sides. In total, the BLASTN hits represented 3,202 unique *P. tricornutum* genes and 3,386 unique diatom BRAKER4 ET genes.

This analysis provides a gene-linked nucleotide similarity table between the diatom genome and *P. tricornutum*. Because nucleotide-level similarity does not fully capture protein-level conservation or gene orthology, the table should be interpreted as a genome-level similarity screen rather than a definitive orthology map.

</details>

---

<details>
<summary><strong>17. Historical pairwise genome comparison with <em>Thalassiosira pseudonana</em></strong> - BLASTN, GFF3, bedtools, Python, and SLURM; superseded by Section 21</summary>

A second whole-genome nucleotide comparison was performed using *Thalassiosira pseudonana* CCMP1335. The analysis follows the same logic as the *Phaeodactylum tricornutum* comparison: the complete reference genome is used as the BLASTN query, the Deer Lake diatom genome is used as the nucleotide database, all reported alignments are retained, and overlapping gene models are assigned on both genomes.

This analysis is a nucleotide-level gene-linked similarity screen. It is not a reciprocal-best-hit analysis and should not be interpreted as confirmed orthology.

### 17.1 Reference assembly and working directories
The NCBI RefSeq assembly used for *T. pseudonana* was:

```text
Species:             Thalassiosira pseudonana
Strain:              CCMP1335
RefSeq accession:    GCF_000149405.2
Assembly name:       ASM14940v2
```

The main analysis directory was:

```text
/work/ebg_lab/eb/diatom_consortia/thalassiosira_to_diatom_blastn_redo
```

To reduce project-directory storage, the downloaded reference genome and GFF3 were stored under the home directory:

```text
$HOME/databases/thalassiosira_pseudonana/GCF_000149405.2_ASM14940v2/
```

Reference files:

```text
Thalassiosira_pseudonana_ASM14940v2_GCF_000149405.2_genomic.fna
thalassiosira_ASM14940v2.gff3
```

Only symbolic links to these files were created in the project analysis directory. The two Python helper scripts used by this section were:

```text
scripts/07_merge_thalassiosira_blast_hits.py
scripts/13_add_thalassiosira_yes_no.py
```

On ARC, `07_merge_thalassiosira_blast_hits.py` was placed at:

```text
/work/ebg_lab/eb/diatom_consortia/thalassiosira_pipeline/07_merge_thalassiosira_blast_hits.py
```

### 17.2 Comparative-genomics environment
The comparison was run in the `diatom_blast` Conda environment. The environment requires BLAST+, bedtools, seqkit, wget, Python 3, and pandas.

```bash
conda activate diatom_blast
conda install -c conda-forge pandas
```

The Python interpreter and pandas installation should be checked before submission:

```bash
which python3
python3 --version
python3 -c "import pandas; print(pandas.__version__)"
```

### 17.3 Complete SLURM workflow
The following block is the complete SLURM workflow used for the comparison. It can be saved on ARC as `run_thalassiosira_comparison_FIXED_PY3.sh`. The shell workflow is documented here directly rather than maintained as a separate script entry in the repository structure.

```bash
#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=120:00:00
#SBATCH --mem=100G
#SBATCH --partition=cpu2025
####### Run your script #########################
set -euo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate diatom_blast

THREADS=${SLURM_CPUS_PER_TASK:-32}

# Whole-genome nucleotide comparison:
# Thalassiosira pseudonana CCMP1335 -> Deer Lake diatom genome
# Mirrors the existing Phaeodactylum tricornutum dc-megablast workflow.

BASE=/work/ebg_lab/eb/diatom_consortia
WORKDIR=${BASE}/thalassiosira_to_diatom_blastn_redo
DIATOM_GENOME=${BASE}/metatranscriptomics/genome_index/18_diatom.fasta
DIATOM_GFF=${BASE}/metatranscriptomics/BRAKER4/final_annotation_ET/DL_diatom.braker4.ET.gff3

REF_ACCESSION=GCF_000149405.2
REF_ASSEMBLY=ASM14940v2
REF_PREFIX=${REF_ACCESSION}_${REF_ASSEMBLY}
REF_FTP=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/405/${REF_PREFIX}

# Store the downloaded Thalassiosira reference under the home directory
# rather than in the project analysis directory.
REF_HOME=${HOME}/databases/thalassiosira_pseudonana/${REF_PREFIX}
REF_GENOME=${REF_HOME}/Thalassiosira_pseudonana_${REF_ASSEMBLY}_${REF_ACCESSION}_genomic.fna
REF_GFF=${REF_HOME}/thalassiosira_${REF_ASSEMBLY}.gff3

# Use an explicit absolute path for the companion Python script.
# This avoids SLURM resolving the script location under /var/spool/slurmd/.
PIPELINE_DIR=/work/ebg_lab/eb/diatom_consortia/thalassiosira_pipeline
MERGE_SCRIPT=${PIPELINE_DIR}/07_merge_thalassiosira_blast_hits.py

echo "Pipeline directory: ${PIPELINE_DIR}"
echo "Merge script: ${MERGE_SCRIPT}"

for program in wget gzip seqkit makeblastdb blastn bedtools python3; do
    command -v "${program}" >/dev/null 2>&1 || {
        echo "ERROR: Required program not found in PATH: ${program}" >&2
        exit 1
    }
done

[[ -s "${DIATOM_GENOME}" ]] || { echo "ERROR: Missing ${DIATOM_GENOME}" >&2; exit 1; }
[[ -s "${DIATOM_GFF}" ]] || { echo "ERROR: Missing ${DIATOM_GFF}" >&2; exit 1; }
[[ -s "${MERGE_SCRIPT}" ]] || { echo "ERROR: Missing ${MERGE_SCRIPT}" >&2; exit 1; }

mkdir -p "${WORKDIR}"/{00_inputs,01_db,02_blast,03_filtered,04_summary,logs,scripts}
mkdir -p "${REF_HOME}"
cd "${WORKDIR}"

ln -sfn "${DIATOM_GENOME}" 00_inputs/diatom_genome.fasta
cp -f "${MERGE_SCRIPT}" scripts/07_merge_thalassiosira_blast_hits.py

# Stream and decompress directly into the home-directory reference folder.
# This avoids storing both .gz and uncompressed copies at the same time.
if [[ ! -s "${REF_GENOME}" ]]; then
    echo "Downloading Thalassiosira genome to ${REF_GENOME}"
    wget -qO- "${REF_FTP}/${REF_PREFIX}_genomic.fna.gz" \
        | gzip -dc > "${REF_GENOME}.tmp"
    mv "${REF_GENOME}.tmp" "${REF_GENOME}"
fi

if [[ ! -s "${REF_GFF}" ]]; then
    echo "Downloading Thalassiosira GFF3 to ${REF_GFF}"
    wget -qO- "${REF_FTP}/${REF_PREFIX}_genomic.gff.gz" \
        | gzip -dc > "${REF_GFF}.tmp"
    mv "${REF_GFF}.tmp" "${REF_GFF}"
fi

# Keep only lightweight symbolic links inside the project analysis directory.
ln -sfn "${REF_GENOME}" 00_inputs/thalassiosira_genome.fna
ln -sfn "${REF_GFF}" 00_inputs/thalassiosira_${REF_ASSEMBLY}.gff3

seqkit stats 00_inputs/diatom_genome.fasta 00_inputs/thalassiosira_genome.fna \
    | tee 04_summary/genome_input_stats.txt

makeblastdb \
    -in 00_inputs/diatom_genome.fasta \
    -dbtype nucl \
    -parse_seqids \
    -out 01_db/diatom_genome_blastdb \
    -title "DL_diatom_genome"

blastn \
    -task dc-megablast \
    -query 00_inputs/thalassiosira_genome.fna \
    -db 01_db/diatom_genome_blastdb \
    -out 02_blast/thalassiosira_vs_diatom_dcmegablast.tsv \
    -evalue 1e-10 \
    -perc_identity 60 \
    -num_threads "${THREADS}" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"

printf "thal_contig\tdiatom_contig\tpident\taln_len\tmismatch\tgapopen\tthal_start\tthal_end\tdiatom_start\tdiatom_end\tevalue\tbitscore\tthal_len\tdiatom_len\tqcovs\n" \
    > 02_blast/thalassiosira_vs_diatom_dcmegablast.header.tsv
cat 02_blast/thalassiosira_vs_diatom_dcmegablast.tsv \
    >> 02_blast/thalassiosira_vs_diatom_dcmegablast.header.tsv

# NCBI Thalassiosira gene models -> BED.
awk -F'\t' '
BEGIN{OFS="\t"}
$3=="gene"{
  id=name=gene=locus="NA"
  n=split($9,a,";")
  for(i=1;i<=n;i++){
    split(a[i],b,"=")
    if(b[1]=="ID") id=b[2]
    else if(b[1]=="Name") name=b[2]
    else if(b[1]=="gene") gene=b[2]
    else if(b[1]=="locus_tag") locus=b[2]
  }
  gsub(/^gene-/,"",id)
  print $1, $4-1, $5, id, name, gene, locus, $7
}' 00_inputs/thalassiosira_${REF_ASSEMBLY}.gff3 \
> 00_inputs/thalassiosira_genes.bed

# Deer Lake BRAKER4 gene models -> BED.
awk -F'\t' '
BEGIN{OFS="\t"}
$3=="gene"{
  id=gene_id="NA"
  n=split($9,a,";")
  for(i=1;i<=n;i++){
    split(a[i],b,"=")
    if(b[1]=="ID") id=b[2]
    else if(b[1]=="gene_id") gene_id=b[2]
  }
  if(gene_id=="NA") gene_id=id
  print $1, $4-1, $5, gene_id, id, $7
}' "${DIATOM_GFF}" \
> 00_inputs/diatom_BRAKER_ET_genes.bed

# Convert each raw BLASTN alignment to BED-like intervals on both genomes.
rm -f 02_blast/thalassiosira_blast_intervals.bed 02_blast/diatom_blast_intervals.bed
awk -F'\t' '
BEGIN{OFS="\t"}
{
  hit=sprintf("hit_%06d", NR)

  qstart=$7; qend=$8
  if(qstart <= qend){qbstart=qstart-1; qbend=qend; qstrand="+"}
  else{qbstart=qend-1; qbend=qstart; qstrand="-"}

  sstart=$9; send=$10
  if(sstart <= send){sbstart=sstart-1; sbend=send; sstrand="+"}
  else{sbstart=send-1; sbend=sstart; sstrand="-"}

  print $1, qbstart, qbend, hit, $3, $4, $11, $12, $13, $2, $9, $10, qstrand > "02_blast/thalassiosira_blast_intervals.bed"
  print $2, sbstart, sbend, hit, $3, $4, $11, $12, $14, $1, $7, $8, sstrand > "02_blast/diatom_blast_intervals.bed"
}' 02_blast/thalassiosira_vs_diatom_dcmegablast.tsv

bedtools intersect \
    -a 02_blast/thalassiosira_blast_intervals.bed \
    -b 00_inputs/thalassiosira_genes.bed \
    -wa -wb -loj \
> 03_filtered/thalassiosira_hits_with_Thalassiosira_genes.tsv

bedtools intersect \
    -a 02_blast/diatom_blast_intervals.bed \
    -b 00_inputs/diatom_BRAKER_ET_genes.bed \
    -wa -wb -loj \
> 03_filtered/diatom_hits_with_BRAKER_ET_genes.tsv

echo "Python executable: $(command -v python3)"
python3 --version
python3 scripts/07_merge_thalassiosira_blast_hits.py

{
    echo -e "raw_blast_hits\t$(wc -l < 02_blast/thalassiosira_vs_diatom_dcmegablast.tsv)"
    echo -e "Thalassiosira_genes\t$(wc -l < 00_inputs/thalassiosira_genes.bed)"
    echo -e "diatom_BRAKER_ET_genes\t$(wc -l < 00_inputs/diatom_BRAKER_ET_genes.bed)"
    echo -e "Thalassiosira_overlap_rows\t$(wc -l < 03_filtered/thalassiosira_hits_with_Thalassiosira_genes.tsv)"
    echo -e "diatom_overlap_rows\t$(wc -l < 03_filtered/diatom_hits_with_BRAKER_ET_genes.tsv)"
} | tee 04_summary/file_counts.tsv

echo "Completed Thalassiosira pseudonana comparison."
echo "Final table: ${WORKDIR}/04_summary/thalassiosira_vs_diatom_BLASTN_with_Thalassiosira_and_BRAKER_ET_genes.tsv"
```

The workflow performs the following operations:

```text
activates the diatom_blast Conda environment
stores the Thalassiosira reference genome and GFF3 in the home directory
creates symbolic links to the reference files in the project directory
builds a Deer Lake nucleotide BLAST database
runs dc-megablast with an e-value threshold of 1e-10 and minimum identity of 60%
retains all reported BLASTN alignments without additional post-BLAST filtering
converts reference and Deer Lake gene models to BED format
assigns overlapping genes with bedtools intersect -loj
collapses overlap results to one row per raw BLASTN hit
writes the final gene-linked table and summary files
```

### 17.4 Submit the SLURM job
The SLURM file and `07_merge_thalassiosira_blast_hits.py` were placed in the same ARC directory:

```text
/work/ebg_lab/eb/diatom_consortia/thalassiosira_pipeline/
```

The job was submitted with:

```bash
cd /work/ebg_lab/eb/diatom_consortia/thalassiosira_pipeline

ls -lh \
    run_thalassiosira_comparison_FIXED_PY3.sh \
    07_merge_thalassiosira_blast_hits.py

sbatch run_thalassiosira_comparison_FIXED_PY3.sh
```

### 17.5 BLASTN comparison settings
The complete *T. pseudonana* genome was used as the query and the Deer Lake assembly was used as the database. The BLASTN settings were identical to the *P. tricornutum* comparison:

```text
BLAST task:             dc-megablast
E-value threshold:      1e-10
Minimum identity:       60%
Additional filtering:   none
Threads:                32
```

No minimum alignment-length, query-coverage, bitscore, or additional identity filter was applied after BLASTN.

### 17.6 Gene-overlap assignment
The NCBI *T. pseudonana* GFF3 was converted to a BED table containing reference gene identifiers, names, symbols, locus tags, coordinates, and strand. The Deer Lake BRAKER4 ET GFF3 was converted to a separate BED table containing BRAKER4 gene identifiers, coordinates, and strand.

Each raw BLASTN alignment received a stable identifier:

```text
hit_000001
hit_000002
hit_000003
...
```

Reference and Deer Lake gene overlaps were assigned using `bedtools intersect -loj`. The `-loj` option retained BLASTN alignments even when no annotated gene overlapped the aligned interval.

The Python merge script then:

```text
retained one output row per raw BLASTN hit
collapsed multiple overlapping genes with semicolons
kept alignments without annotated-gene overlaps
wrote gene-linked counts to a summary file
did not perform additional BLASTN filtering
```

### 17.7 Main comparison outputs
Final gene-linked output:

```text
/work/ebg_lab/eb/diatom_consortia/thalassiosira_to_diatom_blastn_redo/04_summary/thalassiosira_vs_diatom_BLASTN_with_Thalassiosira_and_BRAKER_ET_genes.tsv
```

Summary output:

```text
/work/ebg_lab/eb/diatom_consortia/thalassiosira_to_diatom_blastn_redo/04_summary/thalassiosira_vs_diatom_BLASTN_summary.txt
```

Main columns in the gene-linked table:

```text
blast_hit_id
thal_contig
diatom_contig
pident
aln_len
mismatch
gapopen
thal_start
thal_end
diatom_start
diatom_end
evalue
bitscore
thal_len
diatom_len
qcovs
thal_gene_id
thal_gene_name
thal_gene_symbol
thal_locus_tag
thal_gene_strand
diatom_gene_id
diatom_gene_attr_id
diatom_gene_strand
```

Run-specific hit and gene counts should be added after the SLURM job completes.

### 17.8 Add the *Thalassiosira* yes/no field
The final gene-linked BLASTN table is used as a gene-root yes/no lookup. A Deer Lake BRAKER4 gene is marked `yes` only when the same raw BLASTN hit overlaps both an annotated *T. pseudonana* gene and a Deer Lake BRAKER4 ET gene model.

BRAKER4 isoform IDs are reduced to their gene root during lookup:

```text
g10009.t1 → g10009
```

The existing *Phaeodactylum* field is retained and a separate column is added:

```text
present_in_Thalassiosira_pseudonana
```

Run the integration script from the clean rebuild directory:

```bash
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/transdecoder_to_braker_ID_bridge/CLEAN_REBUILD_FROM_RAW

conda activate diatom_blast

python /work/ebg_lab/eb/diatom_consortia/thalassiosira_pipeline/13_add_thalassiosira_yes_no.py \
    --input-final 09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table.tsv
```

Default outputs:

```text
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table_PT_TP.tsv
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table_PT_TP_sorted_by_Average_TPM.tsv
09_final/DL_diatom_FINAL_gene_table_for_boss_PT_TP.tsv
```

The final comparison columns are:

```text
present_in_Phaeodactylum_tricornutum
present_in_Thalassiosira_pseudonana
```

</details>

---

<details>
<summary><strong>18. Clean BRAKER4 isoform-level gene table construction</strong> - functional annotation, Average_TPM, and reference-diatom comparison fields</summary>

This section describes the clean BRAKER4 isoform-level tables used for pathway curation, manual review, and nucleotide-level comparison with *Phaeodactylum tricornutum* and *Thalassiosira pseudonana*.

The final tables retain one row per BRAKER4 predicted protein isoform. BRAKER4 isoform IDs are not collapsed, GenBank-derived organelle rows are not appended, all-hit TPM summaries are not used, and detailed BLASTN alignment columns are not carried into the final review table.

### 18.1 Clean rebuild directory
The clean rebuild was performed in:

```text
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/transdecoder_to_braker_ID_bridge/CLEAN_REBUILD_FROM_RAW
```

This directory separates the accepted functional annotation, BRAKER4 length information, expression integration, comparative-genomics lookups, and final tables from older intermediate outputs.

### 18.2 Accepted annotation inputs
The rebuild used:

```text
01_input/diatom_predicted_proteins.fa
01_input/DL_diatom.braker4.ET.proteins.faa
01_input/DL_diatom.braker4.ET.gff3
```

Functional annotation layers:

```text
02_diamond/DL_diatom_braker4_ET_vs_swissprot.tsv
02_diamond/DL_diatom_braker4_ET_vs_uniprot_bacillariophyta.tsv
05_interproscan/DL_diatom_braker4_ET_interproscan.tsv
06_combined_annotation/DL_diatom_antifam_flagged_proteins.tsv
```

Expression inputs:

```text
01_input/transcriptome_orfs.transdecoder.clean.pep
07_expression/master_with_custom_broad_categories.csv
```

Comparative-genomics inputs:

```text
/work/ebg_lab/eb/diatom_consortia/phaeodactylum_to_diatom_blastn_redo/04_summary/phaeodactylum_vs_diatom_BLASTN_with_PT_and_BRAKER_ET_genes.tsv

/work/ebg_lab/eb/diatom_consortia/thalassiosira_to_diatom_blastn_redo/04_summary/thalassiosira_vs_diatom_BLASTN_with_Thalassiosira_and_BRAKER_ET_genes.tsv
```

The comparative-genomics files are used only as gene-root yes/no lookups.

### 18.3 Functional annotation integration
Functional annotations were integrated from Swiss-Prot, UniProtKB Bacillariophyta, InterProScan, and AntiFam.

The master table contains one row per BRAKER4 protein isoform:

```text
16,947 BRAKER4 protein isoforms
16,948 lines including the header
```

Master functional annotation table:

```text
06_combined_annotation/DL_diatom_master_functional_annotation.tsv
```

The two AntiFam-flagged proteins were retained as warning flags:

```text
g10893.t1    ANF00012    tRNA
g11404.t1    ANF00005    Antisense to 23S rRNA
```

### 18.4 BRAKER4 coordinates and lengths
BRAKER4 GFF3 coordinates were added to each isoform, including:

```text
contig ID
gene start
gene end
strand
gene length
CDS length
protein length
```

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

### 18.5 TransDecoder ORF to BRAKER4 mapping
TransDecoder peptides were searched against the accepted BRAKER4 ET protein set using DIAMOND BLASTP.

Raw output:

```text
07_expression/transcriptome_ORFs_vs_BRAKER4_ET_proteins.tsv
```

Raw mapping summary:

```text
118,472 DIAMOND hit rows
39,935 unique TransDecoder ORFs with at least one BRAKER4 hit
14,473 unique BRAKER4 proteins hit by at least one reported ORF hit
```

One best BRAKER4 hit was retained per TransDecoder ORF:

```text
07_expression/transdecoder_ORFs_to_BRAKER4_best_hit_per_ORF.tsv
```

Best-hit summary:

```text
39,935 ORF mappings
12,277 unique BRAKER4 proteins represented by best ORF mappings
```

### 18.6 Average_TPM integration
Only the trusted transcriptome fields were used:

```text
orf
Average_TPM
```

If several TransDecoder ORFs mapped best to the same BRAKER4 isoform, the ORF with the highest valid `Average_TPM` was retained.

Output:

```text
07_expression/DL_diatom_master_functional_annotation_lengths_Average_TPM.tsv
```

Final expression counts:

```text
Annotation rows:                         16,947
BRAKER proteins with Average_TPM:        12,276
BRAKER proteins without Average_TPM:      4,671
BRAKER proteins with TransDecoder ORF:   12,276
BRAKER proteins without ORF assignment:   4,671
```

No all-hit TPM sums, means, hit counts, or mapped-ORF lists were retained.

### 18.7 Compartment assignment
Compartment labels were assigned from the BRAKER4 contig ID:

```text
nuclear
plastid_like
mito_like
```

Observed counts:

```text
nuclear        16,873
plastid_like       72
mito_like           2
```

No GenBank-derived organelle gene rows were appended.

### 18.8 *Phaeodactylum tricornutum* yes/no lookup
A BRAKER4 isoform is marked `yes` when its gene root occurs in a BLASTN alignment that overlaps both an annotated *P. tricornutum* gene and a Deer Lake BRAKER4 ET gene.

Column:

```text
present_in_Phaeodactylum_tricornutum
```

Observed counts from the completed redo comparison:

```text
no     13,240
yes     3,707
```

This field represents gene-linked nucleotide similarity and not confirmed orthology.

### 18.9 *Thalassiosira pseudonana* yes/no lookup
The same logic is applied to the *T. pseudonana* comparison. A BRAKER4 isoform is marked `yes` only when its gene root occurs in a BLASTN alignment linked to annotated genes on both genomes.

Column:

```text
present_in_Thalassiosira_pseudonana
```

The lookup preserves full isoform IDs in the final table while matching comparison results by gene root:

```text
g10009.t1 → g10009
```

Run-specific yes/no counts are generated by `13_add_thalassiosira_yes_no.py` and should be added here after the comparison completes.

### 18.10 Final clean isoform-level tables
The base clean table containing the completed *Phaeodactylum* field remains:

```text
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table.tsv
```

After adding the *Thalassiosira* field, the combined outputs are:

```text
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table_PT_TP.tsv
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table_PT_TP_sorted_by_Average_TPM.tsv
```

Each table retains:

```text
16,947 BRAKER4 isoform rows
one header row
one row per predicted protein isoform
```

### 18.11 Simplified review table
The combined review table is:

```text
09_final/DL_diatom_FINAL_gene_table_for_boss_PT_TP.tsv
```

Columns:

```text
gene_id
contig_id
diatom_compartment
diatom_gene_length_bp
functional_annotation
diatom_Average_TPM
present_in_Phaeodactylum_tricornutum
present_in_Thalassiosira_pseudonana
```

Column interpretation:

```text
gene_id
  Original BRAKER4 isoform ID.

contig_id
  Deer Lake diatom contig containing the BRAKER4 gene model.

diatom_compartment
  Nuclear, plastid-like, or mitochondrion-like assignment based on contig identity.

diatom_gene_length_bp
  Gene length calculated from BRAKER4 GFF3 coordinates.

functional_annotation
  Recommended annotation from the integrated Swiss-Prot, Bacillariophyta, InterProScan, and AntiFam evidence layers.

diatom_Average_TPM
  Expression value transferred through the selected TransDecoder ORF-to-BRAKER4 mapping.

present_in_Phaeodactylum_tricornutum
  Gene-linked nucleotide-similarity yes/no field from the Phaeodactylum comparison.

present_in_Thalassiosira_pseudonana
  Gene-linked nucleotide-similarity yes/no field from the Thalassiosira comparison.
```

The review table is intended for manual pathway curation and biological interpretation without carrying forward detailed intermediate BLASTN or annotation fields. The PT and TP nucleotide fields are retained as exploratory similarity screens. The primary final protein comparison and four-species orthogroup-sharing table are documented in Section 19.

</details>

---

<details>
<summary><strong>19. Repeat aware nuclear proteome curation and secondary five species protein orthology</strong> - RepeatModeler, RepeatMasker, OrthoFinder, DIAMOND, FAMSA, FastTree, seqkit, and Python</summary>

This section describes the final protein comparison used to replace the earlier nucleotide similarity screen as the primary comparative analysis. The workflow first generated a nonredundant Deer Lake nuclear representative proteome, screened predicted coding sequences against a de novo repeat library, removed a small set of high confidence transposable element derived models, standardized four reference diatom proteomes, and then inferred orthogroups with OrthoFinder.

The four reference diatoms were:

```text
Phaeodactylum tricornutum       GCF_000150955.2
Thalassiosira pseudonana        GCF_000149405.2
Seminavis robusta               GCA_903772945.1
Nitzschia inconspicua           GCA_019154785.2
```

The *Nitzschia inconspicua* assembly is diploid. Its nuclear protein set was therefore retained without sequence identity based deduplication. OrthoFinder results involving this species were interpreted as orthogroup sharing rather than gene copy number differences.

### 19.1 Working directories

```bash
BASE=/work/ebg_lab/eb/diatom_consortia
COMP=${BASE}/comparative_genomics

mkdir -p ${COMP}/00_DL_reference
mkdir -p ${COMP}/01_SR_reference
mkdir -p ${COMP}/02_PT_reference
mkdir -p ${COMP}/03_TP_reference
mkdir -p ${COMP}/04_repeat_analysis
mkdir -p ${COMP}/orthofinder_input
mkdir -p ${COMP}/orthofinder_results
```

Custom Python scripts used in this section are stored externally in `scripts/`:

```text
scripts/18_make_DL_nuclear_representative_proteome.py
scripts/19_calculate_DL_CDS_repeat_overlap.py
scripts/20_classify_DL_TE_candidates.py
scripts/21_prepare_reference_proteome.py
scripts/22_make_final_orthofinder_gene_table.py
```

### 19.2 Deer Lake nuclear genome used for repeat analysis

The nuclear enriched Deer Lake assembly was:

```text
/work/ebg_lab/eb/diatom_consortia/nuclear_genome_filtering_18_diatom/18_diatom_nuclear_enriched.v1.fasta
```

It was linked into the repeat analysis directory:

```bash
cd ${COMP}/04_repeat_analysis

ln -sfn \
${BASE}/nuclear_genome_filtering_18_diatom/18_diatom_nuclear_enriched.v1.fasta \
DL_nuclear_genome.fasta

seqkit stats DL_nuclear_genome.fasta
```

Observed assembly statistics:

```text
Sequences:       3,007
Total length:   81,911,772 bp
Minimum:               498 bp
Mean:              27,240.4 bp
Maximum:           278,139 bp
```

### 19.3 De novo repeat discovery with RepeatModeler

The ARC software modules were used instead of a Conda RepeatModeler environment because the cluster modules provided a working Perl and RepeatScout dependency stack.

```text
RepeatModeler 2.0.1
RepeatMasker 4.1.1
RepeatScout 1.0.6
RMBlast 2.10.0
TRF 4.09
```

The standard RepeatModeler RECON and RepeatScout workflow was used. `LTRStruct` was not enabled in this run.

The following SLURM script was used directly on ARC:

```bash
#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=120:00:00
#SBATCH --mem=100G
#SBATCH --partition=cpu2025
####### Run your script #########################

set -euo pipefail

module purge
module load repeatmasker/4.1.1
module load repeatmodeler/2.0.1

BASE=/work/ebg_lab/eb/diatom_consortia/comparative_genomics/04_repeat_analysis
GENOME=${BASE}/DL_nuclear_genome.fasta
WORKDIR=${BASE}/01_repeatmodeler
DBNAME=DL_nuclear_repeatmodeler_db

mkdir -p "${WORKDIR}"
mkdir -p "${BASE}/tmp"

export TMPDIR="${BASE}/tmp"

cd "${WORKDIR}"
ln -sfn "${GENOME}" DL_nuclear_genome.fasta

BuildDatabase \
    -name "${DBNAME}" \
    DL_nuclear_genome.fasta

RepeatModeler \
    -database "${DBNAME}" \
    -pa 8 \
    2>&1 | tee repeatmodeler.log

LIB=$(find . -type f -name "consensi.fa.classified" | head -n 1)

if [[ -z "${LIB}" ]]; then
    echo "ERROR: consensi.fa.classified was not produced." >&2
    exit 1
fi

cp "${LIB}" \
    "${BASE}/DL_RepeatModeler_consensi.fa.classified"

grep -c '^>' "${BASE}/DL_RepeatModeler_consensi.fa.classified" \
    > "${BASE}/RepeatModeler_family_count.txt"
```

The final de novo repeat library contained:

```text
Repeat consensus families:   434
Total consensus length:       733,476 bp
Mean consensus length:        1,690 bp
Maximum consensus length:     8,276 bp
```

Family classifications were summarized with:

```bash
grep '^>' DL_RepeatModeler_consensi.fa.classified \
| sed 's/^.*#//' \
| sed 's/ .*//' \
| sort \
| uniq -c \
| sort -nr \
> DL_RepeatModeler_family_classes.txt
```

Observed family counts included:

```text
Unknown                 249
LTR/Copia                89
DNA/PIF-Harbinger        33
LTR/Ngaro                12
LTR/Gypsy                12
DNA/PIF-HarbS            11
LINE/CRE-Ambal            8
DNA/Sola-1                5
DNA/TcMar-Sagan           4
DNA/TcMar-Tc2             3
DNA/TcMar-Ant1            2
DNA/MULE-MuDR             2
DNA/TcMar-Stowaway        1
DNA/TcMar-m44             1
DNA/PiggyBac              1
DNA/Crypton-F             1
```

### 19.4 RepeatMasker annotation of the Deer Lake nuclear genome

The RepeatModeler library was used as a custom RepeatMasker library.

```bash
#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=120:00:00
#SBATCH --mem=100G
#SBATCH --partition=cpu2025
####### Run your script #########################

set -euo pipefail

module purge
module load repeatmasker/4.1.1
module load repeatmodeler/2.0.1

BASE=/work/ebg_lab/eb/diatom_consortia/comparative_genomics/04_repeat_analysis
GENOME=${BASE}/DL_nuclear_genome.fasta
LIB=${BASE}/DL_RepeatModeler_consensi.fa.classified
OUTDIR=${BASE}/02_repeatmasker

mkdir -p "${OUTDIR}"

RepeatMasker \
    -e rmblast \
    -pa 8 \
    -lib "${LIB}" \
    -xsmall \
    -gff \
    -gccalc \
    -dir "${OUTDIR}" \
    "${GENOME}" \
    2>&1 | tee "${OUTDIR}/repeatmasker.log"

TBL=$(find "${OUTDIR}" -maxdepth 1 -type f -name "*.tbl" | head -n 1)
cp "${TBL}" "${BASE}/DL_RepeatMasker_summary.tbl"
```

RepeatMasker identified 45,345,998 bp of masked sequence, corresponding to 55.36% of the 81.91 Mb nuclear assembly. Interspersed repeats accounted for 54.55% of the assembly. LTR elements were the largest classified component, with Ty1/Copia sequence occupying 28.09% of the genome.

```text
All masked sequence                55.36%
Interspersed repeats               54.55%
Retroelements                      37.00%
LTR elements                       36.68%
Ty1/Copia                          28.09%
Gypsy/DIRS1                         3.22%
LINEs                               0.33%
DNA transposons                     7.81%
Tourist/Harbinger                   4.02%
Tc1/IS630/Pogo                      1.65%
Unclassified interspersed repeats   9.74%
Simple repeats                      0.77%
Low complexity                      0.03%
```

This repeat analysis was used as quality control for the predicted protein set. The genome was not reannotated after this step.

### 19.5 Deer Lake representative nuclear proteome

The BRAKER4 ET annotation contained 15,102 genes and 16,947 protein isoforms. A single representative protein was selected per nuclear gene. The two AntiFam flagged gene roots were excluded, and the longest valid protein isoform was retained for each remaining gene. Extremely short proteins were removed only when they were shorter than 50 amino acids, lacked an Average_TPM value, and were unannotated.

Run:

```bash
python scripts/18_make_DL_nuclear_representative_proteome.py \
    --gene-table ${BASE}/metatranscriptomics/transdecoder_to_braker_ID_bridge/CLEAN_REBUILD_FROM_RAW/09_final/DL_diatom_FINAL_gene_table_for_boss_PT_TP_NI.tsv \
    --proteins ${BASE}/metatranscriptomics/BRAKER4/final_annotation_ET/DL_diatom.braker4.ET.proteins.faa \
    --out-fasta ${COMP}/00_DL_reference/DL_diatom_FINAL_nuclear_representative_proteome.faa \
    --out-map ${COMP}/00_DL_reference/DL_diatom_FINAL_nuclear_representative_proteome_map.tsv \
    --antifam-gene-roots g10893,g11404 \
    --min-aa 50
```

Observed counts:

```text
Representative nuclear genes before short protein filter:   15,030
Short unsupported proteins removed:                             13
Representative proteins after this filter:                  15,017
```

### 19.6 CDS overlap with interspersed repeats

Repeat overlap was evaluated across CDS bases rather than complete gene spans so that intronic repeats would not automatically classify an otherwise valid gene model as repeat derived.

The RepeatMasker `.out` file was used because it retains repeat class information. Overlapping RepeatMasker hits were resolved by assigning each genomic segment to the highest scoring hit before CDS overlap was calculated.

```bash
mkdir -p ${COMP}/04_repeat_analysis/03_CDS_repeat_overlap

python scripts/19_calculate_DL_CDS_repeat_overlap.py \
    --map ${COMP}/00_DL_reference/DL_diatom_FINAL_nuclear_representative_proteome_map.tsv \
    --gff ${BASE}/metatranscriptomics/BRAKER4/final_annotation_ET/DL_diatom.braker4.ET.gff3 \
    --repeatmasker-out ${COMP}/04_repeat_analysis/02_repeatmasker/DL_nuclear_genome.fasta.out \
    --outdir ${COMP}/04_repeat_analysis/03_CDS_repeat_overlap \
    --candidate-threshold 50
```

Observed interspersed repeat overlap across the 15,017 representative CDS models:

```text
0% overlap             14,377 genes   95.74%
>0 to 10%                 329 genes    2.19%
>10 to 25%                 63 genes    0.42%
>25 to 50%                 59 genes    0.39%
>50 to 75%                 32 genes    0.21%
>75 to 100%               157 genes    1.05%
```

Thus, 189 genes had at least 50% CDS overlap with an interspersed repeat and were carried forward for targeted review.

### 19.7 Conservative transposable element model filtering

The 189 repeat overlapping models were separated into three review tiers. A model was classified as high confidence transposable element derived only when at least 50% of its CDS overlapped an interspersed repeat and its functional annotation independently indicated transposable element associated activity, including reverse transcriptase, transposase, integrase, gag, Copia, or DDE transposase related functions.

```bash
python scripts/20_classify_DL_TE_candidates.py \
    --candidates ${COMP}/04_repeat_analysis/03_CDS_repeat_overlap/DL_CDS_repeat_overlap_candidates_ge50.tsv \
    --outdir ${COMP}/04_repeat_analysis/03_CDS_repeat_overlap
```

Observed classification:

```text
Tier 1, high confidence TE derived:      38
Tier 2, probable TE, retained for review: 65
Tier 3, retained for review:              86
```

Only the 38 Tier 1 models were removed from the comparative proteome. Tier 2 and Tier 3 models were retained because repeat overlap alone was not considered sufficient evidence for deletion.

```bash
IN=${COMP}/00_DL_reference/DL_diatom_FINAL_nuclear_representative_proteome.faa
REMOVE=${COMP}/04_repeat_analysis/03_CDS_repeat_overlap/DL_TIER1_TE_gene_roots.txt
OUT=${COMP}/00_DL_reference/DL_diatom_FINAL_nuclear_representative_proteome_TEfiltered.faa

seqkit grep \
    -v \
    -f "${REMOVE}" \
    "${IN}" \
    > "${OUT}"

seqkit stats "${OUT}"

# Filter the matching metadata map with the same 38 gene roots.
awk 'NR==FNR {remove[$1]=1; next} FNR==1 || !($1 in remove)' \
    "${REMOVE}" \
    ${COMP}/00_DL_reference/DL_diatom_FINAL_nuclear_representative_proteome_map.tsv \
    > ${COMP}/00_DL_reference/DL_diatom_FINAL_nuclear_representative_proteome_TEfiltered_map.tsv
```

Final Deer Lake comparative proteome:

```text
Proteins:        14,979
Total length:     7,329,909 aa
Minimum:                 67 aa
Mean:                 489.3 aa
Maximum:              7,782 aa
```

The matching metadata map was filtered using the same Tier 1 gene root list so that the final FASTA and table contained the same 14,979 genes.

### 19.8 Reference proteome standardization

Reference proteomes and matching GFF3 files were obtained from NCBI. One representative protein was retained per protein coding locus for *P. tricornutum*, *T. pseudonana*, and *S. robusta*. The longest protein model was used when more than one protein was assigned to the same locus tag.

The generic preparation script can be run as:

```bash
python scripts/21_prepare_reference_proteome.py \
    --proteins <NCBI_PROTEIN_FASTA> \
    --gff <MATCHING_GFF3> \
    --out-fasta <OUTPUT_FASTA> \
    --out-map <OUTPUT_MAP_TSV> \
    --id-mode locus
```

Final reference proteome counts were:

```text
Phaeodactylum tricornutum     10,392 proteins
Thalassiosira pseudonana      11,672 proteins
Seminavis robusta             35,995 proteins
```

#### 19.8.1 Diploid *Nitzschia inconspicua* protein set

The NCBI protein FASTA for *N. inconspicua* contained 38,785 proteins. GFF3 inspection showed that 150 proteins were encoded on plastid accession `MW971520.1` and 34 proteins were encoded on mitochondrial accession `MW971521.1`.

The organelle counts were verified with:

```bash
grep -v '^#' Nitzschia_inconspicua_GCA_019154785.2_genomic.gff3 \
| awk '$1=="MW971520.1" && $3=="CDS"' \
| sed -n 's/.*protein_id=\([^;]*\).*/\1/p' \
| sort -u \
| wc -l

grep -v '^#' Nitzschia_inconspicua_GCA_019154785.2_genomic.gff3 \
| awk '$1=="MW971521.1" && $3=="CDS"' \
| sed -n 's/.*protein_id=\([^;]*\).*/\1/p' \
| sort -u \
| wc -l
```

The observed counts were 150 plastid proteins and 34 mitochondrial proteins. After excluding these organelle contigs, 38,601 nuclear protein coding loci remained.

The generic reference preparation script can reproduce this nuclear set while preserving NCBI protein IDs:

```bash
python scripts/21_prepare_reference_proteome.py \
    --proteins Nitzschia_inconspicua_GCA_019154785.2_protein.faa \
    --gff Nitzschia_inconspicua_GCA_019154785.2_genomic.gff3 \
    --exclude-contigs MW971520.1,MW971521.1 \
    --id-mode protein \
    --out-fasta Nitzschia_inconspicua_nuclear.faa \
    --out-map Nitzschia_inconspicua_nuclear_map.tsv
```

Final *N. inconspicua* nuclear proteome:

```text
Proteins:        38,601
Total length:    19,032,900 aa
Minimum:                 45 aa
Mean:                 493.1 aa
Maximum:              9,933 aa
```

The diploid nuclear gene complement was retained. Sequence identity based deduplication was not used because it could collapse biological paralogs together with allelic copies.

### 19.9 Final five species OrthoFinder input

The final OrthoFinder directory contained exactly five FASTA files:

```text
DeerLake_Nitzschia.faa             14,979 proteins
Nitzschia_inconspicua.faa          38,601 proteins
Phaeodactylum_tricornutum.faa      10,392 proteins
Seminavis_robusta.faa              35,995 proteins
Thalassiosira_pseudonana.faa       11,672 proteins
```

FASTA identifiers were checked for uniqueness:

```bash
cd ${COMP}/orthofinder_input

for f in *.faa; do
    echo "=== $f ==="
    seqkit seq -n "$f" \
    | awk '{print $1}' \
    | sort \
    | uniq -d \
    | wc -l
done
```

No duplicate FASTA IDs were detected. Protein sequences were also checked for stop codons and unusual amino acid characters before analysis.

### 19.10 OrthoFinder installation

A dedicated Miniforge environment was used:

```bash
source ~/miniforge3/etc/profile.d/conda.sh

conda create -n orthofinder \
    --override-channels \
    --strict-channel-priority \
    -c conda-forge \
    -c bioconda \
    orthofinder \
    -y

conda activate orthofinder
```

Versions used:

```text
OrthoFinder 3.1.5
DIAMOND 2.2.6
FAMSA supplied in the OrthoFinder environment
FastTree supplied in the OrthoFinder environment
```

### 19.11 Five species OrthoFinder run

The following SLURM workflow was used:

```bash
#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=120:00:00
#SBATCH --mem=100G
#SBATCH --partition=cpu2025
####### Run your script #########################

set -euo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate /home/ruchita.solanki/miniforge3/envs/orthofinder

INPUT=/work/ebg_lab/eb/diatom_consortia/comparative_genomics/orthofinder_input
OUTPUT=/work/ebg_lab/eb/diatom_consortia/comparative_genomics/orthofinder_results/DL_5species_orthofinder_v3

NFILES=$(find "${INPUT}" -maxdepth 1 -type f -name "*.faa" | wc -l)

if [[ "${NFILES}" -ne 5 ]]; then
    echo "ERROR: Expected exactly 5 proteome FASTA files, found ${NFILES}"
    exit 1
fi

if [[ -e "${OUTPUT}" ]]; then
    echo "ERROR: Output directory already exists: ${OUTPUT}"
    exit 1
fi

orthofinder \
    -f "${INPUT}" \
    -t 32 \
    -a 32 \
    -S diamond \
    -M msa \
    -A famsa \
    -T fasttree \
    -o "${OUTPUT}"

ORTHOGROUP_FILE=$(find "${OUTPUT}" -type f -name "Orthogroups.tsv" | head -n 1)

if [[ -z "${ORTHOGROUP_FILE}" ]]; then
    echo "ERROR: Orthogroups.tsv was not produced."
    exit 1
fi
```

OrthoFinder completed with the following summary:

```text
Total input proteins:                        111,639
Genes assigned to orthogroups:               99,813   89.4%
Orthogroups:                                  16,157
Orthogroups containing all five species:       4,914
All species single copy orthogroups:             276
```

Because the *N. inconspicua* assembly is diploid, family size differences involving that species were not interpreted as gene expansion or contraction. The primary comparison used binary orthogroup sharing for each Deer Lake protein.

### 19.12 Final Deer Lake comparative gene table

The final output was designed to retain the same core fields used in the earlier PT and TP review table while replacing nucleotide similarity calls with OrthoFinder protein orthogroup sharing and adding *N. inconspicua* and *S. robusta*.

The final columns are:

```text
gene_id
contig_id
diatom_compartment
diatom_gene_length_bp
functional_annotation
diatom_Average_TPM
present_in_Thalassiosira_pseudonana
present_in_Phaeodactylum_tricornutum
present_in_Nitzschia_inconspicua
present_in_Seminavis_robusta
```

Run:

```bash
RESULTS=${COMP}/orthofinder_results/DL_5species_orthofinder_v3/Results_Sep09

python scripts/22_make_final_orthofinder_gene_table.py \
    --orthogroups ${RESULTS}/Orthogroups/Orthogroups.tsv \
    --dl-map ${COMP}/00_DL_reference/DL_diatom_FINAL_nuclear_representative_proteome_TEfiltered_map.tsv \
    --metadata ${BASE}/metatranscriptomics/transdecoder_to_braker_ID_bridge/CLEAN_REBUILD_FROM_RAW/09_final/DL_diatom_FINAL_gene_table_for_boss_PT_TP_NI.tsv \
    --out ${RESULTS}/Orthogroups/DL_diatom_FINAL_gene_table_OrthoFinder_PT_TP_NI_SR.tsv \
    --expected-genes 14979
```

Observed final counts:

```text
Final Deer Lake genes:                 14,979
Shared orthogroup with N. inconspicua: 10,163   67.85%
Shared orthogroup with S. robusta:      9,952   66.44%
Shared orthogroup with P. tricornutum:  7,973   53.23%
Shared orthogroup with T. pseudonana:   7,569   50.53%
```

The final file contains 14,979 data rows and one header row:

```text
Orthogroups/DL_diatom_FINAL_gene_table_OrthoFinder_PT_TP_NI_SR.tsv
```

In this table, `yes` indicates that the Deer Lake protein belongs to an OrthoFinder orthogroup containing at least one protein from the corresponding comparator. A `no` value means that no protein from that comparator was present in the Deer Lake protein's orthogroup, or that the Deer Lake protein was unassigned by OrthoFinder. These fields therefore represent orthogroup sharing, not definitive biological gene absence.

### 19.13 Final outputs

```text
comparative_genomics/
├── 00_DL_reference/
│   ├── DL_diatom_FINAL_nuclear_representative_proteome.faa
│   ├── DL_diatom_FINAL_nuclear_representative_proteome_map.tsv
│   ├── DL_diatom_FINAL_nuclear_representative_proteome_TEfiltered.faa
│   └── DL_diatom_FINAL_nuclear_representative_proteome_TEfiltered_map.tsv
├── 04_repeat_analysis/
│   ├── DL_RepeatModeler_consensi.fa.classified
│   ├── DL_RepeatMasker_summary.tbl
│   └── 03_CDS_repeat_overlap/
│       ├── DL_representative_proteins_CDS_repeat_overlap.tsv
│       ├── DL_CDS_repeat_overlap_bins.tsv
│       ├── DL_TE_candidate_classification.tsv
│       ├── DL_TIER1_TE_gene_roots.txt
│       └── TIER1_high_confidence_TE.tsv
├── orthofinder_input/
│   ├── DeerLake_Nitzschia.faa
│   ├── Nitzschia_inconspicua.faa
│   ├── Phaeodactylum_tricornutum.faa
│   ├── Seminavis_robusta.faa
│   └── Thalassiosira_pseudonana.faa
└── orthofinder_results/
    └── DL_5species_orthofinder_v3/
        └── Results_Sep09/
            ├── Orthogroups/
            ├── Orthologues/
            ├── Species_Tree/
            ├── Gene_Duplication_Events/
            └── Comparative_Genomics_Statistics/
```

The OrthoFinder analysis is retained as a secondary protein level comparison. The final four comparator nucleotide comparison used for the current manuscript analysis is documented in Section 21.

### 19.14 Plastid quality control and corrected OrthoFinder rerun

Inspection of photosynthesis related models identified residual plastid derived contigs in the nuclear enriched assembly. Four protein-bearing contigs were removed from the final Deer Lake nuclear query:

```text
contig_475
contig_4813
contig_5686
contig_5702
```

These contigs contained 38 retained representative proteins in total. After removal, the corrected Deer Lake protein set contained:

```text
14,941 proteins
```

Corrected protein FASTA:

```text
/work/ebg_lab/eb/diatom_consortia/comparative_genomics/00_DL_reference/
DL_diatom_FINAL_nuclear_representative_proteome_TEfiltered_plastidQC.faa
```

The corrected OrthoFinder rerun is located at:

```text
/work/ebg_lab/eb/diatom_consortia/comparative_genomics/orthofinder_results/
DL_5species_orthofinder_v3_plastidQC/Results_Sep10
```

Corrected Deer Lake summary:

```text
Input proteins:                    14,941
Assigned to orthogroups:          12,677  (84.8%)
Unassigned:                         2,264  (15.2%)
DL-containing orthogroups:          8,970
Shared with N. inconspicua:        10,154
Shared with S. robusta:             9,952
Shared with P. tricornutum:         7,966
Shared with T. pseudonana:          7,542
```

No mitochondrial based filtering was applied to the nuclear gene set because the candidate mitochondrial contigs were not considered reliable enough to justify gene removal.

</details>

---

<details>
<summary><strong>20. Hi-C read mapping and contig-level proximity-ligation network</strong> - BWA-MEM, samtools, awk, YaHS, and Python</summary>

Hi-C paired-end reads were incorporated after the main assembly, annotation, expression, and comparative-genomics workflow. The goal was to assess how broadly the polished whole assembly was represented in the proximity-ligation dataset, identify contigs connected by Hi-C read pairs, and separately test high-confidence diatom-bacterial read-pair contacts.

The Hi-C analysis was performed on the polished whole assembly rather than the nuclear-enriched subset because the proximity-ligation reads were generated from the complete diatom-associated consortium. This allowed diatom, bacterial, and mixed diatom-bacterial contacts to be evaluated in the same coordinate space.

### 20.1 Input files and working directories
```bash
cd /work/ebg_lab/eb/diatom_consortia

mkdir -p hi-c_diatoms/01_qc
mkdir -p hi-c_diatoms/02_map_to_whole_assembly
mkdir -p hi-c_diatoms/03_yahs_scaffolding
mkdir -p hi-c_diatoms/04_contact_maps
mkdir -p hic_bwa_separate_reads/{00_inputs,01_bwa_index,02_alignments,03_tables,04_logs}
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

Diatom draft genome used for contig-type classification in the separate-read analysis:
```text
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
```

Contigs present in `18_diatom.fasta` were treated as diatom contigs. All remaining contigs in the polished whole assembly were treated as bacterial for the purpose of the diatom-bacterial Hi-C read-pair screen.

### 20.2 Map Hi-C reads to the polished whole assembly as paired-end reads
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
### 20.3 Summarize contig-level Hi-C representation
```bash
awk 'BEGIN {
    OFS="	";
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
        printf "Total contigs: %d
", total;
        printf "Contigs with >=1 Hi-C read mapped: %d
", mapped;
        printf "Contigs with 0 Hi-C reads mapped: %d
", unmapped;
        printf "Percent contigs with Hi-C reads mapped: %.2f%%
", (mapped/total)*100;
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
### 20.4 Exploratory whole-assembly Hi-C scaffolding with YaHS
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
### 20.5 Extract all-primary inter-contig Hi-C contacts
The final contig-contact network used primary mapped Hi-C read pairs from the paired-end BWA-MEM mapping without applying a MAPQ cutoff. Unmapped reads, mate-unmapped reads, secondary alignments, and supplementary alignments were excluded. Each read pair was counted once if the two mates mapped to different contigs.
```bash
cd /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/02_map_to_whole_assembly
BAM=hic_to_whole_assembly.name_sorted.bam
samtools view -@ 8 -f 1 -F 2316 $BAM \
    | awk 'BEGIN{OFS="	"}
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
    | awk 'BEGIN{OFS="	"; print "contig_A","contig_B","HiC_contact_pairs"}
           {print $2,$3,$1}' \
    | sort -k3,3nr \
    > hic_intercontig_contacts_all_primary_pairs.tsv
```
This command extracts inter-contig Hi-C contacts from primary mapped read pairs and counts the number of read pairs supporting each contig-to-contig link.
```bash
awk 'BEGIN{OFS="	"}
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
### 20.6 Summarize connected contigs
```bash
awk 'NR>1 {print $1; print $2}' \
    hic_intercontig_contacts_all_primary_pairs.tsv \
    | sort -u \
    > hic_connected_contigs_all_primary_pairs.txt
```
This command lists every contig involved in at least one inter-contig Hi-C contact.
```bash
{
    echo -e "contig	number_of_connected_contigs	total_intercontig_HiC_pairs	strongest_single_contact_pairs"

    awk 'BEGIN{OFS="	"}
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
### 20.7 Convert the all-primary contact table to network files
The full contig-contact table was converted into GEXF and GraphML network files using a small helper Python script.
The script is saved as:
```text
scripts/14_make_hic_network_files.py
```
Run the script with:
```bash
cd /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/02_map_to_whole_assembly
python scripts/14_make_hic_network_files.py
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

### 20.8 High-confidence separate-read BWA mapping for diatom-bacterial contacts
A second BWA-MEM mapping was performed to keep the two Hi-C read files separate. This made it possible to ask, for each read ID, whether read 1 and read 2 mapped to different biological fractions of the whole assembly.

The separate-read analysis used a stricter read-level filter than the all-primary contact network:
```text
Primary alignments only
MAPQ >= 30
Percent identity >= 95%
```

Set up clean links to the whole assembly and separate Hi-C read files:
```bash
cd /work/ebg_lab/eb/diatom_consortia

mkdir -p hic_bwa_separate_reads/{00_inputs,01_bwa_index,02_alignments,03_tables,04_logs}

ln -sf /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_output/pypolca_corrected.fasta \
    hic_bwa_separate_reads/00_inputs/whole_assembly.fasta

ln -sf /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/1574499_S6_L001_R1_001.fastq.gz \
    hic_bwa_separate_reads/00_inputs/HiC_R1.fastq.gz

ln -sf /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/1574499_S6_L001_R2_001.fastq.gz \
    hic_bwa_separate_reads/00_inputs/HiC_R2.fastq.gz
```
Index the whole assembly for the separate-read analysis:
```bash
cd /work/ebg_lab/eb/diatom_consortia/hic_bwa_separate_reads

conda activate hic_diatom

bwa index \
    -p 01_bwa_index/whole_assembly \
    00_inputs/whole_assembly.fasta \
    2> 04_logs/bwa_index.log
```
The index files created were:
```text
whole_assembly.amb
whole_assembly.ann
whole_assembly.bwt
whole_assembly.pac
whole_assembly.sa
```

The following SLURM script maps read 1 and read 2 independently against the same BWA index and creates two sorted BAM files.
```bash
#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=120:00:00
#SBATCH --mem=100G
#SBATCH --partition=cpu2025
####### Run your script #########################

set -euo pipefail

cd /work/ebg_lab/eb/diatom_consortia/hic_bwa_separate_reads

conda activate hic_diatom

INDEX=01_bwa_index/whole_assembly

mkdir -p 02_alignments 04_logs

rm -f 02_alignments/HiC_R1.sorted.bam 02_alignments/HiC_R1.sorted.bam.bai
rm -f 02_alignments/HiC_R2.sorted.bam 02_alignments/HiC_R2.sorted.bam.bai

echo "Started separate Hi-C BWA mapping at: $(date)"
echo "BWA: $(which bwa)"
echo "samtools: $(which samtools)"

echo "Mapping R1 only..."
bwa mem -t ${THREADS} ${INDEX} 00_inputs/HiC_R1.fastq.gz 2> 04_logs/bwa_mem_R1.log | \
samtools sort -@ ${THREADS} -m2G -T 02_alignments/HiC_R1.tmp -o 02_alignments/HiC_R1.sorted.bam -

samtools index 02_alignments/HiC_R1.sorted.bam
samtools flagstat 02_alignments/HiC_R1.sorted.bam > 04_logs/HiC_R1.flagstat.txt

echo "Mapping R2 only..."
bwa mem -t ${THREADS} ${INDEX} 00_inputs/HiC_R2.fastq.gz 2> 04_logs/bwa_mem_R2.log | \
samtools sort -@ ${THREADS} -m2G -T 02_alignments/HiC_R2.tmp -o 02_alignments/HiC_R2.sorted.bam -

samtools index 02_alignments/HiC_R2.sorted.bam
samtools flagstat 02_alignments/HiC_R2.sorted.bam > 04_logs/HiC_R2.flagstat.txt

echo "Finished separate Hi-C BWA mapping at: $(date)"
ls -lh 02_alignments
```
Submit the separate-read mapping job:
```bash
cd /work/ebg_lab/eb/diatom_consortia/hic_bwa_separate_reads
sbatch 01_run_bwa_mem_separate_reads.slurm
```

Separate-read mapping outputs:
```text
02_alignments/HiC_R1.sorted.bam
02_alignments/HiC_R1.sorted.bam.bai
02_alignments/HiC_R2.sorted.bam
02_alignments/HiC_R2.sorted.bam.bai
04_logs/HiC_R1.flagstat.txt
04_logs/HiC_R2.flagstat.txt
```

Separate-read mapping summary:
```text
R1 primary reads:         445,405
R1 primary mapped reads:  318,884
R1 primary mapping rate:  71.59%

R2 primary reads:         445,405
R2 primary mapped reads:  304,083
R2 primary mapping rate:  68.27%
```

### 20.9 Create high-confidence read-pair tables and classify mixed diatom-bacterial pairs
The high-confidence separate-read tables were generated using custom Python scripts saved outside the markdown file.

The script used to parse the separate BAM files, calculate percent identity from the `NM` tag and aligned CIGAR length, and retain primary MAPQ >= 30 and percent identity >= 95 alignments is saved as:
```text
scripts/15_make_hic_primary_mapq30_pid95_tables.py
```
Run the script with:
```bash
cd /work/ebg_lab/eb/diatom_consortia/hic_bwa_separate_reads
conda activate hic_diatom
python scripts/15_make_hic_primary_mapq30_pid95_tables.py
```
This script generated:
```text
03_tables/HiC_R1.primary_MAPQ30_PID95.tsv
03_tables/HiC_R2.primary_MAPQ30_PID95.tsv
```

The script used to join read 1 and read 2 by read ID, classify contigs as diatom or bacterial, assign pair-type codes, and write full read-pair tables is saved as:
```text
scripts/16_make_hic_pair_type_tables.py
```
Run the script with:
```bash
cd /work/ebg_lab/eb/diatom_consortia/hic_bwa_separate_reads
conda activate hic_diatom
python scripts/16_make_hic_pair_type_tables.py
```
This script used:
```text
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_output/pypolca_corrected.fasta
```
Contigs present in the diatom draft genome were labelled `diatom`, and all remaining whole-assembly contigs were labelled `bacterial`.

The pair-type code system was:
```text
1 = both reads mapped to diatom contigs
2 = both reads mapped to bacterial contigs
3 = read 1 mapped to a diatom contig and read 2 mapped to a bacterial contig
4 = read 1 mapped to a bacterial contig and read 2 mapped to a diatom contig
```

High-confidence pair-type summary:
```text
both diatom:                         39,739 pairs = 64.82%
both bacterial:                      21,224 pairs = 34.62%
read 1 diatom, read 2 bacterial:        164 pairs = 0.27%
read 1 bacterial, read 2 diatom:        177 pairs = 0.29%
```
The total number of high-confidence mixed diatom-bacterial Hi-C read pairs was:
```text
341 mixed read pairs
```

The final simplified mixed-pair table was generated using:
```text
scripts/17_make_hic_simple_mixed_read_table.py
```
Run the script with:
```bash
cd /work/ebg_lab/eb/diatom_consortia/hic_bwa_separate_reads
conda activate hic_diatom
python scripts/17_make_hic_simple_mixed_read_table.py
```
The final simple table is:
```text
03_tables/HiC_DIATOM_BACTERIAL_read_level_SIMPLE_COLUMNS_MAPQ30_PID95.tsv
```

Final simple table columns:
```text
read_id
read_number
contig_id
paired_contig_id
mapq
percent_identity
aligned_length_bp
pair_type_code
pair_type
```
Each mixed Hi-C pair is represented by two rows, one for read 1 and one for read 2. The row count check confirmed:
```text
341 unique mixed Hi-C read pairs
682 read rows
683 lines including header
```
### 20.10 Organize final Hi-C outputs
After generating the final mapping, contact-network, and separate-read mixed-contact outputs, files were organized into final mapping and contact-map folders.
```bash
cd /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms
mkdir -p 02_map_to_whole_assembly/final_mapping
mkdir -p 04_contact_maps/tables
mkdir -p 04_contact_maps/network_files
mkdir -p 04_contact_maps/scripts
```
Final paired-end mapping files:
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
Final all-primary contact-network files:
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
    └── 14_make_hic_network_files.py
```
Final high-confidence separate-read mixed-contact files:
```text
hic_bwa_separate_reads/
├── 02_alignments/
│   ├── HiC_R1.sorted.bam
│   ├── HiC_R1.sorted.bam.bai
│   ├── HiC_R2.sorted.bam
│   └── HiC_R2.sorted.bam.bai
├── 03_tables/
│   ├── whole_assembly_contig_type_map.tsv
│   ├── HiC_R1.primary_MAPQ30_PID95.tsv
│   ├── HiC_R2.primary_MAPQ30_PID95.tsv
│   ├── HiC_read_pairs_MAPQ30_PID95_joined.tsv
│   ├── HiC_read_pairs_MAPQ30_PID95_with_contig_types.tsv
│   ├── HiC_pair_type_summary_MAPQ30_PID95.tsv
│   ├── HiC_DIATOM_BACTERIAL_read_level_MAPQ30_PID95.tsv
│   └── HiC_DIATOM_BACTERIAL_read_level_SIMPLE_COLUMNS_MAPQ30_PID95.tsv
└── 04_logs/
    ├── HiC_R1.flagstat.txt
    ├── HiC_R2.flagstat.txt
    ├── bwa_mem_R1.log
    └── bwa_mem_R2.log
```
### 20.11 Final Hi-C analysis summary
Hi-C reads mapped to 4,010 of 4,925 contigs in the polished whole assembly, corresponding to 81.42% of assembly contigs. At the read level, 622,967 of 890,810 primary reads mapped to the assembly, corresponding to a primary mapping rate of 69.93%.

Inter-contig proximity-ligation contacts were extracted from primary mapped Hi-C read pairs without applying a MAPQ cutoff. The final all-primary contig-contact network contained 3,770 contig nodes and 75,703 inter-contig Hi-C links.

A second high-confidence separate-read analysis was then used to identify mixed diatom-bacterial Hi-C read pairs. Read 1 and read 2 were mapped independently to the polished whole assembly, filtered for primary MAPQ >= 30 and percent identity >= 95 alignments, joined by read ID, and classified using the diatom draft genome as the diatom contig reference. This produced 341 high-confidence mixed diatom-bacterial Hi-C read pairs, represented as 682 read-level rows in the final simplified table.

</details>

---

<details>
<summary><strong>21. Final four comparator BLASTN redo with repeat, plastid QC, annotation, and Average_TPM</strong> - BLAST+, Python, pandas, and SLURM</summary>

This section supersedes the earlier PT only and TP only pairwise nucleotide screens for the current manuscript analysis. The goal was to apply the same nucleotide search procedure to all four reference diatom genomes while retaining every curated Deer Lake query gene in a single master table.

### 21.1 Final Deer Lake nucleotide query

The original clean Deer Lake nucleotide gene FASTA contained:

```text
15,102 BRAKER4 gene sequences
```

Input nucleotide gene FASTA:

```text
/work/ebg_lab/eb/diatom_consortia/thalassiosira_to_diatom_blastn_redo/00_inputs/
diatom_genes.clean.fasta
```

The final query was restricted to gene roots represented in the repeat and plastid quality controlled Deer Lake protein set:

```text
/work/ebg_lab/eb/diatom_consortia/comparative_genomics/00_DL_reference/
DL_diatom_FINAL_nuclear_representative_proteome_TEfiltered_plastidQC.faa
```

Run:

```bash
cd /work/ebg_lab/eb/diatom_consortia/comparative_genomics/blastn_redo

python scripts/23_prepare_DL_BLASTN_query.py \
    --gene-fasta /work/ebg_lab/eb/diatom_consortia/thalassiosira_to_diatom_blastn_redo/00_inputs/diatom_genes.clean.fasta \
    --final-protein-fasta /work/ebg_lab/eb/diatom_consortia/comparative_genomics/00_DL_reference/DL_diatom_FINAL_nuclear_representative_proteome_TEfiltered_plastidQC.faa \
    --out-fasta 00_inputs/DL_final_nuclear_genes_repeat_TE_plastid_clean.fasta \
    --roots-out 00_inputs/DL_final_gene_roots.txt \
    --expected 14941
```

Final query:

```text
00_inputs/DL_final_nuclear_genes_repeat_TE_plastid_clean.fasta
```

Final number of Deer Lake query genes:

```text
14,941
```

The query is repeat/TE filtered and plastid quality controlled. It should not be described as mitochondrial cleaned.

### 21.2 Comparator genome FASTA files

```text
Nitzschia inconspicua
/home/ruchita.solanki/databases/nitzschia_inconspicua/
GCA_019154785.2_GAI293_CANU_175m_combined/
Nitzschia_inconspicua_GCA_019154785.2_genomic.fna

Seminavis robusta
/home/ruchita.solanki/databases/seminavis_robusta/D6/
Seminavis_robusta_GCA_903772945.1_genomic.fna

Phaeodactylum tricornutum
/work/ebg_lab/eb/diatom_consortia/phaeodactylum_to_diatom_blastn_redo/00_inputs/
phaeodactylum_genome.fna

Thalassiosira pseudonana
/home/ruchita.solanki/thalassiosira_pseudonana/
GCF_000149405.2_ASM14940v2/
Thalassiosira_pseudonana_ASM14940v2_GCF_000149405.2_genomic.fna
```

### 21.3 BLAST database construction and four comparator search

The reusable SLURM script is saved as:

```text
slurm/23_run_four_genome_dcmegablast_SLURM.txt
```

The search settings were identical for all four comparators:

```text
BLAST task:        dc-megablast
E value:           1e-10
Threads:           32
Additional fixed identity cutoff: none
Additional fixed coverage cutoff: none
```

BLAST output format:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs
```

Raw output files:

```text
02_blastn_raw/DL_vs_NI.dcmegablast.tsv
02_blastn_raw/DL_vs_SR.dcmegablast.tsv
02_blastn_raw/DL_vs_PT.dcmegablast.tsv
02_blastn_raw/DL_vs_TP.dcmegablast.tsv
```

Observed raw alignment counts:

```text
NI: 27,211
SR: 12,168
PT:  7,915
TP:  4,817
Total: 52,111
```

These values are alignment counts, not unique Deer Lake gene counts.

### 21.4 Master gene table

The metadata map used for Deer Lake annotation and expression was:

```text
/work/ebg_lab/eb/diatom_consortia/comparative_genomics/00_DL_reference/
DL_diatom_FINAL_nuclear_representative_proteome_TEfiltered_map.tsv
```

The map contains gene root, representative transcript ID, contig ID, compartment, protein length, Average_TPM, and functional annotation. Only roots present in the final 14,941 gene query are retained by the master-table script.

Run:

```bash
python scripts/24_build_BLASTN_master_table.py \
    --query-fasta 00_inputs/DL_final_nuclear_genes_repeat_TE_plastid_clean.fasta \
    --metadata-map /work/ebg_lab/eb/diatom_consortia/comparative_genomics/00_DL_reference/DL_diatom_FINAL_nuclear_representative_proteome_TEfiltered_map.tsv \
    --ni 02_blastn_raw/DL_vs_NI.dcmegablast.tsv \
    --sr 02_blastn_raw/DL_vs_SR.dcmegablast.tsv \
    --pt 02_blastn_raw/DL_vs_PT.dcmegablast.tsv \
    --tp 02_blastn_raw/DL_vs_TP.dcmegablast.tsv \
    --out DL_14941_genes_BLASTN_NI_SR_PT_TP_with_TPM_compartment.tsv \
    --expected 14941
```

For each comparator, one best alignment is retained per query gene by highest bitscore, then lower E value, higher query coverage, higher percent identity, and longer alignment as deterministic tie breakers.

The final master table retains all 14,941 query genes, including genes without a BLASTN hit.

Core metadata columns:

```text
gene_root
gene_id
contig_id
diatom_compartment
gene_length_bp
functional_annotation
diatom_Average_TPM
```

Each comparator contributes:

```text
*_hit
*_subject
*_pident
*_alignment_length
*_qcov
*_evalue
*_bitscore
```

A `yes` hit means at least one `dc-megablast` alignment was returned at E <= 1e-10. No additional identity or query coverage threshold is used to delete genes from the master table.

Observed unique Deer Lake gene counts:

```text
NI hit:                         5,767
SR hit:                         4,273
PT hit:                         3,799
TP hit:                         2,636
Hit in at least one comparator: 6,429
No hit in any comparator:       8,512
```

### 21.5 Expression values and unmatched subsets

Average_TPM availability across the 14,941 gene table:

```text
TPM value available: 11,092
TPM > 0:             10,907
TPM = 0:                185
```

Generate the main BLASTN unmatched subsets with:

```bash
python scripts/25_make_BLASTN_subsets.py \
    --master DL_14941_genes_BLASTN_NI_SR_PT_TP_with_TPM_compartment.tsv \
    --query-fasta 00_inputs/DL_final_nuclear_genes_repeat_TE_plastid_clean.fasta \
    --out-prefix DL
```

Expected key outputs:

```text
DL_no_comparator_hit.tsv
DL_no_comparator_hit_gene_roots.txt
DL_no_comparator_hit_genes.fasta
DL_no_comparator_hit_expressed_sorted_by_TPM.tsv
DL_unmatched_expressed_annotated_sorted_TPM.tsv
```

Observed counts:

```text
No comparator hit:                                8,512 genes
No comparator hit + TPM > 0:                     5,332 genes
No comparator hit + TPM > 0 + informative annotation: 3,259 genes
```

No gene length cutoff is applied to these final subsets. Gene length remains a descriptive field only.

Annotation labels that contain taxonomic or organism names should be interpreted as reference or domain based annotation labels unless independent evidence supports that biological identity.

### 21.6 Panel A expression plot

The plotting script is:

```text
scripts/26_plot_BLASTN_TPM_panelA.py
```

Laptop dependencies:

```bash
pip install pandas numpy matplotlib scipy
```

Run:

```bash
python scripts/26_plot_BLASTN_TPM_panelA.py \
    --input DL_14941_genes_BLASTN_NI_SR_PT_TP_with_TPM_compartment.tsv \
    --out-prefix PanelA_BLASTN_TPM
```

The five plotted groups are:

```text
Deer Lake only
Shared with N. inconspicua
Shared with S. robusta
Shared with P. tricornutum
Shared with T. pseudonana
```

`Deer Lake only` is a compact figure label for genes with no detectable BLASTN hit to any of the four comparator genomes under this search. It should not be interpreted as proof of lineage specificity.

The shared comparator groups overlap because a Deer Lake gene may have a BLASTN hit to multiple comparator genomes.

Expression is displayed as:

```text
log10(Average TPM + 1)
```

True TPM values of zero are retained. Blank TPM fields are excluded because they represent missing expression values rather than zero expression.

### 21.7 Interpretation

Preferred manuscript language:

> Genes classified as Deer Lake only had no detectable nucleotide similarity to *N. inconspicua*, *S. robusta*, *P. tricornutum*, or *T. pseudonana* under the dc-megablast search criteria used.

Avoid calling these genes definitively unique, lineage specific, or biologically absent from the four references. Nucleotide search sensitivity, genome assembly quality, annotation completeness, sequence divergence, and query length can affect detection.

</details>

---

<details>
<summary><strong>22. High confidence whole assembly Hi C network visualization</strong> - pandas, NetworkX, NumPy, and Matplotlib</summary>

The final network visualization uses the high confidence separate-read Hi C tables generated in Section 20. The purpose of this figure is to show all contigs in the polished consortium assembly while distinguishing contigs that participate in retained inter-contig Hi C links from contigs that do not.

### 22.1 Input tables

Working directory:

```text
/work/ebg_lab/eb/diatom_consortia/hic_bwa_separate_reads/03_tables
```

Input contact table:

```text
HiC_contig_pair_contacts_MAPQ30_PID95.tsv
```

Columns:

```text
forward_contig_id
reverse_contig_id
forward_contig_type
reverse_contig_type
pair_type_code
pair_type
read_pair_count
```

Whole assembly contig type map:

```text
whole_assembly_contig_type_map.tsv
```

The type map contains 4,925 assembly contigs plus one header line.

### 22.2 Build the undirected high confidence edge table

The original contact table contained self contacts and orientation-specific forward/reverse contig pairs. For the plotting table:

1. self contacts were removed;
2. orientation-specific rows with `read_pair_count >= 2` were retained;
3. each contig pair was converted to a canonical undirected pair; and
4. reciprocal entries were summed.

Run:

```bash
python scripts/27_make_HiC_undirected_min2_edges.py \
    --input HiC_contig_pair_contacts_MAPQ30_PID95.tsv \
    --output HiC_network_edges_min2_undirected.tsv \
    --min-support 2
```

Observed filtering summary:

```text
Non-self contig-pair rows:                    5,146
Non-self rows with >=2 read pairs:              576
Unique undirected edges after collapsing:        471
Unique contigs in retained edges:                487
Minimum retained edge weight:                      2
Maximum retained edge weight:                    194
Median retained edge weight:                       2
```

The threshold is applied before reciprocal orientation collapsing to reproduce the final analysis exactly.

### 22.3 Plot all assembly contigs

The plotting script is:

```text
scripts/28_plot_HiC_whole_assembly_network.py
```

Laptop dependencies:

```bash
pip install pandas numpy networkx matplotlib
```

Run:

```bash
python scripts/28_plot_HiC_whole_assembly_network.py \
    --edges HiC_network_edges_min2_undirected.tsv \
    --types whole_assembly_contig_type_map.tsv \
    --out-prefix HiC_whole_assembly_network
```

The plot contains all 4,925 assembly contigs:

```text
Connected contigs:      487
Isolated contigs:     4,438
Total:                4,925
```

The spring layout is calculated only for the connected contigs. Isolated contigs are placed in concentric rings around the connected network so that their presence is visible without implying spatial relationships among them.

Node colors use a color blind aware palette:

```text
Diatom contigs:       blue
Bacterial contigs:    orange
```

Mixed diatom-bacterial edges are drawn more prominently than within-group edges, but a Hi C edge is interpreted only as proximity ligation support under the filtering criteria used. It is not evidence by itself for mutualism or a direct ecological interaction.

### 22.4 Organelle associated contigs highlighted in the figure

Plastid genome associated contigs:

```text
contig_1443
contig_4315
```

Candidate mitochondrial contigs:

```text
contig_5628
contig_1647
```

The plotting script highlights these four contigs with larger star symbols. The mitochondrial contigs are labelled as candidates because their genome assignment remains uncertain and was not used to filter the final nuclear gene set.

The absence of a retained edge between an organelle contig and the nuclear network does not imply that the organelle never contacts the nucleus. Nuclear, plastid, and mitochondrial DNA are physically separate molecules, and the final network additionally excludes self contacts and weak inter-contig contacts below the selected support threshold.

### 22.5 Inspect raw organelle contact patterns

Before interpreting organelle placement, inspect all raw contact rows involving the highlighted contigs:

```bash
awk -F'\t' '
NR==1 ||
$1=="contig_1443" || $2=="contig_1443" ||
$1=="contig_4315" || $2=="contig_4315" ||
$1=="contig_5628" || $2=="contig_5628" ||
$1=="contig_1647" || $2=="contig_1647"
' HiC_contig_pair_contacts_MAPQ30_PID95.tsv | column -t
```

This allows self contacts, weak contacts removed by the network threshold, diatom contacts, and bacterial contacts to be evaluated separately.

</details>

---

<details>
<summary><strong>23. Comparative metatranscriptomics of expressed Deer Lake ORFs</strong> - nf-core/metatdenovo outputs, BLASTP reciprocal best hits, expression percentiles, KOfam, EggNOG, Pfam/HMMER, and manual functional curation</summary>

This analysis was added after the genome level BLASTN and OrthoFinder comparisons. It addresses a different question: which functional systems are strongly expressed in the Deer Lake diatom metatranscriptome, and how broadly are putative homologs recovered across four reference diatom proteomes?

The analysis starts from the **Deer Lake diatom metatranscriptome generated with nf-core/metatdenovo** in Section 10. It does not reassemble the transcriptome and does not re-map or re-quantify the RNA-seq reads.

The working directory was:

```text
/work/ebg_lab/eb/diatom_consortia/comparative_transcriptomics
```

The main Conda environment used for BLAST and table processing was:

```text
diatom_blast
```

The four reference diatoms were:

```text
NI = Nitzschia inconspicua GAI 293
SR = Seminavis robusta D6
PT = Phaeodactylum tricornutum Phatr3
TP = Thalassiosira pseudonana CCMP1335
```

### 23.1 Deer Lake TransDecoder ORFs and ORF level expression

The ORF level expression table produced from the metatranscriptome analysis was:

```text
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/new_results/
summary_tables/spades.transdecoder.counts.tsv.gz
```

The table contains ORF level counts and TPM values with the columns:

```text
orf  chr  start  end  strand  length  sample  count  tpm
```

The expression table uses identifiers such as `cds.NODE_1.p1`, whereas the TransDecoder peptide and annotation tables use `NODE_1.p1`. The `cds.` prefix was removed during ID matching.

The original TransDecoder peptide set contained **88,924 ORFs**. A total of **87,621 ORFs** had matching TPM values and were retained. A total of 1,303 peptide ORFs did not have a matching TPM value and were not included in this analysis.

The preparation script is:

```text
scripts/29_prepare_DL_expressed_ORFs.py
```

No Salmon re-quantification or new RNA read mapping was performed for this comparison.

### 23.2 Reference proteomes

The four reference proteomes were:

```text
NI: /home/ruchita.solanki/NI_GAI293_NCBI/ncbi_dataset/data/GCA_019154785.2/protein.faa
    38,785 proteins

SR: /home/ruchita.solanki/SR_D6_NCBI/ncbi_dataset/data/GCA_903772945.1/protein.faa
    35,995 proteins

PT: /home/ruchita.solanki/trascriptome/PT_Phatr3_proteins.fasta
    12,178 proteins

TP: /home/ruchita.solanki/TP_CCMP1335_NCBI/ncbi_dataset/data/GCF_000149405.2/protein.faa
    11,673 proteins
```

Project links:

```text
00_inputs/NI_proteins.faa
00_inputs/SR_proteins.faa
00_inputs/PT_proteins.faa
00_inputs/TP_proteins.faa
```

Protein to transcript and gene maps:

```text
03_tables/NI_protein_transcript_gene_map.tsv
03_tables/SR_protein_transcript_gene_map.tsv
03_tables/PT_protein_transcript_gene_map.tsv
03_tables/TP_protein_transcript_gene_map.tsv
```

### 23.3 Forward and reverse BLASTP

The 87,621 expressed Deer Lake ORFs were searched against each reference proteome. The reference proteomes were then searched back against a BLAST database constructed from the expressed Deer Lake ORFs.

SLURM files:

```text
slurm/24_run_DL_vs_reference_proteomes_BLASTP_SLURM.txt
slurm/25_run_reference_proteomes_vs_DL_BLASTP_SLURM.txt
```

BLASTP settings:

```text
E value threshold:    1e-5
max target sequences: 10
SEG filtering:        yes
```

Unique Deer Lake ORFs with at least one forward BLASTP hit:

| Comparator | Deer Lake ORFs |
| --- | ---: |
| NI | 42,369 |
| SR | 42,156 |
| PT | 39,662 |
| TP | 36,231 |

Unique reference proteins with at least one reverse Deer Lake hit:

| Comparator | Reference proteins |
| --- | ---: |
| NI | 27,492 |
| SR | 21,443 |
| PT | 9,740 |
| TP | 8,779 |

These forward and reverse hit counts describe detectable protein similarity only.

### 23.4 Reciprocal best hit extraction

Reciprocal best hits were extracted with:

```text
scripts/30_extract_unique_top_RBH.py
```

For each direction, only a unique maximum bitscore hit was retained. Top bitscore ties to different subjects were treated as ambiguous. An RBH was called only when the unique top relationship was reciprocal.

No additional arbitrary identity or coverage threshold was imposed after the BLASTP search. Percent identity, Deer Lake query coverage, E value, and bitscore were retained for inspection.

Final RBH counts:

| Comparator | RBH pairs |
| --- | ---: |
| NI | 7,110 |
| SR | 8,899 |
| PT | 7,512 |
| TP | 6,121 |

Output files:

```text
03_tables/DL_NI_RBH.tsv
03_tables/DL_SR_RBH.tsv
03_tables/DL_PT_RBH.tsv
03_tables/DL_TP_RBH.tsv
```

RBH is used as a conservative putative homology indicator, not as proof of one to one orthology.

### 23.5 Five species expression and RBH master

The four RBH tables were merged with Deer Lake Average_TPM and reference protein metadata using:

```text
scripts/31_build_5species_RBH_expression_master.py
```

Main output:

```text
03_tables/DL_5species_RBH_expression_master.tsv
```

The RBH pattern is stored in NI, SR, PT, TP order.

Observed pattern counts:

| Pattern | ORFs |
| --- | ---: |
| 0000 | 75,168 |
| 1111 | 2,979 |
| 1000 | 1,306 |
| 0100 | 1,270 |
| 0111 | 1,173 |
| 1110 | 1,051 |
| 0010 | 864 |
| 0110 | 858 |
| 1100 | 743 |
| 0001 | 564 |
| 1101 | 434 |
| 0101 | 391 |
| 1010 | 240 |
| 1001 | 233 |
| 0011 | 223 |
| 1011 | 124 |

```text
RBH in at least one reference: 12,453 ORFs
RBH in all four references:     2,979 ORFs
```

The strict five species core table is:

```text
03_tables/DL_5species_RBH_core1111.tsv
```

A `0000` pattern means that no reciprocal best hit was detected among the four reference proteomes under this analysis. It does not establish that an ORF is unique to Deer Lake.

### 23.6 Within Deer Lake expression percentile

Average_TPM values were ranked across all 87,621 quantified Deer Lake ORFs with:

```text
scripts/32_add_DL_expression_percentile.py
```

Output:

```text
03_tables/DL_5species_RBH_expression_percentile.tsv
```

Percentile definition:

```text
100 x number of ORFs with TPM <= focal ORF TPM / 87,621
```

Example:

```text
NODE_1.p1
DL_Average_TPM = 773.335845
DL_expression_percentile = 99.87
```

The percentile is a within Deer Lake metric and is not a cross species expression comparison.

### 23.7 Functional annotation integration

Three available transcriptome annotation layers were integrated by shared TransDecoder ORF ID.

EggNOG:

```text
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/new_results/
summary_tables/spades.transdecoder.emapper.tsv.gz
```

KOfam:

```text
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/new_results/
summary_tables/spades.transdecoder.kofamscan-uniq.tsv.gz
```

Direct Pfam/HMMER:

```text
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/new_results/
hmmer/spades.transdecoder.Pfam-A.tbl.gz
```

Integration script:

```text
scripts/33_integrate_transcriptome_annotations.py
```

Direct Pfam processing retained only HMMER hits with an included domain count greater than zero.

Observed annotation counts:

```text
EggNOG annotated ORFs:               67,401
KOfam annotated ORFs:                73,877
ORFs with included Pfam hits:        60,667
At least one annotation source:      79,176
Total Deer Lake ORFs in the master:  87,621
```

Combined output:

```text
03_tables/DL_expression_RBH_integrated_annotations.tsv
```

The MMETSP DIAMOND output within the Eukulele workflow was inspected but not used for functional integration because it contained a standard 12 column alignment table with CAMPEP subject identifiers and no functional descriptions.

### 23.8 Broad six system candidate screen

The integrated table was screened for candidates in:

```text
Photosynthesis
Carbon fixation and CCM
Silica metabolism
Ion homeostasis and osmoregulation
Urea and nitrogen metabolism
Oxidative stress
```

Script:

```text
scripts/34_extract_six_system_candidates.py
```

Output:

```text
03_tables/DL_functional_candidates_6systems.tsv
```

Initial broad candidate counts:

| Functional system | Candidate ORFs |
| --- | ---: |
| Photosynthesis | 854 |
| Carbon fixation and CCM | 254 |
| Silica metabolism | 20 |
| Ion homeostasis and osmoregulation | 3,848 |
| Urea and nitrogen metabolism | 832 |
| Oxidative stress | 1,141 |

These counts are not pathway abundance estimates. Broad keywords intentionally recovered false positives and were used only to collect candidates for manual review.

### 23.9 Manual functional review

The broad search was reduced using explicit KOfam, EggNOG, and direct Pfam evidence.

The curated table is generated by:

```text
scripts/35_make_curated_6systems_expression_conservation.py
```

Output:

```text
03_tables/DL_curated_6systems_expression_conservation.tsv
```

The current curated table contains six representative genes per system, for a total of 36 genes.

For silica metabolism, 18 inspected candidates carried both an EggNOG `Silicon transporter` description and a direct `PF03842.19|Silic_transp` Pfam hit. Because some KOfam calls were discordant, these proteins were retained conservatively as **putative silicon transporter family proteins** rather than being assigned SIT1, SIT2, or SIT3 subtype names.

For ion homeostasis, the broad screen was reduced to defined transport machinery such as V type H+ ATPases, chloride channels, potassium channels, SLC9 sodium/hydrogen exchangers, and Na+/H+ antiporters. F type ATP synthase proteins, calcium dependent kinases, and generic sodium coupled nutrient transporters were not treated as core ion homeostasis representatives simply because their descriptions contained ion related terms.

Strong Na+/H+ exchange candidates were supported by SLC9 or Na+/H+ antiporter annotations and, for selected proteins, direct `PF00999` (`Na_H_Exchanger`) support.

### 23.10 Compact 24 gene table

The first compact figure table contains four representatives per system:

```text
scripts/36_make_final_24genes_expression_conservation.py
```

Output:

```text
03_tables/DL_final_24genes_expression_conservation.tsv
```

The compact table retains function, ORF ID, Average_TPM, expression percentile, NI/SR/PT/TP RBH calls, RBH pattern, RBH count, and selected annotation evidence.

Obvious plastid encoded genes are not interpreted from the nuclear reference proteome RBH matrix. Confirmed plastid encoded representatives are marked `NA` rather than `0` in the compact table to avoid implying biological absence from the reference species.

### 23.11 Proposed 30 gene extension

Five genes per functional system were selected as the next figure iteration.

Prepared script:

```text
scripts/37_make_final_30genes_expression_conservation.py
```

Prepared figure script:

```text
scripts/38_plot_DL_expression_RBH_matrix.py
```

At the stopping point documented here, the 30 gene table and final 30 gene figure had **not yet been executed**. The RuBisCO small subunit candidate should have its organelle versus nuclear origin checked before its RBH assessment is finalized.

### 23.12 Exploratory transcript BLASTN

Before the ORF level RBH analysis was adopted, 87,325 assembled Deer Lake transcripts were compared with NI, SR, PT, and TP transcript sets using `dc-megablast`.

SLURM:

```text
slurm/26_run_exploratory_transcript_dcmegablast_SLURM.txt
```

Unique Deer Lake transcripts with at least one hit:

| Comparator | Deer Lake transcripts |
| --- | ---: |
| NI | 12,901 |
| SR | 11,674 |
| PT | 11,360 |
| TP | 8,322 |

Summary:

```text
Total Deer Lake transcripts:            87,325
Hit in at least one comparator:          16,614
No detected hit in the four references: 70,711
```

This was retained as an exploratory nucleotide similarity analysis rather than as the final homology metric.

### 23.13 Public cross species expression was not forced

Public *S. robusta* Salmon output was inspected as a possible cross species expression comparison. The public transcript identifiers used `Sro...` IDs that did not directly map to the NCBI protein/gene identifiers used in the reference proteome, and no simple plain text crosswalk was available from the inspected repository files.

Because independently generated expression studies also differ in experimental design, raw TPM values were not compared directly across species.

The final strategy therefore uses:

```text
Deer Lake Average_TPM
Deer Lake within transcriptome expression percentile
+
NI/SR/PT/TP protein level RBH status
```

### 23.14 Interpretation rules

The analysis can support wording such as:

> Highly expressed Deer Lake ORFs included both broadly conserved proteins with reciprocal best hits across multiple reference diatoms and proteins for which reciprocal best hits were recovered in few or none of the four reference proteomes.

It does not by itself support claims of Deer Lake specificity, definitive biological absence from another species, cross species differential expression, induction, or confirmed one to one orthology.

Preferred wording for `0000`:

> No reciprocal best hit was detected among the four reference proteomes under this analysis.

</details>
