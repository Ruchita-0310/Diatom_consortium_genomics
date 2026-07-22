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
Conda environments, Singularity containers, and local HPC modules were used depending on software availability.

| Step                              | Tools                                                                                                                                     |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| Assembly and polishing              | Flye, Medaka, BWA-MEM, Polypolish, Pypolca                                                                                                |
| Read mapping and coverage           | minimap2, samtools, bedtools, seqkit                                                                                                     |
| Assembly quality                    | BUSCO, QUAST/MetaQUAST                                                                                                                    |
| Binning and bin quality             | MetaBAT2, CheckM2                                                                                                                         |
| Taxonomy and abundance              | GTDB-Tk, MetaEuk, CoverM                                                                                                                  |
| Organelle identification            | MetaQUAST, minimap2, bedtools, seqkit, GeSeq/OGDRAW                                                                                       |
| Phylogenetics                       | Clustal Omega, TrimAl, IQ-TREE 2                                                                                                          |
| Transcriptomics                     | Nextflow, nf-core/metatdenovo, TransDecoder, Barrnap, STAR                                                                                |
| Genome annotation                   | BRAKER4, GeneMark-ET, AUGUSTUS, TSEBRA, STAR, BUSCO/compleasm                                                                             |
| Functional annotation               | DIAMOND, UniProtKB/Swiss-Prot, UniProtKB Bacillariophyta, InterProScan, Pfam, PANTHER, Gene3D, CDD, SMART, SUPERFAMILY, ProSite, Python  |
| Expression integration              | DIAMOND, Python, pandas, TransDecoder ORFs, Average_TPM table                                                                             |
| Comparative genomics                 | NCBI Datasets, BLASTN, bedtools, Python                                                                                                   |
| Hi-C mapping and contact network     | FastQC, MultiQC, BWA-MEM, samtools, seqkit, YaHS, awk, Python                                                                             |

---

## Repository structure for scripts
Custom Python scripts are stored in `scripts/`, numbered in the order they are used.

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
├── 09_add_Average_TPM_to_BRAKER_isoforms.py
├── 10_make_clean_BRAKER_isoform_table.py
├── 11_make_boss_review_gene_table_PTredo.py
├── 12_make_hic_network_files.py
├── 13_make_hic_primary_mapq30_pid95_tables.py
├── 14_make_hic_pair_type_tables.py
└── 15_make_hic_simple_mixed_read_table.py
```

| Script | Purpose |
|---|---|
| `01_classify_metaeuk_contigs.py` | Classifies contigs using MetaEuk ORF-level taxonomy and assigns each contig to a final category. |
| `02_make_swissprot_best_hits.py` | Parses Swiss-Prot DIAMOND output, calculates coverage, assigns confidence classes, and writes best-hit annotation tables. |
| `03_make_bacillariophyta_best_hits.py` | Parses UniProtKB Bacillariophyta DIAMOND output using the same confidence framework as the Swiss-Prot parser. |
| `04_summarize_interproscan.py` | Collapses raw InterProScan TSV output to one row per predicted protein, summarizing domains, GO terms, pathways, and database sources. |
| `05_merge_functional_annotation_layers.py` | Merges Swiss-Prot, Bacillariophyta, InterProScan, and AntiFam evidence into one BRAKER4 functional annotation table. |
| `06_merge_phaeodactylum_blast_hits.py` | Adds overlapping *Phaeodactylum* gene information to the cleaned BLASTN comparison table. |
| `07_add_BRAKER_lengths_clean.py` | Adds BRAKER4 coordinates, contig ID, strand, gene length, CDS length, and protein length to each isoform. |
| `08_make_best_ORF_to_BRAKER_mapping_clean.py` | Parses TransDecoder ORF vs. BRAKER4 DIAMOND output and keeps one best BRAKER4 hit per TransDecoder ORF. |
| `09_add_Average_TPM_to_BRAKER_isoforms.py` | Adds the matching TransDecoder ORF ID and Average_TPM value to BRAKER4 isoforms using the best ORF-to-BRAKER mapping. |
| `10_make_clean_BRAKER_isoform_table.py` | Creates the final clean BRAKER4 isoform-level table, one row per predicted protein isoform. |
| `11_make_boss_review_gene_table_PTredo.py` | Creates the simplified seven-column review table (gene ID, contig, compartment, gene length, functional annotation, Average_TPM, Phaeodactylum yes/no). |
| `12_make_hic_network_files.py` | Converts the Hi-C contig-contact table into GEXF and GraphML network files. |
| `13_make_hic_primary_mapq30_pid95_tables.py` | Parses separate Hi-C R1/R2 BAM files, removes non-primary alignments, keeps MAPQ ≥ 30 and percent identity ≥ 95 alignments. |
| `14_make_hic_pair_type_tables.py` | Joins R1/R2 by read ID, assigns contig types using the diatom draft genome, and creates high-confidence Hi-C pair-type tables. |
| `15_make_hic_simple_mixed_read_table.py` | Creates the final simplified read-level table for mixed diatom-bacterial Hi-C pairs. |

---

# Analysis workflow
Click each section to expand the commands, notes, and outputs.

---

<details>
<summary><strong>1. Genome assembly</strong> — Flye metagenome assembly</summary>

Nanopore reads (basecalled with Guppy) were assembled with Flye in metagenome mode, since the sample was a diatom-associated microbial consortium rather than an isolate genome.

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
This assembly was the starting point for polishing, read mapping, binning, organelle screening, and gene annotation.

</details>

---

<details>
<summary><strong>2. Read mapping and assembly support</strong> — minimap2 and samtools</summary>

Short reads were mapped to the Nanopore assembly to assess read support, mapping rate, and contig-level coverage.

```bash
minimap2 -ax sr \
    guppy_flye_assembly.fasta \
    Diatoms_merged.fastq.gz \
    > sr_alignment.sam

samtools view -S -b sr_alignment.sam > alignment.bam
samtools sort alignment.bam -o alignment_sorted.bam
samtools index alignment_sorted.bam
```
Aligns short reads to the assembly, then converts, sorts, and indexes the BAM file.

```bash
samtools flagstat alignment_sorted.bam > mapping_stats.txt
samtools idxstats alignment_sorted.bam | sort -k3,3rn > sr_all_nanopore_hits.tsv
samtools depth alignment_sorted.bam > sr_depth.txt
awk '{sum[$1]+=$3; count[$1]++} END {for (c in sum) print c, sum[c]/count[c]}' sr_depth.txt \
    | sort -k2,2nr > sr_mean_depth.tsv
```
These summarize overall mapping rate, mapped-read counts per contig, per-base depth, and mean depth per contig (highest to lowest coverage) — used to evaluate short-read support across the assembly.

</details>

---

<details>
<summary><strong>3. Assembly polishing</strong> — Medaka, Polypolish, and Pypolca</summary>

### 3.1 Long-read polishing (Medaka)
```bash
medaka_consensus \
    -i pass_trim.fastq.gz \
    -d guppy_flye_assembly.fasta \
    -o medaka_euk_polished \
    -t 12
```

### 3.2–3.3 Short-read polishing (Polypolish)
```bash
bwa mem -t 16 -a medaka_euk_polished/consensus.fasta \
    Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz > alignments_1.sam
bwa mem -t 16 -a medaka_euk_polished/consensus.fasta \
    Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz > alignments_2.sam

polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam \
    --out1 filtered_1.sam --out2 filtered_2.sam

polypolish polish medaka_euk_polished/consensus.fasta \
    filtered_1.sam filtered_2.sam > sr_poly.fasta
```
Short reads are aligned separately to the Medaka-polished assembly, filtered into Polypolish's expected format, then used to correct the assembly.

### 3.4 Additional polishing (Pypolca)
```bash
pypolca run \
    -a sr_poly.fasta \
    -1 Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz \
    -2 Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz \
    -t 12 \
    -o sr_pypolca_output \
    --careful
```
The final corrected assembly (`pypolca_corrected.fasta`) was used for downstream binning, organelle identification, and gene annotation.

### 3.5 BUSCO assessment
```bash
busco -i 1_sr_pypolca_output/pypolca_corrected.fasta \
      -o stram_pyloca_busco_results \
      -m genome \
      -l stramenopiles_odb12 \
      --metaeuk \
      --cpu 32
```
```text
C:87.5%[S:79.1%,D:8.5%],F:2.2%,M:10.3%,n:697
Complete BUSCOs (C):                 610
Complete and single-copy BUSCOs (S): 551
Complete and duplicated BUSCOs (D):   59
Fragmented BUSCOs (F):                15
Missing BUSCOs (M):                   72
Total BUSCO groups searched:         697

Assembly statistics:
Number of scaffolds: 4,925   Number of contigs: 4,925
Total length:         189,915,395 bp   Percent gaps: 0.000%
Scaffold N50:         72 kbp          Contig N50:    72 kbp

Dependencies: hmmsearch 3.4, metaeuk 7.bba0d80
```
This BUSCO run was performed on the polished whole-assembly (`pypolca_corrected.fasta`, 4,925 contigs), not the nuclear-diatom subset, and used MetaEuk gene prediction for the eukaryotic BUSCO assessment.

</details>

---

<details>
<summary><strong>4. Metagenomic binning</strong> — MetaBAT2</summary>

### 4.1 Map Nanopore reads to the polished assembly
```bash
minimap2 -ax map-ont -t 16 \
    1_sr_pypolca_output/pypolca_corrected.fasta \
    pass_trim.fastq.gz \
    | samtools view -@ 16 -bS - \
    | samtools sort -@ 16 -m 10G -o aligned_reads.sorted.bam

samtools index -@ 16 aligned_reads.sorted.bam
```

### 4.2 Contig depth
```bash
jgi_summarize_bam_contig_depths \
    --outputDepth depth.txt \
    --percentIdentity 85 \
    aligned_reads.sorted.bam
```

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
Bins the polished assembly by contig composition and Nanopore read-depth. Resulting bins were used for quality assessment and taxonomic classification.

</details>

---

<details>
<summary><strong>5. Bin quality assessment</strong> — CheckM2</summary>

```bash
checkm2 predict \
    --threads 16 \
    --input 2_metabat2_bins/ \
    --output_directory 3_checkm2_results
```
Estimates completeness and contamination per bin, used to assess bin quality before taxonomy and interpretation.

</details>

---

<details>
<summary><strong>6. Taxonomic classification</strong> — GTDB-Tk</summary>

```bash
gtdbtk classify_wf \
    --genome_dir 2_metabat2_bins/ \
    --out_dir 4_gtdbtk_output \
    --cpus 16 \
    -x fa
```
GTDB-Tk classifies bacterial/archaeal bins. Because the consortium also contained a dominant eukaryotic diatom, MetaEuk-based ORF taxonomy (Section 7) was added as a contig-level screen for eukaryotic, bacterial, ambiguous, and unclassified fractions.

</details>

---

<details>
<summary><strong>7. MetaEuk-based contig classification</strong> — MetaEuk and custom Python script</summary>

MetaEuk ORF-level taxonomic assignments were used to classify contigs, since bin-level bacterial taxonomy alone doesn't resolve eukaryotic contigs or mixed bins in a diatom-associated consortium.

Because MetaEuk's LCA assignments can place organelle-derived sequences into bacterial lineages, mitochondrial-like (`o_Rickettsiales` / `o__Rickettsiales`) and chloroplast-like (`p_Cyanobacteria`) hits were grouped with direct eukaryotic hits when scoring "eukaryotic" content.

### 7.1 Input files
```text
metaeuk_output_polyp_taxonomy_tax_per_pred.tsv   # ORF-level MetaEuk taxonomy (Contig_ID, Classification)
contig_to_bin.txt                                 # links contig IDs to bin IDs
```

### 7.2 ORF-level labels

| Label | Rule |
| --- | --- |
| Eukaryota | `Classification` contains `d_Eukaryota` |
| Mitochondria-derived | contains `o_Rickettsiales` / `o__Rickettsiales` |
| Chloroplast-derived | contains `p_Cyanobacteria` |
| Bacteria | contains `d_Bacteria`, excluding organelle-derived categories |
| Ambiguous | exactly `_cellular organisms` |
| Other | Archaea, viruses, other biological hits |
| Unclassified | no MetaEuk classification |

### 7.3 Final contig-level classification rule
A contig was classified `Eukaryota` if:
```text
(Eukaryota + Mitochondria-derived + Chloroplast-derived) / Total biological ORFs > 0.30
```
where `Total biological ORFs = Eukaryota + Mitochondria-derived + Chloroplast-derived + Bacteria + Ambiguous + Other`. Otherwise, the remaining labels were assigned by dominant biological category; contigs with no clear dominant category were labeled ambiguous.

### 7.4 Script
`scripts/01_classify_metaeuk_contigs.py` — run from the directory containing the MetaEuk taxonomy table and `contig_to_bin.txt`. Classifies each contig using ORF-level assignments and writes the final classification table.

### 7.5 Output
`contig_classification_final_priority.csv` — bin name, contig ID, ORF-level category counts, and final classification.

</details>

---

<details>
<summary><strong>8. Genome coverage and relative abundance</strong> — CoverM</summary>

```bash
coverm genome \
    --genome-fasta-directory 2_metabat2_bins/bac_bins \
    --genome-fasta-extension fa \
    -1 Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz \
    -2 Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz \
    --mapper bwa-mem \
    -m mean relative_abundance covered_fraction \
    --threads 8 \
    --min-read-percent-identity 95 \
    -o bac_output_coverm.tsv
```
Calculates mean coverage, relative abundance, and covered fraction per bacterial genome bin.

</details>

---

<details>
<summary><strong>9. 18S rRNA phylogenetic analysis</strong> — Clustal Omega, TrimAl, and IQ-TREE 2</summary>

```bash
cat *.fasta > 18S_new.fasta

clustalo -i 18S_new.fasta -o 18S_aligned.fasta

trimal -in 18S_aligned.fasta -out 18S_trimmed.fasta -automated1

iqtree2 \
    -s 18S_trimmed.fasta \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -nt AUTO
```
Combines 18S FASTA files, aligns, trims poorly aligned regions, then infers a tree with ModelFinder best-fit model selection and ultrafast bootstrap/SH-aLRT branch support.

</details>

---

<details>
<summary><strong>10. Transcriptome analysis</strong> — nf-core/metatdenovo</summary>

### 10.1 Java setup
```bash
module purge
module load java/openjdk-23.0.1
export JAVA_HOME=$(dirname $(dirname $(readlink -f $(which java))))
export PATH=$JAVA_HOME/bin:$PATH
```

### 10.2 Run pipeline
```bash
~/nextflow run nf-core/metatdenovo \
    -profile singularity \
    --input samplesheet.csv \
    --outdir new_results \
    -w work \
    --assembler spades \
    --orf_caller transdecoder \
    --eggnog_dbpath eggnog_db \
    --skip_kofam true \
    --hmmfiles Pfam-A.hmm \
    --eukulele_dbpath eukulele_db \
    --eukulele_db mmetsp \
    -resume \
    -with-report report_skipK.html \
    -with-timeline timeline_skipK.html
```
Runs transcript assembly (SPAdes) and ORF prediction (TransDecoder), generating transcript-level annotation files used later for expression integration.

</details>

---

<details>
<summary><strong>11. rRNA gene identification from transcriptome assemblies</strong> — Barrnap</summary>

Barrnap was run for each kingdom to identify eukaryotic, bacterial, and mitochondrial rRNA transcripts from the assembled transcriptome:

```bash
barrnap --kingdom euk  --threads 4 spades.transcripts.fa --outseq euk_transcript_rRNA.fna  > diatom_euk_rRNA.gff
barrnap --kingdom bac  spades.transcripts.fa --outseq bac_transcript_rRNA.fna  > diatom_bac_rRNA.gff
barrnap --kingdom mito spades.transcripts.fa --outseq mito_transcript_rRNA.fna > diatom_mito_rRNA.gff
```
The resulting FASTA/GFF files were used to inspect rRNA transcript origin in the consortium transcriptome.

</details>

---

<details>
<summary><strong>12. Organelle genome identification</strong> — MetaQUAST</summary>

```text
Mitogenome reference:   MT742552
Chloroplast reference:  MT742551
```

Organelle contigs were identified by comparing (a) a candidate diatom bin and (b) the polished whole assembly against the chloroplast/mitochondrial references:

```bash
metaquast.py 8_diatom.fasta \
    -R organelle/ref/ \
    -o ./8_metaquast_output

metaquast.py 1_sr_pypolca_output/pypolca_corrected.fasta \
    -R organelle/ref/ \
    -o ./whole_metaquast_output
```
Used to identify candidate chloroplast and mitochondrial contigs — including ones outside the binned assembly — for downstream organelle genome refinement and annotation.

</details>

---

<details>
<summary><strong>13. Diatom genome annotation with BRAKER4 ET mode</strong> — BRAKER4, STAR, GeneMark-ET, AUGUSTUS, and TSEBRA</summary>

Gene models were generated with BRAKER4 using RNA-seq evidence. The genome was not pre-masked, so repeat masking (RepeatModeler, RepeatMasker, TRF) was performed internally within BRAKER4. The final accepted run used **ET mode** (RNA-seq evidence only).

```text
Genome FASTA → STAR indexing → STAR RNA-seq alignment → coordinate-sorted BAM
→ BRAKER4 internal repeat masking → BRAKER4 ET mode (GeneMark-ET → AUGUSTUS
training/prediction → TSEBRA refinement) → final gene models, proteins, CDS, BUSCO
```

### 13.1 Setup
```bash
git clone https://github.com/Gaius-Augustus/BRAKER4.git
cd BRAKER4
singularity pull braker3.sif docker://teambraker/braker3:latest
```
GeneMark license key stored at `~/.gm_key`.

### 13.2–13.3 STAR indexing and alignment
```bash
STAR --runThreadN 24 --runMode genomeGenerate \
    --genomeDir genome_index \
    --genomeFastaFiles genome_index/18_diatom.fasta \
    --genomeSAindexNbases 10

STAR --runThreadN 24 \
    --genomeDir genome_index \
    --readFilesCommand zcat \
    --readFilesIn R1_rep1.fastq.gz,R1_rep2.fastq.gz,R1_rep3.fastq.gz \
                  R2_rep1.fastq.gz,R2_rep2.fastq.gz,R2_rep3.fastq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix genome_index/Diatoms_Combined_

samtools index genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam
```

### 13.4 Genome input and masking check
`genome_index/18_diatom.fasta` — 3,010 contigs, ~82.17 Mbp.
```bash
seqkit stats genome_index/18_diatom.fasta
grep -v "^>" genome_index/18_diatom.fasta | grep -q '[a-z]' && echo "soft-masked" || echo "not soft-masked"
```
Output: `not soft-masked` → `genome_masked` left empty in `samples.csv`; internal repeat masking enabled in `config.ini`.

### 13.5 BUSCO on the genome input
```bash
busco -i genome_index/18_diatom.fasta \
    -l stramenopiles_odb12 -m genome \
    -o busco_18_diatom_stramenopiles_odb12 -c 24
```
```text
BUSCO 6.0.0, stramenopiles_odb12, euk_genome_met mode (MetaEuk gene predictor)
C:86.1%[S:82.4%,D:3.7%],F:2.2%,M:11.8%,n:697
Scaffolds/Contigs: 3,010   Total length: 82,172,226 bp   N50: 48 kbp
```

### 13.6–13.7 Sample sheet and config
```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage
DL_diatom,genome_index/18_diatom.fasta,,,genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam,,,,,,,,stramenopiles_odb12
```
`protein_fasta` left empty to force ET mode and avoid GeneMark-ETP.

```ini
[paths]
braker_container = BRAKER4/braker3.sif
genemark_key = ~/.gm_key

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
`run_red = True` enables internal repeat masking; `mode = et` selects RNA-seq-only BRAKER4 ET.

### 13.8–13.9 Snakemake dry run and execution
```bash
snakemake -s Snakefile --use-singularity \
    --singularity-args "--bind ...,~/" \
    --cores 24 --latency-wait 120 --printshellcmds --rerun-incomplete -n
```
Dry run confirmed ET-mode rules (`run_stringtie`, `bam2hints`, `run_genemark_et`, `train_augustus`, `run_augustus_hints`, `run_tsebra`, `busco_proteins`, `collect_results`) — presence of `run_genemark_et` and absence of `run_genemark_etp` confirmed ET mode.

```bash
snakemake -s Snakefile --use-singularity \
    --singularity-args "--bind ...,~/" \
    --cores 24 --latency-wait 120 --printshellcmds --rerun-incomplete
```
Full run: internal repeat masking, RNA-seq hint generation, GeneMark-ET prediction, AUGUSTUS training/prediction, TSEBRA refinement, BUSCO, result collection.
```text
StringTie transcripts: 9,000
RNA-seq intron hints:  19,114
```

### 13.10 Final outputs
```text
final_annotation_ET/
├── DL_diatom.braker4.ET.gff3.gz
├── DL_diatom.braker4.ET.gtf.gz
├── DL_diatom.braker4.ET.proteins.faa.gz
├── DL_diatom.braker4.ET.cds.fna.gz
├── DL_diatom.braker4.ET.utr.gtf.gz
├── gene_support.tsv
├── software_versions.tsv
├── braker_report.html
├── braker_citations.bib
└── quality_control/
```
```bash
cd final_annotation_ET
gunzip -c DL_diatom.braker4.ET.gff3.gz     > DL_diatom.braker4.ET.gff3
gunzip -c DL_diatom.braker4.ET.gtf.gz      > DL_diatom.braker4.ET.gtf
gunzip -c DL_diatom.braker4.ET.proteins.faa.gz > DL_diatom.braker4.ET.proteins.faa
gunzip -c DL_diatom.braker4.ET.cds.fna.gz  > DL_diatom.braker4.ET.cds.fna
```
No additional TSEBRA run was needed — refinement was included within BRAKER4.

### 13.11 Annotation statistics
```bash
grep -c $'\tgene\t' DL_diatom.braker4.ET.gff3
grep -c $'\ttranscript\t' DL_diatom.braker4.ET.gff3
grep -c $'\tCDS\t' DL_diatom.braker4.ET.gff3
grep -c "^>" DL_diatom.braker4.ET.proteins.faa
grep -c "^>" DL_diatom.braker4.ET.cds.fna
seqkit stats DL_diatom.braker4.ET.proteins.faa DL_diatom.braker4.ET.cds.fna
```
```text
Genes: 15,102   Transcripts: 16,947   Proteins: 16,947
CDS FASTA: 16,947   CDS features: 31,713   Exons: 31,713   Introns: 14,952
```
```bash
grep -n "\*" DL_diatom.braker4.ET.proteins.faa | head   # check for internal stop codons
```
No internal stop codons were detected.

### 13.12 BUSCO on the final protein set
```bash
busco -i DL_diatom.braker4.ET.proteins.faa \
    -l stramenopiles_odb12 -m proteins \
    -o busco_DL_diatom_braker4_ET_proteins_odb12 -c 24 \
    --download_path databases/busco --offline
```
```text
C:84.8%[S:81.1%,D:3.7%],F:1.9%,M:13.3%,n=697
```

### 13.13 Annotation acceptance
Accepted for downstream analysis: plausible gene set, no internal stop-codon issues, 84.8% BUSCO recovery with low duplication.

Final accepted files: `DL_diatom.braker4.ET.{gff3,gtf,proteins.faa,cds.fna}`

### 13.14 Why ET mode instead of ETP
BRAKER4 was first tested in ETP mode (RNA-seq + protein evidence), but GeneMark-ETP failed during model training:
```text
genes: 0   transcripts: 0   CDS: 1724
ERROR: GeneMark-ETP failed, no genemark.gtf
```
The annotation was rerun in ET mode (RNA-seq only) using GeneMark-ET, AUGUSTUS, and TSEBRA, with `protein_fasta` empty and `mode = et`.

</details>

---

<details>
<summary><strong>14. Functional annotation of BRAKER4-predicted proteins</strong> — DIAMOND, Swiss-Prot, Bacillariophyta, InterProScan, and AntiFam</summary>

The final BRAKER4 ET protein set (**16,947 predicted proteins**) was annotated using three complementary layers: curated Swiss-Prot homology, diatom-focused UniProtKB Bacillariophyta homology, and InterProScan domain/family annotation.

```bash
mkdir -p functional_annotation_swissprot/{00_databases,01_input,02_diamond,03_best_hits,05_interproscan,06_combined_annotation,logs,scripts,slurm}
ln -sfn final_annotation_ET/DL_diatom.braker4.ET.proteins.faa \
    01_input/diatom_predicted_proteins.fa
```

### 14.1 Strategy
Swiss-Prot = conservative curated layer (manually reviewed, fewer hits expected). Bacillariophyta UniProtKB = diatom-specific homologs. InterProScan = conserved domains, families, GO terms, pathways. eggNOG/KEGG/COG/dbCAN were **not** used — the strategy prioritized eukaryote/diatom-focused annotation over broad prokaryotic orthology.

### 14.2–14.3 Swiss-Prot search
```bash
diamond makedb --in uniprot_sprot.fasta.gz --db uniprot_sprot.dmnd
# UniProtKB/Swiss-Prot Release 2026_02 (10-Jun-2026): 575,503 reviewed sequences

diamond blastp \
    --query 01_input/diatom_predicted_proteins.fa \
    --db uniprot_sprot.dmnd \
    --out 02_diamond/DL_diatom_braker4_ET_vs_swissprot.tsv \
    --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle \
    --evalue 1e-5 --max-target-seqs 5 --sensitive --threads 32
```
```text
Total proteins: 16,947   Hit lines: 36,669
With ≥1 hit: 8,078 (47.67%)   Strict hit: 4,284 (25.28%)
Confidence — High: 2,037  Medium: 2,247  Low: 2,826  Weak/fragment: 968
```

### 14.4–14.5 UniProtKB Bacillariophyta search
```bash
curl -L -o uniprotkb_bacillariophyta_taxid2836.fasta.gz \
    "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28taxonomy_id%3A2836%29"
diamond makedb --in uniprotkb_bacillariophyta_taxid2836.fasta.gz --db uniprotkb_bacillariophyta_taxid2836.dmnd

diamond blastp \
    --query 01_input/diatom_predicted_proteins.fa \
    --db uniprotkb_bacillariophyta_taxid2836.dmnd \
    --out 02_diamond/DL_diatom_braker4_ET_vs_uniprot_bacillariophyta.tsv \
    --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle \
    --evalue 1e-5 --max-target-seqs 5 --sensitive --threads 32
```
```text
Total proteins: 16,947   Hit lines: 65,416
With ≥1 hit: 13,574 (80.10%)   Strict hit: 11,294 (66.64%)
Confidence — High: 8,557  Medium: 2,737  Low: 1,668  Weak/fragment: 612
```

### 14.6 Best-hit parsing and confidence filtering
For each query: query/subject coverage calculated, UniProt ID/protein name/organism/gene name/taxon ID/evidence parsed. Best hit chosen by lowest e-value → highest bitscore → highest query coverage → highest percent identity.

```text
High:   e-value ≤ 1e-20, qcov ≥ 70%, pident ≥ 40%
Medium: e-value ≤ 1e-10, qcov ≥ 50%, pident ≥ 30%
Low:    e-value ≤ 1e-5,  qcov ≥ 30%
Weak domain/fragment: everything else

Strict set = e-value ≤ 1e-10, qcov ≥ 50%, pident ≥ 30%
```

### 14.7–14.8 Parsing scripts
`scripts/02_make_swissprot_best_hits.py` and `scripts/03_make_bacillariophyta_best_hits.py` (same coverage/ranking/confidence framework) — run with:
```bash
conda activate swissprot_annot
python scripts/02_make_swissprot_best_hits.py
python scripts/03_make_bacillariophyta_best_hits.py
```
Outputs (per database): `*_all_hits_with_coverage.tsv`, `*_best_hits.tsv`, `*_best_hits_strict.tsv`, `*_all_proteins_with_*_annotation.tsv`, `*_annotation_summary.txt`.

### 14.9 InterProScan
```bash
sbatch slurm/interproscan_diatom_full.sh
```
```bash
interproscan.sh \
    -i 01_input/diatom_predicted_proteins.fa \
    -f TSV -dp -goterms -pa -exclappl MobiDBLite -cpu 32 \
    -o 05_interproscan/DL_diatom_braker4_ET_interproscan.tsv
```
MobiDBLite was excluded after a Python-compatibility failure in its bundled script; it predicts disordered regions and wasn't central to the annotation goals, so it was dropped while keeping the main family/domain/GO databases (AntiFam, CDD, Coils, FunFam, Gene3D, Hamap, NCBIfam, PANTHER, Pfam, PIRSF, PIRSR, PRINTS, ProSitePatterns, ProSiteProfiles, SFLD, SMART, SUPERFAMILY).

```text
Raw rows: 102,153   Proteins with hits: 13,106 / 16,947 (77.34%)
```

### 14.10 Per-protein summary
```bash
conda activate swissprot_annot
python scripts/04_summarize_interproscan.py
```
Outputs: `DL_diatom_interproscan_summary_by_protein.tsv`, `DL_diatom_all_proteins_with_interproscan_summary.tsv`, `DL_diatom_interproscan_analysis_counts.tsv`, `DL_diatom_interproscan_summary_stats.txt` — includes protein_id, protein_length, row/analysis counts, signature accessions/descriptions, InterPro accessions/descriptions, GO terms, and pathway annotations.

Raw-row contribution by database: Pfam 17,891 · Gene3D 17,698 · SUPERFAMILY 13,811 · PANTHER 9,840 · SMART 7,367 · PRINTS 7,354 · ProSiteProfiles 7,018 · Coils 5,536 · CDD 5,075 · FunFam 3,292 · NCBIfam 2,853 · ProSitePatterns 2,748 · Hamap 901 · PIRSF 537 · SFLD 230 · AntiFam 2.

### 14.11–14.12 AntiFam screening and master table
Two AntiFam matches were flagged (not treated as functional annotation, but as spurious/RNA-ORF warnings):
```text
g10893.t1   ANF00012   tRNA
g11404.t1   ANF00005   Antisense to 23S rRNA
```
```bash
awk -F '\t' '$4=="AntiFam"' 05_interproscan/DL_diatom_braker4_ET_interproscan.tsv | column -t
```
A small flag table (`DL_diatom_antifam_flagged_proteins.tsv`) was built, then merged with the Swiss-Prot, Bacillariophyta, and InterProScan tables using:
```bash
conda activate swissprot_annot
python scripts/05_merge_functional_annotation_layers.py
```
```bash
wc -l 06_combined_annotation/DL_diatom_master_functional_annotation.tsv
# 16,948 lines = 16,947 proteins + 1 header
```
Both flagged proteins were confirmed retained and marked in the merged table.

### 14.13 AntiFam-filtered views
```bash
awk -F'\t' 'NR==1 || $9=="yes"'  DL_diatom_master_functional_annotation.tsv > DL_diatom_master_functional_annotation_AntiFam_flagged_only.tsv
awk -F'\t' 'NR==1 || $9!="yes"'  DL_diatom_master_functional_annotation.tsv > DL_diatom_master_functional_annotation_no_AntiFam.tsv
```

### 14.14 Final functional annotation outputs
```text
06_combined_annotation/
├── DL_diatom_interproscan_summary_by_protein.tsv
├── DL_diatom_all_proteins_with_interproscan_summary.tsv
├── DL_diatom_interproscan_analysis_counts.tsv
├── DL_diatom_interproscan_summary_stats.txt
├── DL_diatom_antifam_flagged_proteins.tsv
├── DL_diatom_master_functional_annotation.tsv               # protein_id, recommended annotation/source/confidence,
├── DL_diatom_master_functional_annotation_for_manual_categories.tsv   # Swiss-Prot + Bacillariophyta fields, InterProScan
├── DL_diatom_master_functional_annotation_AntiFam_flagged_only.tsv    # domains/signatures, GO terms, pathways, AntiFam flag
├── DL_diatom_master_functional_annotation_no_AntiFam.tsv
└── DL_diatom_master_functional_annotation_summary.txt
```

### 14.15 Expression integration status
TransDecoder ORFs were aligned to BRAKER4 ET proteins via DIAMOND BLASTP; one best BRAKER4 hit was kept per ORF to avoid duplicating a transcript's `Average_TPM` across multiple BRAKER4 proteins. Only `Average_TPM` was carried forward — no all-hit sums, means, hit counts, or ORF lists. Full workflow: Section 17.

### 14.16 Status checklist
Swiss-Prot DB + DIAMOND search + best-hit parsing ✓ · Bacillariophyta DB + DIAMOND search + best-hit parsing ✓ · InterProScan install/run/summary ✓ · AntiFam screen ✓ · Master table + AntiFam views ✓ · Expression integration (TransDecoder + Average_TPM) ✓ · Final clean BRAKER4 isoform table ✓

</details>

---

<details>
<summary><strong>15. Nuclear-enriched genome generation</strong> — organelle contig removal and BRAKER4 annotation filtering</summary>

The nuclear-enriched genome was generated by removing contigs with strong chloroplast/mitochondrial similarity from the whole diatom assembly.

### 15.1–15.2 Inputs and organelle reference
```bash
WHOLE=genome_index/18_diatom.fasta                                            # 3,010 contigs, 82,172,226 bp
CHLORO=organelle/2_chloro/chloroplast_contig_1443_trimmed.fasta               # 1 contig, 120,429 bp
MITO=organelle/mito/diatom_candidate_mitochondrion_2contigs.fasta             # 2 contigs, 104,526 bp

cat $CHLORO $MITO > organelles_chloro_mito.fasta   # 3 sequences, 224,955 bp
```

### 15.3–15.4 Align genome to organelle references and calculate coverage
```bash
minimap2 -x asm5 -c organelles_chloro_mito.fasta $WHOLE > whole_vs_organelles.paf

awk 'BEGIN{OFS="\t"} {print $1, $3, $4}' whole_vs_organelles.paf \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - > whole_vs_organelles.query_intervals.merged.bed

seqkit fx2tab -n -l $WHOLE > whole_contig_lengths.tsv

awk 'BEGIN{OFS="\t"} {aligned[$1] += ($3 - $2)} END {for (c in aligned) print c, aligned[c]}' \
    whole_vs_organelles.query_intervals.merged.bed > organelle_aligned_length_per_contig.tsv

awk 'BEGIN{OFS="\t"}
FNR==NR {len[$1]=$2; next}
{contig=$1; aligned=$2; pct=(aligned/len[contig])*100; print contig, len[contig], aligned, pct}' \
    whole_contig_lengths.tsv organelle_aligned_length_per_contig.tsv \
    > organelle_coverage_per_contig.tsv
```

### 15.5 Remove organelle-like contigs (≥70% aligned length)
```bash
awk '$4 >= 70 {print $1}' organelle_coverage_per_contig.tsv > organelle_like_contigs.70pct.txt

seqkit grep -v -f organelle_like_contigs.70pct.txt $WHOLE > 18_diatom_nuclear_enriched.v1.fasta
```
Three organelle-like contigs identified: `contig_1443` (chloroplast-like), `contig_5628` and `contig_1647` (mitochondrial-like) — 260,454 bp (0.317%) removed.
```text
18_diatom.fasta                     3,010 contigs   82,172,226 bp
18_diatom_nuclear_enriched.v1.fasta 3,007 contigs   81,911,772 bp
```

### 15.6 Filter BRAKER4 annotation to nuclear contigs
BRAKER4 was not rerun (only 3 contigs removed); the existing annotation was instead filtered to retained nuclear contigs.
```bash
seqkit seq -n 18_diatom_nuclear_enriched.v1.fasta > nuclear_contigs.v1.txt

awk 'BEGIN{FS=OFS="\t"}
FNR==NR {keep[$1]=1; next}
$0 ~ /^#/ {if ($0 !~ /^##FASTA/) print; next}
$1 in keep {print}' nuclear_contigs.v1.txt DL_diatom.braker4.ET.gff3 \
    > braker.18_diatom_nuclear_only.v1.gff3
```
Final nuclear genome files: `18_diatom_nuclear_enriched.v1.fasta`, `braker.18_diatom_nuclear_only.v1.gff3`.

</details>

---

<details>
<summary><strong>16. Pairwise genome comparison with <em>Phaeodactylum tricornutum</em></strong> — BLASTN, GFF3, and bedtools</summary>

A nucleotide-level similarity screen (not a full orthology analysis) between the BRAKER4-annotated diatom genome and *P. tricornutum*. Raw BLASTN hits were retained unfiltered and linked to overlapping gene models in both genomes.

### 16.1–16.3 Setup and reference download
```bash
mkdir -p phaeodactylum_to_diatom_blastn_redo/{00_inputs,01_db,02_blast,03_filtered,04_summary,logs,scripts}
ln -sf genome_index/18_diatom.fasta 00_inputs/diatom_genome.fasta

wget -O 00_inputs/Phaeodactylum_ASM15095v2_genomic.fna.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/955/GCF_000150955.2_ASM15095v2/GCF_000150955.2_ASM15095v2_genomic.fna.gz
gunzip -c 00_inputs/Phaeodactylum_ASM15095v2_genomic.fna.gz > 00_inputs/phaeodactylum_genome.fna

wget -O 00_inputs/phaeodactylum_ASM15095v2.gff3.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/955/GCF_000150955.2_ASM15095v2/GCF_000150955.2_ASM15095v2_genomic.gff.gz
gunzip -c 00_inputs/phaeodactylum_ASM15095v2.gff3.gz > 00_inputs/phaeodactylum_ASM15095v2.gff3
```
```text
diatom_genome.fasta       3,010 seqs   82,172,226 bp
phaeodactylum_genome.fna     88 seqs   27,450,724 bp
```

### 16.4–16.5 BLAST database and search
```bash
makeblastdb -in 00_inputs/diatom_genome.fasta -dbtype nucl -parse_seqids \
    -out 01_db/diatom_genome_blastdb -title "DL_diatom_genome"
# 3,010 sequences added

blastn \
    -task dc-megablast \
    -query 00_inputs/phaeodactylum_genome.fna \
    -db 01_db/diatom_genome_blastdb \
    -out 02_blast/phaeodactylum_vs_diatom_dcmegablast.tsv \
    -evalue 1e-10 -perc_identity 60 -num_threads 16 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"
```
`P. tricornutum` genome as query, diatom genome as subject DB. **28,934 raw BLASTN hits** (retained unfiltered). A headered TSV copy was also created for readability.

### 16.6–16.7 Convert genes to BED
```bash
# Phaeodactylum genes (GFF3 → BED), extracting ID/Name/gene/locus_tag/strand
awk -F'\t' '... $3=="gene" ...' 00_inputs/phaeodactylum_ASM15095v2.gff3 > 00_inputs/phaeodactylum_genes.bed
# 10,392 genes

# Diatom BRAKER4 ET genes (GFF3 → BED)
awk -F'\t' '... $3=="gene" ...' DL_diatom.braker4.ET.gff3 > 00_inputs/diatom_BRAKER_ET_genes.bed
# 15,102 genes
```

### 16.8–16.10 Convert BLASTN hits to BED and intersect with genes
Each raw BLASTN hit was converted into a stable-ID interval on both the *Phaeodactylum* (query) and diatom (subject) side (28,934 intervals each), then intersected against each side's gene BED with `bedtools intersect -wa -wb -loj` (keeping hits with no gene overlap too):
```bash
bedtools intersect -a 02_blast/phaeodactylum_blast_intervals.bed -b 00_inputs/phaeodactylum_genes.bed -wa -wb -loj \
    > 03_filtered/phaeodactylum_hits_with_PT_genes.tsv       # 28,976 rows
bedtools intersect -a 02_blast/diatom_blast_intervals.bed -b 00_inputs/diatom_BRAKER_ET_genes.bed -wa -wb -loj \
    > 03_filtered/diatom_hits_with_BRAKER_ET_genes.tsv       # 29,005 rows
```
(Row counts exceed raw hit counts slightly because some intervals overlap more than one gene.)

### 16.11–16.12 Merge and summarize
```bash
python scripts/06_merge_phaeodactylum_blast_hits.py
```
Keeps one row per raw BLASTN hit, adding overlapping gene info from both genomes (semicolon-collapsed if multiple genes overlap).

Final table: `04_summary/phaeodactylum_vs_diatom_BLASTN_with_PT_and_BRAKER_ET_genes.tsv` — **28,934 hits × 24 columns** (blast_hit_id, pt_contig, diatom_contig, pident, aln_len, mismatch, gapopen, pt/diatom start-end, evalue, bitscore, pt_len, diatom_len, qcovs, pt_gene_id/name/symbol/locus_tag/strand, diatom_gene_id/attr_id/strand).

```text
Total BLASTN hits:                        28,934
Hits with PT gene:                         8,529
Hits with diatom BRAKER ET gene:           6,693
Hits with both:                            6,173
Unique PT genes hit:                       3,202
Unique diatom BRAKER ET genes hit:         3,386
```

### 16.13 Interpretation
This provides a gene-linked nucleotide similarity table rather than a definitive orthology map, since nucleotide-level similarity doesn't fully capture protein-level conservation or gene orthology.

</details>

---

<details>
<summary><strong>17. Clean BRAKER4 isoform-level gene table construction</strong> — functional annotation, Average_TPM, and Phaeodactylum yes/no</summary>

Builds the final clean BRAKER4 isoform-level tables for pathway curation, manual review, and *P. tricornutum* comparison — one row per BRAKER4 predicted protein isoform. Isoform IDs were not collapsed, no GenBank-derived organelle rows were appended, all-hit TPM summaries were not used, and detailed *P. tricornutum* BLASTN hit columns were not retained.

### 17.1–17.2 Clean rebuild directory and inputs
Rebuilt in `CLEAN_REBUILD_FROM_RAW/` to avoid carrying forward columns from older all-hit TPM summaries, earlier *Phaeodactylum* comparisons, or GenBank-based organelle merges. Inputs:
```text
01_input/DL_diatom.braker4.ET.proteins.faa, DL_diatom.braker4.ET.gff3
02_diamond/*_vs_swissprot.tsv, *_vs_uniprot_bacillariophyta.tsv
05_interproscan/DL_diatom_braker4_ET_interproscan.tsv
06_combined_annotation/DL_diatom_antifam_flagged_proteins.tsv
01_input/transcriptome_orfs.transdecoder.clean.pep
07_expression/master_with_custom_broad_categories.csv
phaeodactylum_to_diatom_blastn_redo/04_summary/phaeodactylum_vs_diatom_BLASTN_with_PT_and_BRAKER_ET_genes.tsv
```
The redo BLASTN summary was used **only** as a yes/no lookup — detailed BLASTN fields (PHATRDRAFT IDs, gene names, pident, aln_len, bitscore, evalue) were not retained.

### 17.3 Functional annotation layers
Swiss-Prot (conservative, reviewed), Bacillariophyta (diatom-focused homology), InterProScan (domains/families/GO/pathways), AntiFam (warning flags, not functional annotation). Verified one row per isoform + header (16,948 lines), with both AntiFam-flagged proteins (`g10893.t1`, `g11404.t1`) retained.

### 17.4 BRAKER4 lengths
GFF3 coordinates (contig ID, gene start/end, strand, gene/CDS/protein length) were added to the master table → `07_expression/DL_diatom_master_functional_annotation_with_lengths.tsv`. All 16,947 rows received gene lengths (0 missing).

### 17.5 TransDecoder ORF → BRAKER4 mapping
TransDecoder peptides were searched (DIAMOND BLASTP) against a BRAKER4-protein DIAMOND DB.
```text
Raw DIAMOND hit rows: 118,472
Unique TransDecoder ORFs with ≥1 hit: 39,935
Unique BRAKER4 proteins hit: 14,473
```
One best BRAKER4 hit was retained per ORF (by alignment quality/coverage) → `transdecoder_ORFs_to_BRAKER4_best_hit_per_ORF.tsv` (39,935 mappings + header). This mapping hit 12,277 unique BRAKER4 proteins.

### 17.6 Average_TPM integration
Only `orf` and `Average_TPM` columns were transferred from the trusted expression table via the best-hit mapping; where multiple ORFs mapped to the same protein, the highest valid `Average_TPM` was kept (ORFs lacking a valid numeric TPM were skipped).
```text
Annotation rows: 16,947
BRAKER proteins with Average_TPM: 12,276   without: 4,671
```
Output: `07_expression/DL_diatom_master_functional_annotation_lengths_Average_TPM.tsv` (adds `transdecoder_orf_id`, `Average_TPM` only — no all-hit sums/means/counts).

### 17.7 Phaeodactylum tricornutum yes/no lookup
Rebuilt from the redo BLASTN gene-overlap summary (Section 16). A diatom gene was marked `yes` only if a BLASTN alignment overlapped **both** a BRAKER4 ET gene model and an annotated *P. tricornutum* gene. Lookup was done on the gene root (e.g., BLASTN's `g10009` ↔ isoform `g10009.t1`), while the full isoform ID was preserved in the final table. Only `present_in_Phaeodactylum_tricornutum` was retained — no detailed BLASTN columns. This is a nucleotide-level similarity screen, **not** confirmed orthology, reciprocal-best-hit, or protein-level conservation.

### 17.8 Final clean BRAKER4 isoform-level table
```text
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table.tsv
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table_sorted_by_Average_TPM.tsv
```
Both: 16,947 rows + header (16,948 lines).
```text
Compartment — nuclear: 16,873   plastid_like: 72   mito_like: 2   (assigned from contig ID only; no GenBank rows appended)
Phaeodactylum yes/no  — no: 13,240   yes: 3,707
```

### 17.9 Simplified review table
`09_final/DL_diatom_FINAL_gene_table_for_boss_PTredo.tsv` — 16,947 rows × 7 columns, for manual/pathway review without carrying forward intermediate annotation or BLASTN columns:

| Column | Description |
|---|---|
| `gene_id` | Original (uncollapsed) BRAKER4 isoform ID |
| `contig_id` | Diatom contig containing the gene model |
| `diatom_compartment` | `nuclear` / `plastid_like` / `mito_like`, from contig ID |
| `diatom_gene_length_bp` | Gene length (bp), from BRAKER4 GFF3 |
| `functional_annotation` | Recommended annotation from Swiss-Prot/Bacillariophyta/InterProScan/AntiFam |
| `diatom_Average_TPM` | Expression value via best TransDecoder ORF-to-BRAKER4 mapping |
| `present_in_Phaeodactylum_tricornutum` | Yes/no nucleotide-level gene-linked similarity (redo BLASTN) |

</details>

---

<details>
<summary><strong>18. Hi-C read mapping and contig-level proximity-ligation network</strong> — BWA-MEM, samtools, awk, YaHS, and Python</summary>

Hi-C paired-end reads were added after the main assembly/annotation/expression/comparative-genomics workflow, to assess assembly representation in the proximity-ligation dataset, identify Hi-C-linked contigs, and test high-confidence diatom-bacterial contacts. Mapping was done against the **polished whole assembly** (not the nuclear-enriched subset), since Hi-C reads came from the complete consortium, allowing diatom/bacterial/mixed contacts to be evaluated in one coordinate space.

### 18.1 Inputs
```text
Hi-C reads:    hi-c_diatoms/1574499_S6_L001_R{1,2}_001.fastq.gz
Assembly:      MAGS_guppy/1_sr_pypolca_output/pypolca_corrected.fasta
Diatom ref:    metatranscriptomics/genome_index/18_diatom.fasta   (contigs here = "diatom"; all remaining assembly contigs = "bacterial")
```

### 18.2 Map Hi-C reads to the whole assembly (paired-end)
```bash
ln -sfn pypolca_corrected.fasta pypolca_corrected.hic_input.fasta
bwa index pypolca_corrected.hic_input.fasta
samtools faidx pypolca_corrected.hic_input.fasta
cut -f1,2 pypolca_corrected.hic_input.fasta.fai > pypolca_corrected.chrom.sizes

bwa mem -SP5M -t 24 pypolca_corrected.hic_input.fasta R1.fastq.gz R2.fastq.gz \
    | samtools view -@ 8 -bS - > hic_to_whole_assembly.raw.bam

samtools sort -@ 24 -n -o hic_to_whole_assembly.name_sorted.bam  hic_to_whole_assembly.raw.bam
samtools sort -@ 24    -o hic_to_whole_assembly.coord_sorted.bam hic_to_whole_assembly.raw.bam
samtools index hic_to_whole_assembly.coord_sorted.bam

samtools flagstat hic_to_whole_assembly.coord_sorted.bam > hic_to_whole_assembly.flagstat.txt
samtools idxstats hic_to_whole_assembly.coord_sorted.bam  > hic_to_whole_assembly.idxstats.txt
```
```text
Primary reads: 890,810   Primary mapped: 622,967 (69.93%)
All mapped alignments: 803,799 / 1,071,642 (75.01%)
Read pairs: 445,405   Singletons: 54,019 (6.06%)
```

### 18.3 Contig-level Hi-C representation
```bash
awk 'BEGIN{OFS="\t"; print "contig","length_bp","mapped_HiC_reads","unmapped_HiC_reads","HiC_mapped"}
$1!="*" {status=($3>0?"yes":"no"); print $1,$2,$3,$4,status}' hic_to_whole_assembly.idxstats.txt \
    > hic_contig_mapping_presence.tsv
```
A combined summary file (`hic_contig_mapping_summary.txt`) records the input files plus:
```text
Total contigs: 4,925   With ≥1 Hi-C read mapped: 4,010 (81.42%)   With 0: 915
```

### 18.4 Exploratory whole-assembly scaffolding (YaHS)
```bash
yahs -o DL_diatom_whole_hic_yahs pypolca_corrected.hic_input.fasta hic_to_whole_assembly.name_sorted.bam
seqkit stats pypolca_corrected.hic_input.fasta DL_diatom_whole_hic_yahs_scaffolds_final.fa
```
```text
Input:  4,925 contigs, 189,915,395 bp, max 5,424,378 bp
Output: 5,032 scaffolds, 189,915,395 bp, max 5,424,378 bp
```
YaHS did not increase max scaffold length and increased sequence count, so this was treated as exploratory, not a final scaffolded assembly.

### 18.5–18.6 All-primary inter-contig contacts
Primary mapped read pairs (unmapped/mate-unmapped/secondary/supplementary excluded), no MAPQ cutoff, each pair counted once if mates map to different contigs:
```bash
samtools view -@ 8 -f 1 -F 2316 hic_to_whole_assembly.name_sorted.bam \
    | awk '... pairs mates by read ID, emits distinct contig pairs ...' \
    | sort -S 20G | uniq -c \
    | awk 'BEGIN{OFS="\t"; print "contig_A","contig_B","HiC_contact_pairs"} {print $2,$3,$1}' \
    | sort -k3,3nr > hic_intercontig_contacts_all_primary_pairs.tsv
```
Contig lengths were joined in for a `.with_lengths.tsv` version, and a connectivity summary (partner count, total pairs, strongest single contact per contig) was generated:
```text
Connected contigs: 3,770   Inter-contig Hi-C links: 75,703
```

### 18.7 Contact network files
```bash
python scripts/12_make_hic_network_files.py
```
Produces `hic_contig_network_all_primary_pairs.{gexf,graphml}` — Nodes = contigs (3,770), Edges = Hi-C contact (75,703), edge weight = supporting read-pair count.

### 18.8 High-confidence separate-read mapping (for diatom-bacterial contacts)
R1 and R2 were mapped **independently** (stricter filter: primary only, MAPQ ≥ 30, percent identity ≥ 95%) so that, per read ID, R1/R2 contig assignment could be compared directly.
```bash
ln -sf pypolca_corrected.fasta 00_inputs/whole_assembly.fasta
ln -sf HiC_R1.fastq.gz 00_inputs/HiC_R1.fastq.gz
ln -sf HiC_R2.fastq.gz 00_inputs/HiC_R2.fastq.gz

bwa index -p 01_bwa_index/whole_assembly 00_inputs/whole_assembly.fasta
```
SLURM job maps R1 and R2 separately against the same index:
```bash
bwa mem -t ${THREADS} 01_bwa_index/whole_assembly 00_inputs/HiC_R1.fastq.gz | \
    samtools sort -@ ${THREADS} -m2G -o 02_alignments/HiC_R1.sorted.bam -
samtools index 02_alignments/HiC_R1.sorted.bam
samtools flagstat 02_alignments/HiC_R1.sorted.bam > 04_logs/HiC_R1.flagstat.txt
# (repeated for R2)
```
```text
R1: 445,405 primary reads, 318,884 mapped (71.59%)
R2: 445,405 primary reads, 304,083 mapped (68.27%)
```

### 18.9 High-confidence pair-type tables and mixed-pair classification
```bash
conda activate hic_diatom
python scripts/13_make_hic_primary_mapq30_pid95_tables.py   # parses NM tag + CIGAR → primary, MAPQ≥30, PID≥95
python scripts/14_make_hic_pair_type_tables.py               # joins R1/R2 by read ID; diatom vs. bacterial contig labels
python scripts/15_make_hic_simple_mixed_read_table.py        # final simplified mixed-pair table
```
Pair-type codes: `1` = both diatom, `2` = both bacterial, `3` = R1 diatom/R2 bacterial, `4` = R1 bacterial/R2 diatom.
```text
Both diatom:              39,739 (64.82%)
Both bacterial:           21,224 (34.62%)
R1 diatom, R2 bacterial:     164 (0.27%)
R1 bacterial, R2 diatom:     177 (0.29%)
→ 341 mixed diatom-bacterial pairs
```
Final table: `HiC_DIATOM_BACTERIAL_read_level_SIMPLE_COLUMNS_MAPQ30_PID95.tsv` — columns: `read_id, read_number, contig_id, paired_contig_id, mapq, percent_identity, aligned_length_bp, pair_type_code, pair_type`. Each mixed pair = 2 rows (682 read rows, 683 lines incl. header).

### 18.10 Final output organization
```text
02_map_to_whole_assembly/final_mapping/
├── hic_to_whole_assembly.{flagstat,idxstats}.txt
├── hic_contig_mapping_summary.txt
├── hic_contig_mapping_presence.tsv
├── hic_to_whole_assembly.{name,coord}_sorted.bam(.bai)
└── pypolca_corrected.chrom.sizes

04_contact_maps/
├── tables/  (all-primary contact + connectivity tables)
├── network_files/  (gexf, graphml)
└── scripts/12_make_hic_network_files.py

hic_bwa_separate_reads/
├── 02_alignments/  (HiC_R1/R2 sorted BAM + index)
├── 03_tables/  (contig-type map, MAPQ30/PID95 tables, joined/typed pairs, pair-type summary, final simple table)
└── 04_logs/  (flagstat + bwa logs)
```

### 18.11 Summary
Hi-C reads mapped to 4,010/4,925 contigs (81.42%); 622,967/890,810 primary reads mapped (69.93%). The all-primary contig-contact network contains 3,770 nodes and 75,703 edges. The high-confidence separate-read analysis (primary, MAPQ ≥ 30, PID ≥ 95%) identified 341 mixed diatom-bacterial Hi-C read pairs (682 read-level rows).

</details>
