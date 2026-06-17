# Diatom Consortia: Metagenomic and Metatranscriptomic Pipeline

This repository documents the workflow used to assemble, polish, bin, classify, annotate, and compare genomes and transcriptomes from a diatom-associated microbial consortium. The analysis combines long-read metagenomic assembly, short-read polishing, metagenomic binning, contig-level taxonomic screening, organelle identification, transcriptome analysis, BRAKER4 ET gene prediction, nuclear-enriched genome generation, functional annotation, expression integration, and comparison with *Phaeodactylum tricornutum*.

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
Expression integration using TransDecoder ORFs and Average_TPM only
   ↓
Phaeodactylum tricornutum comparison summarized as yes/no only
   ↓
Final clean BRAKER4 isoform-level gene table
```
## Repository structure
```text
.
├── README.md
├── data_analysis.md
└── scripts/
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
## Main scripts
```text
make_swissprot_best_hits.py
  Parses Swiss-Prot DIAMOND output and writes best-hit annotation tables.

make_bacillariophyta_best_hits.py
  Parses UniProtKB Bacillariophyta DIAMOND output using the same confidence framework.

summarize_interproscan.py
  Collapses InterProScan TSV output to one row per BRAKER4 protein.

merge_functional_annotation_layers.py
  Merges Swiss-Prot, Bacillariophyta, InterProScan, and AntiFam evidence.

07_add_BRAKER_lengths_clean.py
  Adds contig ID, coordinates, strand, gene length, CDS length, and protein length from the BRAKER4 GFF3 file.

08_make_best_ORF_to_BRAKER_mapping_clean.py
  Keeps one best BRAKER4 hit per TransDecoder ORF from the raw DIAMOND output.

09_add_ONLY_Average_TPM_clean.py
  Adds `transdecoder_orf_id` and `Average_TPM` to BRAKER4 proteins. The ORF ID is retained only when the mapped ORF has a valid numeric `Average_TPM`; all-hit TPM sums, means, hit counts, and mapped-ORF lists are not used.

10_make_FINAL_clean_BRAKER_isoform_table.py
  Creates the final clean BRAKER4 isoform table with PT yes/no only.
```
## Detailed workflow
The full command-by-command workflow is documented in [`data_analysis.md`](data_analysis.md). Each major step is stored in a collapsible section with commands, outputs, and interpretation notes.
