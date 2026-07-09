# Diatom Consortia: Metagenomic and Metatranscriptomic Pipeline
This repository documents the analysis workflow used to assemble, polish, classify, annotate, and compare genomes and transcriptomes from a diatom-associated microbial consortium.

The workflow includes long-read metagenomic assembly, short-read polishing, metagenomic binning, contig-level taxonomic screening, organelle identification, BRAKER4 ET gene prediction, functional annotation, expression integration, comparison with *Phaeodactylum tricornutum*, and Hi-C read mapping/contact-network analysis.

---
## Workflow overview

```text
Nanopore reads
   ↓
Flye metagenome assembly
   ↓
Medaka, Polypolish, and Pypolca polishing
   ↓
Assembly assessment and metagenomic binning
   ↓
CheckM2, GTDB-Tk, and MetaEuk classification
   ↓
Organelle identification
   ↓
BRAKER4 ET genome annotation
   ↓
Functional annotation with Swiss-Prot, Bacillariophyta UniProtKB, InterProScan, and AntiFam
   ↓
Expression integration using TransDecoder ORFs and Average_TPM
   ↓
Phaeodactylum tricornutum comparison using redo BLASTN gene-linked overlap table
   ↓
Final clean BRAKER4 isoform-level gene table
   ↓
Simplified seven-column gene table for manual review
   ↓
Hi-C read mapping and contig-contact network
```
---
## Repository structure
```text
.
├── README.md
├── data_analysis.md
└── scripts/
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
    └── 12_make_hic_network_files.py
```

---
## Main outputs
The main output is a clean BRAKER4 isoform-level gene table:
```text
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table.tsv
```
This table contains one row per predicted BRAKER4 protein isoform and integrates:
```text
functional annotation
BRAKER4 coordinate and length fields
TransDecoder ORF mapping
Average_TPM expression values
compartment labels
Phaeodactylum tricornutum yes/no similarity status
AntiFam warning flags
```
A second version sorted by expression is also generated:
```text
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table_sorted_by_Average_TPM.tsv
```
A simplified review table was also generated for manual inspection and discussion:
```text
09_final/DL_diatom_FINAL_gene_table_for_boss_PTredo.tsv
```
This review table keeps seven key columns:
```text
gene_id
contig_id
diatom_compartment
diatom_gene_length_bp
functional_annotation
diatom_Average_TPM
present_in_Phaeodactylum_tricornutum
```
The *Phaeodactylum tricornutum* yes/no column is based on the BLASTN gene-linked overlap table. It should be interpreted as a nucleotide-level similarity screen, not as confirmed orthology.

---
## Scripts
Custom Python scripts are stored in the `scripts/` directory and are numbered in the order they are used in the workflow.
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
  Adds overlapping Phaeodactylum gene information to the cleaned BLASTN comparison table.

07_add_BRAKER_lengths_clean.py
  Adds BRAKER4 coordinates, contig IDs, strand, gene length, CDS length, and protein length to each isoform.

08_make_best_ORF_to_BRAKER_mapping_clean.py
  Parses TransDecoder ORF versus BRAKER4 DIAMOND output and keeps one best BRAKER4 hit per TransDecoder ORF.

09_add_Average_TPM_to_BRAKER_isoforms.py
  Adds the matching TransDecoder ORF ID and Average_TPM value to BRAKER4 isoforms using the best ORF-to-BRAKER mapping.

10_make_clean_BRAKER_isoform_table.py
  Creates the final clean BRAKER4 isoform-level table with one row per predicted protein isoform.

11_make_boss_review_gene_table_PTredo.py
  Creates the simplified seven-column review table with gene ID, contig ID, compartment, gene length, functional annotation, Average_TPM, and redo Phaeodactylum yes/no status.

12_make_hic_network_files.py
  Converts the Hi-C contig-contact table into GEXF and GraphML network files.
```
---
## Hi-C analysis
Hi-C reads were mapped to the polished whole assembly to assess contig-level representation and inter-contig proximity-ligation links.
The YaHS scaffolding test was treated as exploratory. The final Hi-C output is interpreted as contig-level contact evidence rather than as a final chromosome-scale scaffolded assembly.
Final Hi-C network outputs:
```text
hic_contig_network_all_primary_pairs.gexf
hic_contig_network_all_primary_pairs.graphml
```
---
## Detailed workflow
The full workflow, including commands, outputs, and short notes describing what each command does, is documented in [`data_analysis.md`](data_analysis.md).
