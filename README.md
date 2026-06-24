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
Phaeodactylum tricornutum comparison
   ↓
Final clean BRAKER4 isoform-level gene table
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
    ├── 09_add_ONLY_Average_TPM_clean.py
    ├── 10_make_FINAL_clean_BRAKER_isoform_table.py
    └── 11_make_hic_network_files.py
```
---

## Main output
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

---

## Scripts
Custom Python scripts are stored in the `scripts/` directory and are numbered in the order they are used in the workflow.

```text
01_classify_metaeuk_contigs.py
  Classifies contigs using MetaEuk ORF-level taxonomy.

02_make_swissprot_best_hits.py
  Parses Swiss-Prot DIAMOND output and writes best-hit annotation tables.

03_make_bacillariophyta_best_hits.py
  Parses UniProtKB Bacillariophyta DIAMOND output using the same confidence framework as Swiss-Prot.

04_summarize_interproscan.py
  Collapses raw InterProScan TSV output to one row per predicted protein.

05_merge_functional_annotation_layers.py
  Merges Swiss-Prot, Bacillariophyta, InterProScan, and AntiFam evidence.

06_merge_phaeodactylum_blast_hits.py
  Adds overlapping Phaeodactylum gene information to the cleaned BLASTN best-hit table.

07_add_BRAKER_lengths_clean.py
  Adds BRAKER4 coordinates, contig IDs, gene length, CDS length, and protein length.

08_make_best_ORF_to_BRAKER_mapping_clean.py
  Keeps one best BRAKER4 hit per TransDecoder ORF.

09_add_ONLY_Average_TPM_clean.py
  Adds TransDecoder ORF IDs and Average_TPM values to BRAKER4 isoforms.

10_make_FINAL_clean_BRAKER_isoform_table.py
  Creates the final clean BRAKER4 isoform-level table.

11_make_hic_network_files.py
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
The full workflow, including commands, outputs, and short notes describing what each command does, is documented in [`data_analysis.md`](data_analysis.md)
