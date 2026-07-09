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
