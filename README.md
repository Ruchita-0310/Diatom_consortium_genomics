# Diatom Consortia: Metagenomic and Metatranscriptomic Pipeline

This repository documents the analysis workflow used to assemble, polish, classify, annotate, and compare genomes and transcriptomes recovered from a diatom-associated microbial consortium.

The workflow includes long-read metagenomic assembly, short-read polishing, metagenomic binning, contig-level taxonomic screening, organelle identification, BRAKER4 ET gene prediction, functional annotation, expression integration, nucleotide-level comparisons with *Phaeodactylum tricornutum* and *Thalassiosira pseudonana*, and Hi-C read mapping and contact-network analysis.

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
Functional annotation with Swiss-Prot, Bacillariophyta UniProtKB,
InterProScan, and AntiFam
   ↓
Expression integration using TransDecoder ORFs and Average_TPM
   ↓
Phaeodactylum tricornutum nucleotide comparison using BLASTN
and gene-linked overlap analysis
   ↓
Thalassiosira pseudonana nucleotide comparison using BLASTN
and gene-linked overlap analysis
   ↓
Final clean BRAKER4 isoform-level gene table
   ↓
Simplified eight-column gene table for manual review
   ↓
Hi-C read mapping and contig-contact network
```

---
## Main outputs
### Initial clean BRAKER4 isoform table
The initial clean BRAKER4 isoform-level table is:

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
A version sorted by expression is also generated:

```text
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table_sorted_by_Average_TPM.tsv
```
### Final table with both reference comparisons
After the *Thalassiosira pseudonana* comparison is integrated, the updated final table is:

```text
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table_PT_TP.tsv
```
A corresponding version sorted by `Average_TPM` is:
```text
09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table_PT_TP_sorted_by_Average_TPM.tsv
```

The updated table retains all functional annotation, coordinate, expression, compartment, and AntiFam fields and includes two separate nucleotide-similarity columns:
```text
present_in_Phaeodactylum_tricornutum
present_in_Thalassiosira_pseudonana
```
These columns indicate whether a Deer Lake BRAKER4 gene overlapped at least one retained genome-level BLASTN alignment associated with a gene model in the corresponding reference genome.

The comparison fields represent nucleotide-level similarity screens and should not be interpreted as evidence of confirmed orthology.

---
## Simplified review table
A simplified table was generated for manual inspection and pathway review:
```text
09_final/DL_diatom_FINAL_gene_table_for_boss_PT_TP.tsv
```
This table contains eight columns:
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
The two reference-comparison columns are derived independently from the corresponding BLASTN gene-linked overlap tables.

---
## Reference genome comparisons

### *Phaeodactylum tricornutum*
The *Phaeodactylum tricornutum* genome was compared with the Deer Lake diatom assembly using BLASTN.

The reference genome was used as the query, and the Deer Lake diatom assembly was used as the BLAST database. BLASTN alignments were converted to BED coordinates and intersected with gene models from both genomes.

The resulting table links each retained alignment to:

```text
Phaeodactylum tricornutum gene models
Deer Lake BRAKER4 gene models
alignment identity
alignment length
alignment coordinates
e-value
bit score
```

The Deer Lake gene-level result was summarized as:

```text
present_in_Phaeodactylum_tricornutum
```
### *Thalassiosira pseudonana*

The *Thalassiosira pseudonana* CCMP1335 reference genome was compared with the Deer Lake diatom assembly using the same nucleotide-level workflow.

The reference genome and annotation files were stored under the user home directory to avoid duplicating reference files in the project working directory:

```text
$HOME/databases/thalassiosira_pseudonana/GCF_000149405.2_ASM14940v2/
```
Symbolic links to these files were created in the analysis input directory.

The BLASTN alignments were linked to overlapping *Thalassiosira pseudonana* and Deer Lake BRAKER4 gene models. The Deer Lake gene-level result was summarized as:
```text
present_in_Thalassiosira_pseudonana
```
The main gene-linked comparison output is:
```text
thalassiosira_to_diatom_blastn_redo/04_summary/thalassiosira_vs_diatom_BLASTN_with_Thalassiosira_and_BRAKER_ET_genes.tsv
```
---
## Hi-C analysis
Hi-C reads were mapped to the polished whole assembly to assess contig-level representation and inter-contig proximity-ligation links.

The YaHS scaffolding test was treated as exploratory. The final Hi-C output is interpreted as contig-level contact evidence rather than as a chromosome-scale scaffolded assembly.

Final Hi-C network outputs:
```text
hic_contig_network_all_primary_pairs.gexf
hic_contig_network_all_primary_pairs.graphml
```
The network files can be opened in Gephi for visualization and analysis of nuclear, chloroplast-associated, mitochondrial-associated, and mixed contig contacts.

---
## Custom scripts

Custom Python scripts are stored in the `scripts/` directory and numbered according to their order in the overall analysis:

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
└── 17_make_hic_simple_mixed_read_table.py
```
---
## Detailed workflow
The full workflow, including commands, software environments, intermediate files, outputs, and short explanations of each analysis step, is documented in:

[`data_analysis.md`](data_analysis.md)
