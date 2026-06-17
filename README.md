# Diatom Consortia: Metagenomic and Metatranscriptomic Pipeline
This repository contains the workflow used to assemble, polish, classify, annotate, and analyze genomes and transcriptomes from a diatom-associated microbial consortium.
The main workflow is documented in:
[`data_analysis.md`](data_analysis.md)
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
    ├── merge_phaeodactylum_blast_hits.py
    ├── 01_make_clean_transdecoder_to_braker_bridge.py
    ├── 02_merge_TPM_to_BRAKER_ALL_HITS.py
    ├── 03_merge_BRAKER_functional_annotation_with_ALL_hit_TPM.py
    ├── 04_add_BRAKER_gene_lengths.py
    ├── 05_add_ALL_phaeodactylum_hits.py
    ├── 06_make_clean_final_nuclear_gene_table.py
    ├── 07b_make_organelle_tables_and_combine_tolerant_GBK.py
    └── 08_add_organelle_curation_flags.py
```
## Workflow summary
The workflow includes:
```text
1. Long-read metagenomic assembly
2. Long-read and short-read polishing
3. Assembly quality assessment
4. Metagenomic binning and bin quality assessment
5. Taxonomic classification of bacterial, eukaryotic, and organelle-associated contigs
6. Organelle genome identification and annotation
7. Transcriptome assembly and ORF prediction
8. BRAKER4 ET gene prediction using RNA-seq evidence
9. Functional annotation of predicted diatom proteins
10. Expression integration using TransDecoder ORF-to-BRAKER protein mappings
11. Nuclear, plastid, and mitochondrial gene table construction
12. Comparison with Phaeodactylum tricornutum
```
## Main outputs
The final gene-level annotation and expression table is:
```text
DL_diatom_final_gene_table_nuclear_plastid_mito_with_curation_flags.tsv
```
This table combines:
```text
BRAKER4 nuclear gene predictions
Swiss-Prot annotations
UniProtKB Bacillariophyta annotations
InterProScan domain and GO annotations
AntiFam warning flags
Average TPM values from transcriptome ORFs
All Phaeodactylum tricornutum comparison hits
Plastid genes from GenBank annotation
Mitochondrial genes from GenBank annotation
Organelle curation flags
```
A suggested pathway-curation version is also generated:
```text
DL_diatom_final_gene_table_nuclear_plastid_mito_curated_organelle_suggested.tsv
```
## Notes
The repository is intended to document the analysis workflow and custom scripts used for the diatom consortium project. Large input files, intermediate alignment files, databases, and raw sequencing data are not stored in this repository.
