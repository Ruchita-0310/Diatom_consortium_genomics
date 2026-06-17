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
    └── merge_phaeodactylum_blast_hits.py
