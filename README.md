# Diatom Consortia: Metagenomic and Metatranscriptomic Pipeline

This repository documents the analysis workflow used to assemble, polish, classify, annotate, and compare genomic and transcriptomic data recovered from a diatom dominated microbial consortium enriched from Deer Lake, British Columbia.

The workflow includes long read metagenomic assembly, short read polishing, metagenomic binning, contig level taxonomic screening, organelle identification, marker based phylogenetic analyses, BRAKER4 ET gene prediction, functional annotation, transcript expression integration, nuclear genome refinement, repeat analysis, representative proteome curation, comparative protein orthology analysis, and Hi C read mapping and contact network analysis.

The final comparative genomics analysis uses OrthoFinder to compare the Deer Lake diatom proteome with four reference diatoms:

- *Nitzschia inconspicua*
- *Seminavis robusta*
- *Phaeodactylum tricornutum*
- *Thalassiosira pseudonana*

Earlier whole genome BLASTN comparisons with *P. tricornutum* and *T. pseudonana* are retained as exploratory nucleotide similarity analyses but are not used as the primary orthology analysis.

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
18S rRNA, plastid 16S rRNA, and rbcL phylogenetic analyses
   ↓
BRAKER4 ET genome annotation
   ↓
Functional annotation with Swiss-Prot, Bacillariophyta UniProtKB,
InterProScan, and AntiFam
   ↓
Expression integration using TransDecoder ORFs and Average_TPM
   ↓
Nuclear enriched genome generation
   ↓
One representative Deer Lake protein retained per nuclear gene
   ↓
RepeatModeler and RepeatMasker analysis
   ↓
Removal of high confidence TE associated protein models
   ↓
Final Deer Lake nuclear representative proteome
   ↓
Reference proteome preparation for NI, SR, PT, and TP
   ↓
Five proteome comparison with OrthoFinder
   ↓
Deer Lake gene level orthogroup sharing with NI, SR, PT, and TP
   ↓
Final comparative gene table with functional annotation,
Average_TPM, and four reference comparison fields
   ↓
Hi C read mapping and contig contact network
```

---

## Main outputs

### BRAKER4 gene annotation

BRAKER4 ET generated:

```text
15,102 genes
16,947 predicted protein isoforms
```

The accepted BRAKER4 protein set is:

```text
metatranscriptomics/BRAKER4/final_annotation_ET/DL_diatom.braker4.ET.proteins.faa
```

Functional annotation and transcript expression were initially integrated at the BRAKER4 isoform level.

---

### Initial functional annotation and expression table

The original clean BRAKER4 isoform table is:

```text
metatranscriptomics/transdecoder_to_braker_ID_bridge/
CLEAN_REBUILD_FROM_RAW/09_final/
DL_diatom_FINAL_clean_BRAKER_isoform_table.tsv
```

This table contains one row per BRAKER4 predicted protein isoform and integrates:

```text
functional annotation
BRAKER4 coordinate and length information
TransDecoder ORF mapping
Average_TPM
compartment assignment
AntiFam annotation flags
```

The table contains:

```text
16,947 protein isoforms
```

These isoform level tables are retained as intermediate annotation resources.

---

## Nuclear representative proteome

Comparative protein analysis was performed using the nuclear enriched Deer Lake genome rather than the original whole diatom assembly.

The nuclear enriched assembly contains:

```text
3,007 contigs
81,911,772 bp
```

One representative protein was retained per Deer Lake gene. Alternative BRAKER4 isoforms were reduced to a single representative sequence before comparative analysis.

Protein models shorter than 50 amino acids were removed only when they lacked functional annotation and transcript expression support.

Two AntiFam associated gene models were excluded during the initial representative proteome curation:

```text
g10893
g11404
```

The resulting representative proteome initially contained:

```text
15,017 proteins
```

---

## Repeat analysis and TE associated gene filtering

RepeatModeler and RepeatMasker were used to characterize repetitive sequence content in the Deer Lake nuclear assembly.

RepeatMasker identified:

```text
Genome size:                    81,911,772 bp
Total bases masked:            45,345,998 bp
Total masked fraction:             55.36%
Interspersed repeats:              54.55%
LTR retroelements:                 36.68%
Ty1/Copia elements:                28.09%
DNA transposons:                    7.81%
Unclassified repeats:               9.74%
```

Representative CDS coordinates were intersected with RepeatMasker annotations. Repeat hits overlapping the same genomic interval were resolved before calculating CDS overlap to avoid double counting.

High confidence TE associated protein models were defined conservatively using both substantial CDS overlap with interspersed repeats and explicit TE related functional annotation.

A total of:

```text
38 protein models
```

were classified as high confidence TE associated models and removed from the comparative proteome.

The final Deer Lake proteome used for OrthoFinder therefore contains:

```text
14,979 nuclear representative proteins
```

Final FASTA:

```text
comparative_genomics/orthofinder_input/DeerLake_Nitzschia.faa
```

---

## Reference diatom proteomes

Four reference diatom proteomes were included in the OrthoFinder comparison.

| Species | Abbreviation | Proteins used |
| --- | --- | ---: |
| Deer Lake diatom | DL | 14,979 |
| *Nitzschia inconspicua* | NI | 38,601 |
| *Seminavis robusta* | SR | 35,995 |
| *Phaeodactylum tricornutum* | PT | 10,392 |
| *Thalassiosira pseudonana* | TP | 11,672 |

Reference proteomes were standardized before analysis.

For *P. tricornutum* and *T. pseudonana*, one representative protein was retained per annotated locus.

For *N. inconspicua*, organelle encoded proteins were removed before OrthoFinder analysis. This excluded:

```text
150 plastid encoded proteins
34 mitochondrial encoded proteins
```

The final *N. inconspicua* nuclear protein set therefore contained:

```text
38,601 proteins
```

The available *N. inconspicua* annotation represents a diploid genome. Consequently, the analysis uses *N. inconspicua* primarily for orthogroup sharing rather than interpretation of gene copy number or gene family expansion.

Final OrthoFinder input directory:

```text
/work/ebg_lab/eb/diatom_consortia/comparative_genomics/orthofinder_input/
```

containing:

```text
DeerLake_Nitzschia.faa
Nitzschia_inconspicua.faa
Seminavis_robusta.faa
Phaeodactylum_tricornutum.faa
Thalassiosira_pseudonana.faa
```

---

## OrthoFinder comparative genomics

Protein orthology was inferred using OrthoFinder v3.1.5.

The analysis used:

```text
DIAMOND for protein similarity searches
FAMSA for multiple sequence alignment
FastTree for gene tree inference
```

The comparison included:

```text
111,639 proteins
5 diatom proteomes
```

OrthoFinder assigned:

```text
99,813 proteins to orthogroups
16,157 orthogroups
89.4% of proteins assigned
4,914 orthogroups represented in all five proteomes
```

The main OrthoFinder output directory is:

```text
/work/ebg_lab/eb/diatom_consortia/comparative_genomics/
orthofinder_results/DL_5species_orthofinder_v3/Results_Sep09/
```

For each Deer Lake gene, the corresponding orthogroup was examined for the presence of proteins from NI, SR, PT, and TP.

A reference species was recorded as `yes` when the Deer Lake protein belonged to an OrthoFinder orthogroup containing at least one protein from that species.

Therefore, these fields describe **orthogroup sharing** and should not be interpreted as literal genome wide gene presence or absence.

---

## Deer Lake genes shared with reference diatoms

The final Deer Lake proteome contains:

```text
14,979 genes
```

The numbers of Deer Lake genes belonging to orthogroups containing each reference species are:

| Reference diatom | Deer Lake genes in shared orthogroups |
| --- | ---: |
| *Nitzschia inconspicua* | 10,163 |
| *Seminavis robusta* | 9,952 |
| *Phaeodactylum tricornutum* | 7,973 |
| *Thalassiosira pseudonana* | 7,569 |

These categories are not mutually exclusive because a Deer Lake gene can belong to an orthogroup containing several reference species.

---

## Final comparative gene table
The main gene level table for downstream comparative analysis is:
```text
Orthogroups/DL_diatom_FINAL_gene_table_OrthoFinder_PT_TP_NI_SR.tsv
```
The table contains:
```text
14,979 Deer Lake nuclear representative genes
14,980 lines including the header
```
Columns:
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
Column interpretation:
```text
gene_id
  Representative BRAKER4 transcript ID for the Deer Lake gene.

contig_id
  Nuclear Deer Lake contig containing the gene model.

diatom_compartment
  Genomic compartment assignment.

diatom_gene_length_bp
  Gene length calculated from BRAKER4 coordinates.

functional_annotation
  Functional annotation derived from the integrated Swiss-Prot,
  Bacillariophyta UniProtKB, InterProScan, and AntiFam workflow.

diatom_Average_TPM
  Average transcript expression value transferred through the
  selected TransDecoder ORF to BRAKER4 mapping.

present_in_Thalassiosira_pseudonana
  yes when the Deer Lake protein belongs to an orthogroup containing
  at least one T. pseudonana protein.

present_in_Phaeodactylum_tricornutum
  yes when the Deer Lake protein belongs to an orthogroup containing
  at least one P. tricornutum protein.

present_in_Nitzschia_inconspicua
  yes when the Deer Lake protein belongs to an orthogroup containing
  at least one N. inconspicua protein.

present_in_Seminavis_robusta
  yes when the Deer Lake protein belongs to an orthogroup containing
  at least one S. robusta protein.
```
Blank `diatom_Average_TPM` fields are retained as missing expression values rather than being converted to zero.

This table is the primary Deer Lake comparative gene table used for pathway inspection and downstream biological interpretation.

---
## Hi C analysis

Hi C reads were mapped to the polished whole assembly to assess contig representation and proximity ligation links among assembled contigs.

The polished whole metagenomic assembly was used rather than the nuclear enriched subset because the Hi C library was generated from the complete diatom associated consortium.

The YaHS scaffolding analysis was treated as exploratory. The final Hi C results are interpreted primarily as contig level contact evidence rather than as a chromosome scale scaffolded genome.

Final network outputs include:

```text
hic_contig_network_all_primary_pairs.gexf
hic_contig_network_all_primary_pairs.graphml
```

These files can be opened in Gephi for visualization and analysis of contacts among nuclear, organelle associated, bacterial, and mixed contigs.

---

## Custom Python scripts

Custom Python scripts are stored in the `scripts/` directory. Shell commands and SLURM workflows are documented directly in `data_analysis.md`.

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
└── 22_make_final_orthofinder_gene_table.py
```

The scripts relevant to the final comparative analysis are:

```text
18_make_DL_nuclear_representative_proteome.py
  Generates one representative nuclear Deer Lake protein per gene.

19_calculate_DL_CDS_repeat_overlap.py
  Calculates overlap between representative Deer Lake CDS regions
  and RepeatMasker annotations.

20_classify_DL_TE_candidates.py
  Identifies high confidence TE associated protein models using
  repeat overlap and functional annotation.

21_prepare_reference_proteome.py
  Standardizes reference protein sets for comparative analysis.

22_make_final_orthofinder_gene_table.py
  Combines Deer Lake gene metadata, Average_TPM, and OrthoFinder
  orthogroup membership to generate the final 14,979 gene table
  containing PT, TP, NI, and SR comparison fields.
```

---

## Interpretation of comparative fields
The four final comparison columns are based on OrthoFinder protein orthogroups.
For example:

```text
present_in_Nitzschia_inconspicua = yes
```
means that the Deer Lake protein belongs to an orthogroup containing at least one *N. inconspicua* protein.

It does not necessarily indicate a one to one orthologue, identical gene function, or identical copy number.

Similarly:
```text
no
```
means that no protein from that reference proteome was assigned to the same orthogroup. It should not be interpreted as definitive biological absence because genome assembly, annotation quality, sequence divergence, and proteome completeness can affect orthogroup recovery.

---

## Detailed workflow
The complete workflow, including commands, SLURM scripts, software environments, input files, intermediate outputs, quality control steps, and custom Python scripts, is documented in:

[`data_analysis.md`](data_analysis.md)
