# Diatom Consortia: Metagenomic and Metatranscriptomic Pipeline

This repository documents the analysis workflow used to assemble, polish, classify, annotate, and compare genomic and transcriptomic data recovered from a diatom dominated microbial consortium enriched from Deer Lake, British Columbia.

The workflow combines long read metagenomic assembly, short read polishing, metagenomic binning, contig level taxonomic screening, organelle identification, marker based phylogenetic analyses, BRAKER4 ET gene prediction, functional annotation, transcript expression integration, nuclear genome refinement, repeat analysis, representative gene curation, four genome nucleotide comparison, and Hi C contact analysis.

The current comparative nucleotide analysis uses a curated set of **14,941 Deer Lake nuclear genes** and compares them with four reference diatom genomes using `dc-megablast`:

- *Nitzschia inconspicua* (NI)
- *Seminavis robusta* (SR)
- *Phaeodactylum tricornutum* (PT)
- *Thalassiosira pseudonana* (TP)

A protein level OrthoFinder comparison is retained as a secondary exploratory analysis. The BLASTN analysis is the current gene level comparison used for the manuscript figure and expression summaries.

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
Functional annotation and Average_TPM integration
   ↓
Nuclear enriched genome generation
   ↓
One representative Deer Lake gene retained per nuclear locus
   ↓
RepeatModeler and RepeatMasker analysis
   ↓
Removal of high confidence TE associated models
   ↓
Plastid derived contig quality control
   ↓
Final Deer Lake BLASTN query: 14,941 genes
   ↓
Four genome dc-megablast comparison: NI, SR, PT, TP
   ↓
Master 14,941 gene table with BLASTN metrics and Average_TPM
   ↓
Expression comparison and BLASTN unmatched gene subsets
   ↓
Secondary five proteome OrthoFinder analysis
   ↓
Hi C mapping to the polished whole assembly
   ↓
High confidence contig contact table
   ↓
Whole assembly network visualization
```

---

## Main final outputs

### Deer Lake nuclear gene query

The final nucleotide query contains:

```text
14,941 Deer Lake nuclear representative genes
```

Final query FASTA:

```text
comparative_genomics/blastn_redo/00_inputs/
DL_final_nuclear_genes_repeat_TE_plastid_clean.fasta
```

The corresponding final protein set after repeat and plastid quality control is:

```text
comparative_genomics/00_DL_reference/
DL_diatom_FINAL_nuclear_representative_proteome_TEfiltered_plastidQC.faa
```

The query was derived after:

1. retaining one representative nuclear BRAKER4 gene per locus;
2. removing two AntiFam associated models during representative proteome curation;
3. removing 38 high confidence TE associated models using repeat overlap plus annotation evidence; and
4. removing 38 protein models located on four clearly plastid derived contigs (`contig_475`, `contig_4813`, `contig_5686`, and `contig_5702`).

No mitochondrial based nuclear filtering was applied because the candidate mitochondrial assembly was not considered sufficiently reliable for that purpose.

---

### Four genome BLASTN comparison

Each of the 14,941 Deer Lake genes was searched against NI, SR, PT, and TP using `dc-megablast` with:

```text
E value <= 1e-10
```

No additional identity, coverage, or gene length cutoff was used to remove genes from the master table.

Raw alignment counts were:

| Comparator | Raw BLASTN alignments |
| --- | ---: |
| NI | 27,211 |
| SR | 12,168 |
| PT | 7,915 |
| TP | 4,817 |

The final master table is:

```text
DL_14941_genes_BLASTN_NI_SR_PT_TP_with_TPM_compartment.tsv
```

For each comparator, the table retains the highest bitscore BLASTN alignment for each Deer Lake query gene together with percent identity, alignment length, query coverage, E value, and bitscore. A `yes` hit means that at least one BLASTN alignment was returned under the search criterion; it is not a statement of definitive orthology or biological presence.

Observed gene level hit counts:

| Comparison | Deer Lake genes with a BLASTN hit |
| --- | ---: |
| NI | 5,767 |
| SR | 4,273 |
| PT | 3,799 |
| TP | 2,636 |
| Hit in at least one comparator | 6,429 |
| No hit in any comparator | 8,512 |

The comparator categories overlap because one Deer Lake gene can have detectable nucleotide similarity to more than one reference genome.

---

### Expression integration

Average_TPM values were carried into the master BLASTN table from the curated Deer Lake metadata map.

```text
Genes with a TPM value:       11,092
Genes with TPM > 0:           10,907
Genes with TPM = 0:              185
```

Among genes with no BLASTN hit in any of the four reference genomes:

```text
No comparator hit:                         8,512
No comparator hit + TPM > 0:              5,332
No comparator hit + TPM > 0 +
informative functional annotation:         3,259
```

The term **Deer Lake only** is used as a compact figure label for the first group. In text and figure legends it should be defined as genes with no detectable BLASTN hit to NI, SR, PT, or TP under the search criteria used. These genes should not be described as proven lineage specific or unique genes.

---

### Expression figure

`scripts/26_plot_BLASTN_TPM_panelA.py` creates the violin plot comparing:

```text
Deer Lake only
Shared with N. inconspicua
Shared with S. robusta
Shared with P. tricornutum
Shared with T. pseudonana
```

The plot uses `log10(Average TPM + 1)`, retains true TPM values of zero, excludes only missing TPM values, and uses a color blind aware palette. The four shared groups are not mutually exclusive.

---

### Hi C whole assembly network

The polished whole assembly contains:

```text
4,925 contigs
```

The high confidence network plotting workflow starts from:

```text
hic_bwa_separate_reads/03_tables/
HiC_contig_pair_contacts_MAPQ30_PID95.tsv
whole_assembly_contig_type_map.tsv
```

Self contacts were excluded. Oriented contig pair rows supported by at least two Hi C read pairs were retained and then collapsed to unique undirected contig pairs.

Final high confidence network summary:

```text
Unique undirected edges:       471
Connected contigs:             487
Isolated assembly contigs:   4,438
Minimum edge support:            2 read pairs
Maximum edge support:          194 read pairs
Median edge support:             2 read pairs
```

`scripts/28_plot_HiC_whole_assembly_network.py` plots all 4,925 contigs. Connected contigs are shown in the central Hi C network, while contigs without a retained inter contig edge are shown as small isolated nodes around the network.

The figure highlights:

```text
Plastid genome associated contigs:
  contig_1443
  contig_4315

Candidate mitochondrial contigs:
  contig_5628
  contig_1647
```

The mitochondrial contigs are labelled as **candidate mitochondrial contigs** because their assignment requires additional validation. Hi C contact structure is treated as supporting context rather than a standalone taxonomic classifier.

---

## Repeat analysis

RepeatModeler and RepeatMasker were used to characterize the Deer Lake nuclear enriched assembly.

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

Repeat overlap alone was not used to delete genes. Only 38 models with substantial repeat overlap and explicit TE related annotation were removed as high confidence TE associated models.

---

## Secondary OrthoFinder analysis

A five proteome OrthoFinder comparison was also completed using the plastid quality controlled Deer Lake protein set. This analysis is retained as a protein level complement to the BLASTN analysis rather than as the main gene presence or absence result.

Corrected Deer Lake OrthoFinder input:

```text
14,941 proteins
```

Corrected Deer Lake summary:

```text
Assigned to orthogroups:  12,677  (84.8%)
Unassigned:                2,264  (15.2%)
DL containing orthogroups: 8,970
```

Deer Lake genes in orthogroups containing each comparator:

| Comparator | Deer Lake genes |
| --- | ---: |
| NI | 10,154 |
| SR | 9,952 |
| PT | 7,966 |
| TP | 7,542 |

These values represent orthogroup sharing. They should not be interpreted as one to one orthology, definitive biological absence, or gene family expansion.

---

## Repository scripts

The original repository contains scripts `01` through `22`. This update adds:

```text
scripts/
├── 23_prepare_DL_BLASTN_query.py
├── 24_build_BLASTN_master_table.py
├── 25_make_BLASTN_subsets.py
├── 26_plot_BLASTN_TPM_panelA.py
├── 27_make_HiC_undirected_min2_edges.py
└── 28_plot_HiC_whole_assembly_network.py
```

SLURM workflows are stored as plain text files in:

```text
slurm/
└── 23_run_four_genome_dcmegablast_SLURM.txt
```

The complete command history, paths, filtering logic, output counts, and interpretation notes are documented in [`data_analysis.md`](data_analysis.md).

---

## Software used in the added analyses

| Analysis | Software |
| --- | --- |
| Four genome nucleotide comparison | BLAST+, Python, pandas |
| Expression figure | pandas, NumPy, SciPy, Matplotlib |
| Hi C edge processing | Python, pandas |
| Hi C network visualization | pandas, NumPy, NetworkX, Matplotlib |

Laptop plotting dependencies can be installed with:

```bash
pip install -r requirements_laptop.txt
```

---

## Interpretation notes

- `Deer Lake only` is a figure label for genes with no detectable BLASTN hit to the four selected reference genomes under the stated search criteria.
- BLASTN non detection is not proof that a homolog is biologically absent.
- The four comparator hit groups overlap.
- Blank Average_TPM values are missing data and are not converted to zero.
- True TPM values of zero are retained and can be plotted as `log10(TPM + 1) = 0`.
- Hi C edges represent proximity ligation support under the filtering criteria used and should not be interpreted as direct ecological interaction or mutualism.
- Organelle contigs are separate DNA molecules and are not expected to scaffold automatically into the nuclear contact network.
- Candidate mitochondrial contig assignments remain provisional.
