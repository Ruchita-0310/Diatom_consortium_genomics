# Diatom Consortia: Metagenomic and Metatranscriptomic Pipeline

This repository documents the analyses of a diatom dominated microbial consortium enriched from Deer Lake, an alkaline soda lake in British Columbia, Canada. The dominant organism is a pennate diatom of the genus *Nitzschia*.

The workflow reconstructs the diatom genome from a mixed community, annotates and curates its genes, measures their expression, compares them with other diatoms and with cultured eukaryote transcriptomes, maps physical contacts within the consortium with Hi C, and follows the diatom 18S rRNA signal in seasonal samples from the lake.

Full commands, SLURM scripts, settings, file paths and output counts are documented in [`data_analysis.md`](data_analysis.md). Custom Python scripts are in [`scripts/`](scripts/), numbered by their role in the analysis.

---

## Workflow overview

```text
Nanopore + Illumina sequencing of the consortium
   ↓
Flye metagenome assembly
   ↓
Medaka, Polypolish and Pypolca polishing
   ↓
MetaBAT2 binning, CheckM2, GTDB-Tk and MetaEuk classification
   ↓
Organelle identification and nuclear enriched diatom genome
   ↓
18S rRNA, plastid 16S rRNA and rbcL phylogenies
   ↓
BRAKER4 ET gene prediction and functional annotation
   ↓
Metatranscriptome assembly and expression integration
   ↓
Repeat analysis and curation of the nuclear gene set
   ↓
Comparison with four reference diatom genomes
   ↓
Reciprocal best hit analysis of expressed ORFs
   ↓
Expression comparison with MMETSP culture transcriptomes
   ↓
Hi C contact network of the consortium
   ↓
Seasonal 18S rRNA signal and environmental 18S consensus
```

---

## Genome assembly and annotation

Nanopore reads were assembled with Flye in metagenome mode and polished with Medaka (long reads) and Polypolish and Pypolca (short reads). Contigs were binned with MetaBAT2, and bacterial bins were classified with GTDB-Tk and assessed with CheckM2. Because the consortium contains a eukaryote, contigs were also classified with MetaEuk, counting organelle derived hits toward the diatom.

The plastid genome was recovered as one 120 kb contig. Two candidate mitochondrial contigs were identified but remain provisional. After removal of organelle contigs, the nuclear enriched diatom assembly comprised 81.9 Mb in 3,007 contigs and recovered 86.1% of conserved stramenopile BUSCO genes.

Nuclear genes were predicted with BRAKER4 in ET mode using RNA-seq evidence (15,102 genes, 16,947 proteins) and annotated with Swiss-Prot, diatom UniProtKB entries and InterProScan.

---

## Phylogeny

The diatom was placed with three markers: nuclear 18S rRNA from the metatranscriptome, and plastid 16S rRNA and *rbcL* from the plastid genome. Alignments were made with Clustal Omega, trimmed with TrimAl and analysed with IQ-TREE 2.

---

## Repeat analysis and gene set curation

RepeatModeler and RepeatMasker masked 55% of the nuclear assembly, mostly LTR retrotransposons of the Ty1/Copia type. One representative gene was kept per locus, and models removed from the final set were:

- two AntiFam flagged models;
- 38 high confidence transposable element models;
- 38 models on four plastid derived contigs.

The final curated nuclear gene set contains **14,941 genes**.

---

## Comparison with reference diatoms

The 14,941 genes were compared with four reference diatoms:

- *Nitzschia inconspicua* (NI)
- *Seminavis robusta* (SR)
- *Phaeodactylum tricornutum* (PT)
- *Thalassiosira pseudonana* (TP)

**Nucleotide comparison.** A `dc-megablast` search (E ≤ 1e-10) found matches for 6,429 genes in at least one reference. 8,512 genes had no detected match and are labelled "Deer Lake only" in figures.

**Protein orthology.** OrthoFinder assigned 84.8% of Deer Lake proteins to orthogroups. Orthogroup sharing was highest with *N. inconspicua* and *S. robusta*.

**Reciprocal best hits.** The 87,621 expressed ORFs of the metatranscriptome were compared with the four reference proteomes by forward and reverse BLASTP. Hits were integrated with eggNOG, KOfam and Pfam annotation and with the expression percentile of each ORF, and candidates in six functional systems were curated manually.

---

## Expression comparison with MMETSP cultures

Alkaline relevant, fermentation and DUF containing ORFs, plus housekeeping genes, were searched with DIAMOND against 645 MMETSP culture transcriptomes. Each Deer Lake ORF and its best ortholog in each culture library were ranked within their own transcriptome. A transcript was called higher than cultures when it exceeded the 90th percentile of culture libraries both in rank and relative to housekeeping genes.

Enrichment of higher transcripts was tested per group and per module with one sided Fisher exact tests. Results are shown as transcript level and module level dumbbell figures.

---

## Hi C contact network

Hi C reads were mapped to the polished whole consortium assembly (4,925 contigs). High confidence read pairs (MAPQ ≥ 30, ≥ 95% identity) were used to build a contig contact network. After removing self contacts and keeping pairs supported by at least two read pairs, the network contains 471 edges among 487 contigs, including contacts between diatom and bacterial contigs. Plastid and candidate mitochondrial contigs are highlighted in the network figure.

---

## Seasonal 18S signal and environmental consensus

Eukaryotic rRNA reads from 35 seasonal mat and sediment libraries, including sediment incubated in darkness for 3 and 9 months, were mapped to the Deer Lake diatom 18S sequence with BBMap (≥ 98% identity). The diatom signal was expressed as a percentage of eukaryotic rRNA reads per library.

A full length environmental 18S consensus was reconstructed from the best covered untreated library (fall sediment) with `samtools consensus` and added to the 18S phylogeny.

---

## Software

| Analysis | Main tools |
| --- | --- |
| Assembly, polishing and binning | Flye, Medaka, Polypolish, Pypolca, MetaBAT2, CheckM2, CoverM |
| Taxonomy and organelles | GTDB-Tk, MetaEuk, MetaQUAST, fastANI, minimap2 |
| Phylogenetics | Barrnap, Clustal Omega, TrimAl, IQ-TREE 2 |
| Gene prediction and annotation | BRAKER4, STAR, DIAMOND, InterProScan, eggNOG mapper, KOfam, HMMER/Pfam |
| Transcriptomics | nf-core/metatdenovo, SPAdes, TransDecoder |
| Repeats and comparative genomics | RepeatModeler, RepeatMasker, BLAST+, OrthoFinder |
| MMETSP comparison | DIAMOND, seqkit |
| Hi C | BWA-MEM, samtools, YaHS, NetworkX |
| Seasonal 18S | BBMap, samtools, SILVA SSU 138.2 |
| Analysis and figures | Python, pandas, NumPy, SciPy, Matplotlib |
