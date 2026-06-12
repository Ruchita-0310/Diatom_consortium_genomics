# Diatom Consortia: Metagenomic & Metatranscriptomic Pipeline
This repository contains the end-to-end workflow for the assembly, polishing, binning, and annotation of Diatom-associated microbial consortia.                            
🛠 Prerequisites & Environment
The following software environments are required. It is recommended to manage these via Conda or Singularity as noted:
1. Assembly/Polishing: Flye, Medaka, Polypolish, Pypolca
2. Quality & Validation: BUSCO, CheckM2, QUAST/MetaQUAST
3. Binning & Taxonomy: MetaBAT2, GTDB-Tk, CoverM
4. Phylogenetics: ClustalO, TrimAl, IQ-TREE 2, Biopython
5. Annotation: RepeatModeler2, RepeatMasker, STAR, BRAKER4, TSEBRA, StringTie
6. Metatranscriptomics: Nextflow, nf-core/metatdenovo

# 1. Assembly
Basecalled using Guppy
```
flye --nano-raw pass_trim.fastq.gz --meta -g 50m --min-overlap 5000 --out-dir flye_out_new -i 3 --threads 8
```
# 2. Mapping
## Mapping Short Reads (SR) to the assembled reads
```
minimap2 -ax sr /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/guppy_flye_assembly.fasta Diatoms_merged.fastq.gz > sr_alignment.sam
samtools view -S -b sr_alignment.sam > alignment.bam
samtools sort alignment.bam -o alignment_sorted.bam
samtools index alignment_sorted.bam
samtools flagstat alignment_sorted.bam > mapping_stats.txt
samtools idxstats alignment_sorted.bam | sort -k3,3rn > sr_all_nanopore_hits.tsv
samtools depth alignment_sorted.bam > sr_depth.txt
awk '{sum[$1]+=$3; count[$1]++} END {for (c in sum) print c, sum[c]/count[c]}' sr_depth.txt | sort -k2,2nr > sr_mean_depth.tsv
```
# 3. Polishing
## Medaka - LR polising
```
medaka_consensus \
  -i pass_trim.fastq.gz \
  -d guppy_flye_assembly.fasta \
  -o medaka_euk_polished \
  -t 12
```
## Polishing - SR
```
## Mapping ###
bwa mem -t 16 -a /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/medaka_euk_polished/consensus.fasta /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz > alignments_1.sam
bwa mem -t 16 -a /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/medaka_euk_polished/consensus.fasta /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz > alignments_2.sam

### Polypolish filter ###
polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam

### Polypolish filter ###
polypolish polish \
  /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/medaka_euk_polished/consensus.fasta \
  filtered_1.sam filtered_2.sam \
  > sr_poly.fasta

### Pypolca ###
pypolca run -a sr_poly.fasta \
  -1 /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz \
  -2 /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz \
  -t 12 -o sr_pypolca_output --careful
```
Run BUSCO on pyloca_corrected.fasta
```
busco -i pyloca_corrected.fasta -l busco_downloads/lineages/eukaryota_odb10 -o busco_report -m genome
```
# 4. Metabat2
```
conda create -n metabat2_v2 -c conda-forge -c bioconda metabat2 libdeflate=1.10

minimap2 -ax map-ont -t 16 1_sr_pypolca_output/pypolca_corrected.fasta pass_trim.fastq.gz | \
samtools view -@ 16 -bS - | \
samtools sort -@ 16 -m 10G -o aligned_reads.sorted.bam
samtools index -@ 16 aligned_reads.sorted.bam

# Summarize Depth
jgi_summarize_bam_contig_depths --outputDepth depth.txt --percentIdentity 85 aligned_reads.sorted.bam

# Binning
mkdir -p 2_metabat2_bins
metabat2 -i 1_sr_pypolca_output/pypolca_corrected.fasta -a depth.txt -o 2_metabat2_bins/bin -m 1500 -t 16 --unbinned
```
# 5. CheckM2
```
checkm2 predict --threads 16 --input 2_metabat2_bins/ --output_directory 3_checkm2_results
```
# 6. GTDB classification
```
# Assign taxonomy using the Genome Taxonomy Database
gtdbtk classify_wf --genome_dir 2_metabat2_bins/ --out_dir 4_gtdbtk_output --cpus 16 -x fa
```
# 7. CoverM
```
conda create -n coverm -c bioconda -c conda-forge coverm
conda activate coverm

coverm genome \
    --genome-fasta-directory 2_metabat2_bins/bac_bins \
    --genome-fasta-extension fa \
    -1 /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R1_001.fastq.gz \
    -2 /work/ebg_lab/eb/diatom_consortia/sr_diatoms/Li49151-RS-Diatoms-4C_S1_R2_001.fastq.gz \
    --mapper bwa-mem \
    -m mean relative_abundance covered_fraction \
    --threads 8 \
    --min-read-percent-identity 95 \
    -o bac_output_coverm.tsv
```
# 8. Phylogenetic tree
```
cat *.fasta > 18S_new.fasta
clustalo -i 18S_new.fasta -o 18S_aligned.fasta
trimal -in 18S_aligned.fasta -out 18S_trimmed.fasta -automated1
/home/ruchita.solanki/iqtree-2.2.2.7-Linux/bin/iqtree2 -s 18S_trimmed.fasta -m MFP -bb 1000 -alrt 1000 -nt AUTO
```
# 9. Transcriptome analysis 
[Nf core metadenovo](https://github.com/nf-core/metatdenovo)
```
# --- 1. JAVA SETUP ---
module purge
module load java/openjdk-23.0.1
export JAVA_HOME=$(dirname $(dirname $(readlink -f $(which java))))
export PATH=$JAVA_HOME/bin:$PATH

# --- 2. EXECUTION ---
~/nextflow run nf-core/metatdenovo \
    -profile singularity \
    --input samplesheet.csv \
    --outdir /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/new_results \
    -w /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/work \
    --assembler spades \
    --orf_caller transdecoder \
    --eggnog_dbpath /home/ruchita.solanki/eggnog_db \
    --skip_kofam true \
    --hmmfiles /home/ruchita.solanki/Pfam-A.hmm \
    --eukulele_dbpath /home/ruchita.solanki/eukulele_db \
    --eukulele_db mmetsp \
    -resume \
    -with-report report_skipK.html \
    -with-timeline timeline_skipK.html
```

# 10. Identifying rRNA genes from transcriptome
```
barrnap --kingdom euk --threads 4 spades.transcripts.fa --outseq euk_transcript_rRNA.fna > diatom_euk_rRNA.gff
barrnap --kingdom bac spades.transcripts.fa --outseq bac_transcript_rRNA.fna > diatom_bac_rRNA.gff
barrnap --kingdom mito spades.transcripts.fa --outseq mito_transcript_rRNA.fna > diatom_mito_rRNA.gff
```
# 11. Identifying organelle genome
Download mitogenome - MT742552 & chloroplast genome - MT742551
```
conda create -n quast_env quast
metaquast.py 8_diatom.fasta -R /work/ebg_lab/eb/diatom_consortia/organelle/ref/ -o ./8_metaquast_output
metaquast.py /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_output/pypolca_corrected.fasta -R /work/ebg_lab/eb/diatom_consortia/organelle/ref/ -o ./whole_metaquast_output
```
# 12. Diatom Genome Annotation Pipeline
Gene models were generated using the BRAKER4 framework with a genome assembly, RNA-seq alignments, and protein homology evidence. The genome assembly was not soft-masked prior to annotation, so repeat masking was performed internally by BRAKER4 using RED. RNA-seq and protein evidence were supplied through `samples.csv`, which triggered BRAKER4 ETP mode using GeneMark-ETP, AUGUSTUS, and TSEBRA.
```text
Genome FASTA
   ↓
BRAKER4 internal repeat masking with RED
   ↓
RNA-seq alignment evidence: STAR coordinate-sorted BAM
   ↓
Protein homology evidence: OrthoDB Stramenopiles + Phaeodactylum tricornutum
   ↓
BRAKER4 ETP mode: GeneMark-ETP + AUGUSTUS + TSEBRA
   ↓
Gene models: GTF/GFF3
   ↓
Predicted proteins + CDS
   ↓
BUSCO/compleasm assessment
   ↓
Functional annotation and expression quantification
```
## 1. Software Environment and Dependencies
BRAKER4 was used for structural genome annotation. The workflow was run through Snakemake using the official BRAKER4 `Snakefile` and the BRAKER3 Singularity image.
```bash
git clone https://github.com/Gaius-Augustus/BRAKER4.git
cd BRAKER4

singularity pull braker3.sif docker://teambraker/braker3:latest
```
Snakemake was installed in the working environment:
```bash
conda install -c conda-forge -c bioconda snakemake
```
GeneMark requires a license key. The GeneMark key was stored at:
```bash
/home/ruchita.solanki/.gm_key
```
The GeneMark key path was provided in `config.ini`.
## 2. RNA-seq Alignment Evidence
RNA-seq reads were aligned to the diatom genome assembly using STAR. The final alignment file used for BRAKER4 was a coordinate-sorted BAM file:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam
```
The corresponding BAM index was also present:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam.bai
```
The BAM file was checked before BRAKER4 execution:
```bash
ls -lh /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam
ls -lh /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam.bai
```
If the BAM index is missing, it can be generated with:
```bash
samtools index /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam
```
## 3. Protein Evidence Preparation
Protein evidence was prepared from broad stramenopile and diatom-related protein databases rather than from a single reference species. The final evidence set included OrthoDB Stramenopiles proteins and optional *Phaeodactylum tricornutum* proteins.
```bash
mkdir -p /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/proteins
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/proteins
```
### 3.1 Download OrthoDB Stramenopiles Proteins
```bash
wget -c https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Stramenopiles.fa.gz

gzip -t Stramenopiles.fa.gz

gunzip -c Stramenopiles.fa.gz > OrthoDB12_Stramenopiles.fa
```
The OrthoDB Stramenopiles FASTA was checked:
```bash
grep -c "^>" OrthoDB12_Stramenopiles.fa
seqkit stats OrthoDB12_Stramenopiles.fa
```
### 3.2 Download *Phaeodactylum tricornutum* Proteins
```bash
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-63/fasta/phaeodactylum_tricornutum/pep/Phaeodactylum_tricornutum.ASM15095v2.pep.all.fa.gz

gzip -t Phaeodactylum_tricornutum.ASM15095v2.pep.all.fa.gz

gunzip -c Phaeodactylum_tricornutum.ASM15095v2.pep.all.fa.gz \
    > Phaeodactylum_tricornutum.ASM15095v2.pep.all.fa
```
The *Phaeodactylum* FASTA was checked:
```bash
grep -c "^>" Phaeodactylum_tricornutum.ASM15095v2.pep.all.fa
seqkit stats Phaeodactylum_tricornutum.ASM15095v2.pep.all.fa
```
### 3.3 Combine and Clean Protein Evidence
The protein files were combined:
```bash
cat \
    OrthoDB12_Stramenopiles.fa \
    Phaeodactylum_tricornutum.ASM15095v2.pep.all.fa \
    > diatom_protein_evidence.raw.faa
```
The combined protein FASTA was filtered to remove short sequences and exact duplicate protein sequences:
```bash
seqkit seq -m 30 diatom_protein_evidence.raw.faa \
    | seqkit rmdup -s \
    > diatom_protein_evidence.clean.faa
```
The final protein evidence file was checked:
```bash
seqkit stats diatom_protein_evidence.raw.faa diatom_protein_evidence.clean.faa

grep -n "\*" diatom_protein_evidence.clean.faa | head

grep -n -v -E '^>|^[A-Z]+$' diatom_protein_evidence.clean.faa | head
```
No stop codons or malformed sequence lines were detected. The cleaned protein evidence file used for BRAKER4 was:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/proteins/diatom_protein_evidence.clean.faa
```
## 4. Genome Input and Masking Strategy
The genome assembly used for BRAKER4 was:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta
```
The genome was checked for soft masking:
```bash
grep -v "^>" /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta \
    | grep -q '[a-z]' && echo "soft-masked" || echo "not soft-masked"
```
The output was:
```text
not soft-masked
```
Therefore, the `genome_masked` column in `samples.csv` was left empty, and internal repeat masking was enabled in `config.ini` using:
```ini
run_red = True
```
External RepeatModeler and RepeatMasker were not used in this BRAKER4 run.
## 5. BRAKER4 Sample Configuration
A `samples.csv` file was prepared to specify the genome assembly, protein evidence, RNA-seq BAM file, and BUSCO lineage. Both `protein_fasta` and `bam_files` were provided, which triggered ETP mode in BRAKER4.
```bash
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4

nano samples.csv
```
```csv
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage
DL_diatom,/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta,,/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/proteins/diatom_protein_evidence.clean.faa,/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam,,,,,,,,stramenopiles_odb12
```
The `samples.csv` file was checked to confirm the expected number of columns:
```bash
awk -F',' '{print NR, NF}' samples.csv
```
Expected output:
```text
1 13
2 13
```
## 6. BRAKER4 Configuration
The `config.ini` file was edited to specify the Singularity image, GeneMark license key, sample file, repeat masking option, and run parameters.
```bash
nano config.ini
```
```ini
[paths]
braker_container = /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/braker3.sif
genemark_key = /home/ruchita.solanki/.gm_key

[DATA]
samples = samples.csv

[PARAMS]
fungus = False
min_contig = 1000
run_red = True
species = diatom_v1
mode = etp

[SLURM_ARGS]
cpus_per_task = 32
mem_of_node = 350000
max_runtime = 7200
```
The key setting was:

```ini
run_red = True
```
This was required because the genome assembly was not soft-masked.
## 7. Snakemake Dry Run
A Snakemake dry run was performed before launching the full BRAKER4 analysis:
```bash
snakemake \
    -s Snakefile \
    --use-singularity \
    --singularity-args "--bind /work/ebg_lab/eb/diatom_consortia/metatranscriptomics" \
    --cores 24 \
    --latency-wait 120 \
    --printshellcmds \
    -n
```
The dry run successfully built the workflow DAG and included the expected BRAKER4 rules:
```text
run_masking
run_genemark_etp
run_tsebra
run_augustus_hints
busco_genome
busco_proteins
collect_results
```
The presence of `run_genemark_etp` confirmed that BRAKER4 recognized the RNA-seq and protein evidence and configured the analysis in ETP mode. The presence of `run_tsebra` confirmed that TSEBRA refinement would be performed internally by BRAKER4. Therefore, TSEBRA does not need to be run separately after BRAKER4 completion.
## 8. BRAKER4 Execution
After the dry run completed successfully, BRAKER4 was executed using the official BRAKER4 `Snakefile`.
```bash
snakemake \
    -s Snakefile \
    --use-singularity \
    --singularity-args "--bind /work/ebg_lab/eb/diatom_consortia/metatranscriptomics" \
    --cores 24 \
    --latency-wait 120 \
    --printshellcmds \
    --rerun-incomplete
```
If the GeneMark key is not visible inside the container, the home directory can also be bound:
```bash
snakemake \
    -s Snakefile \
    --use-singularity \
    --singularity-args "--bind /work/ebg_lab/eb/diatom_consortia/metatranscriptomics,/home/ruchita.solanki" \
    --cores 24 \
    --latency-wait 120 \
    --printshellcmds \
    --rerun-incomplete
```
The following optional check can be used to confirm whether the GeneMark key is visible inside the container:
```bash
singularity exec \
    --bind /work/ebg_lab/eb/diatom_consortia/metatranscriptomics \
    braker3.sif \
    ls -lh /home/ruchita.solanki/.gm_key
```
BRAKER4 internally performs RED repeat masking, GeneMark-ETP training, AUGUSTUS training and prediction, evidence integration, TSEBRA refinement, BUSCO/compleasm assessment, and final result collection.
## 9. Extraction of Final Annotation Files
After BRAKER4 completed successfully, the final annotation files were extracted from the BRAKER4 results directory.
Output directory:
```bash
/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/output/DL_diatom/results
```
The compressed final files were extracted as follows:
```bash
cd /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/output/DL_diatom/results

gunzip -c braker.gtf.gz > DL_diatom.braker4.gtf
gunzip -c braker.gff3.gz > DL_diatom.braker4.gff3
gunzip -c braker.aa.gz > DL_diatom.braker4.proteins.faa
gunzip -c braker.codingseq.gz > DL_diatom.braker4.cds.fna
```
The final files are:
```text
DL_diatom.braker4.gtf
DL_diatom.braker4.gff3
DL_diatom.braker4.proteins.faa
DL_diatom.braker4.cds.fna
```
No additional TSEBRA run was performed because the BRAKER4 workflow had already included TSEBRA refinement internally.
## 10. Basic Annotation Statistics
Basic annotation statistics were generated from the final GFF3, protein FASTA, and CDS FASTA files.
```bash
grep -c $'\tgene\t' DL_diatom.braker4.gff3
grep -c $'\tmRNA\t' DL_diatom.braker4.gff3
grep -c $'\tCDS\t' DL_diatom.braker4.gff3

grep -c "^>" DL_diatom.braker4.proteins.faa
grep -c "^>" DL_diatom.braker4.cds.fna

seqkit stats DL_diatom.braker4.proteins.faa
seqkit stats DL_diatom.braker4.cds.fna
```
Predicted proteins were checked for stop codons:
```bash
grep -n "\*" DL_diatom.braker4.proteins.faa | head
```
Very short predicted proteins were inspected:
```bash
seqkit fx2tab -n -l DL_diatom.braker4.proteins.faa \
    | awk '$2 < 50' \
    | head
```
## 11. BRAKER4 Reports and Evidence Support
BRAKER4 report files and quality-control summaries were located using:
```bash
find . -type f | grep -Ei "report|summary|busco|compleasm|statistics|support"
```
Relevant outputs include:
```text
braker_report.html
gene_support.tsv
BUSCO/compleasm summaries
gene set statistics
```
These files were used to assess annotation completeness, evidence support, and overall gene model quality.
## 12. Independent BUSCO Assessment
An independent BUSCO assessment was performed on the predicted protein set.
```bash
busco \
    -i DL_diatom.braker4.proteins.faa \
    -l stramenopiles_odb12 \
    -m proteins \
    -o busco_DL_diatom_braker4_proteins \
    -c 24
```
The BUSCO output was used to evaluate the completeness of the predicted protein set.
## 13. Longest Isoform Extraction
If BRAKER4 generated longest-isoform files, these were extracted for downstream functional annotation:
```bash
find . -type f | grep -Ei "longest|aa|codingseq"
```
If present, the longest-isoform protein and CDS files were extracted:
```bash
gunzip -c braker.longest.aa.gz > DL_diatom.braker4.longest.proteins.faa
gunzip -c braker.longest.codingseq.gz > DL_diatom.braker4.longest.cds.fna
```
The longest protein isoform file was used preferentially for downstream functional annotation:
```text
DL_diatom.braker4.longest.proteins.faa
```
If longest-isoform files were not available, the full BRAKER4 protein set was used:
```text
DL_diatom.braker4.proteins.faa
```
## 14. Functional Annotation
The final protein set was used for downstream functional annotation. Annotation tools may include:
```text
eggNOG-mapper
InterProScan
DIAMOND/BLASTP against UniProt or Swiss-Prot
KEGG/KO annotation
```
Functional annotation focused on pathways relevant to the diatom-dominated consortium, including photosynthesis, carbon-concentrating mechanisms, silica and frustule formation, nitrogen assimilation, lipid metabolism, vitamin and cofactor metabolism, stress responses, motility, and extracellular polymeric substance production.
## 15. Expression Quantification Against the Final Gene Models
After the final BRAKER4 annotation was accepted, RNA-seq expression was quantified against the final GTF file.
For gene-level counts, `featureCounts` can be used:
```bash
featureCounts \
    -T 24 \
    -p \
    -t exon \
    -g gene_id \
    -a DL_diatom.braker4.gtf \
    -o DL_diatom.braker4.featureCounts.txt \
    /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam
```
Alternatively, StringTie can be used to estimate transcript abundance:
```bash
stringtie \
    /work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/Diatoms_Combined_Aligned.sortedByCoord.out.bam \
    -G DL_diatom.braker4.gtf \
    -e \
    -B \
    -p 24 \
    -o DL_diatom.stringtie.gtf
```
The resulting expression estimates can be used with the functional annotation to summarize transcriptional activity across major diatom functional categories.
