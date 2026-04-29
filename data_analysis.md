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
# 8. MetaEUK combined with polished assembly
## Mapping
```
# 1. Map the Eukaryotic Coding Sequences to the Polished Assembly
# -x asm5 is for high-identity DNA sequences (>95% match)
minimap2 -x asm5 -t 8 \
/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/sr_pypolca_corrected.fasta \
sr_contigs_metaeuk_output.codon.fas > metaeuk_dna_map.paf

# 2. Extract the IDs of the contigs that hold these genes
# In PAF format, Column 6 is the Target (Assembly) name
awk '($10/$11) >= 0.95 {print $6}' metaeuk_dna_map.paf | sort | uniq > 95_euk_contig_ids.txt

# 3. Create your final Diatom Bin
seqtk subseq /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/sr_pypolca_corrected.fasta \
95_euk_contig_ids.txt > 95_Diatom_Euk_Bin.fasta
```
## Stats
```
echo "Total Contigs:" && grep -c ">" Diatom_Euk_Bin.fasta
#Total Contigs:
#4295

echo "Total Genome Size (bp):" && grep -v ">" Diatom_Euk_Bin.fasta | tr -d '\n' | wc -m
#Total Genome Size (bp):
#184,665,232

stats.sh in=Diatom_Euk_Bin.fasta
# N50 = 73.409 KB

# Corrected Short-Read Mapping 
# 1. Index the bin
minimap2 -d diatom_index.mmi 95_Diatom_Euk_Bin.fasta

# 2. Map & Sort (One efficient pipe)
# Note: Output name changed to '_sorted.bam' here for clarity
minimap2 -ax sr -t 16 diatom_index.mmi \
  /path/to/R1_trimmed.fastq.gz \
  /path/to/R2_trimmed.fastq.gz | \
  samtools view -u - | \
  samtools sort -@ 8 -o 95_Diatom_PE_sorted.bam

# 3. Index the sorted bam
# Use the exact same name from Step 2
samtools index 95_Diatom_PE_sorted.bam

# 4. Generate Flagstat
samtools flagstat 95_Diatom_PE_sorted.bam > 95_diatom_short_read_stats.txt

# 5. Calculate the Mean Depth
samtools depth -a 95_Diatom_PE_sorted.bam | \
  awk '{sum+=$3; cnt++} END {if (cnt > 0) print "Mean Depth = ", sum/cnt; else print "No data"}' > 95_mean_depth_result.txt

# BUSCO
busco -i 95_Diatom_Euk_Bin.fasta \
        -o BUSCO_Diatom_Check \
        -m genome \
        -l /work/ebg_lab/eb/diatom_consortia/MAGS_guppy/busco_downloads/stramenopiles_odb10 \
        --metaeuk \
        --cpu 32 
```
# 9. Phylogenetic tree
```
clustalo -i 18S_new.fasta -o 18S_aligned.fasta
trimal -in 18S_aligned.fasta -out 18S_trimmed.fasta -automated1
/home/ruchita.solanki/iqtree-2.2.2.7-Linux/bin/iqtree2 -s 18S_trimmed.fasta -B 1000 -T 4
```
Once you have the tree file, run this to get names next to the accession IDs
```
python3 -c '
import re
import glob

# Using your specific tree filename
tree_filename = "18S_trimmed.fasta.treefile"
name_map = {}

# 1. Loop through all .fasta files (seqdump, Nitzschia_18S_full, and individual IDs)
fasta_files = glob.glob("*.fasta")
print(f"Reading names from: {len(fasta_files)} files...")

for fasta in fasta_files:
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                parts = line.strip().split()
                acc = parts[0].replace(">", "")
                if len(parts) >= 3:
                    # Extracts Genus_species and cleans special characters
                    species = re.sub(r"[^a-zA-Z0-9_]", "", f"{parts[1]}_{parts[2]}")
                    name_map[acc] = f"{acc}_{species}"

# 2. Perform the swap on the tree
try:
    with open(tree_filename, "r") as t:
        content = t.read()
    
    # Sort keys by length (longest first) to prevent partial matching bugs
    for acc in sorted(name_map.keys(), key=len, reverse=True):
        content = content.replace(acc, name_map[acc])
    
    output_name = "final_species_tree.tre"
    with open(output_name, "w") as out:
        out.write(content)
    
    print(f"\nDone! Mapped {len(name_map)} unique IDs.")
    print(f"Your fixed tree is: {output_name}")

except FileNotFoundError:
    print(f"\nError: Could not find {tree_filename}")
'
```
You can next remove all the "Uncultured", to clean up the tree:
```
from Bio import Phylo

# Load the tree
tree = Phylo.read("new18_species_tree.tre", "newick")

# Identify all terminals (leaves) containing "Uncultured"
to_prune = [leaf for leaf in tree.get_terminals() if "Uncultured" in leaf.name]

# Remove them
for leaf in to_prune:
    tree.prune(leaf)

# Save the cleaned tree
Phylo.write(tree, "cleaned_species_tree.tre", "newick")

print(f"Removed {len(to_prune)} nodes.")
```
# 10. Transcriptome analysis 
[Nf core metadenovo](https://github.com/nf-core/metatdenovo)
```
# --- 1. JAVA SETUP ---
module purge
module load java/openjdk-23.0.1
export JAVA_HOME=$(dirname $(dirname $(readlink -f $(which java))))
export PATH=$JAVA_HOME/bin:$PATH


# --- 3. EXECUTION ---
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

# 11. Identifying rRNA genes from transcriptome
```
barrnap --kingdom euk --threads 4 spades.transcripts.fa --outseq euk_transcript_rRNA.fna > diatom_euk_rRNA.gff
barrnap --kingdom bac spades.transcripts.fa --outseq bac_transcript_rRNA.fna > diatom_bac_rRNA.gff
barrnap --kingdom mito spades.transcripts.fa --outseq mito_transcript_rRNA.fna > diatom_mito_rRNA.gff
```
