#!/bin/bash
set -euo pipefail

# Logic: enter the Hi-C mapping directory containing the samtools idxstats output.
cd /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/02_map_to_whole_assembly

# Logic: convert idxstats into a contig-level yes/no table for Hi-C read representation.
awk 'BEGIN {
    OFS="\t";
    print "contig","length_bp","mapped_HiC_reads","unmapped_HiC_reads","HiC_mapped"
}
$1!="*" {
    status = ($3 > 0 ? "yes" : "no");
    print $1,$2,$3,$4,status
}' hic_to_whole_assembly.idxstats.txt \
> hic_contig_mapping_presence.tsv

# Logic: save the main contig-level and read-level mapping summary in one text file.
{
    echo "Hi-C contig mapping summary"
    echo "Date: $(date)"
    echo
    echo "Input assembly:"
    echo "/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_output/pypolca_corrected.fasta"
    echo
    echo "Hi-C reads:"
    echo "/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/1574499_S6_L001_R1_001.fastq.gz"
    echo "/work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/1574499_S6_L001_R2_001.fastq.gz"
    echo
    awk '
    $1!="*" {
        total++;
        if ($3 > 0) mapped++;
    }
    END {
        unmapped = total - mapped;
        printf "Total contigs: %d\n", total;
        printf "Contigs with >=1 Hi-C read mapped: %d\n", mapped;
        printf "Contigs with 0 Hi-C reads mapped: %d\n", unmapped;
        printf "Percent contigs with Hi-C reads mapped: %.2f%%\n", (mapped/total)*100;
    }
    ' hic_to_whole_assembly.idxstats.txt
    echo
    echo "Read-level mapping summary:"
    cat hic_to_whole_assembly.flagstat.txt
} > hic_contig_mapping_summary.txt

# Logic: report threshold-based contig counts for optional stronger-support summaries.
awk '
$1!="*" {
    total++;
    if ($3 >= 1) c1++;
    if ($3 >= 5) c5++;
    if ($3 >= 10) c10++;
    if ($3 >= 50) c50++;
    if ($3 >= 100) c100++;
}
END {
    printf "Total contigs: %d\n", total;
    printf "Contigs with >=1 mapped Hi-C read: %d (%.2f%%)\n", c1, (c1/total)*100;
    printf "Contigs with >=5 mapped Hi-C reads: %d (%.2f%%)\n", c5, (c5/total)*100;
    printf "Contigs with >=10 mapped Hi-C reads: %d (%.2f%%)\n", c10, (c10/total)*100;
    printf "Contigs with >=50 mapped Hi-C reads: %d (%.2f%%)\n", c50, (c50/total)*100;
    printf "Contigs with >=100 mapped Hi-C reads: %d (%.2f%%)\n", c100, (c100/total)*100;
}
' hic_to_whole_assembly.idxstats.txt \
> hic_contig_mapping_threshold_summary.txt
