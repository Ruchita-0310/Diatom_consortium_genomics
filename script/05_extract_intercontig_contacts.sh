#!/bin/bash
set -euo pipefail

# Logic: enter the Hi-C mapping directory containing the coordinate-sorted BAM file.
cd /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/02_map_to_whole_assembly

# Logic: define the coordinate-sorted BAM used for contact extraction.
BAM=hic_to_whole_assembly.coord_sorted.bam

# Logic: extract MAPQ >= 5 read records whose mates map to different contigs and count contact support for each contig pair.
{
    echo -e "contig_A\tcontig_B\tcontact_read_records\tapprox_contact_pairs"

    samtools view -@ 8 -q 5 -F 2316 ${BAM} \
        | awk 'BEGIN{OFS="\t"}
               $3!="*" && $7!="*" && $7!="=" {
                   a=$3; b=$7;
                   if (a < b) print a,b;
                   else print b,a;
               }' \
        | sort -S 20G \
        | uniq -c \
        | awk 'BEGIN{OFS="\t"} {print $2,$3,$1,$1/2}' \
        | sort -k4,4nr
} > hic_intercontig_contacts_MAPQ5.tsv

# Logic: add contig lengths to each inter-contig contact pair.
awk 'BEGIN{OFS="\t"}
     NR==FNR {len[$1]=$2; next}
     FNR==1 {
         print $0,"len_A","len_B","min_len","records_per_Mb_minlen";
         next
     }
     {
         lA=len[$1];
         lB=len[$2];
         min=(lA<lB ? lA : lB);
         norm=($3/(min/1000000));
         print $0,lA,lB,min,norm
     }' \
     pypolca_corrected.chrom.sizes \
     hic_intercontig_contacts_MAPQ5.tsv \
     > hic_intercontig_contacts_MAPQ5.with_lengths.tsv
