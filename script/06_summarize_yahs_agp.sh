#!/bin/bash
set -euo pipefail

# Logic: enter the YaHS scaffolding output directory.
cd /work/ebg_lab/eb/diatom_consortia/hi-c_diatoms/03_yahs_scaffolding

# Logic: define the YaHS AGP file that records how original contigs were represented in scaffolds.
AGP=DL_diatom_whole_hic_yahs_scaffolds_final.agp

# Logic: count AGP gap lines; gap lines indicate scaffold joins separated by Ns.
echo "Number of gap lines, N rows:"
awk '$5=="N"{n++} END{print n+0}' ${AGP}

# Logic: count scaffolds containing more than one original contig component.
echo "Scaffolds with more than one W component:"
awk '$5=="W"{count[$1]++} END{n=0; for(s in count) if(count[s]>1) n++; print n}' ${AGP}

# Logic: list the scaffolds with the largest number of original contig components.
echo "Top scaffolds by number of W components:"
awk '$5=="W"{count[$1]++} END{for(s in count) if(count[s]>1) print s,count[s]}' ${AGP} \
    | sort -k2,2nr \
    | head -30

# Logic: identify original contigs that appear in multiple AGP components, which indicates splitting.
echo "Original contigs split across multiple AGP components:"
awk '$5=="W"{count[$6]++} END{for(c in count) if(count[c]>1) print c,count[c]}' ${AGP} \
    | sort -k2,2nr \
    | head -50
