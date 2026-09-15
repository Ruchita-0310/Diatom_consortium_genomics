#!/usr/bin/env python3
"""Match TransDecoder ORFs to nf-core/metatdenovo ORF-level TPM values."""

import argparse
import csv
import gzip
from collections import defaultdict

def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")

def fasta_records(path):
    header = None
    seq = []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--peptides", required=True)
    p.add_argument("--counts", required=True)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-expression", required=True)
    args = p.parse_args()

    tpm_sum = defaultdict(float)
    tpm_n = defaultdict(int)

    with open_text(args.counts) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            orf = row["orf"]
            if orf.startswith("cds."):
                orf = orf[4:]
            tpm_sum[orf] += float(row["tpm"])
            tpm_n[orf] += 1

    avg_tpm = {orf: tpm_sum[orf] / tpm_n[orf] for orf in tpm_sum}

    with open(args.out_expression, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["DL_ORF", "DL_Average_TPM"])
        for orf in sorted(avg_tpm):
            writer.writerow([orf, f"{avg_tpm[orf]:.6f}"])

    total = kept = 0
    with open(args.out_fasta, "w") as out:
        for orf, seq in fasta_records(args.peptides):
            total += 1
            if orf not in avg_tpm:
                continue
            kept += 1
            out.write(f">{orf}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

    print("TransDecoder peptide ORFs:", total)
    print("ORFs with TPM support:", kept)
    print("Expression rows:", len(avg_tpm))
    print("Saved FASTA:", args.out_fasta)
    print("Saved expression table:", args.out_expression)

if __name__ == "__main__":
    main()
