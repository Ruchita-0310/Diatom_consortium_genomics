#!/usr/bin/env python3
"""Add within-Deer-Lake Average_TPM percentile to the ORF master table."""

import argparse
import pandas as pd

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    df["DL_expression_percentile"] = (
        df["DL_Average_TPM"].rank(method="max", pct=True) * 100
    ).round(2)

    percentile = df.pop("DL_expression_percentile")
    pos = df.columns.get_loc("DL_Average_TPM") + 1
    df.insert(pos, "DL_expression_percentile", percentile)

    df.to_csv(args.output, sep="\t", index=False)
    print("Rows:", len(df))
    print("Saved:", args.output)

if __name__ == "__main__":
    main()
