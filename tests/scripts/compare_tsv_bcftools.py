#!/usr/bin/env python3
"""Compare vcf-reformatter's TSV annotation columns against `bcftools +split-vep`.

    compare_tsv_bcftools.py --ours X.tsv --bcftools split.tsv --mode first|split \
        --pairs CSQ_SYMBOL=SYMBOL,CSQ_Consequence=Consequence,...

`bcftools` output has no header: CHROM POS REF ALT then the fields, in --pairs order.
Without -d bcftools joins every transcript's value with commas; `first` takes element 0
of each, which is what `-t first` keeps. With -d it is one row per transcript; `split`
compares the two files as multisets. Empty and "." are the same thing on both sides.

Prints one line: `rows=<n> diffs=<d>`; exit 1 if d > 0. Stdlib only.
"""
import argparse
import sys
from collections import Counter


def norm(v):
    return "" if v == "." else v


def read_ours(path, cols):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = [header.index(c) for c in ["CHROM", "POS", "REF", "ALT"] + cols]
        for line in f:
            row = line.rstrip("\n").split("\t")
            yield tuple(norm(row[i]) for i in idx)


def read_bcftools(path, mode, n_fields):
    with open(path) as f:
        for line in f:
            row = line.rstrip("\n").split("\t")
            key = tuple(row[:4])
            vals = row[4 : 4 + n_fields]
            if mode == "first":
                vals = [v.split(",")[0] for v in vals]
            yield key + tuple(norm(v) for v in vals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", required=True)
    ap.add_argument("--bcftools", required=True)
    ap.add_argument("--mode", choices=["first", "split"], required=True)
    ap.add_argument("--pairs", required=True, help="OURCOL=bcftools_field,...")
    a = ap.parse_args()
    ours_cols = [p.split("=")[0] for p in a.pairs.split(",")]

    ours = list(read_ours(a.ours, ours_cols))
    theirs = list(read_bcftools(a.bcftools, a.mode, len(ours_cols)))

    if a.mode == "first":
        # both tools keep input order, one row per VCF line
        diffs = sum(1 for x, y in zip(ours, theirs) if x != y) + abs(len(ours) - len(theirs))
        for x, y in zip(ours, theirs):
            if x != y:
                print(f"first diff: ours={x} bcftools={y}", file=sys.stderr)
                break
    else:
        co, ct = Counter(ours), Counter(theirs)
        diffs = sum((co - ct).values()) + sum((ct - co).values())
        for k in list((co - ct).keys())[:1]:
            print(f"only in ours: {k}", file=sys.stderr)
        for k in list((ct - co).keys())[:1]:
            print(f"only in bcftools: {k}", file=sys.stderr)

    print(f"rows={len(ours)} diffs={diffs}")
    sys.exit(1 if diffs else 0)


if __name__ == "__main__":
    main()
