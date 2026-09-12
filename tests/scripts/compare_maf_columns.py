#!/usr/bin/env python3
"""Diff every shared column between two MAF files, row-matched by variant identity.

The vcf2maf validation script checks a handful of named columns; this checks all of them,
which is how the DNP/TNP/ONP End_Position bug stayed hidden for three validation rounds.

Usage: compare_maf_columns.py <ours.maf> <vcf2maf.maf> [--source-vcf VCF] [--tsv]

Rows are keyed on (chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2) — those
four are the join key, so their agreement is implied by the match count, not by a diff row.
Chromosome is compared with any "chr" prefix stripped, so the key survives either naming.

Every diff lands in a named bucket. Only `unexplained` is a failure; the rest are known,
counted divergences (see BUCKETS). --source-vcf is what makes the multiallelic bucket real
rather than guessed: the positions come out of the input VCF.
"""

import csv
import gzip
import sys
from collections import Counter, defaultdict

BUCKETS = {
    "multiallelic": "row's VCF position carries >1 ALT (from --source-vcf)",
    "entrez": "Entrez_Gene_Id: we emit a real ID, vcf2maf emits 0",
    "filter_tag": "FILTER: vcf2maf appends its own common_variant tag",
    "all_effects": "declared out of scope for 0.8",
    "center": "Center: our --center default is Unknown_Center, vcf2maf writes '.'",
    "matched_normal": "Match_Norm_*: vcf2maf invents NORMAL + REF on a VCF with no normal",
    "unexplained": "",
}

COLUMN_BUCKET = {
    "Entrez_Gene_Id": "entrez",
    "all_effects": "all_effects",
    "Center": "center",
}

# Only counted as matched_normal when *we* left the cell empty, i.e. the VCF has no normal
# sample. On a real tumour/normal file we populate these, and a diff there is a real one.
MATCHED_NORMAL_COLS = {
    "Matched_Norm_Sample_Barcode", "Matched_Norm_Sample_UUID",
    "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
    "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2",
    "n_depth", "n_ref_count", "n_alt_count",
}


def read_maf(path):
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r and not r[0].startswith("#")]
    header, *data = rows
    return header, [dict(zip(header, r)) for r in data]


def key(row):
    chrom = row["Chromosome"].removeprefix("chr")
    return (chrom, row["Start_Position"], row["Reference_Allele"], row["Tumor_Seq_Allele2"])


def multiallelic_positions(vcf_path):
    """Every (chrom, pos) covered by a REF of a multiallelic VCF line.

    The whole REF span, not just the anchor: prefix trimming moves a MAF row downstream of
    the VCF POS it came from (chr3:65439885 in the VCF, 3:65439916 in the MAF).
    """
    covered = set()
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 5)
            if len(f) < 5 or "," not in f[4]:
                continue
            chrom, pos, ref = f[0].removeprefix("chr"), int(f[1]), f[3]
            covered.update((chrom, p) for p in range(pos, pos + len(ref)))
    return covered


def classify(col, ours_val, theirs_val, is_multi):
    if col in COLUMN_BUCKET:
        return COLUMN_BUCKET[col]
    if col == "FILTER" and "common_variant" in theirs_val:
        return "filter_tag"
    if col in MATCHED_NORMAL_COLS and not ours_val:
        return "matched_normal"
    if is_multi:
        return "multiallelic"
    return "unexplained"


def main(argv):
    args = [a for a in argv if not a.startswith("--")]
    as_tsv = "--tsv" in argv
    source_vcf = next((a.split("=", 1)[1] for a in argv if a.startswith("--source-vcf=")), None)
    if "--source-vcf" in argv:
        source_vcf = argv[argv.index("--source-vcf") + 1]
        args = [a for a in args if a != source_vcf]
    ours_path, theirs_path = args[0], args[1]

    ours_header, ours_rows = read_maf(ours_path)
    theirs_header, theirs_rows = read_maf(theirs_path)
    ours = {key(r): r for r in ours_rows}
    theirs = {key(r): r for r in theirs_rows}
    shared_keys = ours.keys() & theirs.keys()
    multi = multiallelic_positions(source_vcf) if source_vcf else set()

    shared_cols = [c for c in ours_header if c in theirs_header]
    buckets = Counter({b: 0 for b in BUCKETS})
    per_col = defaultdict(Counter)
    unexplained_rows = set()
    examples = {}
    for col in shared_cols:
        diffs = [k for k in shared_keys if ours[k][col] != theirs[k][col]]
        if not diffs:
            continue
        for k in diffs:
            b = classify(col, ours[k][col], theirs[k][col], (k[0], int(k[1])) in multi)
            buckets[b] += 1
            per_col[col][b] += 1
            if b == "unexplained":
                unexplained_rows.add(k)
        examples[col] = Counter(
            (ours[k][col] or "<empty>", theirs[k][col] or "<empty>") for k in diffs
        ).most_common(3)

    stats = {
        "ours_rows": len(ours_rows),
        "ours_columns": len(ours_header),
        "vcf2maf_rows": len(theirs_rows),
        "shared_columns": len(shared_cols),
        "matched_rows": len(shared_keys),
        # bucket counts are cells; this is how many *rows* carry at least one unexplained cell
        "unexplained_rows": len(unexplained_rows),
        "ours_only_rows": len(ours.keys() - theirs.keys()),
        "vcf2maf_only_rows": len(theirs.keys() - ours.keys()),
        "source_vcf_multiallelic": 1 if source_vcf else 0,
    }

    if as_tsv:
        for k, v in stats.items():
            print(f"stat\t{k}\t{v}")
        for b in BUCKETS:
            print(f"bucket\t{b}\t{buckets[b]}")
        for col in shared_cols:
            if col in per_col:
                detail = ",".join(f"{b}={n}" for b, n in per_col[col].most_common())
                print(f"col\t{col}\t{sum(per_col[col].values())}\t{detail}")
        return 0

    print(f"ours: {stats['ours_rows']} rows, {stats['ours_columns']} columns")
    print(f"vcf2maf: {stats['vcf2maf_rows']} rows, {len(theirs_header)} columns")
    print(f"matched rows: {stats['matched_rows']}")
    print(f"ours-only rows: {stats['ours_only_rows']}")
    print(f"vcf2maf-only rows: {stats['vcf2maf_only_rows']}")
    print(f"\ncolumns only in ours: {[c for c in ours_header if c not in theirs_header]}")
    print(f"columns only in vcf2maf: {[c for c in theirs_header if c not in ours_header]}")
    if not source_vcf:
        print("\n(no --source-vcf: multiallelic rows cannot be identified, they land in unexplained)")
    print(f"\nper-column diffs on {stats['matched_rows']} matched rows ({len(shared_cols)} shared columns):")
    for col in shared_cols:
        if col not in per_col:
            continue
        n = sum(per_col[col].values())
        pct = 100 * n / max(stats["matched_rows"], 1)
        detail = ",".join(f"{b}={c}" for b, c in per_col[col].most_common())
        shown = "; ".join(f"{o!r} vs {t!r} x{c}" for (o, t), c in examples[col])
        print(f"  {col:<28} {n:>6} ({pct:5.2f}%)  [{detail}]  {shown}")
    print("\nbuckets:")
    for b, desc in BUCKETS.items():
        print(f"  {b:<16} {buckets[b]:>6}  {desc}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
