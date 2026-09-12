#!/usr/bin/env bash
# validate_vep_vs_vcf2maf.sh — compare vcf-reformatter's MAF output against MSKCC vcf2maf,
# on a real, already-VEP-annotated VCF (no re-annotation needed).
#
# vcf2maf normally re-runs VEP itself; --inhibit-vep makes it parse the CSQ field already
# present in the input VCF instead, which is exactly what our real data/ files have. This
# still requires samtools (flanking-bp checks) + tabix (unconditional startup check) + a
# GRCh38 reference FASTA, but skips the VEP binary/cache entirely.
#
# Usage: tests/scripts/validate_vep_vs_vcf2maf.sh <vep-annotated.vcf[.gz]> [output-dir]
#
# One-time setup (once per machine):
#   micromamba create -y -n vcf-reformatter-validation -c bioconda -c conda-forge bcftools samtools
set -euo pipefail

INPUT_VCF="${1:?Usage: $0 <vep-annotated.vcf[.gz]> [output-dir]}"
OUTDIR="${2:-/tmp/vcf2maf_validation}"
ENV_NAME="${VALIDATION_ENV:-vcf-reformatter-validation}"
GENOME_DIR="${GENOME_DIR:-$HOME/genomes/GRCh38}"
VCF2MAF_DIR="${VCF2MAF_DIR:-$HOME/tools/vcf2maf}"
FASTA="$GENOME_DIR/GRCh38.chr.fa"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BIN="$REPO_ROOT/target/release/vcf-reformatter"
MICROMAMBA_BIN="${MICROMAMBA_BIN:-$HOME/.local/bin/micromamba}"
mkdir -p "$OUTDIR"

run() { "$MICROMAMBA_BIN" run -n "$ENV_NAME" "$@"; }

[[ -x "$BIN" ]] || { echo "Build the release binary first: cargo build --release" >&2; exit 1; }

# --- Step 1: GRCh38 reference FASTA (one-time, ~900MB download) ---
if [[ ! -s "$FASTA" ]]; then
    echo "Downloading GRCh38 primary assembly to $GENOME_DIR (one-time, ~900MB)..."
    mkdir -p "$GENOME_DIR"
    curl -L -o "$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
        https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    gunzip "$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    # Ensembl headers are bare ("1", "X"...); add chr prefix to match our chr-prefixed VCFs.
    sed 's/^>\([0-9]\{1,2\}\) />chr\1 /
         s/^>X />chrX /
         s/^>Y />chrY /
         s/^>MT />chrM /' \
        "$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa" > "$FASTA"
    rm "$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    run samtools faidx "$FASTA"
fi

# --- Step 2: vcf2maf.pl itself (one-time, single-file download — core Perl only, no CPAN deps) ---
if [[ ! -f "$VCF2MAF_DIR/vcf2maf.pl" ]]; then
    mkdir -p "$VCF2MAF_DIR"
    curl -L -o "$VCF2MAF_DIR/vcf2maf.pl" \
        https://raw.githubusercontent.com/mskcc/vcf2maf/main/vcf2maf.pl
fi

# --- Step 3: vcf2maf requires an uncompressed input VCF ---
BASENAME="$(basename "$INPUT_VCF")"
BASENAME="${BASENAME%.vcf.gz}"
BASENAME="${BASENAME%.vcf}"
PLAIN_VCF="$OUTDIR/${BASENAME}.plain.vcf"
if [[ "$INPUT_VCF" == *.gz ]]; then
    run bcftools view "$INPUT_VCF" -o "$PLAIN_VCF"
else
    cp "$INPUT_VCF" "$PLAIN_VCF"
fi

# --- Step 4: run vcf2maf (--inhibit-vep parses the existing CSQ, doesn't re-annotate) ---
VCF2MAF_OUT="$OUTDIR/${BASENAME}.vcf2maf.maf"
SAMTOOLS_BIN="$(run which samtools)"
TABIX_BIN="$(run which tabix)"

# Bind vcf2maf to the VCF's first sample column. Without this it defaults to looking for a
# sample literally named TUMOR, finds nothing, and emits empty t_depth/t_ref_count/t_alt_count —
# which silently excludes that whole column family from the comparison.
TUMOR_ID="$(run bcftools query -l "$PLAIN_VCF" | head -1)"
echo "vcf2maf tumor sample: ${TUMOR_ID:-<none found>}"

run perl "$VCF2MAF_DIR/vcf2maf.pl" \
    --inhibit-vep \
    --input-vcf "$PLAIN_VCF" \
    --output-maf "$VCF2MAF_OUT" \
    --samtools-exec "$SAMTOOLS_BIN" \
    --tabix-exec "$TABIX_BIN" \
    --ref-fasta "$FASTA" \
    --ncbi-build GRCh38 \
    ${TUMOR_ID:+--tumor-id "$TUMOR_ID" --vcf-tumor-id "$TUMOR_ID"}

# --- Step 5: run vcf-reformatter on the same VCF ---
OURS_PREFIX="${BASENAME}_ours"
"$BIN" "$INPUT_VCF" -o "$OUTDIR" -p "$OURS_PREFIX" --report none -a vep --output-format maf
OURS_MAF="$OUTDIR/${OURS_PREFIX}_reformatted.maf"

# --- Step 6: compare by variant identity (Chromosome+Start_Position+Reference_Allele+
# Tumor_Seq_Allele2), not by line position. Line-position pairing breaks whenever row counts
# differ, which they will here: vcf-reformatter deliberately splits multiallelic sites into one
# MAF row per ALT allele (see count_multiallelic_sites / from_reformatted_record_multi), while
# vcf2maf run this way (against a non-normalized VCF) does not. The "chr" prefix is also stripped
# for the join key only (vcf-reformatter's Chromosome column omits "chr" per MAF spec; vcf2maf's
# default output keeps it) -- both are reported as informational counts, not failures.
echo
echo "=== $OURS_MAF vs $VCF2MAF_OUT ==="
grep -v "^#" "$OURS_MAF" > "$OUTDIR/.ours_noheader.tsv"
grep -v "^#" "$VCF2MAF_OUT" > "$OUTDIR/.vcf2maf_noheader.tsv"

awk -F'\t' -v OURS="$OUTDIR/.ours_noheader.tsv" -v THEIRS="$OUTDIR/.vcf2maf_noheader.tsv" '
function idx(file,   line, n, i, h) {
    getline line < file
    n = split(line, h, "\t")
    for (i = 1; i <= n; i++) col[file, h[i]] = i
    return n
}
function keyof(a, chrom,   c) {
    c = a[col[F, "Chromosome"]]
    sub(/^chr/, "", c)
    return c "|" a[col[F, "Start_Position"]] "|" a[col[F, "Reference_Allele"]] "|" a[col[F, "Tumor_Seq_Allele2"]]
}
BEGIN {
    idx(OURS); idx(THEIRS)
    while ((getline line < OURS) > 0) {
        n = split(line, a, "\t"); F = OURS; k = keyof(a)
        ours_hs[k] = a[col[F,"Hugo_Symbol"]]; ours_vc[k] = a[col[F,"Variant_Classification"]]
        ours_vt[k] = a[col[F,"Variant_Type"]]; ours_ep[k] = a[col[F,"End_Position"]]
        ours_seen[k] = 1; ours_n++
    }
    while ((getline line < THEIRS) > 0) {
        n = split(line, a, "\t"); F = THEIRS; k = keyof(a)
        their_hs[k] = a[col[F,"Hugo_Symbol"]]; their_vc[k] = a[col[F,"Variant_Classification"]]
        their_vt[k] = a[col[F,"Variant_Type"]]; their_ep[k] = a[col[F,"End_Position"]]
        their_seen[k] = 1; their_n++
    }
    matched = 0; hs_diff = 0; vc_diff = 0; vt_diff = 0; ep_diff = 0; ours_only = 0
    for (k in ours_seen) {
        if (k in their_seen) {
            matched++
            if (ours_hs[k] != their_hs[k]) hs_diff++
            if (ours_vc[k] != their_vc[k]) vc_diff++
            if (ours_vt[k] != their_vt[k]) vt_diff++
            if (ours_ep[k] != their_ep[k]) ep_diff++
        } else ours_only++
    }
    their_only = 0
    for (k in their_seen) if (!(k in ours_seen)) their_only++

    print "ours rows: " ours_n "   vcf2maf rows: " their_n
    print "matched by (chrom,pos,ref,alt): " matched
    print "ours-only (e.g. extra multiallelic ALT rows): " ours_only
    print "vcf2maf-only (variants we produced no row for): " their_only
    print "Hugo_Symbol diffs on matched rows: " hs_diff " / " matched
    print "Variant_Classification diffs on matched rows: " vc_diff " / " matched
    print "Variant_Type diffs on matched rows: " vt_diff " / " matched
    print "End_Position diffs on matched rows: " ep_diff " / " matched
}'
rm -f "$OUTDIR/.ours_noheader.tsv" "$OUTDIR/.vcf2maf_noheader.tsv"