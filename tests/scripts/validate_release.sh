#!/usr/bin/env bash
# validate_release.sh — the release validation matrix. Everything, one file.
# Usage: tests/scripts/validate_release.sh [all|maf|parquet|memory|stdin|reports|bench|edge]
#
# Results land in tests/scripts/results/. Only a FAIL row makes this exit non-zero;
# INFO rows are documented divergences, recorded rather than hidden.
set -euo pipefail

SECTIONS=(maf parquet memory stdin reports flags bench edge)

want="${1:-all}"
if [[ "$want" != all ]] && ! printf '%s\n' "${SECTIONS[@]}" | grep -qx -- "$want"; then
    echo "unknown section: $want (want: all ${SECTIONS[*]})" >&2
    exit 2
fi
run() { [[ "$want" == all || "$want" == "$1" ]]; }
# Every section must leave rows behind. `edge` was declared and empty for a week, so
# `validate_release.sh edge` ran nothing and exited 0 — indistinguishable from a pass.
section_rows() { awk -F'\t' -v s="$1" 'NR>1&&$1==s{n++} END{print n+0}' "$R"; }

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BIN="$ROOT/target/release/vcf-reformatter"
OUT="$ROOT/tests/scripts/results"; WORK="${WORK:-/tmp/vcf_validation}"
mkdir -p "$OUT" "$WORK"
[[ -x "$BIN" ]] || { echo "cargo build --release first" >&2; exit 1; }

R="$OUT/summary.tsv"   # summary.txt and summary.html are written from it at the end
printf 'check\tmetric\tvalue\tstatus\tnote\n' > "$R"
row() { printf '%s\t%s\t%s\t%s\t%s\n' "$@" >> "$R"; }
# eq <section> <metric> <got> <want> <note> — the same four lines forty times over
eq() { local st; [[ "$3" == "$4" ]] && st=PASS || st=FAIL; row "$1" "$2" "$3" "$st" "$5"; }
# exit code of a command that is allowed to fail under set -e
rc() { local c=0; "$@" >/dev/null 2>&1 || c=$?; echo "$c"; }
# peak resident MB of one run; its stdout is noise here.
# BSD time (macOS) prints bytes under -l, GNU time prints kbytes under -f %M, and
# the two flag sets are mutually unintelligible. Detected once, here.
if /usr/bin/time -l true >/dev/null 2>&1; then
    TIME_MODE=bsd
elif /usr/bin/time -f %M true >/dev/null 2>&1; then
    TIME_MODE=gnu
else
    TIME_MODE=none
fi
peak_rss() {
    case "$TIME_MODE" in
        bsd)
            /usr/bin/time -l "$@" >/dev/null 2> "$WORK/rss.err" || true
            awk '/maximum resident set size/{print int($1/1048576); exit}' "$WORK/rss.err"
            ;;
        gnu)
            /usr/bin/time -o "$WORK/rss.err" -f %M "$@" >/dev/null 2>&1 || true
            awk 'NR==1{print int($1/1024); exit}' "$WORK/rss.err"
            ;;
        *) ;;  # unmeasurable: echo nothing, and callers must treat that as a failure
    esac
}

VEP_SMALL="$ROOT/data/B505_V_1/B505_1_V.mutect2.filtered_VEP.ann.vcf.gz"   # 12,239, --pick'd
VEP_BIG="$ROOT/data/B487_1_V_vs_B487_1_cOM.freebayes.annotated.vcf.gz"     # 92,216, T/N
VEP_FB="$ROOT/data/B505_V_1/B505_V_1.freebayes_VEP.ann.vcf.gz"             # 34,415, multi-transcript
SNPEFF="$ROOT/data/variantsNorm3_snpEff_annotated.vcf"                     # 29,589, ANN
MULTISAMPLE="$ROOT/data/multiSamplesMMMT.normalized.filtered.vcf.gz"
UNANNOTATED="$ROOT/data/B487_2_LN_vs_B487_1_cOM/B487_2_LN_vs_B487_1_cOM.mutect2.filtered.vcf.gz"

# --- external tools: bcftools/samtools/tabix live in a micromamba env, vcf2maf is one .pl file
MM="${MICROMAMBA_BIN:-$HOME/.local/bin/micromamba}"
ENV_NAME="${VALIDATION_ENV:-vcf-reformatter-validation}"
VCF2MAF="${VCF2MAF_DIR:-$HOME/tools/vcf2maf}/vcf2maf.pl"
FASTA="${GENOME_DIR:-$HOME/genomes/GRCh38}/GRCh38.chr.fa"
mm() { "$MM" run -n "$ENV_NAME" "$@"; }

# vcf2maf ground truth for one VCF, cached in $WORK (the 92k file takes ~100s).
# Echoes "<maf>\t<tumor sample id>".
ground_truth() {
    local vcf="$1" base plain maf tid
    local -a idargs=()
    base="$(basename "$vcf")"; base="${base%.gz}"; base="${base%.vcf}"
    plain="$WORK/$base.plain.vcf"; maf="$WORK/$base.vcf2maf.maf"
    [[ -s "$plain" ]] || mm bcftools view "$vcf" -o "$plain" >&2
    tid="$(mm bcftools query -l "$plain" | head -1)"
    [[ -z "$tid" ]] || idargs=(--tumor-id "$tid" --vcf-tumor-id "$tid")
    if [[ ! -s "$maf" ]]; then
        # --inhibit-vep parses the CSQ already in the file instead of re-annotating.
        # --tumor-id matters: without it vcf2maf looks for a sample named TUMOR, finds
        # none, and silently empties t_depth/t_ref_count/t_alt_count.
        mm perl "$VCF2MAF" --inhibit-vep --input-vcf "$plain" --output-maf "$maf" \
            --samtools-exec "$(mm which samtools)" --tabix-exec "$(mm which tabix)" \
            --ref-fasta "$FASTA" --ncbi-build GRCh38 "${idargs[@]}" >&2
    fi
    printf '%s\t%s\n' "$maf" "$tid"
}

if run maf; then
    HEADER_EXPECTED="$ROOT/tests/scripts/maf_header.expected"
    # <path>:<picked|nopick> — 'first' can only be held to vcf2maf on --pick'd input.
    MAF_FILES=(
        "$ROOT/data/B505_V_1/B505_1_V.mutect2.filtered_VEP.ann.vcf.gz:picked"
        "$ROOT/data/B505_V_1/B505_V_1.freebayes_VEP.ann.vcf.gz:nopick"
        "$ROOT/data/B487_1_V_vs_B487_1_cOM.freebayes.annotated.vcf.gz:nopick"
    )
    for entry in "${MAF_FILES[@]}"; do
        vcf="${entry%:*}"; picked="${entry##*:}"
        name="$(basename "$vcf")"; name="${name%%.*}.${picked}"
        [[ -s "$vcf" ]] || { row maf "$name" missing FAIL "input VCF not found"; continue; }
        IFS=$'\t' read -r truth tid < <(ground_truth "$vcf")
        oursid=(); [[ -z "$tid" ]] || oursid=(--tumor-id "$tid")

        for mode in first most-severe split; do
            pfx="${name//./_}_${mode//-/_}"
            "$BIN" "$vcf" -o "$WORK" -p "$pfx" --report none -a vep \
                --output-format maf -t "$mode" "${oursid[@]}" >/dev/null
            ours="$WORK/${pfx}_reformatted.maf"

            # --- column presence: a stray tab inside an annotation corrupts the file invisibly
            ncol=$(awk -F'\t' 'NR==1{print NF; exit}' "$ours")
            ragged=$(awk -F'\t' 'NR==1{n=NF} NF!=n{bad++} END{print bad+0}' "$ours")
            [[ "$ncol" == 50 ]] && s=PASS || s=FAIL
            row maf "$name/$mode/columns" "$ncol" "$s" "expected 50"
            [[ "$ragged" == 0 ]] && s=PASS || s=FAIL
            row maf "$name/$mode/ragged_rows" "$ragged" "$s" "rows whose field count != header"
            if [[ -s "$HEADER_EXPECTED" ]]; then
                head -1 "$ours" | tr '\t' '\n' | diff -q - "$HEADER_EXPECTED" >/dev/null \
                    && s=PASS || s=FAIL
                row maf "$name/$mode/header_names" "$s" "$s" "vs maf_header.expected"
            fi

            if [[ "$mode" == split ]]; then
                # vcf2maf emits one row per VCF line, so a whole-file diff is meaningless.
                # Criterion: every vcf2maf variant+transcript appears among our rows.
                python3 - "$ours" "$truth" "$vcf" "$ROOT/tests/scripts" > "$WORK/$pfx.split.txt" <<'SPLITPY'
import sys
sys.path.insert(0, sys.argv[4])
from compare_maf_columns import read_maf, multiallelic_positions
def k(r):
    return (r["Chromosome"].removeprefix("chr"), r["Start_Position"],
            r["Reference_Allele"], r["Tumor_Seq_Allele2"], r["Transcript_ID"])
ours = {k(r) for r in read_maf(sys.argv[1])[1]}
theirs = [k(r) for r in read_maf(sys.argv[2])[1]]
multi = multiallelic_positions(sys.argv[3])
missing = [t for t in theirs if t not in ours]
# vcf2maf emits one row per VCF line, so at a multiallelic site it describes whichever
# allele its single CSQ entry belongs to and skips the rest. A row of its own that we do
# not reproduce there is that known divergence, not a missing variant.
mm = [t for t in missing if (t[0], int(t[1])) in multi]
print(f"ours_rows\t{len(ours)}")
print(f"vcf2maf_rows\t{len(theirs)}")
print(f"missing\t{len(missing)}")
print(f"missing_multiallelic\t{len(mm)}")
print(f"missing_unexplained\t{len(missing) - len(mm)}")
SPLITPY
                v() { awk -F'\t' -v k="$1" '$1==k{print $2}' "$WORK/$pfx.split.txt"; }
                orows=$(v ours_rows)
                [[ "$(v missing_unexplained)" == 0 ]] && s=PASS || s=FAIL
                row maf "$name/split/missing_unexplained" "$(v missing_unexplained)" "$s" \
                    "of $(v missing) vcf2maf (variant,transcript) absent from our $orows rows"
                row maf "$name/split/missing_multiallelic" "$(v missing_multiallelic)" INFO \
                    "vcf2maf emits one row per VCF line at multiallelic sites"
                # annotation count = every CSQ entry in the file, our one-row-per-annotation target
                anns=$(mm bcftools query -f '%INFO/CSQ\n' "$vcf" | awk -F',' '{n+=NF} END{print n}')
                row maf "$name/split/rows_vs_annotations" "$orows/$anns" INFO \
                    "multiallelic ALT expansion accounts for the difference"
                continue
            fi

            python3 "$ROOT/tests/scripts/compare_maf_columns.py" "$ours" "$truth" \
                --source-vcf "$vcf" --tsv > "$OUT/$pfx.buckets.tsv"
            b() { awk -F'\t' -v k="$1" '$1=="bucket"&&$2==k{print $3}' "$OUT/$pfx.buckets.tsv"; }
            matched=$(awk -F'\t' '$1=="stat"&&$2=="matched_rows"{print $3}' "$OUT/$pfx.buckets.tsv")
            unexp=$(b unexplained)
            for bucket in multiallelic entrez filter_tag all_effects center matched_normal; do
                row maf "$name/$mode/$bucket" "$(b "$bucket")" INFO "named divergence, of $matched matched rows"
            done
            if [[ "$mode" == first && "$picked" == nopick ]]; then
                # Documented expected behaviour (2026-08-27): 'first' passes the annotator's
                # own list order through, and this file was not --pick'd.
                urows=$(awk -F'\t' '$1=="stat"&&$2=="unexplained_rows"{print $3}' "$OUT/$pfx.buckets.tsv")
                pct=$(awk -v u="$urows" -v m="$matched" 'BEGIN{printf "%.2f", m?100*u/m:0}')
                row maf "$name/$mode/unexplained" "$unexp" INFO \
                    "cells on $urows rows ($pct%) — -t first on non---pick'd input"
            else
                urows=$(awk -F'\t' '$1=="stat"&&$2=="unexplained_rows"{print $3}' "$OUT/$pfx.buckets.tsv")
                [[ "$unexp" == 0 ]] && s=PASS || s=FAIL
                row maf "$name/$mode/unexplained" "$unexp" "$s" "cells on $urows of $matched matched rows"
            fi
            row maf "$name/$mode/vcf2maf_only_rows" \
                "$(awk -F'\t' '$1=="stat"&&$2=="vcf2maf_only_rows"{print $3}' "$OUT/$pfx.buckets.tsv")" \
                INFO "variants vcf2maf produced and we did not"
        done
    done

    # vcf2maf cannot parse SnpEff's ANN field at all, so there is no ground truth to diff
    # the SnpEff MAF against outside the Fedora VM (Task 14). It is held to structure and
    # internal consistency here instead of to parity — which is still more than the nothing
    # it got before 2026-09-09.
    for mode in first most-severe split; do
        pfx="snpeff_${mode//-/_}"
        "$BIN" "$SNPEFF" -o "$WORK" -p "$pfx" --report none -a snpeff \
            --output-format maf -t "$mode" >/dev/null 2>&1
        ours="$WORK/${pfx}_reformatted.maf"
        [[ -s "$ours" ]] || { row maf "snpeff/$mode" missing FAIL "no MAF produced"; continue; }

        ncol=$(awk -F'\t' 'NR==1{print NF; exit}' "$ours")
        eq maf "snpeff/$mode/columns" "$ncol" 50 "the ANN path must emit the same 50 columns"
        ragged=$(awk -F'\t' 'NR==1{n=NF} NF!=n{bad++} END{print bad+0}' "$ours")
        eq maf "snpeff/$mode/ragged_rows" "$ragged" 0 "rows whose field count != header"
        if [[ -s "$HEADER_EXPECTED" ]]; then
            head -1 "$ours" | tr '\t' '\n' | diff -q - "$HEADER_EXPECTED" >/dev/null \
                && s=PASS || s=FAIL
            row maf "snpeff/$mode/header_names" "$s" "$s" "vs maf_header.expected"
        fi
        # The one cross-column invariant that needs no reference MAF.
        bad=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
             {d=$(h["t_depth"]); r=$(h["t_ref_count"]); a=$(h["t_alt_count"])
              if(d!=""&&d!="."&&r!=""&&r!="."&&a!=""&&a!="."&&r+a>d) n++} END{print n+0}' "$ours")
        eq maf "snpeff/$mode/depth_invariant_violations" "$bad" 0 \
            "t_ref_count + t_alt_count <= t_depth"
        tot=$(( $(wc -l < "$ours") - 1 ))
        for col in Hugo_Symbol Variant_Classification HGVSp_Short Exon_Number; do
            pop=$(awk -F'\t' -v c="$col" 'NR==1{for(i=1;i<=NF;i++) if($i==c) k=i; next}
                  $k!=""&&$k!="."{n++} END{print n+0}' "$ours")
            row maf "snpeff/$mode/populated_$col" "$pop/$tot" INFO \
                "filled from ANN_*; a collapse to 0 is how an ANN regression would show"
        done
    done
fi


if run parquet; then
    # Both annotators: the MAF writer and the parquet writer share code, but the ANN path
    # fills different columns, so losslessness has to be shown on both.
    for src in "vep:$VEP_SMALL" "snpeff:$SNPEFF"; do
      ann="${src%%:*}"; vcf="${src#*:}"
      for fmt in maf tsv; do
        "$BIN" "$vcf" -o "$WORK" -p "pq_${ann}_$fmt" --report none -a "$ann" \
            --output-format "$fmt" --parquet >/dev/null 2>&1
        text="$WORK/pq_${ann}_${fmt}_reformatted.$fmt"; pqf="$text.parquet"
        eq parquet "$ann/$fmt/parquet_written" "$([[ -s "$pqf" ]] && echo yes || echo no)" yes \
            "--parquet writes text AND parquet since 2026-08-19"
        [[ -s "$pqf" ]] || continue

        python3 - "$pqf" "$text" "$fmt" > "$WORK/pq_${ann}_$fmt.txt" <<'PQPY'
import csv, math, sys
import pyarrow.parquet as pq

table = pq.read_table(sys.argv[1])
rows = list(csv.reader(open(sys.argv[2], newline=""), delimiter="\t"))
header, data = rows[0], rows[1:]
# The text writer rounds VAF to 4dp and expands tiny QUALs; parquet keeps the real double.
# Same number, different rendering (2026-08-18), so floats compare with a tolerance.
want = ({"Start_Position": "uint64", "End_Position": "uint64", "t_depth": "uint64",
         "t_ref_count": "uint64", "t_alt_count": "uint64", "QUAL": "double", "VAF": "double"}
        if sys.argv[3] == "maf" else {"POS": "uint64", "QUAL": "double"})

print(f"parquet_rows\t{table.num_rows}")
print(f"text_rows\t{len(data)}")
print(f"columns_match\t{int(table.schema.names == header)}")

idx = {h: i for i, h in enumerate(header)}
bad = nulls = 0
for name in table.schema.names:
    if name not in idx:
        continue
    j = idx[name]
    for i, v in enumerate(table.column(name).to_pylist()):
        txt = data[i][j] if j < len(data[i]) else ""
        if v is None:
            # A numeric column the VCF left as "." is NULL in parquet and "." in text.
            nulls += 1
            bad += txt not in ("", ".")
        elif isinstance(v, float):
            # VAF is written "{:.4}" in text, so the tolerance is absolute, not relative:
            # 0.0150 against parquet's 0.015037594 is the rounding, not a difference.
            try:
                bad += not math.isclose(v, float(txt), rel_tol=1e-6, abs_tol=6e-5)
            except ValueError:
                bad += 1
        else:
            bad += str(v) != txt
print(f"cell_mismatches\t{bad}")
print(f"null_cells\t{nulls}")
violations = [f"{f.name}:{f.type}" for f in table.schema
              if str(f.type) != want.get(f.name, "string")]
print(f"schema_violations\t{len(violations)}")
print(f"schema_detail\t{','.join(violations[:5])}")
PQPY
        g() { awk -F'\t' -v k="$1" '$1==k{print $2}' "$WORK/pq_${ann}_$fmt.txt"; }
        eq parquet "$ann/$fmt/rows" "$(g parquet_rows)" "$(g text_rows)" "parquet vs sibling text"
        eq parquet "$ann/$fmt/columns_match" "$(g columns_match)" 1 "same names, same order"
        eq parquet "$ann/$fmt/cell_mismatches" "$(g cell_mismatches)" 0 \
            "cell-by-cell, keyed by column name"
        eq parquet "$ann/$fmt/schema_violations" "$(g schema_violations)" 0 \
            "positions/depths uint64, QUAL/VAF double, rest string — $(g schema_detail)"
        row parquet "$ann/$fmt/null_cells" "$(g null_cells)" INFO "empty text cell maps to NULL"
        gzip -kf "$text"
        row parquet "$ann/$fmt/size_text_gz_parquet_mb" \
            "$(du -m "$text" "$text.gz" "$pqf" | awk '{printf "%s/", $1}' | sed 's:/$::')" \
            INFO "parquet is ZSTD since 2026-08-18"
      done
    done

    # Parquet streams too since 2026-09-09 — one row group per chunk, not one per file.
    rss=$(peak_rss "$BIN" "$VEP_BIG" -o "$WORK" -p pq_rss --report none -a vep \
        --output-format maf --parquet)
    row parquet peak_rss_mb "${rss:-?}" INFO "92,216 variants, 1.13M annotations, MAF + parquet"
fi

if run memory; then
    # Peak RSS per file, recorded. These three differ in annotation density as well as
    # size, so they are INFO only — a bigger number here can just mean a wider file.
    MEM_FILES=("$VEP_SMALL:12239" "$SNPEFF:29589" "$VEP_FB:34415" "$VEP_BIG:92216")
    for mode in tsv maf tsv+parquet maf+parquet; do
        fmt="${mode%%+*}"; pq=""; [[ "$mode" == *+parquet ]] && pq="--parquet"
        for entry in "${MEM_FILES[@]}"; do
            vcf="${entry%:*}"; n="${entry##*:}"
            r=$(peak_rss "$BIN" "$vcf" -o "$WORK" -p mem_rss --report none -a auto \
                --output-format "$fmt" $pq)   # unquoted: empty means no flag
            row memory "$mode/peak_rss_mb/$n" "${r:-?}" INFO "$n variants"
        done
    done

    # The claim the streaming work actually makes is that memory does not scale with
    # input, so the control is one file against 8 copies of itself — same shape, same
    # density, 8x the bytes. Ratio, not an absolute ceiling, so there is no magic
    # megabyte number to maintain.
    one="$WORK/mem_scale_1x.vcf"; eight="$WORK/mem_scale_8x.vcf"
    [[ -s "$one" ]] || gunzip -c "$VEP_FB" > "$one"
    if [[ ! -s "$eight" ]]; then
        grep '^#' "$one" > "$eight"
        for _ in 1 2 3 4 5 6 7 8; do grep -v '^#' "$one"; done >> "$eight"
    fi
    row memory scale_input_mb "$(du -m "$one" "$eight" | awk '{printf "%s/", $1}' | sed 's:/$::')" \
        INFO "34,415 variants against the same variants 8 times over"
    for mode in tsv maf+parquet; do
        fmt="${mode%%+*}"; pq=""; [[ "$mode" == *+parquet ]] && pq="--parquet"
        small=$(peak_rss "$BIN" "$one" -o "$WORK" -p mem_1x --report none -a vep \
            --output-format "$fmt" $pq)
        big=$(peak_rss "$BIN" "$eight" -o "$WORK" -p mem_8x --report none -a vep \
            --output-format "$fmt" $pq)
        row memory "$mode/peak_rss_mb_1x_8x" "${small:-?}/${big:-?}" INFO "same file, 8x the bytes"
        # Streaming means memory must not track input. TSV is flat to within 1% even at
        # 16x; the MAF path drifts up ~2x over 16x because converting 50 small Strings
        # per record churns the heap — its live footprint stays flat, so this is
        # allocator high-water, not retention. 2x buys room for that without hiding a
        # path that has started holding the file again.
        # An unmeasured run must fail here. Defaulting the two numbers would make
        # "0 <= 2" true and report a pass for a check that measured nothing.
        st=FAIL
        if [[ -z "${small:-}" || -z "${big:-}" ]]; then
            note="no usable /usr/bin/time, nothing measured"
        else
            note="8x input, ${small}MB -> ${big}MB"
            [[ $big -le $((small * 2)) ]] && st=PASS
        fi
        row memory "$mode/sublinear_in_input_size" "${small:-?}/${big:-?}" "$st" "$note"
    done
    rm -f "$eight"
fi

if run stdin; then
    # Both annotators: the gzip sniff is annotator-blind, but -a auto reads the header off
    # a stream it cannot seek back into, which is exactly where a regression would hide.
    for src in "vep:$VEP_SMALL" "snpeff:$SNPEFF"; do
        ann="${src%%:*}"; vcf="${src#*:}"
        gzf="$WORK/sd_${ann}_input.vcf.gz"
        case "$vcf" in *.gz) cp -f "$vcf" "$gzf";; *) gzip -c "$vcf" > "$gzf";; esac
        "$BIN" "$vcf" -o "$WORK" -p "sd_${ann}_path" --report none -a "$ann" >/dev/null 2>&1
        gunzip -c "$gzf" | "$BIN" - -o "$WORK" -p "sd_${ann}_plain" --report none -a "$ann" \
            >/dev/null 2>&1
        "$BIN" - -o "$WORK" -p "sd_${ann}_gz" --report none -a "$ann" < "$gzf" >/dev/null 2>&1
        # -a auto on a pipe, which has no filename and no second pass over the header
        gunzip -c "$gzf" | "$BIN" - -o "$WORK" -p "sd_${ann}_auto" --report none -a auto \
            >/dev/null 2>&1
        for v in plain gz auto; do
            cmp -s "$WORK/sd_${ann}_path_reformatted.tsv" \
                   "$WORK/sd_${ann}_${v}_reformatted.tsv" && st=PASS || st=FAIL
            row stdin "$ann/$v/identical_to_path_arg" "$st" "$st" "piped $v input"
        done
    done

    # Never tested before: a truncated gzip and an empty stream must fail legibly.
    head -c 200000 "$VEP_SMALL" > "$WORK/truncated.vcf.gz"
    code=0; "$BIN" - -o "$WORK" -p sd_trunc --report none -a vep \
        < "$WORK/truncated.vcf.gz" > /dev/null 2> "$WORK/sd_trunc.err" || code=$?
    [[ "$code" != 0 ]] && st=PASS || st=FAIL
    row stdin truncated_gzip/exit_code "$code" "$st" "must be non-zero, not a silent empty file"
    grep -q "panicked" "$WORK/sd_trunc.err" && st=FAIL || st=PASS
    row stdin truncated_gzip/no_panic "$st" "$st" "$(head -c 90 "$WORK/sd_trunc.err" | tr '\n' ' ')"

    code=0; : | "$BIN" - -o "$WORK" -p sd_empty --report none -a vep \
        > /dev/null 2> "$WORK/sd_empty.err" || code=$?
    [[ "$code" != 0 ]] && st=PASS || st=FAIL
    row stdin empty_stdin/exit_code "$code" "$st" "must be non-zero"
    grep -q "panicked" "$WORK/sd_empty.err" && st=FAIL || st=PASS
    row stdin empty_stdin/no_panic "$st" "$st" "$(head -c 90 "$WORK/sd_empty.err" | tr '\n' ' ')"
fi

if run reports; then
    "$BIN" "$VEP_SMALL" -o "$WORK" -p rp_vep --report html -a vep >/dev/null 2>&1
    "$BIN" "$SNPEFF" -o "$WORK" -p rp_snpeff --report html -a snpeff >/dev/null 2>&1
    "$BIN" "$VEP_SMALL" -o "$WORK" -p rp_txt --report txt -a vep >/dev/null 2>&1
    "$BIN" "$VEP_SMALL" -o "$WORK" -p rp_none --report none -a vep >/dev/null 2>&1

    vep_html="$WORK/rp_vep_summary.html"; snpeff_html="$WORK/rp_snpeff_summary.html"
    eq reports vep/html_written "$([[ -s "$vep_html" ]] && echo yes || echo no)" yes ""
    eq reports snpeff/html_written "$([[ -s "$snpeff_html" ]] && echo yes || echo no)" yes ""
    eq reports txt/written "$([[ -s "$WORK/rp_txt_summary.txt" ]] && echo yes || echo no)" yes ""
    eq reports none/no_report_file \
        "$(ls "$WORK"/rp_none_summary.* 2>/dev/null | wc -l | tr -d ' ')" 0 ""

    # Self-contained: nothing fetched from a host at open time.
    eq reports vep/external_refs "$(grep -Eoc '(src|href)="https?://' "$vep_html" || true)" 0 \
        "no script/style/image pulled from a host"
    for term in SIFT PolyPhen; do
        n=$(grep -ci "$term" "$vep_html" || true)
        [[ "$n" -gt 0 ]] && st=PASS || st=FAIL
        row reports "vep/charts_$term" "$n" "$st" "damage breakdown on the VEP page"
        n=$(grep -ci "$term" "$snpeff_html" || true)
        eq reports "snpeff/no_$term" "$n" 0 "SnpEff has neither; the page must not claim them"
    done
    n=$(grep -ci "Impact" "$snpeff_html" || true)
    [[ "$n" -gt 0 ]] && st=PASS || st=FAIL
    row reports snpeff/charts_impact "$n" "$st" "HIGH/MODERATE/LOW/MODIFIER fallback"
    row reports txt/damage_breakdown \
        "$(grep -Eci 'sift|polyphen|impact' "$WORK/rp_txt_summary.txt" || true)" INFO \
        "the txt report carries no damage breakdown at all — Task 12 decides: close or document"
fi

if run flags; then
    # -a auto must reach the same answer as naming the annotator
    "$BIN" "$VEP_SMALL" -o "$WORK" -p fl_vep --report none -a vep >/dev/null 2>&1
    "$BIN" "$VEP_SMALL" -o "$WORK" -p fl_auto_vep --report none -a auto >/dev/null 2>&1
    "$BIN" "$SNPEFF" -o "$WORK" -p fl_snpeff --report none -a snpeff >/dev/null 2>&1
    "$BIN" "$SNPEFF" -o "$WORK" -p fl_auto_snpeff --report none -a auto >/dev/null 2>&1
    for a in vep snpeff; do
        cmp -s "$WORK/fl_${a}_reformatted.tsv" "$WORK/fl_auto_${a}_reformatted.tsv" \
            && st=PASS || st=FAIL
        row flags "annotation_type/auto_matches_$a" "$st" "$st" "-a auto vs -a $a"
    done

    # -t, -j and -c on both annotators. CSQ and ANN go through separate parsing paths, so
    # a VEP-only pass over these flags proves half the tool (gap closed 2026-09-09).
    for src in "vep:$VEP_FB" "snpeff:$SNPEFF"; do
        ann="${src%%:*}"; vcf="${src#*:}"

        # transcript handling: all three modes must run, and split must expand
        declare -a trows=()
        for mode in first most-severe split; do
            tag="fl_${ann}_t_${mode//-/_}"
            eq flags "$ann/transcript/$mode/exit_code" \
                "$(rc "$BIN" "$vcf" -o "$WORK" -p "$tag" --report none -a "$ann" -t "$mode")" \
                0 ""
            n=$(( $(wc -l < "$WORK/${tag}_reformatted.tsv") - 1 ))
            trows+=("$n")
            row flags "$ann/transcript/$mode/rows" "$n" INFO ""
        done
        eq flags "$ann/transcript/first_and_most_severe_same_rows" \
            "$([[ "${trows[0]}" == "${trows[1]}" ]] && echo 1 || echo 0)" 1 \
            "same row count by construction; they differ in content, not shape"
        eq flags "$ann/transcript/split_expands" \
            "$([[ "${trows[2]}" -ge "${trows[0]}" ]] && echo 1 || echo 0)" 1 \
            "one row per annotation: ${trows[2]} vs ${trows[0]}"
        cmp -s "$WORK/fl_${ann}_t_first_reformatted.tsv" \
               "$WORK/fl_${ann}_t_most_severe_reformatted.tsv" && st=INFO || st=PASS
        # On the SnpEff file these two are identical, and that is the input's doing, not a
        # no-op: SnpEff emits ANN already ordered by putative impact — entry 1 carries the
        # best impact on all 26,852 multi-annotation variants (checked 2026-09-09). VEP
        # does not order CSQ that way unless run with --pick.
        row flags "$ann/transcript/first_differs_from_most_severe" \
            "$([[ "$st" == PASS ]] && echo yes || echo no)" "$st" \
            "no is expected for SnpEff, which pre-sorts ANN by impact"

        # threads: 1 and 4 must be byte-identical, 0 (auto-detect) must run
        "$BIN" "$vcf" -o "$WORK" -p "fl_${ann}_j1" --report none -a "$ann" -j 1 >/dev/null 2>&1
        "$BIN" "$vcf" -o "$WORK" -p "fl_${ann}_j4" --report none -a "$ann" -j 4 >/dev/null 2>&1
        cmp -s "$WORK/fl_${ann}_j1_reformatted.tsv" "$WORK/fl_${ann}_j4_reformatted.tsv" \
            && st=PASS || st=FAIL
        row flags "$ann/threads/j1_identical_to_j4" "$st" "$st" "$(basename "$vcf")"
        eq flags "$ann/threads/j0_exit_code" \
            "$(rc "$BIN" "$vcf" -o "$WORK" -p "fl_${ann}_j0" --report none -a "$ann" -j 0)" \
            0 "auto-detect"

        # -c: the gzipped output must decompress to exactly the uncompressed run
        "$BIN" "$vcf" -o "$WORK" -p "fl_${ann}_plain" --report none -a "$ann" >/dev/null 2>&1
        "$BIN" "$vcf" -o "$WORK" -p "fl_${ann}_gz" --report none -a "$ann" -c >/dev/null 2>&1
        gunzip -c "$WORK/fl_${ann}_gz_reformatted.tsv.gz" > "$WORK/fl_${ann}_gz_reformatted.tsv" \
            2>/dev/null || true
        cmp -s "$WORK/fl_${ann}_plain_reformatted.tsv" "$WORK/fl_${ann}_gz_reformatted.tsv" \
            && st=PASS || st=FAIL
        row flags "$ann/compress/gunzip_matches_plain" "$st" "$st" "-c"

        # the four MAF metadata flags reach the MAF on this annotator's path too
        "$BIN" "$vcf" -o "$WORK" -p "fl_${ann}_meta" --report none -a "$ann" \
            --output-format maf --center TESTCENTER --sample-barcode TESTBARCODE \
            --mutation-status Somatic --sequence-source WXS >/dev/null 2>&1
        meta="$WORK/fl_${ann}_meta_reformatted.maf"
        for pair in Center:TESTCENTER Tumor_Sample_Barcode:TESTBARCODE \
                    Mutation_Status:Somatic Sequence_Source:WXS; do
            col="${pair%%:*}"; expect="${pair##*:}"
            got=$(awk -F'\t' -v c="$col" 'NR==1{for(i=1;i<=NF;i++) if($i==c) k=i; next}
                                          NR==2{print $k; exit}' "$meta")
            eq flags "$ann/maf_metadata/$col" "$got" "$expect" "flag value reaches row 1"
        done
    done

    # both stderr notes still fire
    "$BIN" "$VEP_SMALL" -o "$WORK" -p fl_note1 --report none -a vep 2> "$WORK/fl_note1.err" >/dev/null
    grep -q "multiallelic" "$WORK/fl_note1.err" && st=PASS || st=FAIL
    row flags notes/multiallelic "$st" "$st" "242 sites on the mutect2 file"
    "$BIN" "$VEP_FB" -o "$WORK" -p fl_note2 --report none -a vep -t first 2> "$WORK/fl_note2.err" >/dev/null
    grep -q "more than one transcript annotation" "$WORK/fl_note2.err" && st=PASS || st=FAIL
    row flags notes/multi_transcript "$st" "$st" "-t first on non---pick'd input"

    # a VCF with no annotation at all, and the 115-column multi-sample file
    eq flags unannotated/exit_code \
        "$(rc "$BIN" "$UNANNOTATED" -o "$WORK" -p fl_unann --report none)" 0 "no CSQ, no ANN"
    n=$(( $(wc -l < "$WORK/fl_unann_reformatted.tsv") - 1 ))
    [[ "$n" -gt 0 ]] && st=PASS || st=FAIL
    row flags unannotated/rows "$n" "$st" "still converts, just without annotation columns"
    eq flags multisample/exit_code \
        "$(rc "$BIN" "$MULTISAMPLE" -o "$WORK" -p fl_multi --report none)" 0 "72k variants, 119 columns"
    n=$(( $(wc -l < "$WORK/fl_multi_reformatted.tsv") - 1 ))
    [[ "$n" -gt 0 ]] && st=PASS || st=FAIL
    row flags multisample/rows "$n" "$st" ""
fi

if run bench; then
    # Handicaps, stated once so a number is never quoted without them:
    row bench methodology "median of 5" INFO \
        "vcf2maf gets a decompressed copy; we and bcftools read .gz and pay decompression it does not"
    row bench bcftools_split_vep "floor, not peer" INFO \
        "raw CSQ field extraction only — no MAF classification, no position logic"

    # bench_run <tag> <label> <cmd...> — five runs, median wall seconds, worst peak RSS.
    # Leaves the median in BENCH_MED so the speedup can be computed from it.
    bench_run() {
        local tag="$1" label="$2"; shift 2
        local times=() rss=0 t r med
        for _ in 1 2 3 4 5; do
            /usr/bin/time -l "$@" >/dev/null 2> "$WORK/bench.err" || true
            t=$(awk '/ real /{print $1; exit}' "$WORK/bench.err")
            r=$(awk '/maximum resident set size/{print int($1/1048576); exit}' "$WORK/bench.err")
            times+=("${t:-0}")
            [[ "${r:-0}" -gt "$rss" ]] && rss=${r:-0}
        done
        med=$(printf '%s\n' "${times[@]}" | sort -n | sed -n 3p)
        row bench "$tag/$label" "${med}s" INFO "peak RSS ${rss}MB"
        BENCH_MED=$med
    }

    for vcf in "$VEP_SMALL" "$VEP_FB" "$VEP_BIG"; do
        tag="$(basename "$vcf")"; tag="${tag%%.*}"
        base="$(basename "$vcf")"; base="${base%.gz}"; base="${base%.vcf}"
        plain="$WORK/$base.plain.vcf"
        [[ -s "$plain" ]] || mm bcftools view "$vcf" -o "$plain"

        # Annotation density is what actually drives the margin (2026-08-27), so it is
        # recorded next to every speedup rather than left to be inferred from file size.
        variants=$(grep -vc '^#' "$plain")
        anns=$(mm bcftools query -f '%INFO/CSQ\n' "$vcf" | awk -F',' '{n+=NF} END{print n+0}')
        density=$(awk -v a="$anns" -v v="$variants" 'BEGIN{printf "%.2f", v?a/v:0}')
        row bench "$tag/size" "$variants variants, $anns annotations" INFO \
            "$density annotations per variant"

        bench_run "$tag" bcftools_split_vep \
            "$MM" run -n "$ENV_NAME" bcftools +split-vep "$vcf" -f '%CHROM\t%POS\t%Consequence\n' -d >/dev/null
        bench_run "$tag" ours_tsv \
            "$BIN" "$vcf" -o "$WORK" -p bn_tsv --report none -a vep
        bench_run "$tag" ours_maf \
            "$BIN" "$vcf" -o "$WORK" -p bn_maf --report none -a vep --output-format maf -t most-severe
        ours_maf=$BENCH_MED
        bench_run "$tag" ours_maf_parquet \
            "$BIN" "$vcf" -o "$WORK" -p bn_pq --report none -a vep --output-format maf --parquet
        bench_run "$tag" ours_tsv_gzip \
            "$BIN" "$vcf" -o "$WORK" -p bn_gz --report none -a vep -c
        tid="$(mm bcftools query -l "$plain" | head -1)"
        bench_run "$tag" vcf2maf \
            "$MM" run -n "$ENV_NAME" perl "$VCF2MAF" --inhibit-vep --input-vcf "$plain" \
                --output-maf "$WORK/bn_v2m.maf" --samtools-exec "$(mm which samtools)" \
                --tabix-exec "$(mm which tabix)" --ref-fasta "$FASTA" --ncbi-build GRCh38 \
                --tumor-id "$tid" --vcf-tumor-id "$tid"
        v2m=$BENCH_MED

        row bench "$tag/speedup_vs_vcf2maf" \
            "$(awk -v a="$v2m" -v b="$ours_maf" 'BEGIN{printf "%.2fx", b?a/b:0}')" INFO \
            "at $density annotations per variant — quote the two together or not at all"
    done
fi

if run edge; then
    EDGE="$ROOT/tests/fixtures/edge"
    # <fixture>:<expected tsv exit>:<expected maf exit>
    # no_samples is the one non-zero: a sites-only VCF has no Tumor_Sample_Barcode, and MAF
    # requires one. The error names --sample-barcode, and the remedy is asserted below.
    EDGE_CASES=(
        "empty:0:0" "no_samples:0:1" "no_annotation:0:0" "crlf:0:0"
        "star_alt:0:0" "symbolic_alt:0:0" "malformed_csq:0:0"
    )
    for entry in "${EDGE_CASES[@]}"; do
        f="${entry%%:*}"; rest="${entry#*:}"
        want_tsv="${rest%%:*}"; want_maf="${rest##*:}"
        [[ -s "$EDGE/$f.vcf" ]] || { row edge "$f" missing FAIL "fixture not found"; continue; }
        for fmt in tsv maf; do
            [[ "$fmt" == tsv ]] && want="$want_tsv" || want="$want_maf"
            code=0
            "$BIN" "$EDGE/$f.vcf" -o "$WORK" -p "edge_${f}_${fmt}" --report none \
                --output-format "$fmt" > /dev/null 2> "$WORK/edge_$f.$fmt.err" \
                < /dev/null || code=$?
            eq edge "$f/$fmt/exit_code" "$code" "$want" "degenerate but valid input"
            grep -q "panicked" "$WORK/edge_$f.$fmt.err" && st=FAIL || st=PASS
            row edge "$f/$fmt/no_panic" "$st" "$st" "a panic is never an acceptable failure mode"
        done
    done

    # A header-only VCF exits 0 and writes exactly the column header, derived from the
    # VCF's own ##INFO/##FORMAT declarations. It used to write a zero-byte file, which no
    # reader can open and which a pipeline cannot tell apart from a run that never happened.
    eq edge "empty/tsv_lines" "$(wc -l < "$WORK/edge_empty_tsv_reformatted.tsv" | tr -d ' ')" 1 \
        "the column header and no data rows"

    # The documented remedy for the one non-zero exit must actually work.
    code=0
    "$BIN" "$EDGE/no_samples.vcf" -o "$WORK" -p edge_ns_bc --report none \
        --output-format maf --sample-barcode TESTBC > /dev/null 2>&1 < /dev/null || code=$?
    eq edge "no_samples/maf_with_sample_barcode/exit_code" "$code" 0 \
        "the error message names --sample-barcode; it has to be a real remedy"
    n=$(( $(wc -l < "$WORK/edge_ns_bc_reformatted.maf") - 1 ))
    eq edge "no_samples/maf_with_sample_barcode/rows" "$n" 2 "both sites converted"

    # CRLF input must be byte-identical to the LF copy of the same variants.
    cmp -s "$WORK/edge_crlf_tsv_reformatted.tsv" "$WORK/edge_no_annotation_tsv_reformatted.tsv" \
        && st=PASS || st=FAIL
    row edge "crlf/identical_to_lf" "$st" "$st" "same variants, \\r\\n endings"

    # A truncated CSQ must not lose the other variants on the file.
    n=$(( $(wc -l < "$WORK/edge_malformed_csq_maf_reformatted.maf") - 1 ))
    eq edge "malformed_csq/rows" "$n" 3 "all three variants survive a short CSQ on one of them"
    n=$(grep -c "fewer fields than the header declares" "$WORK/edge_malformed_csq.maf.err" || true)
    [[ "$n" -gt 0 ]] && st=PASS || st=FAIL
    row edge "malformed_csq/warns_on_stderr" "$n" "$st" \
        "warns and continues (2026-09-09); one short CSQ must not cost the rest of the file"

    # ALT=* (spanning deletion): kept, an accepted divergence from vcf2maf, which drops those
    # lines. Symbolic <DEL>/<NON_REF>: dropped from MAF since 2026-09-09 — there is no MAF
    # column that can hold one, and they used to come out typed INS.
    star=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
           $(h["Tumor_Seq_Allele2"])=="*"{n++} END{print n+0}' \
           "$WORK/edge_star_alt_maf_reformatted.maf")
    row edge "star_alt/rows_emitted" "$star" INFO \
        "vcf2maf skips ALT=* entirely; we keep them, typed SNP — accepted divergence (user, 2026-09-09)"
    sym=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
          $(h["Tumor_Seq_Allele2"]) ~ /^</{n++} END{print n+0}' \
          "$WORK/edge_symbolic_alt_maf_reformatted.maf")
    eq edge "symbolic_alt/maf_rows_emitted" "$sym" 0 \
        "<DEL>/<NON_REF> have no MAF representation; they were typed INS before 2026-09-09"
    # The real ALT on the same multiallelic line must survive the skip.
    n=$(( $(wc -l < "$WORK/edge_symbolic_alt_maf_reformatted.maf") - 1 ))
    eq edge "symbolic_alt/maf_rows_kept" "$n" 1 "the one literal ALT still converts"
    # The TSV path is a faithful field dump and keeps them, deliberately.
    n=$(grep -c '<' "$WORK/edge_symbolic_alt_tsv_reformatted.tsv" || true)
    [[ "$n" -gt 0 ]] && st=PASS || st=FAIL
    row edge "symbolic_alt/tsv_still_carries" "$n" "$st" "TSV is a field dump; only MAF drops them"

    # --- the large-file path. The 100k chunking threshold is gone (streaming, 2026-09-09),
    # so this is no longer a separate code path — but it is still where a real WGS file lands.
    big="$WORK/edge_184k.vcf"
    if [[ ! -s "$big" ]]; then
        gunzip -c "$VEP_BIG" | grep '^#' > "$big"
        for _ in 1 2; do gunzip -c "$VEP_BIG" | grep -v '^#'; done >> "$big"
    fi
    want_rows=$(( 92216 * 2 ))
    row edge "large/input_variants" "$want_rows" INFO "the 92k file's data lines, twice"

    t0=$(date +%s)
    rss=$(peak_rss "$BIN" "$big" -o "$WORK" -p edge_big --report html -a vep \
        --output-format tsv --parquet)
    row edge "large/tsv_seconds" "$(( $(date +%s) - t0 ))" INFO "TSV + html report + parquet"
    row edge "large/tsv_peak_rss_mb" "${rss:-?}" INFO "$want_rows variants"
    text="$WORK/edge_big_reformatted.tsv"
    n=$(( $(wc -l < "$text") - 1 ))
    eq edge "large/tsv_rows" "$n" "$want_rows" "one row per variant under -t first, no truncation"
    [[ -n "$(tail -c 1 "$text")" ]] && st=FAIL || st=PASS
    row edge "large/tsv_ends_with_newline" "$st" "$st" "a chunk boundary must not eat the last one"
    # 2026-08-18: --parquet above the old threshold wrote nothing and exited 0.
    pqf="$text.parquet"
    eq edge "large/parquet_written" "$([[ -s "$pqf" ]] && echo yes || echo no)" yes \
        "silently wrote nothing above 100k variants before 2026-08-18"
    pqrows=$(python3 -c "import pyarrow.parquet as pq,sys; print(pq.read_metadata(sys.argv[1]).num_rows)" "$pqf")
    eq edge "large/parquet_rows" "$pqrows" "$n" "matches the sibling text file"
    # Reports used to be skipped entirely above 100k variants; they are computed from
    # counters now, so they must exist at any size.
    eq edge "large/report_written_above_100k" \
        "$([[ -s "$WORK/edge_big_summary.html" ]] && echo yes || echo no)" yes \
        "the pre-streaming code skipped reports above 100k"

    rss=$(peak_rss "$BIN" "$big" -o "$WORK" -p edge_bigmaf --report none -a vep \
        --output-format maf -t most-severe)
    row edge "large/maf_peak_rss_mb" "${rss:-?}" INFO "$want_rows variants, MAF"
    m=$(( $(wc -l < "$WORK/edge_bigmaf_reformatted.maf") - 1 ))
    [[ "$m" -ge "$want_rows" ]] && st=PASS || st=FAIL
    row edge "large/maf_rows" "$m" "$st" "at least one row per variant; multiallelic sites expand"
    rm -f "$big"
fi


# Same table three ways: TSV to machine-read, txt to read, html to send someone.
column -t -s $'\t' "$R" | tee "$OUT/summary.txt"

python3 - "$R" "$OUT/summary.html" "$(date '+%Y-%m-%d %H:%M')" "$want" <<'HTMLPY'
import html, sys
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1]) if l.strip()]
head, body = rows[0], rows[1:]
counts = {}
for r in body:
    counts[r[3]] = counts.get(r[3], 0) + 1
cells = lambda r, tag: "".join(f"<{tag}>{html.escape(c)}</{tag}>" for c in r)
trs = "\n".join(
    f'<tr class="{html.escape(r[3].lower())}">{cells(r, "td")}</tr>' for r in body
)
tally = " · ".join(f"{n} {s}" for s, n in sorted(counts.items()))
open(sys.argv[2], "w").write(f"""<!DOCTYPE html>
<meta charset="utf-8"><title>vcf-reformatter release validation</title>
<style>
 body {{ font: 14px/1.5 -apple-system, system-ui, sans-serif; margin: 2rem; color: #222; }}
 table {{ border-collapse: collapse; width: 100%; }}
 th, td {{ padding: .3rem .6rem; border-bottom: 1px solid #e5e5e5; text-align: left;
           font-variant-numeric: tabular-nums; }}
 th {{ background: #f6f6f6; position: sticky; top: 0; }}
 tr.pass td:nth-child(4) {{ color: #137333; font-weight: 600; }}
 tr.fail td {{ background: #fdecea; }}
 tr.fail td:nth-child(4) {{ color: #b3261e; font-weight: 700; }}
 tr.info td:nth-child(4) {{ color: #8a6d00; }}
 td:last-child {{ color: #666; }}
</style>
<h1>vcf-reformatter release validation</h1>
<p>section <b>{html.escape(sys.argv[4])}</b> · {html.escape(sys.argv[3])} · {tally}</p>
<table><thead><tr>{cells(head, "th")}</tr></thead><tbody>
{trs}
</tbody></table>
""")
HTMLPY

# A section that produced no rows ran no checks. `edge` was declared and empty until
# 2026-09-09, and `validate_release.sh edge` exited 0 the whole time — silence has to be
# an error, not a pass.
for sec in "${SECTIONS[@]}"; do
    run "$sec" || continue
    [[ "$(section_rows "$sec")" -gt 0 ]] || {
        echo "section '$sec' produced no rows — it is declared but empty" >&2
        exit 3
    }
done

# Exit non-zero if anything FAILed, so CI and a human read the same verdict.
! grep -q $'\tFAIL\t' "$R"
