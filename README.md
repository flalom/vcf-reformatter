# VCF Reformatter: What is it?

Did it ever happen that you had VCF files and you wanted to have a look at the data as you would do with a normal table? `VCF Reformatter` is here for your rescue!

A Rust command-line tool for parsing and reformatting VCF (Variant Call Format) files, with support for VEP (Variant Effect Predictor) and SnpEff annotations. This tool flattens complex VCF files into tab-separated values (TSV) format for easier downstream analysis.
Also incredibly useful for quick checks to your data!

# VCF Reformatter

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Rust](https://img.shields.io/badge/rust-1.70+-blue.svg)](https://www.rust-lang.org)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)]()
[![Performance](https://img.shields.io/badge/performance-9k--17k%20variants%2Fsec-green.svg)]()
[![Release](https://img.shields.io/github/v/release/flalom/vcf-reformatter)](https://github.com/flalom/vcf-reformatter/releases)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-purple.svg?style=flat)](https://anaconda.org/bioconda/vcf-reformatter)
[![Conda](https://anaconda.org/bioconda/vcf-reformatter/badges/version.svg)](https://anaconda.org/bioconda/vcf-reformatter)

**Transform complex VCF files into clean, analyzable tables with ease**

*A high-performance Rust tool for flattening VCF files with intelligent VEP and SnpEff annotation handling*

</div>

---

## 🚀 Quick Start

```` bash
# Download binary from releases (easiest! You download and use it)
wget https://github.com/flalom/vcf-reformatter/releases/latest/download/vcf-reformatter-v0.7.5-linux-x86_64
chmod +x vcf-reformatter-v0.7.5-linux-x86_64

# Transform your VCF file  
./vcf-reformatter-v0.7.5-linux-x86_64 sample.vcf.gz

# Generate MAF output (validated against vcf2maf for VEP input)
./vcf-reformatter-v0.7.5-linux-x86_64 sample.vcf.gz --output-format maf
````
OR Via Bioconda
```bash
conda install -c bioconda vcf-reformatter
# or
# mamba install vcf-reformatter -c bioconda
```
OR install from [crates.io](https://crates.io/crates/vcf-reformatter):
```bash
cargo install vcf-reformatter
```
OR build from source (you need Rust toolchain):
```` bash
git clone https://github.com/flalom/vcf-reformatter.git
cd vcf-reformatter
cargo build --release
./target/release/vcf-reformatter sample.vcf.gz
````
## 📋 MAF Output Status

**VEP-annotated input: validated.** MAF output from VEP (CSQ) input is checked column-by-column
against [vcf2maf](https://github.com/mskcc/vcf2maf) on three real files (12k, 34k and 92k
variants) in all three transcript modes, with **zero unexplained differences**. The 50 columns and
their order match `vcf2maf.pl`'s own header. The remaining differences are named, counted and
deliberate. See [Known divergences from vcf2maf](#known-divergences-from-vcf2maf).

**SnpEff-annotated input: still beta.** vcf2maf cannot parse SnpEff's `ANN` field at all, so there
is no ground truth to validate the SnpEff MAF path against yet. Its structure is checked (50
columns, correct header, depth invariants hold) but its field-level accuracy is not independently
confirmed. That validation lands in v0.8.0.

**Multi-sample VCFs:** pass `--tumor-id` (and `--normal-id` if you have a matched normal). Without
`--tumor-id` the first sample declaring `DP` is used, which is a guess.


## 🎯 Why VCF Reformatter?

**The Problem:** VCF files are notoriously difficult to analyze. Complex nested annotations, semicolon-separated INFO fields, and multi-transcript VEP annotations make downstream analysis a nightmare.

**The Solution:** VCF Reformatter flattens everything into clean, readable TSV format that works with Excel, R, Python, and any analysis tool (⚠️ beware Excel auto-correction!).

### Before & After

**Before (Raw VCF):**
```
chr1  69511  .  A  G  1294.53  .  DP=65;AF=1;CSQ=G|missense_variant|MODERATE|OR4F5|ENSG00000186092...
```
**After (Reformatted TSV):**
```
CHROM  POS    REF  ALT  QUAL     INFO_DP  INFO_AF  CSQ_Allele  CSQ_Consequence      CSQ_SYMBOL
chr1   69511  A    G    1294.53  65       1        G           missense_variant     OR4F5
```

## ✨ Key Features

| Feature                                 | Description                                      | Benefit                                              |
|-----------------------------------------|--------------------------------------------------|------------------------------------------------------|
| 🧬 **VEP/SnpEff Annotation Parsing**    | Intelligent handling of CSQ/ANN annotations with correct field mapping for both | No more manual parsing of complex VEP/SnpEff output  |
| 👀 **Automatic Annotation Recognition** | Automatic detection of CSQ/ANN annotations       | Saving even more time now for both VEP and SnpEff    |
| 🔀 **Smart Transcript Handling**        | Annotator's first entry, most severe, or split    | Choose the analysis approach that fits your needs    |
| 🚀 **Parallel Processing**              | Multi-threaded processing, 9k-17k variants/sec on real annotated files | Process large cohorts in minutes, not hours          |
| 📁 **Native Compression**               | Direct `.vcf.gz` reading & gzip output           | Seamless workflow with compressed/uncompressed files |
| 🎯 **Production Ready**                 | Comprehensive error handling & logging           | Reliable for automated pipelines                     |
| 📋 **Summary Reports**                  | Per-chromosome variant distribution and processing stats | QC and reproducibility for pipeline logs             |
| 📦 **Parquet Output**                   | Apache Parquet copy alongside the text output, ZSTD-compressed, typed columns | 3-6x smaller than raw text, 38x faster column reads in polars/DuckDB/R |
| 🐳 **Container Support**                | Docker & Singularity ready                       | Deploy anywhere, from laptops to HPC clusters        |

---

## 📦 Installation

### Option 1: Download Pre-compiled Binaries (Easiest!)
**No Rust installation required** - just download and run:

1. **Go to [Releases](https://github.com/flalom/vcf-reformatter/releases/latest)**
2. **Download the binary for your platform:**
    - `vcf-reformatter-v0.7.5-linux-x86_64` → **Linux** (most users)
    - `vcf-reformatter-v0.7.5-linux-x86_64-static` → **HPC clusters** (works everywhere)
    - `vcf-reformatter-v0.7.5-windows-x86_64.exe` → **Windows**
    - `vcf-reformatter-v0.7.5-macos-x86_64` → **Intel Mac**
    - `vcf-reformatter-v0.7.5-macos-arm64` → **Apple Silicon Mac** (all M-series)

3. **Make executable and run:**
````bash
# Linux/Mac
chmod +x vcf-reformatter-*
./vcf-reformatter-* --help

# Windows
# Just double-click or run from command prompt
# C++ might be required, if not already installed
````

### Option 2: **Build from Source**
````bash
git clone https://github.com/flalom/vcf-reformatter.git
cd vcf-reformatter
cargo build --release
````

### Option 3: Docker
```shell script
# Build the container
docker build -t vcf-reformatter .

# Run with your data
docker run --rm -v $(pwd):/data vcf-reformatter /data/sample.vcf.gz
```
### Option 4: Singularity
```shell script
# Build Singularity image
singularity build vcf-reformatter.sif Singularity

# Run on HPC cluster
singularity run --bind $PWD:/data vcf-reformatter.sif /data/sample.vcf.gz -j 16
```

## 🛠️ Usage

### Basic Usage
```shell script
# Simple conversion
vcf-reformatter input.vcf.gz

# Most severe consequence only (recommended for analysis)
vcf-reformatter input.vcf.gz -t most-severe

# All transcripts in separate rows (comprehensive)
vcf-reformatter input.vcf.gz -t split

# Read from stdin: plain or gzipped, auto-detected
bcftools view -f PASS input.vcf.gz | vcf-reformatter - -p filtered
cat input.vcf.gz | vcf-reformatter -
```
### Annotation Type Detection
```shell script
# Auto-detect annotation type (recommended)
vcf-reformatter input.vcf.gz -a auto

# Force VEP processing
vcf-reformatter vep_annotated.vcf.gz -a vep -t most-severe

# Force SnpEff processing  
vcf-reformatter snpeff_annotated.vcf.gz -a snpeff -t most-severe
```
### Summary and Parquet Output
```shell script
# Generate a processing report (HTML by default)
vcf-reformatter input.vcf.gz

# Plain-text report instead of HTML
vcf-reformatter input.vcf.gz --report txt

# Both formats from one pass
vcf-reformatter input.vcf.gz --report html,txt

# Collect reports somewhere else than the data output
vcf-reformatter input.vcf.gz -o results/ --report-dir reports/

# No report file
vcf-reformatter input.vcf.gz --report none

# Write a Parquet copy alongside the TSV
vcf-reformatter input.vcf.gz --parquet

# Same for MAF: writes both the .maf and the .maf.parquet
vcf-reformatter input.vcf.gz --output-format maf --parquet

# Parquet is on by default; build without it with
cargo build --release --no-default-features
```
### MAF Metadata and Sample Selection
```shell script
# Name the tumor sample explicitly. Without this the first sample declaring DP
# is used, which is a guess and wrong on most multi-sample VCFs
vcf-reformatter tumor_only.vcf.gz --output-format maf --tumor-id TUMOR_A

# Tumor/normal pair: t_depth/t_ref_count/t_alt_count come from one sample,
# n_depth/n_ref_count/n_alt_count and Matched_Norm_Sample_Barcode from the other
vcf-reformatter paired.vcf.gz --output-format maf \
  --tumor-id B487_1_V --normal-id B487_1_cOM

# Sample names must match the #CHROM line exactly. Check them first:
bcftools query -l paired.vcf.gz

# Metadata a VCF cannot supply. These columns stay empty unless you set them
vcf-reformatter sample.vcf.gz --output-format maf \
  --center MySeqCenter \
  --sample-barcode TCGA-AB-1234-01A \
  --mutation-status Somatic \
  --sequence-source WXS
```
A `--tumor-id` naming a sample that is not in the VCF yields empty depth columns rather than
falling back to the pooled `INFO` counts, so a typo shows up as missing data rather than as the
wrong sample's read counts.

### Advanced Usage
```shell script
# High-performance processing with compression
vcf-reformatter large_cohort.vcf.gz \
  --transcript-handling most-severe \
  --threads 0 \
  --compress \
  --output-dir results/ \
  --prefix my_analysis \
  --verbose

# Optimized for HPC environments
vcf-reformatter huge_dataset.vcf.gz -t most-severe -j 32 -o /scratch/results/ -c -v
```
### Complete Options
```
Usage: vcf-reformatter [OPTIONS] <INPUT_FILE>

Arguments:
  <INPUT_FILE>  Input VCF file (supports .vcf.gz)

Options:
  --output-format <FORMAT>     Output format [default: tsv] 
                               [values: tsv, maf]
  --center <CENTER>            Sequencing center for MAF output  
  --ncbi-build <BUILD>         Genome build 
                               [default: GRCh38]
  --sample-barcode <BARCODE>   Sample identifier for MAF output
  -t, --transcript-handling <MODE>  How to handle multiple transcripts
                                   [default: first]
                                   [values: most-severe, first, split]
                                   first: the annotator's own first entry, unranked
                                   most-severe: ranked by consequence severity
                                   split: one row per transcript
  -a, --annotation-type <N>        Which annotations to parse VEP/SnpEff
                                   [default: auto]
                                   [values: snpeff, vep, auto]
  -j, --threads <N>                Thread count (0 = auto-detect) [default: 1]
  -o, --output-dir <DIR>           Output directory [default: current]
  -p, --prefix <PREFIX>            Output file prefix [default: input filename]
  -c, --compress                   Compress output with gzip
      --report <FORMAT>            Report format [default: html]
                                    [values: html, txt, none]
                                    html: stat cards + SIFT/PolyPhen/Impact charts
                                    txt:  plain text, no damage breakdown
                                    comma-separate for both: --report html,txt
                                    none anywhere in the list wins
      --report-dir <DIR>           Directory for the report [default: --output-dir]
      --parquet                    Also write an Apache Parquet copy alongside
                                   the text output (mutually exclusive with -c)
      --tumor-id <NAME>            Tumor sample name as it appears in #CHROM.
                                   Source of t_depth/t_ref_count/t_alt_count
      --normal-id <NAME>           Matched normal sample name. Populates
                                   n_depth/n_ref_count/n_alt_count
      --mutation-status <STATUS>   MAF Mutation_Status, e.g. Somatic. Empty
                                   unless set. A VCF does not state this
      --sequence-source <SOURCE>   MAF Sequence_Source, e.g. WXS. Empty unless
                                   set. A VCF does not state this either
  -v, --verbose                    Detailed performance statistics
  -h, --help                       Show help
  -V, --version                    Show version

Pass `-` as <INPUT_FILE> to read the VCF from stdin (plain or gzipped, auto-detected).
```

## 🧬 Transcript Handling Modes

VCF files with VEP annotations often contain multiple transcript annotations per variant. Choose the strategy that fits your analysis:

### 🎯 Most Severe (`--transcript-handling most-severe`)
**Best for:** Clinical analysis, variant prioritization
```shell script
vcf-reformatter input.vcf.gz -t most-severe

# for maf output
vcf-reformatter input.vcf.gz -t most-severe --output-format maf
```
Selects the transcript with the most severe consequence (stop_gained > missense_variant > synonymous, etc.)

### ⚡ First Only (`--transcript-handling first`) *[Default]*
**Best for:** Annotator-picked input, reproducibility, performance-critical workflows
```shell script
vcf-reformatter input.vcf.gz  # Uses first transcript by default
```

Keeps whichever annotation the annotator listed first, without re-ranking it.

When your VCF carries **one transcript annotation per variant** (VEP run with `--pick`, or a
caller that annotates a single transcript, which is what many Mutect2 pipelines produce), that
entry *is* the annotator's own selection, and `first` reports exactly what the annotator decided.
This is why it is the default: it passes the upstream choice through unchanged, and the same input
always gives the same output.

When a variant carries **several transcript annotations**, `first` takes literal list order, which
is not severity order. The consequence reported is then the first one listed, which need not be the
most damaging one present. That is expected behaviour, and vcf-reformatter prints a note on
stderr when it sees such sites so the choice is never silent. Use `-t most-severe` if you want the
annotations ranked by consequence severity, or `-t split` to keep every transcript.

### 📊 Split All (`--transcript-handling split`)
**Best for:** Comprehensive analysis, transcript-level studies
```shell script
vcf-reformatter input.vcf.gz -t split
```
Creates separate rows for each transcript (most detailed output)

## 🧾 MAF Output

```shell script
# Single-sample tumor-only
vcf-reformatter input.vcf.gz --output-format maf -t most-severe

# Multi-sample: name the samples, do not let the tool guess
vcf-reformatter tn.vcf.gz --output-format maf --tumor-id TUMOR --normal-id NORMAL

# Add the metadata a VCF cannot supply
vcf-reformatter input.vcf.gz --output-format maf \
  --center MySeqCenter --sequence-source WXS --mutation-status Somatic
```

### Columns

50 columns: `vcf2maf.pl`'s 46 core columns in its exact order and naming, plus four extras this
tool adds at the end: `FILTER`, `QUAL`, `VAF`, `Protein_Position`.

Columns with no source in a single VCF (validation and lab metadata, UUIDs, and the matched-normal
columns when no `--normal-id` is given) are written **empty**, which is how vcf2maf itself renders
a tumor-only MAF. A `.` in the output is a `.` that came from the VCF. Absence and content are
not the same thing.

### Known divergences from vcf2maf

| Column | Difference | Why |
|---|---|---|
| Multiallelic rows | We annotate the allele the annotation actually describes; vcf2maf hangs the first CSQ entry on whichever ALT it emitted, and emits only one row per VCF line | vcf2maf's own guard keys on `ALLELE_NUM`, which VEP writes only under `--allele_number`.
| `Entrez_Gene_Id` | We emit the real ID; vcf2maf emits `0` | vcf2maf does not look it up |
| `all_effects` | Empty | Deliberate gap: it needs the full per-transcript consequence list. Use `-t split` to get the same information as rows instead of one `;`-joined cell |
| `Center` | Defaults to `Unknown_Center`; vcf2maf leaves it empty | Set it with `--center` |
| `Strand` | Always `+` | Per the MAF spec the column is genomic, not the transcript strand (which is still in the TSV's `CSQ_STRAND`) |
| `ALT=*` | Kept as a row typed `SNP`; vcf2maf drops the line | Accepted divergence |
| Symbolic ALTs (`<DEL>`, breakends) | Skipped in MAF, kept verbatim in TSV | MAF has no column that can hold one |

### Reproducing the validation

```shell script
tests/scripts/validate_release.sh          # all sections
tests/scripts/validate_release.sh maf      # MAF parity only
```

Results land in `tests/scripts/results/summary.{tsv,txt,html}`; the script exits non-zero if any
check fails.

## 📈 Performance

### Benchmarks

Median of 5 runs, real VEP-annotated files, MAF conversion against
[vcf2maf](https://github.com/mskcc/vcf2maf) `--inhibit-vep`:

| variants | annotations/variant | vcf-reformatter (MAF) | vcf2maf (Perl) | speedup |
|---|---|---|---|---|
| 12,239 | 1.00 | 0.75s | 1.24s | **1.65x** |
| 34,415 | 2.01 | 2.58s | 5.47s | **2.12x** |
| 92,216 | 12.05 | 12.48s | 97.79s | **7.84x** |


### Memory

Both output paths stream: the VCF is read, converted, written and dropped in chunks, so peak
memory stays roughly constant as the input grows.

| input | TSV | MAF | MAF + `--parquet` |
|---|---|---|---|
| 12,239 variants | 149 MB | 158 MB | 181 MB |
| 34,415 variants | 235 MB | 253 MB | 288 MB |
| 92,216 variants | 330 MB | 397 MB | 442 MB |

Those three files differ in annotation density, so they are not a scaling test. The honest control
is one file against eight copies of itself, same shape, 8x the bytes (46 MB to 362 MB): **TSV 232
MB → 238 MB**, MAF + parquet 290 MB → 485 MB. 

### Internal Optimizations
- **Streaming I/O**: the reader hands back an iterator; both output paths convert and write in 10,000-line chunks and drop them, including `--parquet` (one row group per 65,536 rows)
- **Zero-copy output**: TSV and MAF output use borrowed references (`Cow<str>`) instead of cloning strings, reducing heap allocations per variant
- **Move semantics**: Single-transcript variants (the common case) avoid all string cloning during record construction

### Optimization Tips
```shell script
# Auto-detect optimal thread count
vcf-reformatter input.vcf.gz -j 0

# For files > 10K variants, use parallel processing
vcf-reformatter input.vcf.gz -t most-severe -j 0 -v

# Combine with compression for large outputs
vcf-reformatter input.vcf.gz -t split -j 0 -c -v
```

## 📊 Output Format

### File Structure
- `{prefix}_header.txt`: original VCF header and metadata (TSV mode only)
- `{prefix}_reformatted.tsv`: flattened tabular data (or `_reformatted.maf` with `--output-format maf`)
- `{prefix}_summary.html`: processing report, unless `--report txt|none` (`--report html,txt` writes both `_summary.html` and `_summary.txt`)
- `{prefix}_reformatted.tsv.parquet`: with `--parquet`, written *alongside* the text file, never instead of it

Reading from stdin (`-`) with no `--prefix` names the outputs `stdin_*`.

### Parquet column types
`POS`, `Start_Position`, `End_Position` and the depth columns are `uint64`; `QUAL` and `VAF` are
`double`; everything else is a string. Unpopulated cells are real NULLs, ZSTD-compressed, one row
group per 65,536 rows. Verified lossless against the sibling text output: 0 cell mismatches
across 35.9M cells.

### Column Types
1. **Standard VCF**: `CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`
2. **INFO Fields**: `INFO_DP`, `INFO_AF`, `INFO_AC`, etc.
3. **VEP Annotations**: `CSQ_Allele`, `CSQ_Consequence`, `CSQ_SYMBOL`, `CSQ_Gene`, etc.
3. **SnpEff Annotations**: `ANN_Allele`, `ANN_Annotation_Impact`, `ANN_Gene_Name`, `ANN_Distance`, etc.
4. **Sample Data**: `SAMPLE1_GT`, `SAMPLE1_DP`, `SAMPLE1_AD`, etc.

The INFO and FORMAT columns come from the VCF header's own `##INFO` / `##FORMAT` declarations, so a
field that only appears on later variants still gets a column — a caller that writes `LOF` on 251
of 29,589 variants, or `PGT`/`PID`/`PS` on a third of them, is not silently dropped. A declared
field that never occurs costs one column of `.`.

### Example Output VEP
```
CHROM  POS    ID     REF  ALT  QUAL     FILTER  INFO_DP  CSQ_Consequence      CSQ_SYMBOL  SAMPLE1_GT
chr1   69511  .      A    G    1294.53  PASS    65       missense_variant     OR4F5       1/1
chr1   69761  rs123  C    T    892.15   PASS    42       synonymous_variant   OR4F5       0/1
```

### Example Output SnpEff
```
CHROM  POS    ID     REF  ALT  QUAL     FILTER  INFO_DP  ANN_Annotation          ANN_Gene_Name  SAMPLE1_GT
chr1   69761  rs587   C    T  730  PASS   .     214      synonymous_variant      OR4F5          0/1
chr1   924024  .      A    G  53   PASS   .     409      5_prime_UTR_variant     SAMD11         1/1
```

## 🔧 Integration Examples

### With R
```textmate
# Read compressed output directly
library(data.table)
data <- fread("output_reformatted.tsv.gz")

# Quick variant summary
summary(data$CSQ_Consequence)
```

### With Python
```textmate
import pandas as pd

# Load and analyze
df = pd.read_csv("output_reformatted.tsv.gz", sep="\t", compression="gzip")
df['CSQ_Consequence'].value_counts()
```

### In Workflows
```shell script
# Nextflow pipeline
vcf-reformatter ${vcf} -t most-severe -j ${task.cpus} -o results/ -c

# Snakemake rule
shell: "vcf-reformatter {input.vcf} -t most-severe -j {threads} -o {params.outdir} -c"
```

## 🐳 Container Usage

### Docker
```shell script
# Build once
docker build -t vcf-reformatter .

# Run anywhere
docker run --rm \
  -v $(pwd):/data \
  vcf-reformatter \
  /data/input.vcf.gz \
  -t most-severe -j 4 -o /data/results/ -c
```

### Singularity (HPC)
```shell script
# On HPC cluster
singularity run \
  --bind $PWD:/data \
  --bind /scratch:/scratch \
  vcf-reformatter.sif \
  /data/large_cohort.vcf.gz \
  -t most-severe -j 16 -o /scratch/results/ -c -v
```
## 🧪 Use Cases

| Use Case | Command | Why It Works |
|----------|---------|--------------|
| **Clinical Variant Review** | `vcf-reformatter variants.vcf.gz -t most-severe` | Prioritizes clinically relevant consequences |
| **Population Analysis** | `vcf-reformatter cohort.vcf.gz -t first -j 0 -c` | Fast processing of large cohorts |
| **Transcript Studies** | `vcf-reformatter genes.vcf.gz -t split -v` | Comprehensive transcript-level analysis |
| **Quick Data Exploration** | `vcf-reformatter sample.vcf.gz` | Simple, fast conversion for immediate analysis |
| **HPC Batch Processing** | `vcf-reformatter huge.vcf.gz -t most-severe -j 32 -c` | Optimized for high-performance computing |

## 🚀 What's New in v0.7.5
- ✅ **MAF output validated against vcf2maf.** 50 columns in vcf2maf's own order, 0 unexplained differences across three real VEP files and all three transcript modes
- ✅ **Streaming I/O.** Memory is now flat in file size: 333 MB on a 92k-variant file, down from 2.6 GB
- ✅ **Matched-normal and multi-sample support.** `--tumor-id` / `--normal-id` populate the tumor and normal depth columns from named samples instead of pooled INFO counts
- ✅ **stdin support.** Pass `-` to pipe from `bcftools`; plain or gzipped, auto-detected
- ✅ **HTML report.** `--report html|txt|none`, a self-contained page with per-chromosome SIFT / PolyPhen / IMPACT charts
- ✅ **Parquet alongside text.** ZSTD, typed columns, bounded row groups, streamed
- ✅ **MAF correctness fixes.** HGVS accession stripping, splice-variant protein coordinates, shared REF/ALT prefix trimming, per-allele annotation on multiallelic sites, `dbSNP_RS` from VEP's `Existing_variation`
- ✅ **195 tests** (97 unit + 98 integration) plus a release validation harness in `tests/scripts/`

## Previous Releases
### 🚀 What's New in v0.4.0
- ✅ **Performance: Reduced memory allocations** - Replaced ~30 unnecessary `.clone()` calls with zero-cost borrows using `Cow<str>`, `as_str()`, and move semantics
- ✅ **Comprehensive Testing** - 86 test cases ensure reliability across VEP and SnpEff pipelines

### 🚀 What's New in v0.3.0
- ✅ **MAF Output Support (in Beta⚠️)** - Direct conversion to Mutation Annotation Format
- ✅ **Auto-metadata Detection (in Beta⚠️)** - Extracts center/sample info from VCF headers for MAF
- ✅ **Memory-Efficient Processing (streaming)** - Chunked streaming for large files (>>100K variants)
- ✅ **Enhanced Error Handling** - Better processing of malformed files
- ✅ **Comprehensive Testing** - 70+ test cases ensure reliability

### 🚀 What's New in v0.2.0
- ✅ **SnpEff Support** - Full ANN field parsing with intelligent detection
- ✅ **Smart Auto-Detection** - Automatically identifies VEP vs SnpEff annotations
- ✅ **Enhanced Error Handling** - Better processing of malformed or headerless files

## TODOs
- ~~Add SnpEff support✅~~
- ~~Output MAF format option✅~~
- ~~Reduce `.clone()` allocations for better performance✅~~
- ~~Add `stdin` to combine with other tools, such as `bcftools`✅~~
- ~~Support for multi-sample VCF files in MAF output✅~~
- ~~Streaming MAF output for large files✅~~
- ~~Unified annotation parsing (deduplicate CSQ/ANN code paths)✅~~
- Validate the SnpEff MAF path against a ground truth (v0.8.0)
- Populate `all_effects` (today: use `-t split`)
- Enriched reports: Ti/Tv ratio, variant-type breakdown, filter distribution

## 🤝 Contributing

We welcome contributions! Here's how to get started:

1. **Fork** the repository
2. **Create** a feature branch: `git checkout -b feature-name`
3. **Add tests** for new functionality
4. **Commit** your changes: `git commit -am 'Add feature'`
5. **Push** to the branch: `git push origin feature-name`
6. **Submit** a pull request

### Development Setup
```shell script
git clone https://github.com/flalom/vcf-reformatter.git
cd vcf-reformatter
cargo test  # Run the test suite
cargo run -- data/sample.vcf.gz -v  # Test with sample data
```

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- **VCF Format Contributors** - For the standard that enables genomic data sharing
- **VEP Team** - For the powerful variant annotation framework
- **Rust Community** - For the incredible ecosystem that makes this possible
- **Bioinformatics Community** - For feedback and feature requests

---

## Frequently Asked Questions

### Q: Which transcript handling mode should I use?
- **Clinical analysis**: `--transcript-handling most-severe`
- **Quick exploration**: `--transcript-handling first`
- **Comprehensive analysis**: `--transcript-handling split`

### Q: How does this compare to other VCF tools?
VCF Reformatter is specifically designed for:
- Converting complex VEP/SnpEff annotations to tabular format
- Handling multiple transcripts intelligently
- High-performance parallel processing
- Easy integration with R/Python workflows

### Q: Can I use this in production pipelines?
Yes! VCF Reformatter is designed for production use with:
- Comprehensive error handling
- Docker/Singularity support
- Automated testing
- Stable CLI interface

### Q: What's the difference between TSV and MAF output?
- **TSV**: Direct flattening of VCF fields (default)
- **MAF**: Standardized cancer genomics format for downstream tools. Validated against vcf2maf for VEP input; still beta for SnpEff input

### Q: What if I get out-of-memory errors?
Both paths stream since v0.7.5, so a bigger file does not cost more memory: an input 8x larger
costs the TSV path 3% more. If you still hit a limit, drop `--parquet` (it buffers a
row group) and monitor with `-v`.

___

## 📞 Support

- **📋 Issues**: [GitHub Issues](https://github.com/flalom/vcf-reformatter/issues)
- **📧 Email**: [fl@flaviolombardo.site](mailto:fl@flaviolombardo.site)

---

<div align="center">

**⭐ Star this repo if VCF Reformatter helps your research!**

Made with ❤️ by [Flavio Lombardo](https://github.com/flalom)

</div>


