# VCF Reformatter

A Rust command-line tool for parsing and reformatting VCF (Variant Call Format) files, with special support for VEP (Variant Effect Predictor) annotations. This tool flattens complex VCF files into tab-separated values (TSV) format for easier downstream analysis. 
Also incredibly useful for quick checks to your data!

## Features

- 📁 **Compressed file support**: Reads `.vcf.gz` files directly
- 🧬 **VEP annotation parsing**: Extracts and separates CSQ (Consequence) annotations into individual columns
- 📊 **INFO field expansion**: Converts semicolon-separated INFO fields into individual columns
- 👥 **Sample data handling**: Properly parses FORMAT and sample columns
- 📈 **TSV output**: Generates clean, tabular output suitable for deeper analysis
- 🔍 **Header extraction**: Separates header information for reference

## Installation

### Prerequisites

- Rust 1.70+
- Cargo package manager

### Dependencies

Add these to your `Cargo.toml`:

```toml
[dependencies]
flate2 = "1.0"
regex = "1.0"
```

### Build from source

```bash
git clone https://github.com/flalom/vcf-reformatter.git
cd vcf-reformatter
cargo build --release
```

## Usage

### Basic usage

```bash
./target/release/vcf-reformatter input_file.vcf.gz
```

### Example

```bash
./target/release/vcf-reformatter data/sample.freebayes_VEP.ann.vcf.gz
```

### Output files

The tool generates two output files based on the input filename vcf.gz file:

- `{basename}_header.txt`: Contains VCF header information and metadata
- `{basename}_reformatted.tsv`: Main output with flattened VCF data

For input `sample.vcf.gz`, you'll get:
- `sample_header.txt`
- `sample_reformatted.tsv`

The `header` contains the whole header of the VCF file; while the reformatted samples contains the vcf data, in tab separated format `.tsv`.

## Output Format

### Reformatted TSV Structure

The reformatted file contains these column types:

1. **Standard VCF columns**:
    - `CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`

2. **INFO fields** (prefixed with `INFO_`):
    - `INFO_DP`, `INFO_AF`, `INFO_AC`, etc.

3. **CSQ annotations** (prefixed with `CSQ_`):
    - `CSQ_Allele`, `CSQ_Consequence`, `CSQ_SYMBOL`, `CSQ_Gene`, etc.

4. **Sample data** (format: `{SAMPLE_NAME}_{FORMAT_KEY}`):
    - `SAMPLE1_GT`, `SAMPLE1_DP`, `SAMPLE1_AD`, etc.

### Example output

```tsv
CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO_DP	INFO_AF	CSQ_Allele	CSQ_Consequence	CSQ_SYMBOL	SAMPLE1_GT	SAMPLE1_DP
chr1	69511	.	A	G	1294.53	.	65	1	G	missense_variant	OR4F5	1/1	65
```

## Project Structure

```
src/
├── lib.rs                      # Library entry point
├── main.rs                     # CLI application entry point
├── read_vcf_gz.rs             # Compressed VCF file reader
├── essentials_fields.rs       # Core VCF data structures
├── extract_csq_and_csq_names.rs # CSQ annotation extraction
├── get_info_from_header.rs    # Header parsing utilities
├── reformat_vcf.rs            # Main reformatting logic
└── extract_sample_info.rs     # Sample/FORMAT field parsing
```

## Key Components

### VCF Reading (`read_vcf_gz.rs`)
- Handles gzip-compressed VCF files
- Separates header, column names, and data lines
- Memory-efficient line-by-line processing

### Data Structures (`essentials_fields.rs`)
- `VcfVariant`: Represents a single variant record
- `MafRecord`: For potential MAF format support

### CSQ Processing (`extract_csq_and_csq_names.rs`)
- Extracts VEP CSQ annotations
- Removes CSQ from original INFO field
- Supports multiple transcript annotations

### Reformatting (`reformat_vcf.rs`)
- `ReformattedVcfRecord`: Flattened variant representation
- Separates INFO fields into individual key-value pairs
- Handles CSQ field expansion using header format information

### Sample Data (`extract_sample_info.rs`)
- `ParsedFormatSample`: Structured sample data
- Maps FORMAT keys to sample values
- Supports multiple samples per variant

## Supported VCF Features

- ✅ Standard VCF format (v4.x)
- ✅ Gzip compressed files
- ✅ VEP CSQ annotations
- ✅ Multiple samples
- ✅ Complex INFO fields
- ✅ FORMAT/genotype data
- ⚠️ Limited support for structural variants

## Performance

- **Memory efficient**: Processes files line-by-line
- **Fast parsing**: Optimized for large genomics files
- **Concurrent safe**: No shared state between records

## Use Cases

- 🔬 **Genomics research**: Converting VCF files for statistical analysis
- 📊 **Data analysis**: Preparing variant data for R, Python, or Excel
- 🧪 **Pipeline integration**: Standardizing VCF output format
- 📈 **Visualization**: Creating data suitable for plotting tools

## Contributing

1. Fork the repository
2. Create a feature branch: `git checkout -b feature-name`
3. Make your changes
4. Add tests for new functionality
5. Commit your changes: `git commit -am 'Add feature'`
6. Push to the branch: `git push origin feature-name`
7. Submit a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this tool in your research, please cite:

```
VCF Reformatter: A Rust tool for flattening VCF files with VEP annotations
```

## Acknowledgments

- VCF format specification contributors
- Rust community