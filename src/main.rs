//! # VCF Reformatter
//!
//! A fast, parallel VCF (Variant Call Format) parser and reformatter for bioinformatics.
//! Supports VEP and SnpEff annotations with configurable transcript handling.
//!
//! ## Quick Start
//!
//! ```bash
//! # Basic usage with auto-detection
//! vcf-reformatter sample.vcf.gz
//!
//! # Use VEP annotations with most severe consequence
//! vcf-reformatter sample.vcf.gz -a vep -t most-severe -j 4
//!
//! # Use SnpEff annotations with all transcripts
//! vcf-reformatter sample.vcf.gz -a snpeff -t split -o results/
//! ```
use clap::{Parser, ValueEnum};
use essentials_fields::MafRecord;
use reformat_vcf::{
    reformat_vcf_data_with_header, reformat_vcf_data_with_header_parallel, AnnotationType,
    TranscriptHandling,
};
use std::path::Path;
use std::time::Instant;

use std::fs::File;
use std::io::Write;

mod essentials_fields;
mod extract_ann_and_ann_names;
mod extract_csq_and_csq_names;
mod extract_sample_info;
mod get_info_from_header;
mod html_report;
mod read_vcf_gz;
mod reformat_vcf;
mod summary;

#[cfg(feature = "parquet_out")]
mod parquet_writer;

use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::BufWriter;

// Add these helper structs for MAF metadata extraction
#[derive(Debug, Default)]
struct MafMetadata {
    ncbi_build: Option<String>,
    center: Option<String>,
    sample_names: Vec<String>,
    primary_sample: Option<String>,
}
#[derive(Debug, Default)]
struct MafConfig {
    center: String,
    sample_barcode: String,
    ncbi_build: String,
}
#[derive(Debug)]
struct ProcessingTiming {
    read_time: std::time::Duration,
    process_time: std::time::Duration,
    /// `None` when writing is streamed inline with processing and has no
    /// separately measurable duration (e.g. TSV chunked/parallel output).
    write_time: Option<std::time::Duration>,
    /// `None` for MAF output, which has no separate header file.
    header_write_time: Option<std::time::Duration>,
    total_time: std::time::Duration,
}

#[derive(Debug)]
struct ProcessingStats {
    variants_per_sec: f64,
    thread_count: usize,
    use_parallel: bool,
}

#[derive(Parser)]
#[command(
    name = "vcf-reformatter",
    version,
    about = "🧬 Fast VCF file parser and reformatter with VEP and SnpEff annotation support",
    long_about = "A Rust command-line tool for parsing and reformatting VCF (Variant Call Format) files, with support for VEP (Variant Effect Predictor) and SnpEff annotations. This tool flattens complex VCF files into tab-separated values (TSV) or Mutation Annotation Format (MAF) for easier downstream analysis.",
    after_help = "EXAMPLES:
    Basic usage (auto-detect annotation type):
      vcf-reformatter sample.vcf.gz

    Generate MAF output (auto-detects metadata from header):
      vcf-reformatter sample.vcf.gz --output-format maf

    Generate MAF output with manual parameters:
      vcf-reformatter sample.vcf.gz --output-format maf --center TCGA --sample-barcode TCGA-01

    Use VEP annotations with most severe consequence:
      vcf-reformatter sample.vcf.gz -a vep -t most-severe -j 4

    Use SnpEff annotations with all transcripts:
      vcf-reformatter sample.vcf.gz -a snpeff -t split -o results/ -p analysis

    Auto-detect annotation type with parallel processing:
      vcf-reformatter sample.vcf.gz -a auto -j 0 -v

    Complete example with SnpEff and compression:
      vcf-reformatter sample.vcf.gz -a snpeff -t most-severe -j 4 -o results/ -p my_analysis -v --compress"
)]
struct Cli {
    /// Input VCF file (supports .vcf.gz compressed files), or "-" to read from stdin
    #[arg(value_name = "INPUT_FILE")]
    input_file: String,

    /// Annotation type to parse
    #[arg(short = 'a', long = "annotation-type", value_enum, default_value_t = AnnotationTypeCli::Auto)]
    annotation_type: AnnotationTypeCli,

    /// Transcript handling mode
    ///
    /// 'first' (the default) keeps whichever annotation the annotator listed first, without
    /// re-ranking it. On VEP output run with --pick, and on callers that emit a single
    /// transcript per variant, that is the annotator's own selection and there is nothing to
    /// choose between. Where a variant carries several transcript annotations it is literal
    /// list order, not severity order, so the consequence reported need not be the most
    /// damaging one present. Use 'most-severe' to rank the annotations by consequence
    /// severity instead, or 'split' to keep every transcript as its own row.
    #[arg(short = 't', long = "transcript-handling", value_enum, default_value_t = TranscriptHandlingCli::FirstOnly)]
    transcript_handling: TranscriptHandlingCli,

    /// Number of threads to use for parallel processing (0 = auto-detect)
    #[arg(short = 'j', long = "threads", default_value_t = 1)]
    threads: usize,

    /// Output directory (default: current directory)
    #[arg(short = 'o', long = "output-dir")]
    output_dir: Option<String>,

    /// Prefix for output files (default: input filename)
    #[arg(short = 'p', long = "prefix")]
    prefix: Option<String>,

    /// Verbose output
    #[arg(short = 'v', long = "verbose")]
    verbose: bool,

    /// Compress output files with gzip
    #[arg(short = 'c', long = "compress")]
    compress: bool,

    /// Output format: tsv or maf
    #[arg(long, value_enum, default_value_t = OutputFormatCli::Tsv)]
    output_format: OutputFormatCli,

    /// Center name for MAF output (auto-detected from header if not provided)
    #[arg(long)]
    center: Option<String>,

    /// NCBI build for MAF output (auto-detected from header, defaults to GRCh38)
    #[arg(long, default_value = "GRCh38")]
    ncbi_build: String,

    /// Sample barcode for MAF output (auto-detected from header if not provided)
    #[arg(long)]
    sample_barcode: Option<String>,

    /// Mutation_Status for MAF output, e.g. Somatic or Germline. A VCF does not state this,
    /// so the column is left empty unless you set it.
    #[arg(long)]
    mutation_status: Option<String>,

    /// Sequence_Source for MAF output, e.g. WXS or WGS. A VCF does not state this, so the
    /// column is left empty unless you set it.
    #[arg(long)]
    sequence_source: Option<String>,

    /// Name of the tumor sample, exactly as it appears in the VCF's #CHROM line.
    ///
    /// The t_depth / t_ref_count / t_alt_count columns are read from this sample. Without it
    /// the first sample declaring DP is used, which is a guess — set this on any multi-sample
    /// VCF, where guessing can report the wrong sample's read counts.
    #[arg(long)]
    tumor_id: Option<String>,

    /// Name of the matched normal sample, exactly as it appears in the VCF's #CHROM line.
    /// Populates the n_depth / n_ref_count / n_alt_count and Matched_Norm_Sample_Barcode
    /// columns, which are otherwise left empty.
    #[arg(long)]
    normal_id: Option<String>,

    /// Report format: html (default), txt, or none (no report generated)
    #[arg(long, value_enum, default_value_t = ReportFormatCli::Html)]
    report: ReportFormatCli,

    /// Also write an Apache Parquet copy of the output, alongside the text file
    #[arg(long)]
    parquet: bool,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum AnnotationTypeCli {
    /// VEP annotations (CSQ field)
    #[value(name = "vep")]
    Vep,
    /// SnpEff annotations (ANN field)
    #[value(name = "snpeff")]
    SnpEff,
    /// Auto-detect from header
    #[value(name = "auto")]
    Auto,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum TranscriptHandlingCli {
    /// Extract only the most severe consequence for each variant
    #[value(name = "most-severe")]
    MostSevere,
    /// Keep the annotation the annotator listed first, unranked (default, fastest)
    #[value(name = "first")]
    FirstOnly,
    /// Split every transcript into separate rows
    #[value(name = "split")]
    SplitRows,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum OutputFormatCli {
    /// Tab-separated values (default)
    #[value(name = "tsv")]
    Tsv,
    /// Mutation Annotation Format
    #[value(name = "maf")]
    Maf,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum ReportFormatCli {
    /// Self-contained HTML report with stat cards and damage-metric charts (default)
    #[value(name = "html")]
    Html,
    /// Plain-text report (legacy format; no damage breakdown, that is HTML-only)
    #[value(name = "txt")]
    Txt,
    /// No report file
    #[value(name = "none")]
    None,
}

impl From<AnnotationTypeCli> for AnnotationType {
    fn from(cli: AnnotationTypeCli) -> Self {
        match cli {
            AnnotationTypeCli::Vep => AnnotationType::Vep,
            AnnotationTypeCli::SnpEff => AnnotationType::SnpEff,
            AnnotationTypeCli::Auto => AnnotationType::Auto,
        }
    }
}

impl From<TranscriptHandlingCli> for TranscriptHandling {
    fn from(cli: TranscriptHandlingCli) -> Self {
        match cli {
            TranscriptHandlingCli::MostSevere => TranscriptHandling::MostSevere,
            TranscriptHandlingCli::FirstOnly => TranscriptHandling::FirstOnly,
            TranscriptHandlingCli::SplitRows => TranscriptHandling::SplitRows,
        }
    }
}

// Helper function to extract MAF metadata from VCF header
fn extract_maf_metadata_from_header(header: &str, column_line: &str) -> MafMetadata {
    let mut metadata = MafMetadata {
        center: Some("Unknown_Center".to_string()),
        ncbi_build: Some("GRCh38".to_string()),
        sample_names: Vec::new(),
        primary_sample: None,
    };

    for line in header.lines() {
        if line.starts_with("##reference=") {
            // Extract reference genome info
            if line.contains("GRCh37") || line.contains("hg19") {
                metadata.ncbi_build = Some("GRCh37".to_string());
            } else if line.contains("GRCh38") || line.contains("hg38") {
                metadata.ncbi_build = Some("GRCh38".to_string());
            }
        } else if line.starts_with("##center=") {
            // Look for center info (rare but possible)
            metadata.center = Some(line.trim_start_matches("##center=").to_string());
        } else if line.starts_with("##source=") {
            // Don't put tool info in center - we'll handle this separately
            // This could be used for sequencer field instead
        }
    }

    // Extract sample names from column line
    if !column_line.is_empty() {
        let columns: Vec<&str> = column_line.split('\t').collect();
        if columns.len() > 9 {
            metadata.sample_names = columns[9..].iter().map(|s| s.to_string()).collect();
            metadata.primary_sample = metadata.sample_names.first().cloned();
        }
    }

    metadata
}

// Helper function to validate MAF arguments with auto-detection
fn validate_maf_arguments(cli: &Cli, metadata: &MafMetadata) -> Result<MafConfig, String> {
    if cli.output_format != OutputFormatCli::Maf {
        return Ok(MafConfig::default()); // Not needed for TSV output
    }

    // Determine center
    let center = if let Some(ref user_center) = cli.center {
        user_center.clone()
    } else if let Some(ref detected_center) = metadata.center {
        println!(
            "🔍 Auto-detected center from VCF header: {}",
            detected_center
        );
        detected_center.clone()
    } else {
        return Err("Center is required for MAF output. Please use --center <CENTER_NAME> or ensure your VCF has center information in the header".to_string());
    };
    // Determine sample barcode
    let sample_barcode = if let Some(ref user_sample) = cli.sample_barcode {
        user_sample.clone()
    } else if let Some(ref detected_sample) = metadata.primary_sample {
        println!(
            "🔍 Auto-detected sample barcode from VCF header: {}",
            detected_sample
        );
        detected_sample.clone()
    } else {
        return Err("Sample barcode is required for MAF output. Use --sample-barcode or ensure your VCF has sample columns".to_string());
    };

    // Determine NCBI build (prefer user input, then detected, then CLI default)
    let ncbi_build = if cli.ncbi_build != "GRCh38" {
        // User explicitly set it to something else
        cli.ncbi_build.clone()
    } else if let Some(ref detected_build) = metadata.ncbi_build {
        println!(
            "🔍 Auto-detected NCBI build from VCF header: {}",
            detected_build
        );
        detected_build.clone()
    } else {
        cli.ncbi_build.clone() // Use CLI default
    };

    Ok(MafConfig {
        center,
        sample_barcode,
        ncbi_build,
    })
}

/// Lines pulled from the reader before each parse/convert/write cycle. The knob that decides
/// peak memory: ~200MB on a dense file at 10k, ~500MB at 25k.
const STREAM_CHUNK: usize = 10_000;

/// Take up to `n` lines. An empty result means the input is exhausted.
fn next_chunk(
    lines: &mut dyn Iterator<Item = std::io::Result<String>>,
    n: usize,
) -> std::io::Result<Vec<String>> {
    let mut chunk = Vec::with_capacity(n);
    while chunk.len() < n {
        match lines.next() {
            Some(line) => chunk.push(line?),
            None => break,
        }
    }
    Ok(chunk)
}

fn warn_about_multiallelic(count: usize) {
    if count > 0 {
        eprintln!(
            "⚠️  Warning: {count} multiallelic site(s) detected (ALT field lists more than one allele)."
        );
        eprintln!("   Each alternate allele is expanded into its own output record.");
        eprintln!("   To get one ALT per line instead, normalize first: bcftools norm -m- <input>");
    }
}

fn note_multi_transcript(count: usize) {
    if count > 0 {
        eprintln!(
            "ℹ️  Note: {count} site(s) carry more than one transcript annotation; -t first is in use."
        );
        eprintln!("   'first' reports the annotation listed first by the annotator, without re-ranking it,");
        eprintln!(
            "   so at those sites the consequence shown need not be the most damaging one present."
        );
        eprintln!(
            "   -t most-severe ranks by consequence severity; -t split keeps every transcript."
        );
    }
}

fn warn_about_malformed_annotations(count: usize, key: &str) {
    if count > 0 {
        eprintln!(
            "⚠️  Warning: {count} {key} annotation entr(ies) carry fewer fields than the header declares."
        );
        eprintln!(
            "   Those entries are read as far as they go; the missing fields come out empty,"
        );
        eprintln!("   so affected variants are still converted — expect Unknown/blank annotation columns.");
    }
}

/// The annotation key in use and how many `|`-separated fields its header declares, for the
/// malformed-entry note. CSQ wins if both are present, matching the parser's own order.
fn annotation_format(header: &str) -> Option<(&'static str, usize)> {
    get_info_from_header::extract_csq_format_from_header(header)
        .filter(|f| !f.is_empty())
        .map(|f| ("CSQ", f.len()))
        .or_else(|| {
            get_info_from_header::extract_ann_format_from_header(header)
                .filter(|f| !f.is_empty())
                .map(|f| ("ANN", f.len()))
        })
}

fn main() {
    let cli = Cli::parse();

    // Validate --parquet + --compress conflict
    if cli.parquet && cli.compress {
        eprintln!(
            "Error: Parquet has built-in compression; --compress is not needed with --parquet"
        );
        std::process::exit(1);
    }

    #[cfg(not(feature = "parquet_out"))]
    if cli.parquet {
        eprintln!("Error: Parquet support not compiled. Rebuild with: cargo build --release --features parquet_out");
        std::process::exit(1);
    }

    // Validate input file ("-" means read from stdin, so it has no path to check)
    if cli.input_file != "-" && !Path::new(&cli.input_file).exists() {
        eprintln!("❌ Error: File '{}' not found", cli.input_file);
        std::process::exit(1);
    }

    // Parse thread count
    let thread_count = if cli.threads == 0 {
        num_cpus::get()
    } else {
        cli.threads
    };

    let transcript_handling = TranscriptHandling::from(cli.transcript_handling);
    let annotation_type = AnnotationType::from(cli.annotation_type);
    let use_parallel = thread_count > 1;

    // Start total timing
    let total_start = Instant::now();

    // Read VCF file first to extract metadata
    println!("📖 Reading VCF file...");
    let read_start = Instant::now();
    let stream = match read_vcf_gz::open_vcf(&cli.input_file) {
        Ok(stream) => stream,
        Err(e) => {
            eprintln!("❌ Error reading VCF file: {e}");
            std::process::exit(1);
        }
    };
    let read_vcf_gz::VcfStream {
        header,
        columns_title,
        mut lines,
    } = stream;

    // Extract metadata from header for MAF
    let maf_metadata = extract_maf_metadata_from_header(&header, &columns_title);

    // Validate MAF-specific arguments with auto-detection
    let maf_config = match validate_maf_arguments(&cli, &maf_metadata) {
        Ok(config) => config,
        Err(e) => {
            eprintln!("❌ Error: {}", e);
            std::process::exit(1);
        }
    };

    // Print startup information
    print_startup_info(
        &cli,
        thread_count,
        use_parallel,
        transcript_handling,
        annotation_type,
    );

    // Set thread pool size if using parallel processing
    if use_parallel {
        if let Err(e) = rayon::ThreadPoolBuilder::new()
            .num_threads(thread_count)
            .build_global()
        {
            eprintln!("⚠️  Warning: Could not set thread pool size: {e}");
            eprintln!("   Continuing with default thread pool...");
        }
    }

    // Use the function to tell the user if VEP or SNPEFF were detected
    if let Err(e) = detect_and_print_annotation_type(&header, cli.annotation_type) {
        eprintln!("❌ {}", e);
        std::process::exit(1);
    }

    let read_time = read_start.elapsed();
    println!("✅ Header read in {read_time:.2?}");
    println!("   📑 Header lines: {}", header.matches('\n').count());
    println!();

    // Branch based on output format
    match cli.output_format {
        OutputFormatCli::Tsv => {
            let (header_file, reformatted_file) = generate_output_filenames(&cli);

            println!("🔄 Processing VCF data...");
            let process_start = Instant::now();

            let output_file = match File::create(&reformatted_file) {
                Ok(file) => file,
                Err(e) => {
                    eprintln!("❌ Error creating output file: {e}");
                    std::process::exit(1);
                }
            };
            let mut writer: Box<dyn Write> = if cli.compress {
                Box::new(BufWriter::new(GzEncoder::new(
                    output_file,
                    Compression::default(),
                )))
            } else {
                Box::new(BufWriter::new(output_file))
            };

            if cli.verbose {
                if use_parallel {
                    println!("   🚀 Using parallel processing with {thread_count} threads...");
                } else {
                    println!("   🐌 Using sequential processing...");
                }
            }

            let mut input_variants = 0usize;
            let mut output_record_count = 0usize;
            let mut input_chrom_counts: indexmap::IndexMap<String, usize> =
                indexmap::IndexMap::new();
            let mut output_chrom_counts: indexmap::IndexMap<String, usize> =
                indexmap::IndexMap::new();
            let mut damage_breakdowns: Vec<summary::DamageBreakdown> = Vec::new();
            let mut multiallelic_count = 0usize;
            let mut multi_transcript_count = 0usize;
            let mut malformed_count = 0usize;
            let ann_format = annotation_format(&header);
            let mut progress_marks = 0usize;
            // Columns come from the first record, exactly as the buffered path always has.
            let mut headers: Vec<String> = Vec::new();
            // Parquet is written a chunk at a time too; the sink opens once the first
            // chunk has settled the column list.
            #[cfg(feature = "parquet_out")]
            let parquet_file = format!("{reformatted_file}.parquet");
            #[cfg(feature = "parquet_out")]
            let mut parquet_sink: Option<parquet_writer::ParquetSink> = None;

            loop {
                let chunk = match next_chunk(lines.as_mut(), STREAM_CHUNK) {
                    Ok(chunk) => chunk,
                    Err(e) => {
                        eprintln!("❌ Error reading VCF file: {e}");
                        std::process::exit(1);
                    }
                };
                if chunk.is_empty() {
                    break;
                }
                input_variants += chunk.len();
                for line in &chunk {
                    if let Some(chrom) = line.split('\t').next() {
                        *input_chrom_counts.entry(chrom.to_string()).or_insert(0) += 1;
                    }
                }

                let chunk_multiallelic = summary::count_multiallelic_sites(&chunk);
                if chunk_multiallelic > 0 && multiallelic_count == 0 {
                    warn_about_multiallelic(chunk_multiallelic);
                }
                multiallelic_count += chunk_multiallelic;
                if let Some((key, n_fields)) = ann_format {
                    let chunk_malformed =
                        summary::count_malformed_annotation_entries(&chunk, key, n_fields);
                    if chunk_malformed > 0 && malformed_count == 0 {
                        warn_about_malformed_annotations(chunk_malformed, key);
                    }
                    malformed_count += chunk_malformed;
                }
                if transcript_handling == TranscriptHandling::FirstOnly {
                    let chunk_multi_transcript = summary::count_multi_transcript_sites(&chunk);
                    if chunk_multi_transcript > 0 && multi_transcript_count == 0 {
                        note_multi_transcript(chunk_multi_transcript);
                    }
                    multi_transcript_count += chunk_multi_transcript;
                }

                let converted = if use_parallel {
                    reformat_vcf_data_with_header_parallel(
                        &header,
                        &columns_title,
                        &chunk,
                        transcript_handling,
                    )
                } else {
                    reformat_vcf_data_with_header(
                        &header,
                        &columns_title,
                        &chunk,
                        transcript_handling,
                    )
                };
                let (chunk_headers, records) = match converted {
                    Ok(result) => result,
                    Err(e) => {
                        eprintln!("❌ Error reformatting VCF data: {e}");
                        std::process::exit(1);
                    }
                };

                if headers.is_empty() && !records.is_empty() {
                    headers = chunk_headers;
                    if let Err(e) = reformat_vcf::write_tsv_header(&mut writer, &headers) {
                        eprintln!("❌ Error writing file: {e}");
                        std::process::exit(1);
                    }
                }
                if let Err(e) = reformat_vcf::write_tsv_rows(&mut writer, &headers, &records) {
                    eprintln!("❌ Error writing file: {e}");
                    std::process::exit(1);
                }

                output_record_count += records.len();
                for record in &records {
                    *output_chrom_counts
                        .entry(record.chromosome.clone())
                        .or_insert(0) += 1;
                }
                if cli.report == ReportFormatCli::Html {
                    summary::merge_damage_breakdowns(
                        &mut damage_breakdowns,
                        summary::compute_damage_breakdowns(&records),
                    );
                }
                #[cfg(feature = "parquet_out")]
                if cli.parquet && !headers.is_empty() {
                    let sink = match parquet_sink.as_mut() {
                        Some(sink) => sink,
                        None => {
                            match parquet_writer::ParquetSink::create_tsv(&parquet_file, &headers) {
                                Ok(sink) => parquet_sink.insert(sink),
                                Err(e) => {
                                    eprintln!("❌ Error writing parquet file: {e}");
                                    std::process::exit(1);
                                }
                            }
                        }
                    };
                    if let Err(e) = sink.write_tsv(&records) {
                        eprintln!("❌ Error writing parquet file: {e}");
                        std::process::exit(1);
                    }
                }

                if input_variants / 200_000 > progress_marks {
                    progress_marks = input_variants / 200_000;
                    println!("   … {input_variants} variants processed");
                }
            }

            if let Err(e) = writer.flush() {
                eprintln!("❌ Error writing file: {e}");
                std::process::exit(1);
            }
            drop(writer);

            let process_time = process_start.elapsed();
            let variants_per_sec = input_variants as f64 / process_time.as_secs_f64();
            println!("✅ Data processing completed in {process_time:.2?}");
            println!("   📊 Total variants: {input_variants}");
            println!("   📈 Processing rate: {variants_per_sec:.0} variants/sec");
            println!(
                "   📊 Output records: {} ({:.2}x expansion)",
                output_record_count,
                output_record_count as f64 / input_variants.max(1) as f64
            );
            if use_parallel {
                println!(
                    "   🚀 Parallel efficiency: {:.1}x speedup potential",
                    thread_count as f64
                );
            }
            println!();

            // The header file records the variant count, so it is written once that count is
            // known — which, streaming, is after the data rather than before it.
            println!("📝 Writing header file...");
            let header_write_start = Instant::now();
            if let Err(e) = write_header_file(&header_file, &header, input_variants) {
                eprintln!("❌ Error writing header file: {e}");
                std::process::exit(1);
            }
            let header_write_time = header_write_start.elapsed();
            println!("✅ Header file written in {header_write_time:.2?}");
            println!();

            // The text file's own name plus .parquet, so both output formats read
            // X_reformatted.<tsv|maf>.parquet. --parquet and --compress are mutually
            // exclusive (checked at startup), so there is no .gz case.
            #[cfg(feature = "parquet_out")]
            if let Some(sink) = parquet_sink.take() {
                match sink.close() {
                    Ok(()) => println!("✅ Parquet file written: {}", parquet_file),
                    Err(e) => {
                        eprintln!("❌ Error writing parquet file: {e}");
                        std::process::exit(1);
                    }
                }
            }

            write_summary_from_counts(
                &cli,
                input_variants,
                input_chrom_counts,
                output_record_count,
                output_chrom_counts,
                process_time.as_secs_f64(),
                variants_per_sec,
                damage_breakdowns,
            );

            let total_time = total_start.elapsed();
            println!(
                "✅ Output file ready: {}{}",
                reformatted_file,
                if cli.compress { " (compressed)" } else { "" }
            );
            println!();

            let timing = ProcessingTiming {
                read_time,
                process_time,
                write_time: None,
                header_write_time: Some(header_write_time),
                total_time,
            };
            let stats = ProcessingStats {
                variants_per_sec,
                thread_count,
                use_parallel,
            };

            print_final_summary(
                &cli,
                &header_file,
                &reformatted_file,
                input_variants,
                output_record_count,
                &timing,
                &stats,
            );
        }
        OutputFormatCli::Maf => {
            println!("🔄 Processing VCF data for MAF format...");
            let maf_output_file = generate_maf_output_filename(&cli);
            let process_start = Instant::now();

            let mut writer = match reformat_vcf::MafWriter::create(&maf_output_file, cli.compress) {
                Ok(writer) => writer,
                Err(e) => {
                    eprintln!("❌ Error writing MAF file: {e}");
                    std::process::exit(1);
                }
            };

            let mut input_variants = 0usize;
            let mut maf_record_count = 0usize;
            let mut input_chrom_counts: indexmap::IndexMap<String, usize> =
                indexmap::IndexMap::new();
            let mut output_chrom_counts: indexmap::IndexMap<String, usize> =
                indexmap::IndexMap::new();
            let mut damage_breakdowns: Vec<summary::DamageBreakdown> = Vec::new();
            let mut multiallelic_count = 0usize;
            let mut multi_transcript_count = 0usize;
            let mut malformed_count = 0usize;
            let ann_format = annotation_format(&header);
            let mut progress_marks = 0usize;
            #[cfg(feature = "parquet_out")]
            let parquet_file = format!("{maf_output_file}.parquet");
            #[cfg(feature = "parquet_out")]
            let mut parquet_sink = if cli.parquet {
                match parquet_writer::ParquetSink::create_maf(&parquet_file) {
                    Ok(sink) => Some(sink),
                    Err(e) => {
                        eprintln!("❌ Error writing MAF parquet file: {e}");
                        std::process::exit(1);
                    }
                }
            } else {
                None
            };

            loop {
                let chunk = match next_chunk(lines.as_mut(), STREAM_CHUNK) {
                    Ok(chunk) => chunk,
                    Err(e) => {
                        eprintln!("❌ Error reading VCF file: {e}");
                        std::process::exit(1);
                    }
                };
                if chunk.is_empty() {
                    break;
                }
                input_variants += chunk.len();
                for line in &chunk {
                    if let Some(chrom) = line.split('\t').next() {
                        *input_chrom_counts.entry(chrom.to_string()).or_insert(0) += 1;
                    }
                }

                // Warn the first time one is seen rather than after a whole extra pass over
                // the file; the totals go out with the final summary.
                let chunk_multiallelic = summary::count_multiallelic_sites(&chunk);
                if chunk_multiallelic > 0 && multiallelic_count == 0 {
                    warn_about_multiallelic(chunk_multiallelic);
                }
                multiallelic_count += chunk_multiallelic;
                if let Some((key, n_fields)) = ann_format {
                    let chunk_malformed =
                        summary::count_malformed_annotation_entries(&chunk, key, n_fields);
                    if chunk_malformed > 0 && malformed_count == 0 {
                        warn_about_malformed_annotations(chunk_malformed, key);
                    }
                    malformed_count += chunk_malformed;
                }
                if transcript_handling == TranscriptHandling::FirstOnly {
                    let chunk_multi_transcript = summary::count_multi_transcript_sites(&chunk);
                    if chunk_multi_transcript > 0 && multi_transcript_count == 0 {
                        note_multi_transcript(chunk_multi_transcript);
                    }
                    multi_transcript_count += chunk_multi_transcript;
                }

                let params = MafConversionParams {
                    header: &header,
                    columns_title: &columns_title,
                    data_lines: &chunk,
                    center: &maf_config.center,
                    ncbi_build: &maf_config.ncbi_build,
                    sample_barcode: &maf_config.sample_barcode,
                    tumor_id: cli.tumor_id.as_deref(),
                    normal_id: cli.normal_id.as_deref(),
                    transcript_handling,
                    verbose: cli.verbose && input_variants <= STREAM_CHUNK,
                    use_parallel,
                    compute_damage: cli.report == ReportFormatCli::Html,
                };
                let (mut records, breakdowns) = match convert_to_maf_records(&params) {
                    Ok(result) => result,
                    Err(e) => {
                        eprintln!("❌ Error converting to MAF: {e}");
                        std::process::exit(1);
                    }
                };
                summary::merge_damage_breakdowns(&mut damage_breakdowns, breakdowns);

                // Study metadata no VCF carries; "." unless the user states it.
                for record in &mut records {
                    if let Some(status) = &cli.mutation_status {
                        record.mutation_status = status.clone();
                    }
                    if let Some(source) = &cli.sequence_source {
                        record.sequence_source = source.clone();
                    }
                    *output_chrom_counts
                        .entry(record.chromosome.clone())
                        .or_insert(0) += 1;
                }
                maf_record_count += records.len();

                if let Err(e) = writer.write_rows(&records) {
                    eprintln!("❌ Error writing MAF file: {e}");
                    std::process::exit(1);
                }
                #[cfg(feature = "parquet_out")]
                if let Some(sink) = parquet_sink.as_mut() {
                    if let Err(e) = sink.write_maf(&records) {
                        eprintln!("❌ Error writing MAF parquet file: {e}");
                        std::process::exit(1);
                    }
                }

                // With streaming the total is unknowable, so report what has been done.
                if input_variants / 200_000 > progress_marks {
                    progress_marks = input_variants / 200_000;
                    println!("   … {input_variants} variants processed");
                }
            }

            if let Err(e) = writer.finish() {
                eprintln!("❌ Error writing MAF file: {e}");
                std::process::exit(1);
            }

            let process_time = process_start.elapsed();
            let variants_per_sec = input_variants as f64 / process_time.as_secs_f64();
            println!("✅ MAF processing completed in {process_time:.2?}");
            println!("   📊 Total variants: {input_variants}");
            println!("   📈 Processing rate: {variants_per_sec:.0} variants/sec");
            println!(
                "   📊 MAF records: {} ({:.2}x expansion)",
                maf_record_count,
                maf_record_count as f64 / input_variants.max(1) as f64
            );
            println!();

            #[cfg(feature = "parquet_out")]
            if let Some(sink) = parquet_sink.take() {
                match sink.close() {
                    Ok(()) => println!("✅ MAF Parquet file written: {}", parquet_file),
                    Err(e) => {
                        eprintln!("❌ Error writing MAF parquet file: {e}");
                        std::process::exit(1);
                    }
                }
            }

            write_summary_from_counts(
                &cli,
                input_variants,
                input_chrom_counts,
                maf_record_count,
                output_chrom_counts,
                process_time.as_secs_f64(),
                variants_per_sec,
                damage_breakdowns,
            );

            let total_time = total_start.elapsed();
            println!(
                "✅ MAF file written{}",
                if cli.compress { " (compressed)" } else { "" }
            );
            println!();

            let timing = ProcessingTiming {
                read_time,
                process_time,
                write_time: None,
                header_write_time: None,
                total_time,
            };
            let stats = ProcessingStats {
                variants_per_sec,
                thread_count,
                use_parallel,
            };

            print_maf_summary(&cli, &maf_output_file, maf_record_count, &timing, &stats);
        }
    }
}

/// Build an output path under `cli.output_dir`, using `cli.prefix` (or the
/// input file's stem) as the base name, creating the directory if needed.
fn output_path(cli: &Cli, suffix: &str, extension: &str) -> String {
    let base_name = if cli.input_file == "-" {
        "stdin"
    } else {
        let input_path = Path::new(&cli.input_file);
        let base_name = input_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output");
        base_name.strip_suffix(".vcf").unwrap_or(base_name)
    };

    let output_dir = cli.output_dir.as_deref().unwrap_or(".");
    let prefix = cli.prefix.as_deref().unwrap_or(base_name);

    if let Err(e) = std::fs::create_dir_all(output_dir) {
        eprintln!(
            "⚠️  Warning: Could not create output directory '{}': {}",
            output_dir, e
        );
    }

    format!("{}/{}{}.{}", output_dir, prefix, suffix, extension)
}

fn generate_maf_output_filename(cli: &Cli) -> String {
    let extension = if cli.compress { "maf.gz" } else { "maf" };
    output_path(cli, "_reformatted", extension)
}

fn print_maf_summary(
    cli: &Cli,
    maf_file: &str,
    maf_record_count: usize,
    timing: &ProcessingTiming,
    stats: &ProcessingStats,
) {
    println!("🎉 MAF CONVERSION COMPLETED SUCCESSFULLY!");
    println!("════════════════════════════════════════");
    println!("📁 Input file:      {}", cli.input_file);
    println!("📁 MAF output:      {}", maf_file);
    println!();
    println!("📊 PROCESSING STATISTICS:");
    println!("─────────────────────────");
    println!("📖 File reading:    {:.2?}", timing.read_time);
    println!("🔄 Data processing:  {:.2?}", timing.process_time);
    if let Some(write_time) = timing.write_time {
        println!("💾 File writing:     {:.2?}", write_time);
    }
    println!("⏱️  Total time:       {:.2?}", timing.total_time);
    println!();
    println!("📈 PERFORMANCE METRICS:");
    println!("─────────────────────────");
    println!(
        "🚀 Processing rate:  {:.0} variants/sec",
        stats.variants_per_sec
    );
    println!("📊 MAF records:      {maf_record_count}");
    if stats.use_parallel {
        println!("🧵 Threads used:     {}", stats.thread_count);
    }
    println!();
}

fn print_startup_info(
    cli: &Cli,
    thread_count: usize,
    use_parallel: bool,
    transcript_handling: TranscriptHandling,
    _annotation_type: AnnotationType,
) {
    // Welcome messages
    println!("🧬 VCF REFORMATTER v0.4.0");
    println!("═══════════════════════════");
    println!("📁 Input file: {}", cli.input_file);
    println!("🧵 Transcript handling: {:?}", transcript_handling);
    println!("⚡ Threads: {}", thread_count);
    if cli.output_format == OutputFormatCli::Maf {
        println!("📋 Output format: MAF");
    } else {
        println!("📋 Output format: TSV");
    }
    if use_parallel {
        println!("🚀 Parallel processing: enabled");
    }
    println!();
}

fn detect_and_print_annotation_type(
    header: &str,
    annotation_type: AnnotationTypeCli,
) -> Result<(), String> {
    let has_csq = header.contains("##INFO=<ID=CSQ");
    let has_ann = header.contains("##INFO=<ID=ANN");

    match (has_csq, has_ann) {
        (true, true) => {
            println!("🔍 Detected: Both VEP (CSQ) and SnpEff (ANN) annotations");

            // Check if user specified auto - this is an error
            if annotation_type == AnnotationTypeCli::Auto {
                return Err(format!(
                    "Ambiguous annotation types detected!\n\
                    🔍 Found: Both VEP (CSQ) and SnpEff (ANN) annotations in the header.\n\
                    \n\
                    Please specify which annotation type to use:\n\
                      --annotation-type vep     (to use VEP/CSQ annotations)\n\
                      --annotation-type snpeff  (to use SnpEff/ANN annotations)\n\
                    \n\
                    Example: {} --annotation-type vep [other options...]",
                    std::env::args()
                        .next()
                        .unwrap_or_else(|| "vcf-reformatter".to_string())
                ));
            }
        }
        (true, false) => println!("🔍 Detected: VEP (CSQ) annotations"),
        (false, true) => println!("🔍 Detected: SnpEff (ANN) annotations"),
        (false, false) => println!("⚠️  No standard annotations detected"),
    }

    Ok(())
}

fn generate_output_filenames(cli: &Cli) -> (String, String) {
    let extension = if cli.compress { "tsv.gz" } else { "tsv" };
    (
        output_path(cli, "_header", "txt"),
        output_path(cli, "_reformatted", extension),
    )
}

fn generate_summary_filename(cli: &Cli, extension: &str) -> String {
    output_path(cli, "_summary", extension)
}

/// The same report, built from counters instead of the retained input lines — which is what a
/// streaming path has. `write_summary_if_requested` is the collect-everything caller.
#[allow(clippy::too_many_arguments)]
fn write_summary_from_counts(
    cli: &Cli,
    input_variant_count: usize,
    input_chrom_counts: indexmap::IndexMap<String, usize>,
    output_record_count: usize,
    output_chrom_counts: indexmap::IndexMap<String, usize>,
    process_time_secs: f64,
    variants_per_sec: f64,
    damage_breakdowns: Vec<summary::DamageBreakdown>,
) {
    if cli.report == ReportFormatCli::None {
        return;
    }
    let input_chrom_counts = summary::sort_chromosomes(input_chrom_counts);
    let output_chrom_counts = summary::sort_chromosomes(output_chrom_counts);

    let output_format = match cli.output_format {
        OutputFormatCli::Tsv => "TSV",
        OutputFormatCli::Maf => "MAF",
    };
    let transcript_mode = match cli.transcript_handling {
        TranscriptHandlingCli::FirstOnly => "first",
        TranscriptHandlingCli::MostSevere => "most-severe",
        TranscriptHandlingCli::SplitRows => "split",
    };

    let summary_stats = summary::SummaryStats {
        input_file: cli.input_file.clone(),
        output_format: output_format.to_string(),
        transcript_handling: transcript_mode.to_string(),
        input_variant_count,
        output_record_count,
        input_chrom_counts,
        output_chrom_counts,
        processing_time_secs: process_time_secs,
        variants_per_sec,
    };

    match cli.report {
        ReportFormatCli::None => {}
        ReportFormatCli::Txt => {
            let summary_file = generate_summary_filename(cli, "txt");
            if let Err(e) = summary_stats.write_to_file(&summary_file) {
                eprintln!("Warning: Could not write summary file: {}", e);
            } else {
                println!("Summary written to: {}", summary_file);
            }
        }
        ReportFormatCli::Html => {
            let summary_file = generate_summary_filename(cli, "html");
            if let Err(e) =
                html_report::write_html_report(&summary_stats, &damage_breakdowns, &summary_file)
            {
                eprintln!("Warning: Could not write summary report: {}", e);
            } else {
                println!("Summary report written to: {}", summary_file);
            }
        }
    }
}

fn write_header_file(filename: &str, header: &str, variant_count: usize) -> std::io::Result<()> {
    let mut file = File::create(filename)?;
    writeln!(file, "VCF Header Information")?;
    writeln!(file, "======================")?;
    writeln!(file, "Total variants: {}", variant_count)?;
    writeln!(file)?;
    writeln!(file, "Header content:")?;
    writeln!(file, "{}", header)?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn print_final_summary(
    cli: &Cli,
    header_file: &str,
    reformatted_file: &str,
    input_variant_count: usize,
    output_record_count: usize,
    timing: &ProcessingTiming,
    stats: &ProcessingStats,
) {
    println!("🎉 PROCESSING COMPLETED SUCCESSFULLY!");
    println!("════════════════════════════════════════");
    println!("📁 Input file:       {}", cli.input_file);
    println!("📁 Header file:      {}", header_file);
    println!("📁 Reformatted file: {}", reformatted_file);
    println!();
    println!("📊 STATISTICS:");
    println!("──────────────");
    println!("📖 File reading:     {:.2?}", timing.read_time);
    if let Some(header_write_time) = timing.header_write_time {
        println!("📝 Header writing:   {:.2?}", header_write_time);
    }
    println!("🔄 Data processing:  {:.2?}", timing.process_time);
    if let Some(write_time) = timing.write_time {
        println!("💾 File writing:     {:.2?}", write_time);
    }
    println!("⏱️  Total time:       {:.2?}", timing.total_time);
    println!();
    println!("📈 PERFORMANCE:");
    println!("──────────────");
    println!(
        "🚀 Processing rate:  {:.0} variants/sec",
        stats.variants_per_sec
    );
    println!("📊 Input variants:   {input_variant_count}");
    println!("📊 Output records:   {output_record_count}");
    if stats.use_parallel {
        println!("🧵 Threads used:     {}", stats.thread_count);
    }
}

#[derive(Debug)]
struct MafConversionParams<'a> {
    header: &'a str,
    columns_title: &'a str,
    data_lines: &'a [String],
    center: &'a str,
    ncbi_build: &'a str,
    sample_barcode: &'a str,
    tumor_id: Option<&'a str>,
    normal_id: Option<&'a str>,
    transcript_handling: TranscriptHandling,
    verbose: bool,
    use_parallel: bool,
    /// Only the HTML report reads the damage counts; with no report they are pure waste.
    compute_damage: bool,
}

fn convert_to_maf_records(
    params: &MafConversionParams,
) -> Result<(Vec<MafRecord>, Vec<summary::DamageBreakdown>), String> {
    // Lines are converted a chunk at a time and each chunk's ReformattedVcfRecords are dropped
    // as soon as its MAF rows exist. Holding the whole intermediate vec alongside the MAF vec
    // cost ~1GB on a 92k-variant file, purely so the HTML report could read it afterwards; the
    // report only needs the damage counts, which are folded in per chunk instead.
    const CHUNK: usize = 25_000;

    let mut maf_records = Vec::new();
    let mut breakdowns: Vec<summary::DamageBreakdown> = Vec::new();

    if params.verbose && params.use_parallel {
        println!("   🚀 Using parallel processing for MAF conversion...");
    }

    for chunk in params.data_lines.chunks(CHUNK) {
        let (_, reformatted) = if params.use_parallel {
            reformat_vcf_data_with_header_parallel(
                params.header,
                params.columns_title,
                chunk,
                params.transcript_handling,
            )
        } else {
            reformat_vcf_data_with_header(
                params.header,
                params.columns_title,
                chunk,
                params.transcript_handling,
            )
        }
        .map_err(|e| format!("Failed to process VCF data: {}", e))?;

        if params.compute_damage {
            summary::merge_damage_breakdowns(
                &mut breakdowns,
                summary::compute_damage_breakdowns(&reformatted),
            );
        }

        for record in &reformatted {
            let converted = MafRecord::from_reformatted_record_multi_for_samples(
                record,
                params.center,
                params.ncbi_build,
                params.sample_barcode,
                params.tumor_id,
                params.normal_id,
            )
            .map_err(|e| format!("Failed to convert to MAF format: {}", e))?;
            maf_records.extend(converted);
        }
    }

    if params.verbose {
        println!("   ✅ Generated {} MAF records", maf_records.len());
    }

    Ok((maf_records, breakdowns))
}
