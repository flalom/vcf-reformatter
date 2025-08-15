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
use read_vcf_gz::read_vcf_gz;
use reformat_vcf::{
    reformat_vcf_data_with_header,
    reformat_vcf_data_with_header_parallel,
    reformat_vcf_data_with_header_parallel_chunked, // Add this line
    write_maf_file,
    write_reformatted_vcf,
    AnnotationType,
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
mod read_vcf_gz;
mod reformat_vcf;

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
    write_time: std::time::Duration,
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
    version = "0.3.0",
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
    /// Input VCF file (supports .vcf.gz compressed files)
    #[arg(value_name = "INPUT_FILE")]
    input_file: String,

    /// Annotation type to parse
    #[arg(short = 'a', long = "annotation-type", value_enum, default_value_t = AnnotationTypeCli::Auto)]
    annotation_type: AnnotationTypeCli,

    /// Transcript handling mode
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
    /// Keep first transcript only (fastest)
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
    // Set better default values
    metadata.center = Some("Unknown_Center".to_string());
    metadata.ncbi_build = Some("GRCh38".to_string()); // reasonable default

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

fn main() {
    let cli = Cli::parse();

    // Validate input file
    if !Path::new(&cli.input_file).exists() {
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
    let data = match read_vcf_gz(&cli.input_file) {
        Ok(data) => data,
        Err(e) => {
            eprintln!("❌ Error reading VCF file: {e}");
            std::process::exit(1);
        }
    };

    // Extract metadata from header for MAF
    let maf_metadata = extract_maf_metadata_from_header(&data.0, &data.1);

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
    if let Err(e) = detect_and_print_annotation_type(&data.0, cli.annotation_type) {
        eprintln!("❌ {}", e);
        std::process::exit(1);
    }

    let read_time = read_start.elapsed();
    println!("✅ File read completed in {read_time:.2?}");
    println!("   📊 Total variants: {}", data.2.len());
    println!("   📑 Header lines: {}", data.0.matches('\n').count());

    if cli.verbose && data.2.len() > 50_000 {
        println!(
            "   📊 Large dataset detected - processing {} variants in optimized chunks...",
            data.2.len()
        );
    }

    println!();

    // Branch based on output format
    match cli.output_format {
        OutputFormatCli::Tsv => {
            // Generate output filenames for TSV
            let (header_file, reformatted_file) = generate_output_filenames(&cli);

            // Write header file
            println!("📝 Writing header file...");
            let header_write_start = Instant::now();
            if let Err(e) = write_header_file(&header_file, &data.0, data.2.len()) {
                eprintln!("❌ Error writing header file: {e}");
                std::process::exit(1);
            }
            let header_write_time = header_write_start.elapsed();
            println!("✅ Header file written in {header_write_time:.2?}");
            println!();

            // Process VCF data for TSV
            println!("🔄 Processing VCF data...");
            let process_start = Instant::now();

            // Create output writer for streaming
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

            let result = if use_parallel {
                if cli.verbose {
                    println!("   🚀 Using parallel processing with {thread_count} threads...");
                }

                // Use smart threshold logic that accounts for transcript expansion
                let estimated_output_records =
                    if transcript_handling == TranscriptHandling::SplitRows {
                        // With split transcripts, estimate 10x expansion for typical VEP/SnpEff data
                        data.2.len() * 10
                    } else {
                        // FirstOnly or MostSevere produce 1:1 record mapping
                        data.2.len()
                    };

                // Use chunked processing if:
                // 1. Large input (>100K variants) OR
                // 2. Large estimated output (>500K records)
                if data.2.len() > 100_000 || estimated_output_records > 500_000 {
                    if cli.verbose {
                        if estimated_output_records > data.2.len() {
                            println!("   📦 Large output expected ({} estimated records) - using chunked processing to prevent memory issues", estimated_output_records);
                        } else {
                            println!("   📦 Large dataset detected ({} variants) - using chunked processing to prevent memory issues", data.2.len());
                        }
                    }

                    // Use streaming chunked function
                    reformat_vcf_data_with_header_parallel_chunked(
                        &data.0,
                        &data.1,
                        &data.2,
                        transcript_handling,
                        &mut *writer,
                    )
                    .map(|headers| {
                        // Flush and close writer
                        writer.flush().expect("Failed to flush writer");
                        drop(writer);
                        (headers, Vec::new()) // Return empty records since we streamed them
                    })
                } else {
                    // Use regular parallel function for small files
                    let result = reformat_vcf_data_with_header_parallel(
                        &data.0,
                        &data.1,
                        &data.2,
                        transcript_handling,
                    );

                    // Write results using traditional method for small files
                    if let Ok((headers, records)) = &result {
                        // Drop the streaming writer
                        drop(writer);
                        // Use the existing write function for small files
                        if let Err(e) =
                            write_reformatted_vcf(&reformatted_file, headers, records, cli.compress)
                        {
                            eprintln!("❌ Error writing file: {e}");
                            std::process::exit(1);
                        }
                    }
                    result
                }
            } else {
                if cli.verbose {
                    println!("   🐌 Using sequential processing...");
                }
                let result =
                    reformat_vcf_data_with_header(&data.0, &data.1, &data.2, transcript_handling);

                // Write results using traditional method for sequential processing
                if let Ok((headers, records)) = &result {
                    drop(writer);
                    if let Err(e) =
                        write_reformatted_vcf(&reformatted_file, headers, records, cli.compress)
                    {
                        eprintln!("❌ Error writing file: {e}");
                        std::process::exit(1);
                    }
                }
                result
            };

            match result {
                Ok((_headers, records)) => {
                    let process_time = process_start.elapsed();
                    let variants_per_sec = data.2.len() as f64 / process_time.as_secs_f64();
                    println!("✅ Data processing completed in {process_time:.2?}");
                    println!("   📈 Processing rate: {variants_per_sec:.0} variants/sec");

                    // For large files (streamed), records will be empty
                    if records.is_empty() && data.2.len() > 100_000 {
                        println!("   📊 Output: Streamed directly to file (memory-efficient)");
                    } else {
                        println!(
                            "   📊 Output records: {} ({:.2}x expansion)",
                            records.len(),
                            records.len() as f64 / data.2.len() as f64
                        );
                    }

                    if use_parallel {
                        println!(
                            "   🚀 Parallel efficiency: {:.1}x speedup potential",
                            thread_count as f64
                        );
                    }
                    println!();

                    // File writing is already handled above for streaming case
                    let write_time = std::time::Duration::from_millis(10); // Nominal time for streamed files
                    let total_time = total_start.elapsed();
                    println!(
                        "✅ Output file ready: {}{}",
                        reformatted_file,
                        if cli.compress { " (compressed)" } else { "" }
                    );
                    println!();

                    // Print final summary - create the structs with correct data
                    let timing = ProcessingTiming {
                        read_time,
                        process_time,
                        write_time,
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
                        &data,
                        &records,
                        &timing,
                        header_write_time,
                        &stats,
                    );
                }
                Err(e) => {
                    eprintln!("❌ Error reformatting VCF data: {e}");
                    std::process::exit(1);
                }
            }
        }
        OutputFormatCli::Maf => {
            // Process VCF data for MAF
            println!("🔄 Processing VCF data for MAF format...");

            if data.2.len() > 1_000_000 {
                println!("⚠️  WARNING: Large file with MAF output detected!");
                println!("   📁 File size: {} variants", data.2.len());
                println!(
                    "   💾 Estimated memory: ~{} GB",
                    (data.2.len() as f64 * 6.0 / 1_000_000.0) as u32
                );
                println!("   💡 Consider using TSV format for better memory efficiency:");
                println!(
                    "      vcf-reformatter {} -j {} -c -v",
                    cli.input_file, cli.threads
                );
                println!("   🤔 Continue with MAF anyway? This may cause out-of-memory errors.");
                println!();
            }

            let process_start = Instant::now();

            // Create MafConversionParams struct
            let params = MafConversionParams {
                header: &data.0,
                columns_title: &data.1,
                data_lines: &data.2,
                center: &maf_config.center,
                ncbi_build: &maf_config.ncbi_build,
                sample_barcode: &maf_config.sample_barcode,
                transcript_handling,
                verbose: cli.verbose,
                use_parallel,
            };

            // Convert to MAF records using the fixed function
            let maf_records = match convert_to_maf_records(&params) {
                Ok(records) => records,
                Err(e) => {
                    eprintln!("❌ Error converting to MAF: {e}");
                    std::process::exit(1);
                }
            };

            let process_time = process_start.elapsed();
            let variants_per_sec = data.2.len() as f64 / process_time.as_secs_f64();
            println!("✅ MAF processing completed in {process_time:.2?}");
            println!("   📈 Processing rate: {variants_per_sec:.0} variants/sec");
            println!(
                "   📊 MAF records: {} ({:.2}x expansion)",
                maf_records.len(),
                maf_records.len() as f64 / data.2.len() as f64
            );
            println!();

            // Generate MAF output filename
            let maf_output_file = generate_maf_output_filename(&cli);

            // Write MAF file
            println!(
                "💾 Writing MAF file{}...",
                if cli.compress { " (compressed)" } else { "" }
            );
            let write_start = Instant::now();
            if let Err(e) = write_maf_file(&maf_output_file, &maf_records, cli.compress) {
                eprintln!("❌ Error writing MAF file: {e}");
                std::process::exit(1);
            }
            let write_time = write_start.elapsed();
            let total_time = total_start.elapsed();
            println!(
                "✅ MAF file written in {:.2?}{}",
                write_time,
                if cli.compress { " (compressed)" } else { "" }
            );
            println!();

            // Print MAF summary - create the structs with correct data
            let timing = ProcessingTiming {
                read_time,
                process_time,
                write_time,
                total_time,
            };

            let stats = ProcessingStats {
                variants_per_sec,
                thread_count,
                use_parallel,
            };

            print_maf_summary(&cli, &maf_output_file, &maf_records, &timing, &stats);
        }
    }
}

fn generate_maf_output_filename(cli: &Cli) -> String {
    let input_path = Path::new(&cli.input_file);
    let base_name = input_path
        .file_stem()
        .and_then(|s| s.to_str())
        .ok_or("Invalid input filename")
        .unwrap();

    let base_name = if let Some(_stripped) = base_name.strip_suffix(".vcf") {
        &base_name[..base_name.len() - 4]
    } else {
        base_name
    };

    let output_dir = cli.output_dir.as_deref().unwrap_or(".");
    let prefix = cli.prefix.as_deref().unwrap_or(base_name);

    // Create output directory if it doesn't exist
    if let Err(e) = std::fs::create_dir_all(output_dir) {
        eprintln!(
            "⚠️  Warning: Could not create output directory '{}': {}",
            output_dir, e
        );
    }

    let extension = if cli.compress { "maf.gz" } else { "maf" };
    format!("{}/{}_reformatted.{}", output_dir, prefix, extension)
}

fn print_maf_summary(
    cli: &Cli,
    maf_file: &str,
    maf_records: &[MafRecord],
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
    println!("💾 File writing:     {:.2?}", timing.write_time);
    println!("⏱️  Total time:       {:.2?}", timing.total_time);
    println!();
    println!("📈 PERFORMANCE METRICS:");
    println!("─────────────────────────");
    println!(
        "🚀 Processing rate:  {:.0} variants/sec",
        stats.variants_per_sec
    );
    println!("📊 MAF records:      {}", maf_records.len());
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
    println!("🧬 VCF REFORMATTER v0.3.0");
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
    let input_path = Path::new(&cli.input_file);
    let base_name = input_path.file_stem().unwrap().to_str().unwrap();
    let base_name = if let Some(_stripped) = base_name.strip_suffix(".vcf") {
        &base_name[..base_name.len() - 4]
    } else {
        base_name
    };

    let output_dir = cli.output_dir.as_deref().unwrap_or(".");
    let prefix = cli.prefix.as_deref().unwrap_or(base_name);

    // Create output directory if it doesn't exist
    if let Err(e) = std::fs::create_dir_all(output_dir) {
        eprintln!(
            "⚠️  Warning: Could not create output directory '{}': {}",
            output_dir, e
        );
    }

    let header_file = format!("{}/{}_header.txt", output_dir, prefix);
    let extension = if cli.compress { "tsv.gz" } else { "tsv" };
    let reformatted_file = format!("{}/{}_reformatted.{}", output_dir, prefix, extension);

    (header_file, reformatted_file)
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
    data: &(String, String, Vec<String>),
    records: &[reformat_vcf::ReformattedVcfRecord],
    timing: &ProcessingTiming,
    header_write_time: std::time::Duration,
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
    println!("📝 Header writing:   {:.2?}", header_write_time);
    println!("🔄 Data processing:  {:.2?}", timing.process_time);
    println!("💾 File writing:     {:.2?}", timing.write_time);
    println!("⏱️  Total time:       {:.2?}", timing.total_time);
    println!();
    println!("📈 PERFORMANCE:");
    println!("──────────────");
    println!(
        "🚀 Processing rate:  {:.0} variants/sec",
        stats.variants_per_sec
    );
    println!("📊 Input variants:   {}", data.2.len());
    println!("📊 Output records:   {}", records.len());
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
    transcript_handling: TranscriptHandling,
    verbose: bool,
    use_parallel: bool,
}

fn convert_to_maf_records(params: &MafConversionParams) -> Result<Vec<MafRecord>, String> {
    // For MAF conversion, we need the actual records, so we can't use the streaming function
    // We'll use the regular functions that return (headers, records) tuples
    let reformatted_records = if params.use_parallel {
        if params.verbose {
            println!("   🚀 Using parallel processing for MAF conversion...");
            if params.data_lines.len() > 100_000 {
                println!("   📦 Large dataset detected ({} variants) - using parallel processing for MAF conversion", params.data_lines.len());
            }
        }

        // For MAF, we always use the regular parallel function (not streaming)
        // because we need access to all records to convert them to MAF
        reformat_vcf_data_with_header_parallel(
            params.header,
            params.columns_title,
            params.data_lines,
            params.transcript_handling,
        )
    } else {
        reformat_vcf_data_with_header(
            params.header,
            params.columns_title,
            params.data_lines,
            params.transcript_handling,
        )
    }
        .map_err(|e| format!("Failed to process VCF data: {}", e))?;

    if params.verbose {
        println!(
            "   ✅ Processed {} VCF records",
            reformatted_records.1.len()
        );
        println!("🔄 Converting to MAF format...");
    }

    // Convert to MAF records
    let maf_records: Result<Vec<MafRecord>, _> = reformatted_records
        .1 // Access the records part of the tuple (headers, records)
        .iter()
        .flat_map(|record| {
            match MafRecord::from_reformatted_record_multi(
                record,
                params.center,
                params.ncbi_build,
                params.sample_barcode,
            ) {
                Ok(records) => records.into_iter().map(Ok).collect::<Vec<_>>(),
                Err(e) => vec![Err(e)],
            }
        })
        .collect();

    let maf_records = maf_records.map_err(|e| format!("Failed to convert to MAF format: {}", e))?;

    if params.verbose {
        println!("   ✅ Generated {} MAF records", maf_records.len());
    }

    Ok(maf_records)
}
