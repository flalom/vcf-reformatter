#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]

mod essentials_fields;
mod extract_csq_and_csq_names;
mod extract_sample_info;
mod get_info_from_header;
mod read_vcf_gz;
mod reformat_vcf;

use clap::{Parser, ValueEnum};
use std::fs::{create_dir_all, File};
use std::io::{Result, Write};
use std::path::Path;
use std::time::Instant;

use crate::essentials_fields::*;
use crate::essentials_fields::*;
use crate::extract_csq_and_csq_names::extract_csq_regex;
use crate::extract_sample_info::*;
use crate::extract_sample_info::*;
use crate::read_vcf_gz::read_vcf_gz;
use crate::reformat_vcf::{
    reformat_vcf_data_with_header, reformat_vcf_data_with_header_parallel, write_reformatted_vcf,
    TranscriptHandling,
};

#[derive(Parser)]
#[command(
    name = "vcf-reformatter",
    version = "0.1.0",
    about = "🧬 Fast VCF file parser and reformatter with VEP annotation support",
    long_about = "A Rust command-line tool for parsing and reformatting VCF (Variant Call Format) files, with special support for VEP (Variant Effect Predictor) annotations. This tool flattens complex VCF files into tab-separated values (TSV) format for easier downstream analysis.",
    after_help = "EXAMPLES:
    Basic usage:
      vcf-reformatter sample.vcf.gz

    Use most severe consequence with parallel processing:
      vcf-reformatter sample.vcf.gz -t most-severe -j 4

    Split all transcripts with custom output:
      vcf-reformatter sample.vcf.gz -t split -o results/ -p analysis

    Auto-detect threads with verbose output:
      vcf-reformatter sample.vcf.gz -j 0 -v

    Compress output with gzip:
      vcf-reformatter sample.vcf.gz -t most-severe -j 4 -o results/ -p pippo -v --compress

    Complete example with all options:
      vcf-reformatter sample.vcf.gz -t most-severe -j 4 -o results/ -p my_analysis -v --compress

    Container usage:
      docker run --rm -v $(pwd):/data vcf-reformatter /data/sample.vcf.gz -t split --compress

    For more examples and documentation, visit:
    https://github.com/flalom/vcf-reformatter"
)]
struct Cli {
    /// Input VCF file (supports .vcf.gz compressed files)
    #[arg(value_name = "INPUT_FILE")]
    input_file: String,

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

impl From<TranscriptHandlingCli> for TranscriptHandling {
    fn from(cli: TranscriptHandlingCli) -> Self {
        match cli {
            TranscriptHandlingCli::MostSevere => TranscriptHandling::MostSevere,
            TranscriptHandlingCli::FirstOnly => TranscriptHandling::FirstOnly,
            TranscriptHandlingCli::SplitRows => TranscriptHandling::SplitRows,
        }
    }
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
    let use_parallel = thread_count > 1;

    // Start total timing
    let total_start = Instant::now();

    // Print startup information
    print_startup_info(&cli, thread_count, use_parallel, transcript_handling);

    // Set thread pool size if using parallel processing
    if use_parallel {
        rayon::ThreadPoolBuilder::new()
            .num_threads(thread_count)
            .build_global()
            .unwrap();
    }

    // Read VCF file
    println!("📖 Reading VCF file...");
    let read_start = Instant::now();
    let data = match read_vcf_gz(&cli.input_file) {
        Ok(data) => data,
        Err(e) => {
            eprintln!("❌ Error reading VCF file: {e}");
            std::process::exit(1);
        }
    };
    let read_time = read_start.elapsed();

    println!("✅ File read completed in {read_time:.2?}");
    println!("   📊 Total variants: {}", data.2.len());
    println!("   📑 Header lines: {}", data.0.matches('\n').count());
    println!();

    // Generate output filenames
    let (header_file, reformatted_file) = generate_output_filenames(&cli);

    // Write heade file
    println!("📝 Writing header file...");
    let header_write_start = Instant::now();
    if let Err(e) = write_header_file(&header_file, &data.0, data.2.len()) {
        eprintln!("❌ Error writing header file: {e}");
        std::process::exit(1);
    }
    let header_write_time = header_write_start.elapsed();
    println!("✅ Header file written in {header_write_time:.2?}");
    println!();

    // Process VCF data
    println!("🔄 Processing VCF data...");
    let process_start = Instant::now();

    let result = if use_parallel {
        if cli.verbose {
            println!("   🚀 Using parallel processing with {thread_count} threads...");
        }
        reformat_vcf_data_with_header_parallel(&data.0, &data.1, &data.2, transcript_handling)
    } else {
        if cli.verbose {
            println!("   🐌 Using sequential processing...");
        }
        reformat_vcf_data_with_header(&data.0, &data.1, &data.2, transcript_handling)
    };

    match result {
        Ok((headers, records)) => {
            let process_time = process_start.elapsed();
            let variants_per_sec = data.2.len() as f64 / process_time.as_secs_f64();

            println!("✅ Data processing completed in {process_time:.2?}");
            println!("   📈 Processing rate: {variants_per_sec:.0} variants/sec");
            println!(
                "   📊 Output records: {} ({:.2}x expansion)",
                records.len(),
                records.len() as f64 / data.2.len() as f64
            );
            if use_parallel {
                println!(
                    "   🚀 Parallel efficiency: {:.1}x speedup potential",
                    thread_count as f64
                );
            }
            println!();

            // Write reformatted file
            println!(
                "💾 Writing reformatted file{}...",
                if cli.compress { " (compressed)" } else { "" }
            );
            let write_start = Instant::now();
            if let Err(e) =
                write_reformatted_vcf(&reformatted_file, &headers, &records, cli.compress)
            {
                eprintln!("❌ Error writing reformatted file: {e}");
                std::process::exit(1);
            }
            let write_time = write_start.elapsed();
            let total_time = total_start.elapsed();

            println!(
                "✅ Reformatted file written in {:.2?}{}",
                write_time,
                if cli.compress { " (compressed)" } else { "" }
            );
            println!();

            // Print final summary
            print_final_summary(
                &cli,
                &header_file,
                &reformatted_file,
                &data,
                &records,
                read_time,
                process_time,
                write_time,
                header_write_time,
                total_time,
                variants_per_sec,
                thread_count,
                use_parallel,
            );
        }
        Err(e) => {
            eprintln!("❌ Error reformatting VCF data: {e}");
            std::process::exit(1);
        }
    }

    /*
        //println!("{:?}", data.1);
        let first_char = data.1.chars().next().unwrap_or('?');
        //println!("First character: {}", first_char);

        let column_names: Vec<&str> = data.1.split('\t').collect();
        //println!("Column names: {:?}", column_names);

        // Parse first variant (correct syntax ✅)
        let vcf_variant = VcfVariant::from_line(&data.2[0], &column_names).unwrap();
        //println!("First variant: {:?}", vcf_variant);
        //println!("Reference: {:?}", vcf_variant.reference);
        //println!("Alternative: {:?}", vcf_variant.alternate);

        // Split ALL data lines by tabs
        let mut all_parsed_lines: Vec<Vec<&str>> = data.2
            .iter()
            .map(|line| line.split('\t').collect())
            .collect();

        //println!("Parsed lines: {:?}", all_parsed_lines[1]);
        //let mut csq = extract_and_remove_csq(extract_and_remove_csq[1..2]);
        println!("Parsed lines: {:?}", all_parsed_lines[1]);
        // Convert Vec<Vec<&str>> to Vec<String> by joining each line's parts with tabs
        let mut string_lines: Vec<String> = all_parsed_lines
            .into_iter()
            .map(|line_parts| line_parts.join("\t"))
            .collect();

        println!("Parsed lines 2: {:?}", extract_csq_regex(&mut string_lines).unwrap().split('|'));



        //println!("Header lines count: {}", data.0.len());
        //println!("Example header {}", data.0);  // First element (header)
        //println!("Example data (first few): {:?}", data.1.iter().take(3).collect::<Vec<_>>());

        /*
        for (i, line) in data.1.iter().take(5).enumerate() {
            println!("  {}: {}", i + 1, line);
        }
        */
        //println!("Total data lines: {}", data.1.len());
    */
}

fn print_startup_info(
    cli: &Cli,
    thread_count: usize,
    use_parallel: bool,
    transcript_handling: TranscriptHandling,
) {
    println!("🧬 VCF Reformatter Starting...");
    println!("📁 Input file: {}", cli.input_file);
    println!(
        "🔧 Transcript handling: {}",
        match transcript_handling {
            TranscriptHandling::MostSevere => "Most severe consequence only",
            TranscriptHandling::FirstOnly => "First transcript only",
            TranscriptHandling::SplitRows => "All transcripts in separate rows",
        }
    );
    println!(
        "⚡ Processing mode: {} (using {} thread{})",
        if use_parallel {
            "Parallel"
        } else {
            "Sequential"
        },
        thread_count,
        if thread_count == 1 { "" } else { "s" }
    );

    if cli.compress {
        println!("🗜️  Output compression: Enabled (gzip)");
    }

    if let Some(ref output_dir) = cli.output_dir {
        println!("📂 Output directory: {output_dir}");
    }

    if let Some(ref prefix) = cli.prefix {
        println!("🏷️  Output prefix: {prefix}");
    }

    println!();
}

fn generate_output_filenames(cli: &Cli) -> (String, String) {
    let base_name = if let Some(ref prefix) = cli.prefix {
        prefix.clone()
    } else {
        get_base_filename(&cli.input_file)
    };

    let output_dir = cli.output_dir.as_deref().unwrap_or(".");

    let header_file = format!("{output_dir}/{base_name}_header.txt");
    let reformatted_file = if cli.compress {
        format!("{output_dir}/{base_name}_reformatted.tsv.gz")
    } else {
        format!("{output_dir}/{base_name}_reformatted.tsv")
    };

    (header_file, reformatted_file)
}

#[allow(clippy::too_many_arguments)]
fn print_final_summary(
    cli: &Cli,
    header_file: &str,
    reformatted_file: &str,
    data: &(String, String, Vec<String>),
    records: &[reformat_vcf::ReformattedVcfRecord],
    read_time: std::time::Duration,
    process_time: std::time::Duration,
    write_time: std::time::Duration,
    header_write_time: std::time::Duration,
    total_time: std::time::Duration,
    variants_per_sec: f64,
    thread_count: usize,
    use_parallel: bool,
) {
    // Final summary
    println!("🎉 Processing completed successfully!");
    println!("📋 Summary:");
    println!("   📁 Files created:");
    println!("     • Header: {header_file}");
    println!(
        "     • Reformatted: {}{}",
        reformatted_file,
        if cli.compress {
            " (gzip compressed)"
        } else {
            ""
        }
    );
    println!();

    if cli.verbose {
        println!("⏱️  Performance Statistics:");
        println!(
            "   📖 File reading:     {:.2?} ({:.1}% of total)",
            read_time,
            (read_time.as_secs_f64() / total_time.as_secs_f64()) * 100.0
        );
        println!(
            "   🔄 Data processing:  {:.2?} ({:.1}% of total)",
            process_time,
            (process_time.as_secs_f64() / total_time.as_secs_f64()) * 100.0
        );
        println!(
            "   💾 File writing:     {:.2?} ({:.1}% of total)",
            write_time,
            (write_time.as_secs_f64() / total_time.as_secs_f64()) * 100.0
        );
        println!(
            "   📝 Header writing:   {:.2?} ({:.1}% of total)",
            header_write_time,
            (header_write_time.as_secs_f64() / total_time.as_secs_f64()) * 100.0
        );
        println!("   🕐 Total time:       {total_time:.2?}");
        println!();
        println!("📊 Processing Rates:");
        println!(
            "   🧬 Overall rate:     {:.0} variants/sec",
            data.2.len() as f64 / total_time.as_secs_f64()
        );
        println!(
            "   📝 Output rate:      {:.0} records/sec",
            records.len() as f64 / total_time.as_secs_f64()
        );
        println!("   💽 Processing only:  {variants_per_sec:.0} variants/sec");

        // Memory and efficiency stats
        println!();
        println!("📈 Efficiency Metrics:");
        println!(
            "   📊 Expansion ratio:  {:.2}x (output/input records)",
            records.len() as f64 / data.2.len() as f64
        );
        println!(
            "   🏃 Speed vs baseline: {}x faster than 1,000 variants/sec",
            (variants_per_sec / 1000.0).round() as u32
        );

        // Size estimation
        let estimated_input_size = estimate_file_size(&data.2);
        println!(
            "   💾 Estimated input:  {}",
            format_file_size(estimated_input_size)
        );

        if cli.compress {
            println!("   🗜️  Compression:     Enabled (estimated 60-80% size reduction)");
        }

        // Performance analysis
        print_performance_analysis(
            read_time,
            process_time,
            write_time,
            total_time,
            variants_per_sec,
            thread_count,
            use_parallel,
            data.2.len(),
        );
    } else {
        println!("⏱️  Total time: {total_time:.2?}");
        println!("📈 Processing rate: {variants_per_sec:.0} variants/sec");
    }
}
#[allow(clippy::too_many_arguments)]
fn print_performance_analysis(
    read_time: std::time::Duration,
    process_time: std::time::Duration,
    write_time: std::time::Duration,
    total_time: std::time::Duration,
    variants_per_sec: f64,
    thread_count: usize,
    use_parallel: bool,
    variant_count: usize,
) {
    println!();
    println!("🔍 Performance Analysis:");

    if use_parallel && thread_count > 1 {
        println!(
            "   ⚡ Parallel efficiency: {:.1}%",
            (variants_per_sec / (1000.0 * thread_count as f64)) * 100.0
        );
    }

    let bottleneck = if read_time > process_time && read_time > write_time {
        "File I/O (reading)"
    } else if process_time > write_time {
        "Data processing"
    } else {
        "File I/O (writing)"
    };

    println!("   🚧 Primary bottleneck: {bottleneck}");

    if variant_count > 10000 && !use_parallel {
        println!("   💡 Suggestion: Try parallel processing with -j 0 for better performance");
    }

    if variants_per_sec < 1000.0 {
        println!("   ⚠️  Low throughput detected - consider optimizing or using more threads");
    }
}

pub fn write_header_file(filename: &str, header: &str, content_length: usize) -> Result<()> {
    // Create the directory structure if it doesn't exist
    if let Some(parent) = Path::new(filename).parent() {
        create_dir_all(parent)?;
    }

    let mut file = File::create(filename)?;
    writeln!(file, "Content Length: {content_length}")?;
    writeln!(file, "Header:")?;
    file.write_all(header.as_bytes())?;
    Ok(())
}

fn get_base_filename(file_path: &str) -> String {
    let path = std::path::Path::new(file_path);
    let filename = path
        .file_name()
        .unwrap_or_default()
        .to_str()
        .unwrap_or("output");

    // Remove common VCF extensions
    let base = filename
        .strip_suffix(".vcf.gz")
        .unwrap_or(filename)
        .strip_suffix(".vcf")
        .unwrap_or(filename);

    base.to_string()
}

fn estimate_file_size(data_lines: &[String]) -> u64 {
    data_lines.iter().map(|line| line.len() as u64).sum()
}

fn format_file_size(size: u64) -> String {
    const UNITS: &[&str] = &["B", "KB", "MB", "GB", "TB"];
    let mut size = size as f64;
    let mut unit_index = 0;

    while size >= 1024.0 && unit_index < UNITS.len() - 1 {
        size /= 1024.0;
        unit_index += 1;
    }

    format!("{:.2} {}", size, UNITS[unit_index])
}
