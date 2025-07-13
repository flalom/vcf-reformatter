mod read_vcf_gz;
mod essentials_fields;
mod extract_csq_and_csq_names;
mod get_info_from_header;
mod reformat_vcf;
mod extract_sample_info;

use flate2::read::{GzDecoder, MultiGzDecoder};
use std::io::{BufReader, Write, Result};
use std::fs::File;
use std::path::Path;
use std::env;
use crate::read_vcf_gz::read_vcf_gz;
use crate::essentials_fields::*;
use crate::extract_csq_and_csq_names::extract_csq_regex;
use crate::reformat_vcf::{reformat_vcf_data, write_reformatted_vcf,*};
use crate::extract_sample_info::*;

fn main() {
    //println!("Hello, world!");
    //let data = read_vcf_gz("data/B505_V_1/B505_V_1.freebayes_VEP.ann.vcf.gz").unwrap();
    //println!("{:?}", data.0.len());
    //println!("{:?}", data.1.len());
    //println!("{:?}", data.2.len());
    let args: Vec<String> = env::args().collect();

    // Check if input file is provided
    if args.len() != 2 {
        eprintln!("Usage: {} <input_vcf_file>", args[0]);
        eprintln!("Example: {} data/sample.vcf.gz", args[0]);
        std::process::exit(1);
    }

    let input_path = &args[1];

    // Check if file exists
    if !Path::new(input_path).exists() {
        eprintln!("Error: File '{}' not found", input_path);
        std::process::exit(1);
    }

    // Read VCF file
    let data = match read_vcf_gz(input_path) {
        Ok(data) => data,
        Err(e) => {
            eprintln!("Error reading VCF file: {}", e);
            std::process::exit(1);
        }
    };

    // Generate output filenames based on input filename
    let base_name = get_base_filename(input_path);
    let header_file = format!("{}_header.txt", base_name);
    let content_file = format!("{}_content.tsv", base_name);
    let reformatted_file = format!("{}_reformatted.tsv", base_name);

    // Write original files
    if let Err(e) = write_header_file(&header_file, &data.0, data.2.len()) {
        eprintln!("Error writing header file: {}", e);
        std::process::exit(1);
    }

    /*
    if let Err(e) = write_vcf_content(&content_file, &data.1, &data.2) {
        eprintln!("Error writing content file: {}", e);
        std::process::exit(1);
    }
    */

    // Create reformatted version with separated INFO fields and FORMAT/sample fields
    match reformat_vcf_data_with_header(&data.0, &data.1, &data.2) {
        Ok((headers, records)) => {
            if let Err(e) = write_reformatted_vcf(&reformatted_file, &headers, &records) {
                eprintln!("Error writing reformatted file: {}", e);
                std::process::exit(1);
            }
            println!("Files created successfully:");
            println!("  Header: {}", header_file);
            //println!("  Content: {}", content_file);
            println!("  Reformatted: {}", reformatted_file);
        }
        Err(e) => {
            eprintln!("Error reformatting VCF data: {}", e);
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


pub fn write_header_file(filename: &str, header: &str, content_length: usize) -> Result<()> {
    let mut file = File::create(filename)?;
    writeln!(file, "Content Length: {}", content_length)?;
    writeln!(file, "Header:")?;
    file.write_all(header.as_bytes())?;
    Ok(())
}

pub fn write_vcf_content(filename: &str, column_names: &str, data_lines: &[String]) -> Result<()> {
    let mut file = File::create(filename)?;

    // Write column headers first (remove the # prefix for cleaner TSV)
    let clean_column_names = column_names.trim_start_matches('#');
    writeln!(file, "{}", clean_column_names)?;

    // Write data lines
    for line in data_lines {
        writeln!(file, "{}", line)?;
    }
    Ok(())
}

fn get_base_filename(file_path: &str) -> String {
    Path::new(file_path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("output")
        .to_string()
}


