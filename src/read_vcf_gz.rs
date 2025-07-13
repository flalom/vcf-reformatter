use std::fs::File;
use std::io::{BufRead, BufReader};
use flate2::read::MultiGzDecoder;


pub fn read_vcf_gz(file_path: &str) -> Result<(String, String, Vec<String>), Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
    println!("File {} opened successfully", file_path);

    let reader: Box<dyn BufRead> = if file_path.ends_with(".gz") {
        // Handle gzip compressed files
        let decoder = MultiGzDecoder::new(file);
        Box::new(BufReader::new(decoder))
    } else {
        // Handle uncompressed files
        Box::new(BufReader::new(file))
    };

    let mut header = String::new();
    let mut data_lines = Vec::new();
    let mut line_count = 0;
    let mut columns_title = String::new();

    for line in reader.lines() {
        match line {
            Ok(line_content) => {
                line_count += 1;

                if line_content.starts_with("##") {
                    header.push_str(&line_content);
                    header.push('\n');
                } else if line_content.starts_with("#CHROM") {
                    columns_title = line_content.clone();
                } else if !line_content.trim().is_empty() {
                    data_lines.push(line_content);
                }
            }
            Err(e) => {
                println!("Error reading line {}: {}", line_count + 1, e);
                return Err(e.into());
            }
        }
    }

    println!("Total lines read: {}", line_count);
    println!("Header lines: {}", header.matches('\n').count());
    println!("Data lines: {}", data_lines.len());

    Ok((header, columns_title, data_lines))
}

