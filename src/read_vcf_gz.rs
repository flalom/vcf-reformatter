use std::fs::File;
use std::io::{BufRead, BufReader};

pub fn read_vcf_gz(file_path: &str) -> Result<(String, String, Vec<String>), Box<dyn std::error::Error>> {
    // Open the file
    let mut file = File::open(file_path)?;
    println!("File {} opened successfully", file_path);
    // Make a buffer from the compressed file
    let mut compressed_file = flate2::read::MultiGzDecoder::new(&mut file);
    let buffer = BufReader::new(&mut compressed_file);
    // I could have used read_to_string() but then all is in memory and that is bad

    let mut header = String::new();
    let mut data_lines = Vec::new();
    let mut line_count = 0;
    let mut columns_title: String = String::new();

    for line in buffer.lines() {
        match line {
            Ok(line_content) => {
                line_count += 1;
                //println!("Line {}: {}", line_count, line_content); // Debug: print each line

                if line_content.starts_with("##") {
                    header.push_str(&line_content);
                    header.push('\n');
                } else if line_content.starts_with("#CHROM") {
                    columns_title = line_content.clone();
                    //println!("Columns title: {}", &columns_title);
                } else {
                    //let processed_line = line_content.replace("\t", " ");
                    //data_lines.push(processed_line.clone());
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
