use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read};

pub fn read_vcf_gz(
    file_path: &str,
) -> Result<(String, String, Vec<String>), Box<dyn std::error::Error>> {
    let reader: Box<dyn BufRead> = if file_path == "-" {
        println!("Reading VCF from stdin");
        sniff_gzip_reader(io::stdin())?
    } else {
        let file = File::open(file_path)?;
        println!("File {file_path} opened successfully");

        if file_path.ends_with(".gz") {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        }
    };

    parse_vcf_reader(reader)
}

/// Peek the first two bytes of `source` to detect the gzip magic number
/// (0x1f 0x8b), since a stream like stdin has no filename extension to check.
/// The peeked bytes are preserved via a `Chain`, so no data is lost.
fn sniff_gzip_reader<R: Read + 'static>(mut source: R) -> io::Result<Box<dyn BufRead>> {
    let mut magic = [0u8; 2];
    let bytes_read = source.read(&mut magic)?;
    let chained = Cursor::new(magic[..bytes_read].to_vec()).chain(source);

    if bytes_read == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(chained))))
    } else {
        Ok(Box::new(BufReader::new(chained)))
    }
}

fn parse_vcf_reader(
    reader: Box<dyn BufRead>,
) -> Result<(String, String, Vec<String>), Box<dyn std::error::Error>> {
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

    println!("Total lines read: {line_count}");
    println!("Header lines: {}", header.matches('\n').count());
    println!("Data lines: {}", data_lines.len());

    Ok((header, columns_title, data_lines))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plain_vcf() -> &'static str {
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t60\tPASS\tDP=10\n"
    }

    #[test]
    fn test_sniff_gzip_reader_plain_text() {
        let reader = sniff_gzip_reader(Cursor::new(plain_vcf().as_bytes().to_vec())).unwrap();
        let (header, columns, data) = parse_vcf_reader(reader).unwrap();
        assert!(header.contains("fileformat"));
        assert!(columns.contains("CHROM"));
        assert_eq!(data.len(), 1);
    }

    #[test]
    fn test_sniff_gzip_reader_gzip_bytes() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;

        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(plain_vcf().as_bytes()).unwrap();
        let compressed = encoder.finish().unwrap();

        let reader = sniff_gzip_reader(Cursor::new(compressed)).unwrap();
        let (header, columns, data) = parse_vcf_reader(reader).unwrap();
        assert!(header.contains("fileformat"));
        assert!(columns.contains("CHROM"));
        assert_eq!(data.len(), 1);
    }

    #[test]
    fn test_sniff_gzip_reader_short_input() {
        // Fewer than 2 bytes total: must not panic on the magic-number slice.
        let reader = sniff_gzip_reader(Cursor::new(b"a".to_vec())).unwrap();
        let (_, _, data) = parse_vcf_reader(reader).unwrap();
        assert_eq!(data.len(), 1); // treated as a single data line "a"
    }

    #[test]
    fn test_sniff_gzip_reader_empty_input() {
        let reader = sniff_gzip_reader(Cursor::new(Vec::new())).unwrap();
        let (header, columns, data) = parse_vcf_reader(reader).unwrap();
        assert!(header.is_empty());
        assert!(columns.is_empty());
        assert!(data.is_empty());
    }
}