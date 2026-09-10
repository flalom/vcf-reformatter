use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read};

/// A VCF opened for reading: header and column line in hand, data lines still to come.
///
/// The data lines are an iterator, not a `Vec`, so a caller can process the file in chunks and
/// never hold it whole. `read_vcf_gz` keeps the old collect-everything behaviour for callers
/// that have not been converted yet.
pub struct VcfStream {
    pub header: String,
    pub columns_title: String,
    pub lines: Box<dyn Iterator<Item = io::Result<String>>>,
}

/// Open a VCF (plain, gzipped, or `-` for stdin) and read only as far as the `#CHROM` line.
pub fn open_vcf(file_path: &str) -> Result<VcfStream, Box<dyn std::error::Error>> {
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

    split_header(reader)
}

// Kept for the library API — `tests/test.rs` exercises it. The binary now streams instead.
#[allow(dead_code)]
pub fn read_vcf_gz(
    file_path: &str,
) -> Result<(String, String, Vec<String>), Box<dyn std::error::Error>> {
    let stream = open_vcf(file_path)?;
    let data_lines = stream.lines.collect::<io::Result<Vec<String>>>()?;

    println!(
        "Total lines read: {}",
        stream.header.matches('\n').count() + 1 + data_lines.len()
    );
    println!("Header lines: {}", stream.header.matches('\n').count());
    println!("Data lines: {}", data_lines.len());

    Ok((stream.header, stream.columns_title, data_lines))
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

/// Consume the `##` lines and the `#CHROM` line, leaving the data lines unread.
fn split_header(mut reader: Box<dyn BufRead>) -> Result<VcfStream, Box<dyn std::error::Error>> {
    let mut header = String::new();
    let mut columns_title = String::new();
    let mut line = String::new();
    let mut lines_read = 0usize;

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        lines_read += 1;
        let trimmed = line.trim_end_matches(['\n', '\r']);
        if trimmed.starts_with("##") {
            header.push_str(trimmed);
            header.push('\n');
        } else if trimmed.starts_with("#CHROM") {
            columns_title = trimmed.to_string();
            break;
        }
    }

    // No #CHROM line means this was not a VCF at all — an empty stream, a truncated one, or
    // the wrong file. Reporting "completed successfully" over a 1-byte output file is worse
    // than failing. A real VCF carrying zero variants still has the line, and still succeeds.
    if columns_title.is_empty() {
        return Err(format!(
            "no #CHROM header line found after {lines_read} line(s) — is this a VCF?"
        )
        .into());
    }

    let lines = reader.lines().filter(|line| match line {
        Ok(text) => !text.trim().is_empty() && !text.starts_with('#'),
        Err(_) => true,
    });

    Ok(VcfStream {
        header,
        columns_title,
        lines: Box::new(lines),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The old collect-everything entry point, kept for the tests that predate `VcfStream`.
    fn parse_vcf_reader(
        reader: Box<dyn BufRead>,
    ) -> Result<(String, String, Vec<String>), Box<dyn std::error::Error>> {
        let stream = split_header(reader)?;
        let data = stream.lines.collect::<io::Result<Vec<String>>>()?;
        Ok((stream.header, stream.columns_title, data))
    }

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
        // Fewer than 2 bytes total: must not panic on the magic-number slice. It is still not
        // a VCF, so it is rejected rather than read as one data line.
        let reader = sniff_gzip_reader(Cursor::new(b"a".to_vec())).unwrap();
        assert!(parse_vcf_reader(reader).is_err());
    }

    #[test]
    fn test_sniff_gzip_reader_empty_input() {
        // An empty stream used to "complete successfully" over a 1-byte output file.
        let reader = sniff_gzip_reader(Cursor::new(Vec::new())).unwrap();
        let err = parse_vcf_reader(reader).unwrap_err().to_string();
        assert!(err.contains("#CHROM"), "unhelpful message: {err}");
    }

    #[test]
    fn test_vcf_with_header_but_no_variants_is_accepted() {
        // Zero variants is a legitimate VCF; only a missing #CHROM line is an error.
        let text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let reader = sniff_gzip_reader(Cursor::new(text.as_bytes().to_vec())).unwrap();
        let (_, columns, data) = parse_vcf_reader(reader).unwrap();
        assert!(columns.contains("CHROM"));
        assert!(data.is_empty());
    }
}
