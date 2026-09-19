//! `--parquet`: the sink writes the same rows and columns as the text output, and the CLI
//! flags around it (`--parquet -c`, `--parquet-only`) do what the help text says.
#![cfg(feature = "parquet_out")]
mod common;

use parquet::file::reader::{FileReader, SerializedFileReader};
use tempfile::tempdir;

use common::{bin, maf, write_vcf, TINY_VCF};
use vcf_reformatter::essentials_fields::MafRecord;
use vcf_reformatter::parquet_writer::ParquetSink;
use vcf_reformatter::reformat_vcf::{reformat_vcf_data_with_header, TranscriptHandling};

const SNPEFF_HEADER: &str = "##fileformat=VCFv4.2\n##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'\">";
const COLUMNS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

fn two_snpeff_lines() -> Vec<String> {
    [
        "chr1\t100\trs1\tA\tG\t60\tPASS\tDP=50;ANN=G|missense_variant|MODERATE|BRCA1|ENSG1|transcript|ENST1|protein_coding|5/24|c.1T>C|p.M1T|1/100|1/100|1/33||",
        "chr2\t200\t.\tC\tT\t40\tPASS\tDP=30;ANN=T|stop_gained|HIGH|TP53|ENSG2|transcript|ENST2|protein_coding|7/11|c.2C>T|p.Q2*|2/200|2/200|2/66||",
    ]
    .map(String::from)
    .to_vec()
}

fn open(path: &str) -> SerializedFileReader<std::fs::File> {
    SerializedFileReader::new(std::fs::File::open(path).unwrap()).unwrap()
}

fn row_count(path: &str) -> i64 {
    open(path).metadata().file_metadata().num_rows()
}

fn column_names(path: &str) -> Vec<String> {
    open(path)
        .metadata()
        .file_metadata()
        .schema_descr()
        .columns()
        .iter()
        .map(|c| c.name().to_string())
        .collect()
}

// --- the sink -------------------------------------------------------------------------------

#[test]
fn tsv_sink_rows_and_columns_match_the_text_layout() {
    let (headers, records) = reformat_vcf_data_with_header(
        SNPEFF_HEADER,
        COLUMNS,
        &two_snpeff_lines(),
        TranscriptHandling::FirstOnly,
    )
    .unwrap();
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().to_str().unwrap();

    let mut sink = ParquetSink::create_tsv(path, &headers).unwrap();
    sink.write_tsv(&records).unwrap();
    sink.close().unwrap();

    assert_eq!(row_count(path), records.len() as i64);
    assert_eq!(column_names(path), headers);
}

#[test]
fn maf_sink_rows_and_columns_match_get_maf_headers() {
    let (_, records) = reformat_vcf_data_with_header(
        SNPEFF_HEADER,
        COLUMNS,
        &two_snpeff_lines(),
        TranscriptHandling::FirstOnly,
    )
    .unwrap();
    let rows: Vec<MafRecord> = records.iter().map(maf).collect();
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().to_str().unwrap();

    let mut sink = ParquetSink::create_maf(path).unwrap();
    sink.write_maf(&rows).unwrap();
    sink.close().unwrap();

    assert_eq!(row_count(path), rows.len() as i64);
    assert_eq!(column_names(path), MafRecord::get_maf_headers());
}

// --- CLI flags ------------------------------------------------------------------------------

#[test]
fn parquet_filename_keeps_the_format_extension_on_both_paths() {
    let dir = tempdir().unwrap();
    let vcf = write_vcf(dir.path(), "in.vcf", TINY_VCF);
    let out = dir.path().to_str().unwrap();
    for (format, expected) in [
        ("tsv", "pq_reformatted.tsv.parquet"),
        ("maf", "pq_reformatted.maf.parquet"),
    ] {
        let o = bin()
            .args([
                &vcf,
                "-o",
                out,
                "-p",
                "pq",
                "--output-format",
                format,
                "--report",
                "none",
                "--parquet",
            ])
            // no sample columns: MAF refuses to run without a barcode
            .args(["--sample-barcode", "TUMOR"])
            .output()
            .unwrap();
        assert!(
            o.status.success(),
            "stderr: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        assert!(
            dir.path().join(expected).exists(),
            "--output-format {format} should write {expected}"
        );
    }
}

#[test]
fn parquet_with_compress_gzips_the_text_and_names_parquet_without_gz() {
    // Issue #6: --parquet used to reject -c outright. Now the text file is gzipped and the
    // parquet sits next to it under the un-gzipped name, not "X.tsv.gz.parquet".
    let dir = tempdir().unwrap();
    let vcf = write_vcf(dir.path(), "in.vcf", TINY_VCF);
    let out = dir.path().to_str().unwrap();
    for (format, text, parquet) in [
        ("tsv", "pq_reformatted.tsv.gz", "pq_reformatted.tsv.parquet"),
        ("maf", "pq_reformatted.maf.gz", "pq_reformatted.maf.parquet"),
    ] {
        let o = bin()
            .args([
                &vcf,
                "-o",
                out,
                "-p",
                "pq",
                "--output-format",
                format,
                "--report",
                "none",
            ])
            .args(["--parquet", "-c", "--sample-barcode", "TUMOR"])
            .output()
            .unwrap();
        assert!(
            o.status.success(),
            "stderr: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        let bytes = std::fs::read(dir.path().join(text)).expect(text);
        assert_eq!(&bytes[..2], &[0x1f, 0x8b], "{text} should be gzip");
        assert!(dir.path().join(parquet).exists(), "{parquet} missing");
    }
}

#[test]
fn parquet_only_writes_no_text_file() {
    let dir = tempdir().unwrap();
    let vcf = write_vcf(dir.path(), "in.vcf", TINY_VCF);
    let out = dir.path().to_str().unwrap();
    for (format, text, parquet) in [
        ("tsv", "pq_reformatted.tsv", "pq_reformatted.tsv.parquet"),
        ("maf", "pq_reformatted.maf", "pq_reformatted.maf.parquet"),
    ] {
        let o = bin()
            .args([
                &vcf,
                "-o",
                out,
                "-p",
                "pq",
                "--output-format",
                format,
                "--report",
                "none",
            ])
            .args(["--parquet-only", "--sample-barcode", "TUMOR"])
            .output()
            .unwrap();
        assert!(
            o.status.success(),
            "stderr: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        assert!(!dir.path().join(text).exists(), "{text} should not exist");
        let p = dir.path().join(parquet);
        assert!(p.exists(), "{parquet} missing");
        assert_eq!(row_count(p.to_str().unwrap()), 1, "{parquet}");
    }
}

#[test]
fn parquet_only_rejects_compress() {
    let dir = tempdir().unwrap();
    let vcf = write_vcf(dir.path(), "in.vcf", TINY_VCF);
    let o = bin()
        .args([
            &vcf,
            "-o",
            dir.path().to_str().unwrap(),
            "--report",
            "none",
            "--parquet-only",
            "-c",
        ])
        .output()
        .unwrap();
    assert!(!o.status.success(), "--parquet-only -c should fail");
    let stderr = String::from_utf8_lossy(&o.stderr);
    assert!(
        stderr.contains("--parquet-only") && !stderr.contains("unexpected argument"),
        "stderr: {stderr}"
    );
}
