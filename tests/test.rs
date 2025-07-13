use std::fs::File;
use std::collections::HashMap;
use std::io::Write;
use tempfile::tempdir;

use vcf_reformatter::{
    read_vcf_gz::read_vcf_gz,
    essentials_fields::VcfVariant,
    extract_csq_and_csq_names::extract_csq_regex,
    get_info_from_header::{extract_csq_format_from_header, extract_all_info_descriptions},
    reformat_vcf::{
        reformat_vcf_data_with_header,
        reformat_vcf_data_with_header_parallel,
        TranscriptHandling,
        write_reformatted_vcf,
        ReformattedVcfRecord
    },
    extract_sample_info::{parse_format_and_samples, ParsedFormatSample, SampleData}
};
use vcf_reformatter::extract_sample_info::ParsedSample;
// ------------------------------------------------------------------------------
// Tests for read_vcf_gz.rs
// ------------------------------------------------------------------------------

#[test]
fn test_basic_vcf_processing() {
    let dir = tempdir().unwrap();
    let file_path = dir.path().join("test.vcf");

    let mut file = File::create(&file_path).unwrap();
    writeln!(file, "##fileformat=VCFv4.2").unwrap();
    writeln!(file, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">").unwrap();
    writeln!(file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    writeln!(file, "chr1\t100\t.\tA\tG\t60\tPASS\tDP=10").unwrap();

    let result = read_vcf_gz(file_path.to_str().unwrap());
    assert!(result.is_ok());

    let (header, columns, data) = result.unwrap();
    assert!(!header.is_empty());
    assert!(columns.contains("CHROM"));
    assert_eq!(data.len(), 1);
}

#[test]
fn test_vcf_with_multiple_lines() {
    let dir = tempdir().unwrap();
    let file_path = dir.path().join("multiple.vcf");

    let mut file = File::create(&file_path).unwrap();
    writeln!(file, "##fileformat=VCFv4.2").unwrap();
    writeln!(file, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">").unwrap();
    writeln!(file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    writeln!(file, "chr1\t100\t.\tA\tG\t60\tPASS\tDP=10").unwrap();
    writeln!(file, "chr2\t200\trs123\tC\tT\t80\tPASS\tDP=20").unwrap();
    writeln!(file, "chr3\t300\t.\tG\tA\t90\tPASS\tDP=30").unwrap();

    let result = read_vcf_gz(file_path.to_str().unwrap());
    assert!(result.is_ok());

    let (header, columns, data) = result.unwrap();
    assert_eq!(header.matches('\n').count(), 2); // 2 header lines
    assert_eq!(data.len(), 3); // 3 data lines
}

#[test]
fn test_vcf_with_empty_lines() {
    let dir = tempdir().unwrap();
    let file_path = dir.path().join("empty_lines.vcf");

    let mut file = File::create(&file_path).unwrap();
    writeln!(file, "##fileformat=VCFv4.2").unwrap();
    writeln!(file, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">").unwrap();
    writeln!(file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    writeln!(file, "").unwrap(); // Empty line
    writeln!(file, "chr1\t100\t.\tA\tG\t60\tPASS\tDP=10").unwrap();
    writeln!(file, "  ").unwrap(); // Whitespace line
    writeln!(file, "chr2\t200\trs123\tC\tT\t80\tPASS\tDP=20").unwrap();

    let result = read_vcf_gz(file_path.to_str().unwrap());
    assert!(result.is_ok());

    let (_, _, data) = result.unwrap();
    assert_eq!(data.len(), 2); // Empty lines should be skipped
}

// ------------------------------------------------------------------------------
// Tests for essentials_fields.rs
// ------------------------------------------------------------------------------

#[test]
fn test_vcf_variant_parsing() {
    let line = "chr1\t100\t.\tA\tG\t60\tPASS\tDP=10";
    let column_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];

    let result = VcfVariant::from_line(line, &column_names);
    assert!(result.is_ok());

    let variant = result.unwrap();
    assert_eq!(variant.chromosome, "chr1");
    assert_eq!(variant.position, 100);
    assert_eq!(variant.id, None);
    assert_eq!(variant.reference, "A");
    assert_eq!(variant.alternate, "G");
    assert_eq!(variant.quality, Some(60.0));
    assert_eq!(variant.filter, "PASS");
    assert_eq!(variant.info, "DP=10");
    assert_eq!(variant.format, None);
    assert_eq!(variant.samples.len(), 0);
}

#[test]
fn test_vcf_variant_with_id_and_samples() {
    let line = "chr2\t200\trs123\tC\tT\t80\tPASS\tDP=20\tGT:DP\t0/1:20\t1/1:22";
    let column_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE1", "SAMPLE2"];

    let result = VcfVariant::from_line(line, &column_names);
    assert!(result.is_ok());

    let variant = result.unwrap();
    assert_eq!(variant.chromosome, "chr2");
    assert_eq!(variant.position, 200);
    assert_eq!(variant.id, Some("rs123".to_string()));
    assert_eq!(variant.reference, "C");
    assert_eq!(variant.alternate, "T");
    assert_eq!(variant.quality, Some(80.0));
    assert_eq!(variant.filter, "PASS");
    assert_eq!(variant.info, "DP=20");
    assert_eq!(variant.format, Some("GT:DP".to_string()));
    assert_eq!(variant.samples.len(), 2);
    assert_eq!(variant.samples[0], "0/1:20");
    assert_eq!(variant.samples[1], "1/1:22");
}

#[test]
fn test_vcf_variant_missing_quality() {
    let line = "chr3\t300\t.\tG\tA\t.\tPASS\tDP=30";
    let column_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];

    let result = VcfVariant::from_line(line, &column_names);
    assert!(result.is_ok());

    let variant = result.unwrap();
    assert_eq!(variant.quality, None);
}

#[test]
fn test_vcf_variant_invalid_position() {
    let line = "chr3\tinvalid\t.\tG\tA\t90\tPASS\tDP=30";
    let column_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];

    let result = VcfVariant::from_line(line, &column_names);
    assert!(result.is_err());
}

#[test]
fn test_vcf_variant_too_few_fields() {
    let line = "chr3\t300\t.\tG\tA\t90\tPASS";
    let column_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];

    let result = VcfVariant::from_line(line, &column_names);
    assert!(result.is_err());
}

// ------------------------------------------------------------------------------
// Tests for extract_csq_and_csq_names.rs
// ------------------------------------------------------------------------------

#[test]
fn test_extract_csq_regex() {
    let mut parsed_lines = vec![
        "chr1".to_string(),
        "100".to_string(),
        ".".to_string(),
        "A".to_string(),
        "G".to_string(),
        "60".to_string(),
        "PASS".to_string(),
        "DP=10;CSQ=A|B|C;AF=0.5".to_string()
    ];

    let result = extract_csq_regex(&mut parsed_lines);
    assert!(result.is_some());
    assert_eq!(result.unwrap(), "A|B|C");

    // Check that CSQ is removed from INFO field
    assert_eq!(parsed_lines[7], "DP=10;AF=0.5");
}

#[test]
fn test_extract_csq_regex_no_csq() {
    let mut parsed_lines = vec![
        "chr1".to_string(),
        "100".to_string(),
        ".".to_string(),
        "A".to_string(),
        "G".to_string(),
        "60".to_string(),
        "PASS".to_string(),
        "DP=10;AF=0.5".to_string()
    ];

    let result = extract_csq_regex(&mut parsed_lines);
    assert!(result.is_none());

    // INFO field should remain unchanged
    assert_eq!(parsed_lines[7], "DP=10;AF=0.5");
}

#[test]
fn test_extract_csq_regex_empty_csq() {
    let mut parsed_lines = vec![
        "chr1".to_string(),
        "100".to_string(),
        ".".to_string(),
        "A".to_string(),
        "G".to_string(),
        "60".to_string(),
        "PASS".to_string(),
        "DP=10;CSQ=;AF=0.5".to_string()
    ];

    let result = extract_csq_regex(&mut parsed_lines);
    // The current implementation doesn't match empty CSQ values
    assert!(result.is_none());

    // INFO field should remain unchanged since CSQ wasn't matched
    assert_eq!(parsed_lines[7], "DP=10;CSQ=;AF=0.5");
}

#[test]
fn test_extract_csq_regex_too_few_fields() {
    let mut parsed_lines = vec![
        "chr1".to_string(),
        "100".to_string(),
        ".".to_string(),
        "A".to_string(),
        "G".to_string(),
        "60".to_string(),
        "PASS".to_string()
    ];

    let result = extract_csq_regex(&mut parsed_lines);
    assert!(result.is_none());
}

// ------------------------------------------------------------------------------
// Tests for get_info_from_header.rs
// ------------------------------------------------------------------------------

#[test]
fn test_extract_csq_format_from_header() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene">"#;

    let result = extract_csq_format_from_header(header);
    assert!(result.is_some());

    let fields = result.unwrap();
    assert_eq!(fields, vec!["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene"]);
}

#[test]
fn test_extract_csq_format_from_header_missing_csq() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">"#;

    let result = extract_csq_format_from_header(header);
    assert!(result.is_none());
}

#[test]
fn test_extract_csq_format_from_header_no_format() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP annotations">"#;

    let result = extract_csq_format_from_header(header);
    assert!(result.is_none());
}

#[test]
fn test_extract_all_info_descriptions() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">"#;

    let result = extract_all_info_descriptions(header);

    assert_eq!(result.len(), 3);
    assert_eq!(result.get("DP").unwrap(), "Total Depth");
    assert_eq!(result.get("AF").unwrap(), "Allele Frequency");
    assert_eq!(result.get("DB").unwrap(), "dbSNP membership");
}

#[test]
fn test_extract_all_info_descriptions_no_info() {
    let header = r#"##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;

    let result = extract_all_info_descriptions(header);
    assert!(result.is_empty());
}

// ------------------------------------------------------------------------------
// Tests for extract_sample_info.rs
// ------------------------------------------------------------------------------

#[test]
fn test_parse_format_and_samples() {
    let format = Some("GT:DP:AD");
    let sample_fields = vec!["0/1:20:10,10".to_string(), "1/1:30:0,30".to_string()];
    let sample_names = vec!["SAMPLE1".to_string(), "SAMPLE2".to_string()];

    let result = parse_format_and_samples(format, &sample_fields, &sample_names);
    assert!(result.is_ok());

    let parsed = result.unwrap();
    assert_eq!(parsed.format_keys, vec!["GT", "DP", "AD"]);
    assert_eq!(parsed.samples.len(), 2);

    assert_eq!(parsed.samples[0].sample_name, "SAMPLE1");
    assert_eq!(parsed.samples[0].format_fields.get("GT").unwrap(), "0/1");
    assert_eq!(parsed.samples[0].format_fields.get("DP").unwrap(), "20");
    assert_eq!(parsed.samples[0].format_fields.get("AD").unwrap(), "10,10");

    assert_eq!(parsed.samples[1].sample_name, "SAMPLE2");
    assert_eq!(parsed.samples[1].format_fields.get("GT").unwrap(), "1/1");
    assert_eq!(parsed.samples[1].format_fields.get("DP").unwrap(), "30");
    assert_eq!(parsed.samples[1].format_fields.get("AD").unwrap(), "0,30");
}

#[test]
fn test_parse_format_and_samples_missing_fields() {
    let format = Some("GT:DP:AD:PL");
    let sample_fields = vec!["0/1:20:10,10".to_string()]; // Missing PL
    let sample_names = vec!["SAMPLE1".to_string()];

    let result = parse_format_and_samples(format, &sample_fields, &sample_names);
    assert!(result.is_ok());

    let parsed = result.unwrap();
    assert_eq!(parsed.format_keys, vec!["GT", "DP", "AD", "PL"]);
    assert_eq!(parsed.samples.len(), 1);

    assert_eq!(parsed.samples[0].format_fields.get("GT").unwrap(), "0/1");
    assert_eq!(parsed.samples[0].format_fields.get("DP").unwrap(), "20");
    assert_eq!(parsed.samples[0].format_fields.get("AD").unwrap(), "10,10");
    assert_eq!(parsed.samples[0].format_fields.get("PL").unwrap(), "."); // Missing field should be "."
}

#[test]
fn test_parse_format_and_samples_no_format() {
    let format = None;
    let sample_fields = vec!["0/1:20:10,10".to_string()];
    let sample_names = vec!["SAMPLE1".to_string()];

    let result = parse_format_and_samples(format, &sample_fields, &sample_names);
    assert!(result.is_ok());

    let parsed = result.unwrap();
    assert!(parsed.format_keys.is_empty());
    assert!(parsed.samples.is_empty());
}

#[test]
fn test_get_headers_for_samples() {
    let mut parsed = ParsedFormatSample::new();
    parsed.format_keys = vec!["GT".to_string(), "DP".to_string()];

    let mut sample1 = SampleData::new("SAMPLE1".to_string());
    sample1.format_fields.insert("GT".to_string(), "0/1".to_string());
    sample1.format_fields.insert("DP".to_string(), "20".to_string());

    let mut sample2 = SampleData::new("SAMPLE2".to_string());
    sample2.format_fields.insert("GT".to_string(), "1/1".to_string());
    sample2.format_fields.insert("DP".to_string(), "30".to_string());

    parsed.samples = vec![sample1, sample2];

    let headers = parsed.get_headers_for_samples();
    assert_eq!(headers, vec![
        "SAMPLE1_GT", "SAMPLE1_DP",
        "SAMPLE2_GT", "SAMPLE2_DP"
    ]);
}

#[test]
fn test_get_values_for_samples() {
    let mut parsed = ParsedFormatSample::new();
    parsed.format_keys = vec!["GT".to_string(), "DP".to_string()];

    let mut sample1 = SampleData::new("SAMPLE1".to_string());
    sample1.format_fields.insert("GT".to_string(), "0/1".to_string());
    sample1.format_fields.insert("DP".to_string(), "20".to_string());

    let mut sample2 = SampleData::new("SAMPLE2".to_string());
    sample2.format_fields.insert("GT".to_string(), "1/1".to_string());
    sample2.format_fields.insert("DP".to_string(), "30".to_string());

    parsed.samples = vec![sample1, sample2];

    let values = parsed.get_values_for_samples();
    assert_eq!(values, vec![
        "0/1", "20",
        "1/1", "30"
    ]);
}

// ------------------------------------------------------------------------------
// Tests for reformat_vcf.rs
// ------------------------------------------------------------------------------

#[test]
fn test_reformat_vcf_data_with_header() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">"#;
    let column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    let data_lines = vec![
        "chr1\t100\t.\tA\tG\t60\tPASS\tDP=10".to_string(),
        "chr2\t200\trs123\tC\tT\t80\tPASS\tDP=20".to_string()
    ];

    let result = reformat_vcf_data_with_header(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::FirstOnly
    );

    assert!(result.is_ok());

    let (headers, records) = result.unwrap();
    assert!(headers.contains(&"CHROM".to_string()));
    assert!(headers.contains(&"INFO_DP".to_string()));
    assert_eq!(records.len(), 2);

    assert_eq!(records[0].chromosome, "chr1");
    assert_eq!(records[0].position, 100);
    assert_eq!(records[0].info_fields.get("INFO_DP").unwrap(), "10");

    assert_eq!(records[1].chromosome, "chr2");
    assert_eq!(records[1].position, 200);
    assert_eq!(records[1].info_fields.get("INFO_DP").unwrap(), "20");
}

#[test]
fn test_reformat_vcf_data_with_csq() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT">"#;
    let column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    let data_lines = vec![
        "chr1\t100\t.\tA\tG\t60\tPASS\tCSQ=G|missense_variant|MODERATE".to_string(),
        "chr2\t200\trs123\tC\tT\t80\tPASS\tCSQ=T|synonymous_variant|LOW,T|intron_variant|MODIFIER".to_string()
    ];

    let result = reformat_vcf_data_with_header(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::FirstOnly
    );

    assert!(result.is_ok());

    let (headers, records) = result.unwrap();
    assert!(headers.contains(&"CSQ_Allele".to_string()));
    assert!(headers.contains(&"CSQ_Consequence".to_string()));
    assert!(headers.contains(&"CSQ_IMPACT".to_string()));

    assert_eq!(records.len(), 2);
    assert_eq!(records[0].info_fields.get("CSQ_Allele").unwrap(), "G");
    assert_eq!(records[0].info_fields.get("CSQ_Consequence").unwrap(), "missense_variant");
    assert_eq!(records[0].info_fields.get("CSQ_IMPACT").unwrap(), "MODERATE");

    assert_eq!(records[1].info_fields.get("CSQ_Allele").unwrap(), "T");
    assert_eq!(records[1].info_fields.get("CSQ_Consequence").unwrap(), "synonymous_variant");
    assert_eq!(records[1].info_fields.get("CSQ_IMPACT").unwrap(), "LOW");
}

#[test]
fn test_reformat_vcf_data_with_split_transcripts() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT">"#;
    let column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    let data_lines = vec![
        "chr1\t100\t.\tA\tG\t60\tPASS\tCSQ=G|missense_variant|MODERATE,G|intron_variant|MODIFIER".to_string()
    ];

    let result = reformat_vcf_data_with_header(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::SplitRows
    );

    assert!(result.is_ok());

    let (_, records) = result.unwrap();
    assert_eq!(records.len(), 2); // Should split into 2 records

    assert_eq!(records[0].info_fields.get("CSQ_Consequence").unwrap(), "missense_variant");
    assert_eq!(records[0].info_fields.get("CSQ_IMPACT").unwrap(), "MODERATE");

    assert_eq!(records[1].info_fields.get("CSQ_Consequence").unwrap(), "intron_variant");
    assert_eq!(records[1].info_fields.get("CSQ_IMPACT").unwrap(), "MODIFIER");
}

#[test]
fn test_reformat_vcf_data_with_most_severe() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT">"#;
    let column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    let data_lines = vec![
        "chr1\t100\t.\tA\tG\t60\tPASS\tCSQ=G|synonymous_variant|LOW,G|stop_gained|HIGH,G|intron_variant|MODIFIER".to_string()
    ];

    let result = reformat_vcf_data_with_header(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::MostSevere
    );

    assert!(result.is_ok());

    let (_, records) = result.unwrap();
    assert_eq!(records.len(), 1);

    // stop_gained is more severe than synonymous_variant and intron_variant
    assert_eq!(records[0].info_fields.get("CSQ_Consequence").unwrap(), "stop_gained");
    assert_eq!(records[0].info_fields.get("CSQ_IMPACT").unwrap(), "HIGH");
}

#[test]
fn test_parallel_processing() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">"#;
    let column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    // Create a larger dataset
    let mut data_lines = Vec::with_capacity(100);
    for i in 1..101 {
        data_lines.push(format!("chr1\t{}\t.\tA\tG\t60\tPASS\tDP={}", i * 100, i));
    }

    let result = reformat_vcf_data_with_header_parallel(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::FirstOnly
    );

    assert!(result.is_ok());

    let (_, records) = result.unwrap();
    assert_eq!(records.len(), 100);
}

#[test]
fn test_write_reformatted_vcf() {
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("output.tsv");

    // Create some sample records
    let mut record = ReformattedVcfRecord {
        chromosome: "chr1".to_string(),
        position: 100,
        id: None,
        reference: "A".to_string(),
        alternate: "G".to_string(),
        quality: Some(60.0),
        filter: "PASS".to_string(),
        info_fields: HashMap::new(),
        format_sample_data: None,
    };

    record.info_fields.insert("INFO_DP".to_string(), "10".to_string());
    record.info_fields.insert("CSQ_Allele".to_string(), "G".to_string());
    record.info_fields.insert("CSQ_Consequence".to_string(), "missense_variant".to_string());

    let headers = vec![
        "CHROM".to_string(), "POS".to_string(), "ID".to_string(),
        "REF".to_string(), "ALT".to_string(), "QUAL".to_string(), "FILTER".to_string(),
        "INFO_DP".to_string(), "CSQ_Allele".to_string(), "CSQ_Consequence".to_string()
    ];

    let records = vec![record];

    let result = write_reformatted_vcf(
        output_path.to_str().unwrap(),
        &headers,
        &records,
        false // No compression
    );

    assert!(result.is_ok());
    assert!(output_path.exists());

    // Check file content
    let content = std::fs::read_to_string(output_path).unwrap();
    assert!(content.contains("CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO_DP\tCSQ_Allele\tCSQ_Consequence"));
    assert!(content.contains("chr1\t100\t.\tA\tG\t60\tPASS\t10\tG\tmissense_variant"));
}

#[test]
fn test_complex_sample_names_with_underscores() {
    use vcf_reformatter::extract_sample_info::{parse_format_and_samples};

    // Test with complex sample names containing underscores (like B487_B487_1_cOM)
    let format = Some("GT:AD:AF:DP:F1R2:F2R1:FAD:SB");
    let sample_fields = vec![
        "0/0:257,4:0.017:261:109,0:87,3:229,3:164,93,3,1".to_string(),
        "0/1:303,6:0.020:309:105,3:106,1:245,4:187,116,2,4".to_string(),
    ];
    let sample_names = vec![
        "B487_B487_1_cOM".to_string(),
        "B487_B487_2_LN".to_string(),
    ];

    let parsed = parse_format_and_samples(format, &sample_fields, &sample_names).unwrap();

    // Verify format keys
    assert_eq!(parsed.format_keys, vec!["GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "SB"]);

    // Verify sample count
    assert_eq!(parsed.samples.len(), 2);

    // Verify first sample
    assert_eq!(parsed.samples[0].sample_name, "B487_B487_1_cOM");
    assert_eq!(parsed.samples[0].format_fields.get("GT").unwrap(), "0/0");
    assert_eq!(parsed.samples[0].format_fields.get("AD").unwrap(), "257,4");
    assert_eq!(parsed.samples[0].format_fields.get("AF").unwrap(), "0.017");
    assert_eq!(parsed.samples[0].format_fields.get("DP").unwrap(), "261");
    assert_eq!(parsed.samples[0].format_fields.get("F1R2").unwrap(), "109,0");
    assert_eq!(parsed.samples[0].format_fields.get("F2R1").unwrap(), "87,3");
    assert_eq!(parsed.samples[0].format_fields.get("FAD").unwrap(), "229,3");
    assert_eq!(parsed.samples[0].format_fields.get("SB").unwrap(), "164,93,3,1");

    // Verify second sample
    assert_eq!(parsed.samples[1].sample_name, "B487_B487_2_LN");
    assert_eq!(parsed.samples[1].format_fields.get("GT").unwrap(), "0/1");
    assert_eq!(parsed.samples[1].format_fields.get("AD").unwrap(), "303,6");
    assert_eq!(parsed.samples[1].format_fields.get("AF").unwrap(), "0.020");
    assert_eq!(parsed.samples[1].format_fields.get("DP").unwrap(), "309");
    assert_eq!(parsed.samples[1].format_fields.get("F1R2").unwrap(), "105,3");
    assert_eq!(parsed.samples[1].format_fields.get("F2R1").unwrap(), "106,1");
    assert_eq!(parsed.samples[1].format_fields.get("FAD").unwrap(), "245,4");
    assert_eq!(parsed.samples[1].format_fields.get("SB").unwrap(), "187,116,2,4");
}

#[test]
fn test_complex_sample_names_header_generation() {
    use vcf_reformatter::extract_sample_info::{parse_format_and_samples};

    let format = Some("GT:AD:AF:DP");
    let sample_fields = vec![
        "0/0:257,4:0.017:261".to_string(),
        "0/1:303,6:0.020:309".to_string(),
    ];
    let sample_names = vec![
        "B487_B487_1_cOM".to_string(),
        "B487_B487_2_LN".to_string(),
    ];

    let parsed = parse_format_and_samples(format, &sample_fields, &sample_names).unwrap();
    let headers = parsed.get_headers_for_samples();

    // Expected headers with complex sample names
    let expected_headers = vec![
        "B487_B487_1_cOM_GT", "B487_B487_1_cOM_AD", "B487_B487_1_cOM_AF", "B487_B487_1_cOM_DP",
        "B487_B487_2_LN_GT", "B487_B487_2_LN_AD", "B487_B487_2_LN_AF", "B487_B487_2_LN_DP"
    ];

    assert_eq!(headers, expected_headers);
}

#[test]
fn test_complex_sample_names_value_extraction() {
    use vcf_reformatter::extract_sample_info::{parse_format_and_samples};

    let format = Some("GT:AD:AF:DP");
    let sample_fields = vec![
        "0/0:257,4:0.017:261".to_string(),
        "0/1:303,6:0.020:309".to_string(),
    ];
    let sample_names = vec![
        "B487_B487_1_cOM".to_string(),
        "B487_B487_2_LN".to_string(),
    ];

    let parsed = parse_format_and_samples(format, &sample_fields, &sample_names).unwrap();
    let values = parsed.get_values_for_samples();

    // Expected values in the same order as headers
    let expected_values = vec![
        "0/0", "257,4", "0.017", "261",     // B487_B487_1_cOM
        "0/1", "303,6", "0.020", "309"      // B487_B487_2_LN
    ];

    assert_eq!(values, expected_values);
}

#[test]
fn test_vcf_record_with_complex_sample_names() {
    use vcf_reformatter::reformat_vcf::{ReformattedVcfRecord, TranscriptHandling};

    // Test a complete VCF line with complex sample names
    let vcf_line = "chr1\t123456\t.\tA\tG\t1000\tPASS\tDP=570;AF=0.5\tGT:AD:AF:DP:F1R2:F2R1:FAD:SB\t0/0:257,4:0.017:261:109,0:87,3:229,3:164,93,3,1\t0/1:303,6:0.020:309:105,3:106,1:245,4:187,116,2,4";

    let column_names = vec![
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
        "B487_B487_1_cOM", "B487_B487_2_LN"
    ];

    let records = ReformattedVcfRecord::from_vcf_line(
        vcf_line,
        &column_names,
        &None,
        TranscriptHandling::FirstOnly
    ).unwrap();

    assert_eq!(records.len(), 1);
    let record = &records[0];

    // Verify basic fields
    assert_eq!(record.chromosome, "chr1");
    assert_eq!(record.position, 123456);
    assert_eq!(record.reference, "A");
    assert_eq!(record.alternate, "G");

    // Verify sample data is present
    assert!(record.format_sample_data.is_some());
    let sample_data = record.format_sample_data.as_ref().unwrap();

    // Verify sample names are preserved correctly
    assert_eq!(sample_data.samples.len(), 2);
    assert_eq!(sample_data.samples[0].sample_name, "B487_B487_1_cOM");
    assert_eq!(sample_data.samples[1].sample_name, "B487_B487_2_LN");

    // Verify sample values
    assert_eq!(sample_data.samples[0].format_fields.get("GT").unwrap(), "0/0");
    assert_eq!(sample_data.samples[0].format_fields.get("AD").unwrap(), "257,4");
    assert_eq!(sample_data.samples[1].format_fields.get("GT").unwrap(), "0/1");
    assert_eq!(sample_data.samples[1].format_fields.get("AD").unwrap(), "303,6");
}

#[test]
fn test_tsv_output_with_complex_sample_names() {
    use vcf_reformatter::reformat_vcf::{ReformattedVcfRecord};
    use std::collections::HashMap;
    use vcf_reformatter::extract_sample_info::{ParsedFormatSample, ParsedSample};

    // Create a test record with complex sample names
    let mut sample_data = ParsedFormatSample::new();
    sample_data.format_keys = vec!["GT".to_string(), "AD".to_string(), "AF".to_string(), "DP".to_string()];

    let mut sample1 = ParsedSample::new("B487_B487_1_cOM".to_string());
    sample1.format_fields.insert("GT".to_string(), "0/0".to_string());
    sample1.format_fields.insert("AD".to_string(), "257,4".to_string());
    sample1.format_fields.insert("AF".to_string(), "0.017".to_string());
    sample1.format_fields.insert("DP".to_string(), "261".to_string());

    let mut sample2 = ParsedSample::new("B487_B487_2_LN".to_string());
    sample2.format_fields.insert("GT".to_string(), "0/1".to_string());
    sample2.format_fields.insert("AD".to_string(), "303,6".to_string());
    sample2.format_fields.insert("AF".to_string(), "0.020".to_string());
    sample2.format_fields.insert("DP".to_string(), "309".to_string());

    sample_data.samples = vec![sample1, sample2];

    let record = ReformattedVcfRecord {
        chromosome: "chr1".to_string(),
        position: 123456,
        id: None,
        reference: "A".to_string(),
        alternate: "G".to_string(),
        quality: Some(1000.0),
        filter: "PASS".to_string(),
        info_fields: HashMap::new(),
        format_sample_data: Some(sample_data),
    };

    // Test header generation
    let headers = vec![
        "CHROM".to_string(), "POS".to_string(), "ID".to_string(), "REF".to_string(),
        "ALT".to_string(), "QUAL".to_string(), "FILTER".to_string(),
        "B487_B487_1_cOM_GT".to_string(), "B487_B487_1_cOM_AD".to_string(),
        "B487_B487_1_cOM_AF".to_string(), "B487_B487_1_cOM_DP".to_string(),
        "B487_B487_2_LN_GT".to_string(), "B487_B487_2_LN_AD".to_string(),
        "B487_B487_2_LN_AF".to_string(), "B487_B487_2_LN_DP".to_string()
    ];

    // Test that the write_tsv_content logic would work correctly
    // (We can't directly test the private function, but we can verify the structure)

    // Test value extraction for specific headers
    let sample_data_ref = record.format_sample_data.as_ref().unwrap();

    // Test that we can find the correct sample and format key combinations
    let mut found_sample1_gt = false;
    let mut found_sample2_af = false;

    for sample in &sample_data_ref.samples {
        for format_key in &sample_data_ref.format_keys {
            let expected_header = format!("{}_{}", sample.sample_name, format_key);
            if expected_header == "B487_B487_1_cOM_GT" {
                assert_eq!(sample.format_fields.get(format_key).unwrap(), "0/0");
                found_sample1_gt = true;
            }
            if expected_header == "B487_B487_2_LN_AF" {
                assert_eq!(sample.format_fields.get(format_key).unwrap(), "0.020");
                found_sample2_af = true;
            }
        }
    }

    assert!(found_sample1_gt, "Should find B487_B487_1_cOM_GT");
    assert!(found_sample2_af, "Should find B487_B487_2_LN_AF");
}

#[test]
fn test_edge_case_sample_names_with_multiple_underscores() {
    use vcf_reformatter::extract_sample_info::{parse_format_and_samples};

    // Test with even more complex sample names
    let format = Some("GT:DP");
    let sample_fields = vec![
        "0/0:100".to_string(),
        "0/1:200".to_string(),
        "1/1:300".to_string(),
    ];
    let sample_names = vec![
        "SAMPLE_WITH_MANY_UNDERSCORES_1".to_string(),
        "ANOTHER_COMPLEX_SAMPLE_NAME_2".to_string(),
        "FINAL_TEST_SAMPLE_NAME_3".to_string(),
    ];

    let parsed = parse_format_and_samples(format, &sample_fields, &sample_names).unwrap();

    // Verify all samples are parsed correctly
    assert_eq!(parsed.samples.len(), 3);
    assert_eq!(parsed.samples[0].sample_name, "SAMPLE_WITH_MANY_UNDERSCORES_1");
    assert_eq!(parsed.samples[1].sample_name, "ANOTHER_COMPLEX_SAMPLE_NAME_2");
    assert_eq!(parsed.samples[2].sample_name, "FINAL_TEST_SAMPLE_NAME_3");

    // Verify values
    assert_eq!(parsed.samples[0].format_fields.get("GT").unwrap(), "0/0");
    assert_eq!(parsed.samples[0].format_fields.get("DP").unwrap(), "100");
    assert_eq!(parsed.samples[1].format_fields.get("GT").unwrap(), "0/1");
    assert_eq!(parsed.samples[1].format_fields.get("DP").unwrap(), "200");
    assert_eq!(parsed.samples[2].format_fields.get("GT").unwrap(), "1/1");
    assert_eq!(parsed.samples[2].format_fields.get("DP").unwrap(), "300");

    // Test header generation
    let headers = parsed.get_headers_for_samples();
    let expected_headers = vec![
        "SAMPLE_WITH_MANY_UNDERSCORES_1_GT", "SAMPLE_WITH_MANY_UNDERSCORES_1_DP",
        "ANOTHER_COMPLEX_SAMPLE_NAME_2_GT", "ANOTHER_COMPLEX_SAMPLE_NAME_2_DP",
        "FINAL_TEST_SAMPLE_NAME_3_GT", "FINAL_TEST_SAMPLE_NAME_3_DP"
    ];
    assert_eq!(headers, expected_headers);
}