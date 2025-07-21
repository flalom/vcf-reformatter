use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use tempfile::tempdir;

use vcf_reformatter::{
    essentials_fields::VcfVariant,
    extract_csq_and_csq_names::extract_csq_regex,
    //extract_ann_and_ann_names::extract_ann_regex,
    extract_sample_info::{parse_format_and_samples, ParsedFormatSample, _SampleData},
    get_info_from_header::{_extract_all_info_descriptions, extract_csq_format_from_header},
    read_vcf_gz::read_vcf_gz,
    reformat_vcf::{
        get_ann_impact_severity, parse_info_field, reformat_vcf_data_with_header,
        reformat_vcf_data_with_header_parallel, sanitize_field_name, write_reformatted_vcf,
        ReformattedVcfRecord, TranscriptHandling,
    },
};

// ------------------------------------------------------------------------------
// Tests for read_vcf_gz.rs
// ------------------------------------------------------------------------------

#[test]
fn test_basic_vcf_processing() {
    let dir = tempdir().unwrap();
    let file_path = dir.path().join("test.vcf");

    let mut file = File::create(&file_path).unwrap();
    writeln!(file, "##fileformat=VCFv4.2").unwrap();
    writeln!(
        file,
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
    )
    .unwrap();
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
    writeln!(
        file,
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
    )
    .unwrap();
    writeln!(file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    writeln!(file, "chr1\t100\t.\tA\tG\t60\tPASS\tDP=10").unwrap();
    writeln!(file, "chr2\t200\trs123\tC\tT\t80\tPASS\tDP=20").unwrap();
    writeln!(file, "chr3\t300\t.\tG\tA\t90\tPASS\tDP=30").unwrap();

    let result = read_vcf_gz(file_path.to_str().unwrap());
    assert!(result.is_ok());

    let (header, _columns, data) = result.unwrap();
    assert_eq!(header.matches('\n').count(), 2); // 2 header lines
    assert_eq!(data.len(), 3); // 3 data lines
}

#[test]
fn test_vcf_with_empty_lines() {
    let dir = tempdir().unwrap();
    let file_path = dir.path().join("empty_lines.vcf");

    let mut file = File::create(&file_path).unwrap();
    writeln!(file, "##fileformat=VCFv4.2").unwrap();
    writeln!(
        file,
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
    )
    .unwrap();
    writeln!(file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    writeln!(file).unwrap(); // Empty line
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
    let column_names = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE1",
        "SAMPLE2",
    ];

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
        "DP=10;CSQ=A|B|C;AF=0.5".to_string(),
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
        "DP=10;AF=0.5".to_string(),
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
        "DP=10;CSQ=;AF=0.5".to_string(),
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
        "PASS".to_string(),
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
    assert_eq!(
        fields,
        vec!["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene"]
    );
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

    let result = _extract_all_info_descriptions(header);

    assert_eq!(result.len(), 3);
    assert_eq!(result.get("DP").unwrap(), "Total Depth");
    assert_eq!(result.get("AF").unwrap(), "Allele Frequency");
    assert_eq!(result.get("DB").unwrap(), "dbSNP membership");
}

#[test]
fn test_extract_all_info_descriptions_no_info() {
    let header = r#"##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;

    let result = _extract_all_info_descriptions(header);
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

    let mut sample1 = _SampleData::_new("SAMPLE1".to_string());
    sample1
        .format_fields
        .insert("GT".to_string(), "0/1".to_string());
    sample1
        .format_fields
        .insert("DP".to_string(), "20".to_string());

    let mut sample2 = _SampleData::_new("SAMPLE2".to_string());
    sample2
        .format_fields
        .insert("GT".to_string(), "1/1".to_string());
    sample2
        .format_fields
        .insert("DP".to_string(), "30".to_string());

    parsed.samples = vec![sample1, sample2];

    let headers = parsed.get_headers_for_samples();
    assert_eq!(
        headers,
        vec!["SAMPLE1_GT", "SAMPLE1_DP", "SAMPLE2_GT", "SAMPLE2_DP"]
    );
}

#[test]
fn test_get_values_for_samples() {
    let mut parsed = ParsedFormatSample::new();
    parsed.format_keys = vec!["GT".to_string(), "DP".to_string()];

    let mut sample1 = _SampleData::_new("SAMPLE1".to_string());
    sample1
        .format_fields
        .insert("GT".to_string(), "0/1".to_string());
    sample1
        .format_fields
        .insert("DP".to_string(), "20".to_string());

    let mut sample2 = _SampleData::_new("SAMPLE2".to_string());
    sample2
        .format_fields
        .insert("GT".to_string(), "1/1".to_string());
    sample2
        .format_fields
        .insert("DP".to_string(), "30".to_string());

    parsed.samples = vec![sample1, sample2];

    let values = parsed._get_values_for_samples();
    assert_eq!(values, vec!["0/1", "20", "1/1", "30"]);
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
        "chr2\t200\trs123\tC\tT\t80\tPASS\tDP=20".to_string(),
    ];

    let result = reformat_vcf_data_with_header(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::FirstOnly,
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
        "chr2\t200\trs123\tC\tT\t80\tPASS\tCSQ=T|synonymous_variant|LOW,T|intron_variant|MODIFIER"
            .to_string(),
    ];

    let result = reformat_vcf_data_with_header(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::FirstOnly,
    );

    assert!(result.is_ok());

    let (headers, records) = result.unwrap();
    assert!(headers.contains(&"CSQ_Allele".to_string()));
    assert!(headers.contains(&"CSQ_Consequence".to_string()));
    assert!(headers.contains(&"CSQ_IMPACT".to_string()));

    assert_eq!(records.len(), 2);
    assert_eq!(records[0].info_fields.get("CSQ_Allele").unwrap(), "G");
    assert_eq!(
        records[0].info_fields.get("CSQ_Consequence").unwrap(),
        "missense_variant"
    );
    assert_eq!(
        records[0].info_fields.get("CSQ_IMPACT").unwrap(),
        "MODERATE"
    );

    assert_eq!(records[1].info_fields.get("CSQ_Allele").unwrap(), "T");
    assert_eq!(
        records[1].info_fields.get("CSQ_Consequence").unwrap(),
        "synonymous_variant"
    );
    assert_eq!(records[1].info_fields.get("CSQ_IMPACT").unwrap(), "LOW");
}

#[test]
fn test_reformat_vcf_data_with_split_transcripts() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT">"#;
    let column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    let data_lines = vec![
        "chr1\t100\t.\tA\tG\t60\tPASS\tCSQ=G|missense_variant|MODERATE,G|intron_variant|MODIFIER"
            .to_string(),
    ];

    let result = reformat_vcf_data_with_header(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::SplitRows,
    );

    assert!(result.is_ok());

    let (_, records) = result.unwrap();
    assert_eq!(records.len(), 2); // Should split into 2 records

    assert_eq!(
        records[0].info_fields.get("CSQ_Consequence").unwrap(),
        "missense_variant"
    );
    assert_eq!(
        records[0].info_fields.get("CSQ_IMPACT").unwrap(),
        "MODERATE"
    );

    assert_eq!(
        records[1].info_fields.get("CSQ_Consequence").unwrap(),
        "intron_variant"
    );
    assert_eq!(
        records[1].info_fields.get("CSQ_IMPACT").unwrap(),
        "MODIFIER"
    );
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
        TranscriptHandling::MostSevere,
    );

    assert!(result.is_ok());

    let (_, records) = result.unwrap();
    assert_eq!(records.len(), 1);

    // stop_gained is more severe than synonymous_variant and intron_variant
    assert_eq!(
        records[0].info_fields.get("CSQ_Consequence").unwrap(),
        "stop_gained"
    );
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
        TranscriptHandling::FirstOnly,
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

    record
        .info_fields
        .insert("INFO_DP".to_string(), "10".to_string());
    record
        .info_fields
        .insert("CSQ_Allele".to_string(), "G".to_string());
    record.info_fields.insert(
        "CSQ_Consequence".to_string(),
        "missense_variant".to_string(),
    );

    let headers = vec![
        "CHROM".to_string(),
        "POS".to_string(),
        "ID".to_string(),
        "REF".to_string(),
        "ALT".to_string(),
        "QUAL".to_string(),
        "FILTER".to_string(),
        "INFO_DP".to_string(),
        "CSQ_Allele".to_string(),
        "CSQ_Consequence".to_string(),
    ];

    let records = vec![record];

    let result = write_reformatted_vcf(
        output_path.to_str().unwrap(),
        &headers,
        &records,
        false, // No compression
    );

    assert!(result.is_ok());
    assert!(output_path.exists());

    // Check file content
    let content = std::fs::read_to_string(output_path).unwrap();
    assert!(content
        .contains("CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO_DP\tCSQ_Allele\tCSQ_Consequence"));
    assert!(content.contains("chr1\t100\t.\tA\tG\t60\tPASS\t10\tG\tmissense_variant"));
}

#[test]
fn test_complex_sample_names_with_underscores() {
    use vcf_reformatter::extract_sample_info::parse_format_and_samples;

    // Test with complex sample names containing underscores (like B487_B487_1_cOM)
    let format = Some("GT:AD:AF:DP:F1R2:F2R1:FAD:SB");
    let sample_fields = vec![
        "0/0:257,4:0.017:261:109,0:87,3:229,3:164,93,3,1".to_string(),
        "0/1:303,6:0.020:309:105,3:106,1:245,4:187,116,2,4".to_string(),
    ];
    let sample_names = vec!["B487_B487_1_cOM".to_string(), "B487_B487_2_LN".to_string()];

    let parsed = parse_format_and_samples(format, &sample_fields, &sample_names).unwrap();

    // Verify format keys
    assert_eq!(
        parsed.format_keys,
        vec!["GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "SB"]
    );

    // Verify sample count
    assert_eq!(parsed.samples.len(), 2);

    // Verify first sample
    assert_eq!(parsed.samples[0].sample_name, "B487_B487_1_cOM");
    assert_eq!(parsed.samples[0].format_fields.get("GT").unwrap(), "0/0");
    assert_eq!(parsed.samples[0].format_fields.get("AD").unwrap(), "257,4");
    assert_eq!(parsed.samples[0].format_fields.get("AF").unwrap(), "0.017");
    assert_eq!(parsed.samples[0].format_fields.get("DP").unwrap(), "261");
    assert_eq!(
        parsed.samples[0].format_fields.get("F1R2").unwrap(),
        "109,0"
    );
    assert_eq!(parsed.samples[0].format_fields.get("F2R1").unwrap(), "87,3");
    assert_eq!(parsed.samples[0].format_fields.get("FAD").unwrap(), "229,3");
    assert_eq!(
        parsed.samples[0].format_fields.get("SB").unwrap(),
        "164,93,3,1"
    );

    // Verify second sample
    assert_eq!(parsed.samples[1].sample_name, "B487_B487_2_LN");
    assert_eq!(parsed.samples[1].format_fields.get("GT").unwrap(), "0/1");
    assert_eq!(parsed.samples[1].format_fields.get("AD").unwrap(), "303,6");
    assert_eq!(parsed.samples[1].format_fields.get("AF").unwrap(), "0.020");
    assert_eq!(parsed.samples[1].format_fields.get("DP").unwrap(), "309");
    assert_eq!(
        parsed.samples[1].format_fields.get("F1R2").unwrap(),
        "105,3"
    );
    assert_eq!(
        parsed.samples[1].format_fields.get("F2R1").unwrap(),
        "106,1"
    );
    assert_eq!(parsed.samples[1].format_fields.get("FAD").unwrap(), "245,4");
    assert_eq!(
        parsed.samples[1].format_fields.get("SB").unwrap(),
        "187,116,2,4"
    );
}

#[test]
fn test_complex_sample_names_header_generation() {
    use vcf_reformatter::extract_sample_info::parse_format_and_samples;

    let format = Some("GT:AD:AF:DP");
    let sample_fields = vec![
        "0/0:257,4:0.017:261".to_string(),
        "0/1:303,6:0.020:309".to_string(),
    ];
    let sample_names = vec!["B487_B487_1_cOM".to_string(), "B487_B487_2_LN".to_string()];

    let parsed = parse_format_and_samples(format, &sample_fields, &sample_names).unwrap();
    let headers = parsed.get_headers_for_samples();

    // Expected headers with complex sample names
    let expected_headers = vec![
        "B487_B487_1_cOM_GT",
        "B487_B487_1_cOM_AD",
        "B487_B487_1_cOM_AF",
        "B487_B487_1_cOM_DP",
        "B487_B487_2_LN_GT",
        "B487_B487_2_LN_AD",
        "B487_B487_2_LN_AF",
        "B487_B487_2_LN_DP",
    ];

    assert_eq!(headers, expected_headers);
}

#[test]
fn test_complex_sample_names_value_extraction() {
    use vcf_reformatter::extract_sample_info::parse_format_and_samples;

    let format = Some("GT:AD:AF:DP");
    let sample_fields = vec![
        "0/0:257,4:0.017:261".to_string(),
        "0/1:303,6:0.020:309".to_string(),
    ];
    let sample_names = vec!["B487_B487_1_cOM".to_string(), "B487_B487_2_LN".to_string()];

    let parsed = parse_format_and_samples(format, &sample_fields, &sample_names).unwrap();
    let values = parsed._get_values_for_samples();

    // Expected values in the same order as headers
    let expected_values = vec![
        "0/0", "257,4", "0.017", "261", // B487_B487_1_cOM
        "0/1", "303,6", "0.020", "309", // B487_B487_2_LN
    ];

    assert_eq!(values, expected_values);
}

#[test]
fn test_vcf_record_with_complex_sample_names() {
    use vcf_reformatter::reformat_vcf::{ReformattedVcfRecord, TranscriptHandling};

    // Test a complete VCF line with complex sample names
    let vcf_line = "chr1\t123456\t.\tA\tG\t1000\tPASS\tDP=570;AF=0.5\tGT:AD:AF:DP:F1R2:F2R1:FAD:SB\t0/0:257,4:0.017:261:109,0:87,3:229,3:164,93,3,1\t0/1:303,6:0.020:309:105,3:106,1:245,4:187,116,2,4";

    let column_names = vec![
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "B487_B487_1_cOM",
        "B487_B487_2_LN",
    ];

    let records = ReformattedVcfRecord::from_vcf_line(
        vcf_line,
        &column_names,
        &None,                         // csq_field_names
        &None,                         // ann_field_names
        TranscriptHandling::FirstOnly, // transcript_handling
    )
    .unwrap();

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
    assert_eq!(
        sample_data.samples[0].format_fields.get("GT").unwrap(),
        "0/0"
    );
    assert_eq!(
        sample_data.samples[0].format_fields.get("AD").unwrap(),
        "257,4"
    );
    assert_eq!(
        sample_data.samples[1].format_fields.get("GT").unwrap(),
        "0/1"
    );
    assert_eq!(
        sample_data.samples[1].format_fields.get("AD").unwrap(),
        "303,6"
    );
}

#[test]
fn test_tsv_output_with_complex_sample_names() {
    use std::collections::HashMap;
    use vcf_reformatter::extract_sample_info::{ParsedFormatSample, ParsedSample};
    use vcf_reformatter::reformat_vcf::ReformattedVcfRecord;

    // Create a test record with complex sample names
    let mut sample_data = ParsedFormatSample::new();
    sample_data.format_keys = vec![
        "GT".to_string(),
        "AD".to_string(),
        "AF".to_string(),
        "DP".to_string(),
    ];

    let mut sample1 = ParsedSample::_new("B487_B487_1_cOM".to_string());
    sample1
        .format_fields
        .insert("GT".to_string(), "0/0".to_string());
    sample1
        .format_fields
        .insert("AD".to_string(), "257,4".to_string());
    sample1
        .format_fields
        .insert("AF".to_string(), "0.017".to_string());
    sample1
        .format_fields
        .insert("DP".to_string(), "261".to_string());

    let mut sample2 = ParsedSample::_new("B487_B487_2_LN".to_string());
    sample2
        .format_fields
        .insert("GT".to_string(), "0/1".to_string());
    sample2
        .format_fields
        .insert("AD".to_string(), "303,6".to_string());
    sample2
        .format_fields
        .insert("AF".to_string(), "0.020".to_string());
    sample2
        .format_fields
        .insert("DP".to_string(), "309".to_string());

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
    let _headers = vec![
        "CHROM".to_string(),
        "POS".to_string(),
        "ID".to_string(),
        "REF".to_string(),
        "ALT".to_string(),
        "QUAL".to_string(),
        "FILTER".to_string(),
        "B487_B487_1_cOM_GT".to_string(),
        "B487_B487_1_cOM_AD".to_string(),
        "B487_B487_1_cOM_AF".to_string(),
        "B487_B487_1_cOM_DP".to_string(),
        "B487_B487_2_LN_GT".to_string(),
        "B487_B487_2_LN_AD".to_string(),
        "B487_B487_2_LN_AF".to_string(),
        "B487_B487_2_LN_DP".to_string(),
    ];

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
    use vcf_reformatter::extract_sample_info::parse_format_and_samples;

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
    assert_eq!(
        parsed.samples[0].sample_name,
        "SAMPLE_WITH_MANY_UNDERSCORES_1"
    );
    assert_eq!(
        parsed.samples[1].sample_name,
        "ANOTHER_COMPLEX_SAMPLE_NAME_2"
    );
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
        "SAMPLE_WITH_MANY_UNDERSCORES_1_GT",
        "SAMPLE_WITH_MANY_UNDERSCORES_1_DP",
        "ANOTHER_COMPLEX_SAMPLE_NAME_2_GT",
        "ANOTHER_COMPLEX_SAMPLE_NAME_2_DP",
        "FINAL_TEST_SAMPLE_NAME_3_GT",
        "FINAL_TEST_SAMPLE_NAME_3_DP",
    ];
    assert_eq!(headers, expected_headers);
}

#[test]
fn test_parse_info_field_empty() {
    let result = parse_info_field("", &None, &None, TranscriptHandling::FirstOnly);
    assert!(result.is_ok());
    assert_eq!(result.unwrap().len(), 1);
}

#[test]
fn test_parse_info_field_no_annotations() {
    let csq_fields = None;
    let ann_fields = None;
    let result = parse_info_field(
        "DP=10;AF=0.5",
        &csq_fields,
        &ann_fields,
        TranscriptHandling::FirstOnly,
    );

    assert!(result.is_ok());
    let records = result.unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].get("INFO_DP"), Some(&"10".to_string()));
    assert_eq!(records[0].get("INFO_AF"), Some(&"0.5".to_string()));
}

#[test]
fn test_parse_info_field_with_csq() {
    let csq_fields = Some(vec!["Allele".to_string(), "Consequence".to_string()]);
    let ann_fields = None;
    let result = parse_info_field(
        "DP=10;CSQ=A|missense_variant,T|synonymous_variant;AF=0.5",
        &csq_fields,
        &ann_fields,
        TranscriptHandling::SplitRows,
    );

    assert!(result.is_ok());
    let records = result.unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].get("CSQ_Allele"), Some(&"A".to_string()));
    assert_eq!(
        records[0].get("CSQ_Consequence"),
        Some(&"missense_variant".to_string())
    );
    assert_eq!(records[0].get("INFO_DP"), Some(&"10".to_string()));
}

#[test]
fn test_sanitize_field_name() {
    assert_eq!(sanitize_field_name("HGVS.c"), "HGVS_c");
    assert_eq!(
        sanitize_field_name("cDNA.pos / cDNA.length"),
        "cDNA_pos___cDNA_length"
    );
    assert_eq!(sanitize_field_name("normal_name"), "normal_name");
    assert_eq!(sanitize_field_name(""), "");
}

#[test]
fn test_get_ann_impact_severity() {
    assert_eq!(get_ann_impact_severity("HIGH"), 4);
    assert_eq!(get_ann_impact_severity("MODERATE"), 3);
    assert_eq!(get_ann_impact_severity("LOW"), 2);
    assert_eq!(get_ann_impact_severity("MODIFIER"), 1);
    assert_eq!(get_ann_impact_severity("UNKNOWN"), 0);
}

#[test]
fn test_parse_info_field_with_ann() {
    let csq_fields = None;
    let ann_fields = Some(vec![
        "Allele".to_string(),
        "Annotation".to_string(),
        "Annotation_Impact".to_string(),
        "Gene_Name".to_string(),
    ]);
    let result = parse_info_field(
        "DP=15;ANN=T|missense_variant|MODERATE|BRCA1,G|synonymous_variant|LOW|TP53;AF=0.3",
        &csq_fields,
        &ann_fields,
        TranscriptHandling::SplitRows,
    );

    assert!(result.is_ok());
    let records = result.unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].get("ANN_Allele"), Some(&"T".to_string()));
    assert_eq!(
        records[0].get("ANN_Annotation"),
        Some(&"missense_variant".to_string())
    );
    assert_eq!(
        records[0].get("ANN_Annotation_Impact"),
        Some(&"MODERATE".to_string())
    );
    assert_eq!(records[0].get("ANN_Gene_Name"), Some(&"BRCA1".to_string()));
    assert_eq!(records[0].get("INFO_DP"), Some(&"15".to_string()));
    assert_eq!(records[0].get("INFO_AF"), Some(&"0.3".to_string()));
}

//############# add tests for snpeff
#[test]
fn test_snpeff_ann_parsing_basic() {
    use vcf_reformatter::reformat_vcf::{ReformattedVcfRecord, TranscriptHandling};

    let vcf_line = "chr1\t100\t.\tA\tG\t60\tPASS\tANN=G|missense_variant|MODERATE|BRCA1|ENSG00000012048|transcript|ENST00000357654|protein_coding|5/24|c.181T>C|p.Cys61Arg|181/5592|181/4863|61/1620||\tGT\t0/1";

    let column_names = vec![
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE1",
    ];

    // Define ANN field names as they would be extracted from header
    let ann_field_names = Some(vec![
        "Allele".to_string(),
        "Annotation".to_string(),
        "Annotation_Impact".to_string(),
        "Gene_Name".to_string(),
        "Gene_ID".to_string(),
        "Feature_Type".to_string(),
        "Feature_ID".to_string(),
        "Transcript_BioType".to_string(),
        "Rank".to_string(),
        "HGVS.c".to_string(),
        "HGVS.p".to_string(),
        "cDNA.pos / cDNA.length".to_string(),
        "CDS.pos / CDS.length".to_string(),
        "AA.pos / AA.length".to_string(),
        "Distance".to_string(),
    ]);

    let records = ReformattedVcfRecord::from_vcf_line(
        vcf_line,
        &column_names,
        &None, // csq_field_names
        &ann_field_names,
        TranscriptHandling::FirstOnly,
    )
    .unwrap();

    assert_eq!(records.len(), 1);
    let record = &records[0];

    // Verify ANN fields are properly parsed
    assert_eq!(record.info_fields.get("ANN_Allele"), Some(&"G".to_string()));
    assert_eq!(
        record.info_fields.get("ANN_Annotation"),
        Some(&"missense_variant".to_string())
    );
    assert_eq!(
        record.info_fields.get("ANN_Annotation_Impact"),
        Some(&"MODERATE".to_string())
    );
    assert_eq!(
        record.info_fields.get("ANN_Gene_Name"),
        Some(&"BRCA1".to_string())
    );
    assert_eq!(
        record.info_fields.get("ANN_Gene_ID"),
        Some(&"ENSG00000012048".to_string())
    );
    assert_eq!(
        record.info_fields.get("ANN_Feature_ID"),
        Some(&"ENST00000357654".to_string())
    );
    assert_eq!(
        record.info_fields.get("ANN_HGVS_c"),
        Some(&"c.181T>C".to_string())
    );
    assert_eq!(
        record.info_fields.get("ANN_HGVS_p"),
        Some(&"p.Cys61Arg".to_string())
    );
}

#[test]
fn test_snpeff_ann_multiple_transcripts_split_rows() {
    use vcf_reformatter::reformat_vcf::{ReformattedVcfRecord, TranscriptHandling};

    let vcf_line = "chr1\t100\t.\tA\tG\t60\tPASS\tANN=G|missense_variant|MODERATE|BRCA1|ENSG00000012048|transcript|ENST00000357654|protein_coding|5/24|c.181T>C|p.Cys61Arg||||,G|synonymous_variant|LOW|BRCA1|ENSG00000012048|transcript|ENST00000123456|protein_coding|6/25|c.200A>G|p.Leu67Leu||||";

    let column_names = vec!["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];

    let ann_field_names = Some(vec![
        "Allele".to_string(),
        "Annotation".to_string(),
        "Annotation_Impact".to_string(),
        "Gene_Name".to_string(),
        "Gene_ID".to_string(),
        "Feature_Type".to_string(),
        "Feature_ID".to_string(),
        "Transcript_BioType".to_string(),
        "Rank".to_string(),
        "HGVS.c".to_string(),
        "HGVS.p".to_string(),
        "cDNA.pos / cDNA.length".to_string(),
        "CDS.pos / CDS.length".to_string(),
        "AA.pos / AA.length".to_string(),
        "Distance".to_string(),
    ]);

    let records = ReformattedVcfRecord::from_vcf_line(
        vcf_line,
        &column_names,
        &None,
        &ann_field_names,
        TranscriptHandling::SplitRows, // create separate records?
    )
    .unwrap();

    assert_eq!(records.len(), 2);

    // First record - missense variant
    assert_eq!(
        records[0].info_fields.get("ANN_Annotation"),
        Some(&"missense_variant".to_string())
    );
    assert_eq!(
        records[0].info_fields.get("ANN_Annotation_Impact"),
        Some(&"MODERATE".to_string())
    );
    assert_eq!(
        records[0].info_fields.get("ANN_Feature_ID"),
        Some(&"ENST00000357654".to_string())
    );

    // Second record - synonymous variant
    assert_eq!(
        records[1].info_fields.get("ANN_Annotation"),
        Some(&"synonymous_variant".to_string())
    );
    assert_eq!(
        records[1].info_fields.get("ANN_Annotation_Impact"),
        Some(&"LOW".to_string())
    );
    assert_eq!(
        records[1].info_fields.get("ANN_Feature_ID"),
        Some(&"ENST00000123456".to_string())
    );
}

#[test]
fn test_snpeff_ann_most_severe_selection() {
    use vcf_reformatter::reformat_vcf::{ReformattedVcfRecord, TranscriptHandling};

    let vcf_line = "chr1\t100\t.\tA\tG\t60\tPASS\tANN=G|synonymous_variant|LOW|BRCA1||||||||||||,G|missense_variant|MODERATE|BRCA1||||||||||||,G|stop_gained|HIGH|BRCA1||||||||||||";

    let column_names = vec!["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];

    let ann_field_names = Some(vec![
        "Allele".to_string(),
        "Annotation".to_string(),
        "Annotation_Impact".to_string(),
        "Gene_Name".to_string(),
        "Gene_ID".to_string(),
        "Feature_Type".to_string(),
        "Feature_ID".to_string(),
        "Transcript_BioType".to_string(),
        "Rank".to_string(),
        "HGVS.c".to_string(),
        "HGVS.p".to_string(),
        "cDNA.pos / cDNA.length".to_string(),
        "CDS.pos / CDS.length".to_string(),
        "AA.pos / AA.length".to_string(),
        "Distance".to_string(),
    ]);

    let records = ReformattedVcfRecord::from_vcf_line(
        vcf_line,
        &column_names,
        &None,
        &ann_field_names,
        TranscriptHandling::MostSevere, // should select the HIGH impact variant
    )
    .unwrap();

    assert_eq!(records.len(), 1);

    // select the stop_gained variant (HIGH impact)
    assert_eq!(
        records[0].info_fields.get("ANN_Annotation"),
        Some(&"stop_gained".to_string())
    );
    assert_eq!(
        records[0].info_fields.get("ANN_Annotation_Impact"),
        Some(&"HIGH".to_string())
    );
}

#[test]
fn test_snpeff_ann_impact_severity_ranking() {
    use vcf_reformatter::reformat_vcf::get_ann_impact_severity;

    // the severity function
    assert_eq!(get_ann_impact_severity("HIGH"), 4);
    assert_eq!(get_ann_impact_severity("MODERATE"), 3);
    assert_eq!(get_ann_impact_severity("LOW"), 2);
    assert_eq!(get_ann_impact_severity("MODIFIER"), 1);
    assert_eq!(get_ann_impact_severity("UNKNOWN"), 0);

    // Test case insensitivity
    assert_eq!(get_ann_impact_severity("high"), 4);
    assert_eq!(get_ann_impact_severity("moderate"), 3);

    // Test with whitespace
    assert_eq!(get_ann_impact_severity(" HIGH "), 4);
}

#[test]
fn test_vep_and_snpeff_field_name_sanitization() {
    use vcf_reformatter::reformat_vcf::sanitize_field_name;

    // Test VEP field name sanitization
    assert_eq!(sanitize_field_name("HGVS.c"), "HGVS_c");
    assert_eq!(sanitize_field_name("HGVS.p"), "HGVS_p");

    // Test SnpEff field name sanitization
    assert_eq!(
        sanitize_field_name("cDNA.pos / cDNA.length"),
        "cDNA_pos___cDNA_length"
    );
    assert_eq!(
        sanitize_field_name("CDS.pos / CDS.length"),
        "CDS_pos___CDS_length"
    );
    assert_eq!(
        sanitize_field_name("AA.pos / AA.length"),
        "AA_pos___AA_length"
    );

    // Test other special characters
    assert_eq!(
        sanitize_field_name("ERRORS / WARNINGS / INFO"),
        "ERRORS___WARNINGS___INFO"
    );

    // Test normal names
    assert_eq!(sanitize_field_name("Gene_Name"), "Gene_Name");
    assert_eq!(sanitize_field_name("Allele"), "Allele");
}

#[test]
fn test_mixed_info_fields_with_ann() {
    use vcf_reformatter::reformat_vcf::{ReformattedVcfRecord, TranscriptHandling};

    let vcf_line = "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50;AF=0.25;ANN=G|missense_variant|MODERATE|BRCA1||||||||||||;AC=2;AN=4\tGT\t0/1";

    let column_names = vec![
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE1",
    ];

    let ann_field_names = Some(vec![
        "Allele".to_string(),
        "Annotation".to_string(),
        "Annotation_Impact".to_string(),
        "Gene_Name".to_string(),
    ]);

    let records = ReformattedVcfRecord::from_vcf_line(
        vcf_line,
        &column_names,
        &None,
        &ann_field_names,
        TranscriptHandling::FirstOnly,
    )
    .unwrap();

    assert_eq!(records.len(), 1);
    let record = &records[0];

    // Verify both INFO and ANN fields are present
    assert_eq!(record.info_fields.get("INFO_DP"), Some(&"50".to_string()));
    assert_eq!(record.info_fields.get("INFO_AF"), Some(&"0.25".to_string()));
    assert_eq!(record.info_fields.get("INFO_AC"), Some(&"2".to_string()));
    assert_eq!(record.info_fields.get("INFO_AN"), Some(&"4".to_string()));

    assert_eq!(record.info_fields.get("ANN_Allele"), Some(&"G".to_string()));
    assert_eq!(
        record.info_fields.get("ANN_Annotation"),
        Some(&"missense_variant".to_string())
    );
    assert_eq!(
        record.info_fields.get("ANN_Annotation_Impact"),
        Some(&"MODERATE".to_string())
    );
    assert_eq!(
        record.info_fields.get("ANN_Gene_Name"),
        Some(&"BRCA1".to_string())
    );
}

#[test]
fn test_real_world_snpeff_header_extraction() {
    use vcf_reformatter::get_info_from_header::extract_ann_format_from_header;

    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">"#;

    let ann_format = extract_ann_format_from_header(header);
    assert!(ann_format.is_some());

    let fields = ann_format.unwrap();
    assert_eq!(fields.len(), 16);
    assert_eq!(fields[0], "Allele");
    assert_eq!(fields[1], "Annotation");
    assert_eq!(fields[2], "Annotation_Impact");
    assert_eq!(fields[3], "Gene_Name");
    assert_eq!(fields[9], "HGVS.c");
    assert_eq!(fields[10], "HGVS.p");
    assert_eq!(fields[11], "cDNA.pos / cDNA.length");
    assert_eq!(fields[15], "ERRORS / WARNINGS / INFO");
}

#[test]
fn test_full_pipeline_with_snpeff_data() {
    use vcf_reformatter::reformat_vcf::{reformat_vcf_data_with_header, TranscriptHandling};

    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">"#;

    let column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    let data_lines = vec![
        "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50;ANN=G|missense_variant|MODERATE|BRCA1|ENSG00000012048|transcript|ENST00000357654|protein_coding|5/24|c.181T>C|p.Cys61Arg|181/5592|181/4863|61/1620||".to_string(),
        "chr2\t200\t.\tC\tT\t40\tPASS\tDP=30;ANN=T|stop_gained|HIGH|TP53|ENSG00000141510|transcript|ENST00000269305|protein_coding|7/11|c.916C>T|p.Arg306Ter|916/1182|916/1182|306/393||".to_string()
    ];

    let result = reformat_vcf_data_with_header(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::FirstOnly,
    );

    assert!(result.is_ok());

    let (headers, records) = result.unwrap();
    assert_eq!(records.len(), 2);

    // Check that headers include both INFO and ANN fields
    assert!(headers.contains(&"INFO_DP".to_string()));
    assert!(headers.contains(&"ANN_Allele".to_string()));
    assert!(headers.contains(&"ANN_Annotation".to_string()));
    assert!(headers.contains(&"ANN_Gene_Name".to_string()));
    assert!(headers.contains(&"ANN_HGVS_c".to_string()));

    // Verify first record (BRCA1 missense)
    assert_eq!(records[0].chromosome, "chr1");
    assert_eq!(
        records[0].info_fields.get("ANN_Gene_Name"),
        Some(&"BRCA1".to_string())
    );
    assert_eq!(
        records[0].info_fields.get("ANN_Annotation"),
        Some(&"missense_variant".to_string())
    );

    // Verify second record (TP53 stop_gained)
    assert_eq!(records[1].chromosome, "chr2");
    assert_eq!(
        records[1].info_fields.get("ANN_Gene_Name"),
        Some(&"TP53".to_string())
    );
    assert_eq!(
        records[1].info_fields.get("ANN_Annotation"),
        Some(&"stop_gained".to_string())
    );
}

#[test]
fn test_empty_ann_field_handling() {
    use vcf_reformatter::reformat_vcf::{ReformattedVcfRecord, TranscriptHandling};

    let vcf_line = "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50;ANN=;AF=0.25";

    let column_names = vec!["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];

    let ann_field_names = Some(vec!["Allele".to_string(), "Annotation".to_string()]);

    let records = ReformattedVcfRecord::from_vcf_line(
        vcf_line,
        &column_names,
        &None,
        &ann_field_names,
        TranscriptHandling::FirstOnly,
    )
    .unwrap();

    assert_eq!(records.len(), 1);
    let record = &records[0];

    // Should still have INFO fields
    assert_eq!(record.info_fields.get("INFO_DP"), Some(&"50".to_string()));
    assert_eq!(record.info_fields.get("INFO_AF"), Some(&"0.25".to_string()));

    // Should not have ANN fields since ANN was empty
    assert!(!record.info_fields.contains_key("ANN_Allele"));
    assert!(!record.info_fields.contains_key("ANN_Annotation"));
}
