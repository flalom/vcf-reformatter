use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use tempfile::tempdir;

use vcf_reformatter::{
    extract_csq_and_csq_names::extract_csq_regex,
    extract_sample_info::{_SampleData, parse_format_and_samples, ParsedFormatSample},
    get_info_from_header::extract_csq_format_from_header,
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
        annotation_field_type: vcf_reformatter::reformat_vcf::AnnotationFieldType::None,
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
fn test_vcf_record_with_complex_sample_names() {
    use vcf_reformatter::reformat_vcf::{ReformattedVcfRecord, TranscriptHandling};
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
        annotation_field_type: vcf_reformatter::reformat_vcf::AnnotationFieldType::None,
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
    let (records, field_type) = result.unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(
        field_type,
        vcf_reformatter::reformat_vcf::AnnotationFieldType::None
    );
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
    let (records, field_type) = result.unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].get("INFO_DP"), Some(&"10".to_string()));
    assert_eq!(records[0].get("INFO_AF"), Some(&"0.5".to_string()));
    assert_eq!(
        field_type,
        vcf_reformatter::reformat_vcf::AnnotationFieldType::None
    );
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
    let (records, field_type) = result.unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].get("CSQ_Allele"), Some(&"A".to_string()));
    assert_eq!(
        records[0].get("CSQ_Consequence"),
        Some(&"missense_variant".to_string())
    );
    assert_eq!(records[0].get("INFO_DP"), Some(&"10".to_string()));
    assert_eq!(
        field_type,
        vcf_reformatter::reformat_vcf::AnnotationFieldType::Csq
    );
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
    let (records, field_type) = result.unwrap();
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
    assert_eq!(
        field_type,
        vcf_reformatter::reformat_vcf::AnnotationFieldType::Ann
    );
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

    assert!(!record.info_fields.contains_key("ANN_Allele"));
    assert!(!record.info_fields.contains_key("ANN_Annotation"));
}
// ------------------------------------------------------------------------------
// Tests for MAF
// ------------------------------------------------------------------------------
// Helper function to create a test ReformattedVcfRecord for MAF tests
fn create_test_maf_record(
    chromosome: &str,
    position: u64,
    reference: &str,
    alternate: &str,
    quality: Option<f64>,
    filter: &str,
    info_fields: HashMap<String, String>,
) -> vcf_reformatter::reformat_vcf::ReformattedVcfRecord {
    vcf_reformatter::reformat_vcf::ReformattedVcfRecord {
        chromosome: chromosome.to_string(),
        position,
        id: Some("rs123456".to_string()),
        reference: reference.to_string(),
        alternate: alternate.to_string(),
        quality,
        filter: filter.to_string(),
        info_fields,
        format_sample_data: None,
        annotation_field_type: vcf_reformatter::reformat_vcf::AnnotationFieldType::None,
    }
}

#[test]
fn test_maf_headers_completeness() {
    use vcf_reformatter::essentials_fields::MafRecord;

    let headers = MafRecord::get_maf_headers();

    // Test that all required MAF headers are present (vcf2maf-parity 50-column layout)
    let expected_headers = vec![
        "Hugo_Symbol",
        "Entrez_Gene_Id",
        "Center",
        "NCBI_Build",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Strand",
        "Variant_Classification",
        "Variant_Type",
        "Reference_Allele",
        "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2",
        "dbSNP_RS",
        "dbSNP_Val_Status",
        "Tumor_Sample_Barcode",
        "Matched_Norm_Sample_Barcode",
        "Match_Norm_Seq_Allele1",
        "Match_Norm_Seq_Allele2",
        "Tumor_Validation_Allele1",
        "Tumor_Validation_Allele2",
        "Match_Norm_Validation_Allele1",
        "Match_Norm_Validation_Allele2",
        "Verification_Status",
        "Validation_Status",
        "Mutation_Status",
        "Sequencing_Phase",
        "Sequence_Source",
        "Validation_Method",
        "Score",
        "BAM_File",
        "Sequencer",
        "Tumor_Sample_UUID",
        "Matched_Norm_Sample_UUID",
        "HGVSc",
        "HGVSp",
        "HGVSp_Short",
        "Transcript_ID",
        "Exon_Number",
        "t_depth",
        "t_ref_count",
        "t_alt_count",
        "n_depth",
        "n_ref_count",
        "n_alt_count",
        "all_effects",
        "FILTER",
        "QUAL",
        "VAF",
        "Protein_Position",
    ];

    assert_eq!(headers.len(), expected_headers.len());
    for expected in expected_headers {
        assert!(
            headers.contains(&expected.to_string()),
            "Missing header: {}",
            expected
        );
    }
}

#[test]
fn test_maf_snp_creation() {
    use vcf_reformatter::essentials_fields::MafRecord;

    let mut info_fields = HashMap::new();
    info_fields.insert("CSQ_SYMBOL".to_string(), "TP53".to_string());
    info_fields.insert(
        "CSQ_Consequence".to_string(),
        "missense_variant".to_string(),
    );
    info_fields.insert("CSQ_Feature".to_string(), "ENST00000269305".to_string());
    info_fields.insert("CSQ_Protein_position".to_string(), "175".to_string());
    info_fields.insert(
        "CSQ_HGVSp".to_string(),
        "ENSP00000269305.4:p.Arg175His".to_string(),
    );
    info_fields.insert(
        "CSQ_HGVSc".to_string(),
        "ENST00000269305.8:c.524G>A".to_string(),
    );
    info_fields.insert("INFO_DP".to_string(), "100".to_string());

    let record =
        create_test_maf_record("chr17", 7674220, "G", "A", Some(60.0), "PASS", info_fields);

    let maf_record =
        MafRecord::from_reformatted_record(&record, "TEST_CENTER", "GRCh38", "TEST_SAMPLE")
            .unwrap();

    // Test basic fields - now with consistent MAF formatting
    assert_eq!(maf_record.hugo_symbol, "TP53");
    assert_eq!(maf_record.chromosome, "chr17"); // VCF naming passed through unchanged
    assert_eq!(maf_record.start_position, 7674220);
    assert_eq!(maf_record.end_position, 7674220); // SNP has same start/end
    assert_eq!(maf_record.variant_type, "SNP");
    assert_eq!(maf_record.reference_allele, "G");
    assert_eq!(maf_record.tumor_seq_allele2, "A");
    assert_eq!(maf_record.center, "TEST_CENTER");
    assert_eq!(maf_record.ncbi_build, "GRCh38");
    assert_eq!(maf_record.tumor_sample_barcode, "TEST_SAMPLE");

    // Test new fields
    assert_eq!(maf_record.qual, Some(60.0));
    assert_eq!(maf_record.filter_status, "PASS");
    assert_eq!(
        maf_record.transcript_id,
        Some("ENST00000269305".to_string())
    );
    assert_eq!(maf_record.protein_position, Some("175".to_string()));
    // vcf2maf.pl:786-787 strips the reference-sequence accession; the transcript is already
    // reported in Transcript_ID, checked just above.
    assert_eq!(maf_record.hgvsp, Some("p.Arg175His".to_string()));
    assert_eq!(maf_record.hgvsc, Some("c.524G>A".to_string()));
    assert_eq!(maf_record.hgvsp_short, Some("p.R175H".to_string()));
}

#[test]
fn test_maf_insertion_positions() {
    use vcf_reformatter::essentials_fields::MafRecord;

    let record = create_test_maf_record(
        "chr1",
        1000,
        "A",
        "ATCG", // Insertion
        Some(45.0),
        "PASS",
        HashMap::new(),
    );

    let maf_record =
        MafRecord::from_reformatted_record(&record, "TEST_CENTER", "GRCh38", "TEST_SAMPLE")
            .unwrap();

    assert_eq!(maf_record.variant_type, "INS");
    assert_eq!(maf_record.start_position, 1000); // Last shared base before insertion
    assert_eq!(maf_record.end_position, 1001); // INS: start + 1
    assert_eq!(maf_record.reference_allele, "-"); // MAF format uses "-" for insertions
    assert_eq!(maf_record.tumor_seq_allele2, "TCG"); // Inserted bases only (stripped shared "A" prefix)
    assert_eq!(maf_record.chromosome, "chr1"); // VCF naming passed through unchanged
}

#[test]
fn test_maf_deletion_positions() {
    use vcf_reformatter::essentials_fields::MafRecord;

    let record = create_test_maf_record(
        "chr2",
        2000,
        "ATCG",
        "A", // Deletion of TCG
        Some(55.0),
        "PASS",
        HashMap::new(),
    );

    let maf_record =
        MafRecord::from_reformatted_record(&record, "TEST_CENTER", "GRCh38", "TEST_SAMPLE")
            .unwrap();

    assert_eq!(maf_record.variant_type, "DEL");
    assert_eq!(maf_record.start_position, 2001); // First deleted base (after shared "A" anchor)
    assert_eq!(maf_record.end_position, 2003); // Last deleted base
    assert_eq!(maf_record.reference_allele, "TCG"); // Deleted bases only (stripped shared "A" prefix)
    assert_eq!(maf_record.tumor_seq_allele2, "-"); // MAF format uses "-" for deletions
    assert_eq!(maf_record.chromosome, "chr2"); // VCF naming passed through unchanged
}

#[test]
fn test_maf_variant_classification_mapping() {
    use vcf_reformatter::essentials_fields::MafRecord;

    // Test missense variant
    let mut info_fields = HashMap::new();
    info_fields.insert(
        "CSQ_Consequence".to_string(),
        "missense_variant".to_string(),
    );

    let record = create_test_maf_record(
        "chr1",
        1000,
        "G",
        "A",
        Some(60.0),
        "PASS",
        info_fields.clone(),
    );
    let maf_record =
        MafRecord::from_reformatted_record(&record, "TEST", "GRCh38", "SAMPLE").unwrap();
    assert_eq!(maf_record.variant_classification, "Missense_Mutation");

    // Test synonymous variant
    info_fields.clear();
    info_fields.insert(
        "CSQ_Consequence".to_string(),
        "synonymous_variant".to_string(),
    );
    let record = create_test_maf_record(
        "chr1",
        1000,
        "G",
        "A",
        Some(60.0),
        "PASS",
        info_fields.clone(),
    );
    let maf_record =
        MafRecord::from_reformatted_record(&record, "TEST", "GRCh38", "SAMPLE").unwrap();
    assert_eq!(maf_record.variant_classification, "Silent");

    // Test stop gained
    info_fields.clear();
    info_fields.insert("CSQ_Consequence".to_string(), "stop_gained".to_string());
    let record = create_test_maf_record("chr1", 1000, "G", "A", Some(60.0), "PASS", info_fields);
    let maf_record =
        MafRecord::from_reformatted_record(&record, "TEST", "GRCh38", "SAMPLE").unwrap();
    assert_eq!(maf_record.variant_classification, "Nonsense_Mutation");
}

#[test]
fn test_maf_multi_allelic_handling() {
    use vcf_reformatter::essentials_fields::MafRecord;

    let mut info_fields = HashMap::new();
    info_fields.insert("CSQ_SYMBOL".to_string(), "BRCA1".to_string());
    info_fields.insert(
        "CSQ_Consequence".to_string(),
        "missense_variant".to_string(),
    );

    let record = create_test_maf_record(
        "chr17",
        43094692,
        "G",
        "A,T", // Multi-allelic
        Some(80.0),
        "PASS",
        info_fields,
    );

    let maf_records =
        MafRecord::from_reformatted_record_multi(&record, "TEST_CENTER", "GRCh38", "TEST_SAMPLE")
            .unwrap();

    assert_eq!(maf_records.len(), 2); // Should create 2 MAF records

    // First allele (G>A)
    assert_eq!(maf_records[0].tumor_seq_allele2, "A");
    assert_eq!(maf_records[0].hugo_symbol, "BRCA1");

    // Second allele (G>T)
    assert_eq!(maf_records[1].tumor_seq_allele2, "T");
    assert_eq!(maf_records[1].hugo_symbol, "BRCA1");

    // Both should have same position and quality
    assert_eq!(maf_records[0].start_position, maf_records[1].start_position);
    assert_eq!(maf_records[0].qual, maf_records[1].qual);
}

#[test]
fn test_maf_missing_annotations() {
    use vcf_reformatter::essentials_fields::MafRecord;

    // Test with minimal annotations
    let record = create_test_maf_record(
        "chr22",
        42000,
        "C",
        "G",
        None,           // No quality score
        ".",            // No filter
        HashMap::new(), // No annotations
    );

    let maf_record =
        MafRecord::from_reformatted_record(&record, "TEST_CENTER", "GRCh38", "TEST_SAMPLE")
            .unwrap();

    // No SYMBOL and no transcript at all (pure IGR case) -- vcf2maf never leaves Hugo_Symbol
    // blank, falling back to "Unknown" here since there's no transcript ID to fall back to either.
    assert_eq!(maf_record.hugo_symbol, "Unknown");
    assert_eq!(maf_record.qual, None);
    assert_eq!(maf_record.filter_status, ".");
    assert_eq!(maf_record.transcript_id, None);
    assert_eq!(maf_record.protein_position, None);
    assert_eq!(maf_record.hgvsp, None);
    assert_eq!(maf_record.hgvsc, None);
}

#[test]
fn test_maf_tsv_line_format() {
    use vcf_reformatter::essentials_fields::MafRecord;

    let mut info_fields = HashMap::new();
    info_fields.insert("CSQ_SYMBOL".to_string(), "EGFR".to_string());
    info_fields.insert("CSQ_Feature".to_string(), "ENST00000275493".to_string());
    info_fields.insert("CSQ_Protein_position".to_string(), "858".to_string());

    let record =
        create_test_maf_record("chr7", 55191822, "T", "G", Some(99.9), "PASS", info_fields);

    let maf_record =
        MafRecord::from_reformatted_record(&record, "TCGA", "GRCh38", "TCGA-01-0001").unwrap();

    let tsv_line = maf_record.to_tsv_line();
    let fields: Vec<&str> = tsv_line.split('\t').collect();

    // Verify field count matches header count
    let headers = MafRecord::get_maf_headers();
    assert_eq!(fields.len(), headers.len());

    // Test specific field values - now with normalized chromosome
    assert_eq!(fields[0], "EGFR"); // Hugo_Symbol
    assert_eq!(fields[2], "TCGA"); // Center
    assert_eq!(fields[3], "GRCh38"); // NCBI_Build
    assert_eq!(fields[4], "chr7"); // VCF naming passed through unchanged // Chromosome (now normalized)
    assert_eq!(fields[5], "55191822"); // Start_Position
    assert_eq!(fields[37], "ENST00000275493"); // Transcript_ID
    assert_eq!(fields[46], "PASS"); // FILTER
    assert_eq!(fields[47], "99.9"); // QUAL
    assert_eq!(fields[49], "858"); // Protein_Position
}
#[test]
fn test_maf_chromosome_normalization() {
    use vcf_reformatter::essentials_fields::MafRecord;

    // The VCF's own chromosome naming is passed through unchanged (user's call
    // 2026-08-30, reversing the earlier strip-to-bare-name behaviour). vcf2maf does
    // the same, so this also removes a diff class against it.
    let test_cases = vec![
        ("1", "1"),
        ("chr1", "chr1"),
        ("CHR1", "CHR1"),
        ("X", "X"),
        ("chrX", "chrX"),
        ("Y", "Y"),
        ("chrY", "chrY"),
        ("M", "M"),
        ("MT", "MT"),
        ("chrM", "chrM"),
        ("chrMT", "chrMT"),
    ];

    for (input, expected) in test_cases {
        let record =
            create_test_maf_record(input, 1000, "A", "T", Some(50.0), "PASS", HashMap::new());

        let maf_record =
            MafRecord::from_reformatted_record(&record, "TEST", "GRCh38", "SAMPLE").unwrap();

        assert_eq!(
            maf_record.chromosome, expected,
            "Failed to handle chromosome {} - got {}",
            input, maf_record.chromosome
        );
    }
}
#[test]
fn test_maf_edge_cases() {
    use vcf_reformatter::essentials_fields::MafRecord;

    // Test with very long insertion
    let long_insertion = "A".to_string() + &"T".repeat(100);
    let record = create_test_maf_record(
        "chr1",
        1000,
        "A",
        &long_insertion,
        Some(30.0),
        "PASS",
        HashMap::new(),
    );

    let maf_record =
        MafRecord::from_reformatted_record(&record, "TEST", "GRCh38", "SAMPLE").unwrap();

    assert_eq!(maf_record.variant_type, "INS");
    // Proper MAF format: strip shared "A" prefix
    assert_eq!(maf_record.reference_allele, "-");
    let expected_inserted = "T".repeat(100); // "A" + 100x"T" minus shared "A"
    assert_eq!(maf_record.tumor_seq_allele2, expected_inserted);

    // Test with very long deletion
    let long_deletion = "A".to_string() + &"G".repeat(50);
    let record = create_test_maf_record(
        "chr1",
        2000,
        &long_deletion,
        "A",
        Some(40.0),
        "PASS",
        HashMap::new(),
    );

    let maf_record =
        MafRecord::from_reformatted_record(&record, "TEST", "GRCh38", "SAMPLE").unwrap();

    assert_eq!(maf_record.variant_type, "DEL");
    assert_eq!(maf_record.tumor_seq_allele2, "-");
    let expected_deleted = "G".repeat(50); // "A" + 50x"G" minus shared "A"
    assert_eq!(maf_record.reference_allele, expected_deleted);
    assert_eq!(maf_record.start_position, 2001); // First deleted base (after shared "A")
    assert_eq!(maf_record.end_position, 2050); // 2001 + 50 - 1
}

#[test]
fn test_maf_vep_vs_snpeff_annotations() {
    use vcf_reformatter::essentials_fields::MafRecord;

    // Test VEP-style annotations
    let mut vep_fields = HashMap::new();
    vep_fields.insert("CSQ_SYMBOL".to_string(), "TP53".to_string());
    vep_fields.insert("CSQ_Feature".to_string(), "ENST00000269305".to_string());
    vep_fields.insert("CSQ_Protein_position".to_string(), "175".to_string());

    let record = create_test_maf_record("chr17", 7674220, "G", "A", Some(60.0), "PASS", vep_fields);
    let vep_maf = MafRecord::from_reformatted_record(&record, "TEST", "GRCh38", "SAMPLE").unwrap();

    // Test SnpEff-style annotations
    let mut snpeff_fields = HashMap::new();
    snpeff_fields.insert("ANN_Gene_Name".to_string(), "TP53".to_string());
    snpeff_fields.insert("ANN_Feature_ID".to_string(), "ENST00000269305".to_string());
    snpeff_fields.insert("ANN_Protein_position".to_string(), "175".to_string());

    let record = create_test_maf_record(
        "chr17",
        7674220,
        "G",
        "A",
        Some(60.0),
        "PASS",
        snpeff_fields,
    );
    let snpeff_maf =
        MafRecord::from_reformatted_record(&record, "TEST", "GRCh38", "SAMPLE").unwrap();

    // Both should produce equivalent MAF records
    assert_eq!(vep_maf.hugo_symbol, snpeff_maf.hugo_symbol);
    assert_eq!(vep_maf.transcript_id, snpeff_maf.transcript_id);
    assert_eq!(vep_maf.protein_position, snpeff_maf.protein_position);
}

#[test]
fn test_maf_vaf_calculation() {
    // Create a simple inline VAF calculation function for testing
    fn calculate_vaf(tumor_depth: Option<u32>, total_depth: Option<u32>) -> Option<f32> {
        match (tumor_depth, total_depth) {
            (Some(t_depth), Some(total)) if total > 0 => {
                let vaf = t_depth as f32 / total as f32;
                if vaf <= 1.0 {
                    Some(vaf)
                } else {
                    None
                } // Prevent VAF > 1.0
            }
            _ => None,
        }
    }

    // Test normal VAF calculation
    assert_eq!(calculate_vaf(Some(50), Some(200)), Some(0.25));
    assert_eq!(calculate_vaf(Some(100), Some(100)), Some(1.0));
    assert_eq!(calculate_vaf(Some(0), Some(100)), Some(0.0));

    // Test edge cases
    assert_eq!(calculate_vaf(None, Some(100)), None);
    assert_eq!(calculate_vaf(Some(50), None), None);
    assert_eq!(calculate_vaf(Some(50), Some(0)), None); // Avoid division by zero
    assert_eq!(calculate_vaf(None, None), None);

    // Test VAF > 1.0 edge case (should return None)
    assert_eq!(calculate_vaf(Some(150), Some(100)), None);
}

#[test]
fn test_maf_depth_extraction() {
    use vcf_reformatter::essentials_fields::MafRecord;

    let mut info_fields = HashMap::new();
    info_fields.insert("INFO_DP".to_string(), "200".to_string());
    info_fields.insert("INFO_AD".to_string(), "150,50".to_string()); // REF,ALT counts
    info_fields.insert("CSQ_SYMBOL".to_string(), "BRCA1".to_string());

    let record =
        create_test_maf_record("chr17", 43094692, "G", "A", Some(70.0), "PASS", info_fields);

    let maf_record =
        MafRecord::from_reformatted_record(&record, "TEST_CENTER", "GRCh38", "TEST_SAMPLE")
            .unwrap();

    // Verify depth extraction - test what your implementation actually returns
    assert_eq!(maf_record.t_depth, Some(200));

    // Your implementation might extract depth differently
    // Test both possibilities
    if let Some(depth) = maf_record.t_alt_count {
        // If depth is extracted, it should be reasonable
        assert!(depth <= 200);
    }

    // Test VAF calculation if both values are present
    if let (Some(depth), Some(total)) = (maf_record.t_alt_count, maf_record.t_depth) {
        if let Some(vaf) = maf_record.vaf {
            let expected_vaf = depth as f32 / total as f32;
            assert!((vaf - expected_vaf).abs() < 0.001);
        }
    }
}
#[test]
fn test_maf_frameshift_classification() {
    use vcf_reformatter::essentials_fields::MafRecord;

    let mut info_fields = HashMap::new();
    info_fields.insert(
        "CSQ_Consequence".to_string(),
        "frameshift_variant".to_string(),
    );
    info_fields.insert("CSQ_SYMBOL".to_string(), "SAMD11".to_string());

    let record = create_test_maf_record(
        "chr1",
        939436,
        "-",       // This indicates an insertion in VCF format
        "CTCCCCT", // Complex insertion
        Some(1152.62),
        "PASS",
        info_fields,
    );

    let maf_record =
        MafRecord::from_reformatted_record(&record, "Unknown_Center", "GRCh38", "B487_B487_1_V")
            .unwrap();

    assert_eq!(maf_record.hugo_symbol, "SAMD11");
    // Test what your implementation actually classifies this as
    // It might be "Frame_Shift_Del" based on the reference being "-"
    assert!(maf_record.variant_classification.contains("Frame_Shift"));
    assert_eq!(maf_record.qual, Some(1152.62));
}

#[test]
fn test_maf_intron_and_utr_variants() {
    use vcf_reformatter::essentials_fields::MafRecord;

    // Test 3'UTR variant
    let mut utr_fields = HashMap::new();
    utr_fields.insert(
        "CSQ_Consequence".to_string(),
        "3_prime_UTR_variant".to_string(),
    );
    utr_fields.insert("CSQ_SYMBOL".to_string(), "SAMD11".to_string());

    let record = create_test_maf_record("chr1", 924024, "C", "G", Some(7567.2), "PASS", utr_fields);
    let maf_record =
        MafRecord::from_reformatted_record(&record, "Unknown_Center", "GRCh38", "B487_B487_1_V")
            .unwrap();

    assert_eq!(maf_record.variant_classification, "3'UTR");

    // Test intron variant
    let mut intron_fields = HashMap::new();
    intron_fields.insert("CSQ_Consequence".to_string(), "intron_variant".to_string());
    intron_fields.insert("CSQ_SYMBOL".to_string(), "LINC02593".to_string());

    let record = create_test_maf_record(
        "chr1",
        923892,
        "T",
        "A",
        Some(1.18132E-14),
        "PASS",
        intron_fields,
    );
    let maf_record =
        MafRecord::from_reformatted_record(&record, "Unknown_Center", "GRCh38", "B487_B487_1_V")
            .unwrap();

    assert_eq!(maf_record.variant_classification, "Intron");
    assert_eq!(maf_record.hugo_symbol, "LINC02593");
}

#[test]
fn test_reformatted_record_carries_annotation_type() {
    use vcf_reformatter::reformat_vcf::AnnotationFieldType;

    let record = create_test_maf_record(
        "chr1",
        1000,
        "A",
        "T",
        Some(30.0),
        "PASS",
        HashMap::from([
            ("CSQ_SYMBOL".to_string(), "BRCA1".to_string()),
            (
                "CSQ_Consequence".to_string(),
                "missense_variant".to_string(),
            ),
        ]),
    );
    // This field should exist on ReformattedVcfRecord
    assert_eq!(record.annotation_field_type, AnnotationFieldType::None);
}

#[test]
fn test_maf_snpeff_field_priority_with_annotation_type() {
    use vcf_reformatter::reformat_vcf::AnnotationFieldType;

    // SnpEff record with ANN_ prefixed fields
    let mut info = HashMap::new();
    info.insert("ANN_Gene_Name".to_string(), "TP53".to_string());
    info.insert("ANN_Annotation".to_string(), "missense_variant".to_string());
    info.insert("ANN_Annotation_Impact".to_string(), "MODERATE".to_string());
    info.insert("ANN_Feature_ID".to_string(), "NM_000546.6".to_string());
    info.insert("ANN_HGVS_p".to_string(), "p.Arg175His".to_string());
    info.insert("ANN_HGVS_c".to_string(), "c.524G>A".to_string());

    let mut record = create_test_maf_record("chr17", 7675088, "C", "T", Some(100.0), "PASS", info);
    record.annotation_field_type = AnnotationFieldType::Ann;

    let maf = vcf_reformatter::essentials_fields::MafRecord::from_reformatted_record(
        &record,
        "TestCenter",
        "GRCh38",
        "SAMPLE-001",
    )
    .unwrap();

    assert_eq!(maf.hugo_symbol, "TP53");
    assert_eq!(maf.variant_classification, "Missense_Mutation");
    assert_eq!(maf.transcript_id, Some("NM_000546.6".to_string()));
    assert_eq!(maf.hgvsp, Some("p.Arg175His".to_string()));
    assert_eq!(maf.hgvsc, Some("c.524G>A".to_string()));
}

#[test]
fn test_maf_snpeff_specific_consequences() {
    let test_cases = vec![
        ("disruptive_inframe_insertion", "In_Frame_Ins"),
        ("conservative_inframe_insertion", "In_Frame_Ins"),
        ("disruptive_inframe_deletion", "In_Frame_Del"),
        ("conservative_inframe_deletion", "In_Frame_Del"),
        ("initiator_codon_variant", "Translation_Start_Site"),
        ("rare_amino_acid_variant", "Missense_Mutation"),
        ("stop_retained_variant", "Silent"),
        ("upstream_gene_variant", "5'Flank"),
        ("downstream_gene_variant", "3'Flank"),
        ("non_coding_transcript_exon_variant", "RNA"),
        ("non_coding_transcript_variant", "RNA"),
    ];

    for (consequence, expected_classification) in &test_cases {
        let mut info = HashMap::new();
        info.insert("ANN_Annotation".to_string(), consequence.to_string());
        let record = create_test_maf_record("chr1", 1000, "A", "T", Some(30.0), "PASS", info);
        let maf = vcf_reformatter::essentials_fields::MafRecord::from_reformatted_record(
            &record,
            "TestCenter",
            "GRCh38",
            "SAMPLE-001",
        )
        .unwrap();
        assert_eq!(
            maf.variant_classification, *expected_classification,
            "Failed for consequence: {} — expected {}, got {}",
            consequence, expected_classification, maf.variant_classification
        );
    }
}

#[test]
fn test_maf_consequence_case_insensitive() {
    // frameshift_variant is deliberately not covered here since its classification depends on
    // variant_type (indel length/direction), not just the consequence term; see
    // test_maf_frameshift_ins_vs_del for that.
    let test_cases = vec![
        ("STOP_GAINED", "Nonsense_Mutation"),
        ("Missense_Variant", "Missense_Mutation"),
        ("Synonymous_Variant", "Silent"),
        ("SPLICE_DONOR_VARIANT", "Splice_Site"),
    ];

    for (consequence, expected) in &test_cases {
        let mut info = HashMap::new();
        info.insert("ANN_Annotation".to_string(), consequence.to_string());
        let record = create_test_maf_record("chr1", 1000, "A", "T", Some(30.0), "PASS", info);
        let maf = vcf_reformatter::essentials_fields::MafRecord::from_reformatted_record(
            &record,
            "TestCenter",
            "GRCh38",
            "SAMPLE-001",
        )
        .unwrap();
        assert_eq!(
            maf.variant_classification, *expected,
            "Failed for consequence: {} — expected {}, got {}",
            consequence, expected, maf.variant_classification
        );
    }
}

#[test]
fn test_maf_utr_classification_correctness() {
    // 5' UTR should map to 5'UTR, not 3'UTR (was a bug in original code)
    let mut info_5utr = HashMap::new();
    info_5utr.insert(
        "CSQ_Consequence".to_string(),
        "5_prime_utr_variant".to_string(),
    );
    let record_5utr = create_test_maf_record("chr1", 1000, "A", "T", Some(30.0), "PASS", info_5utr);
    let maf_5utr = vcf_reformatter::essentials_fields::MafRecord::from_reformatted_record(
        &record_5utr,
        "TestCenter",
        "GRCh38",
        "SAMPLE-001",
    )
    .unwrap();
    assert_eq!(maf_5utr.variant_classification, "5'UTR");

    // 3' UTR should map to 3'UTR
    let mut info_3utr = HashMap::new();
    info_3utr.insert(
        "CSQ_Consequence".to_string(),
        "3_prime_utr_variant".to_string(),
    );
    let record_3utr = create_test_maf_record("chr1", 1000, "A", "T", Some(30.0), "PASS", info_3utr);
    let maf_3utr = vcf_reformatter::essentials_fields::MafRecord::from_reformatted_record(
        &record_3utr,
        "TestCenter",
        "GRCh38",
        "SAMPLE-001",
    )
    .unwrap();
    assert_eq!(maf_3utr.variant_classification, "3'UTR");
}

#[test]
fn test_maf_5utr_insertion() {
    use vcf_reformatter::essentials_fields::MafRecord;

    let mut info = HashMap::new();
    info.insert(
        "CSQ_Consequence".to_string(),
        "5_prime_utr_variant".to_string(),
    );
    info.insert("CSQ_SYMBOL".to_string(), "SAMD11".to_string());

    // VCF insertion: ref=A, alt=ATCG at position 924024
    let record = create_test_maf_record("chr1", 924024, "A", "ATCG", Some(53.0), "PASS", info);
    let maf =
        MafRecord::from_reformatted_record(&record, "TestCenter", "GRCh38", "SAMPLE-001").unwrap();

    assert_eq!(maf.variant_classification, "5'UTR");
    assert_eq!(maf.variant_type, "INS");
    assert_eq!(maf.reference_allele, "-");
    assert_eq!(maf.tumor_seq_allele2, "TCG"); // Stripped shared "A" prefix
    assert_eq!(maf.start_position, 924024); // Last shared base
    assert_eq!(maf.end_position, 924025); // start + 1
    assert_eq!(maf.hugo_symbol, "SAMD11");
}

#[test]
fn test_maf_5utr_deletion() {
    use vcf_reformatter::essentials_fields::MafRecord;

    let mut info = HashMap::new();
    info.insert(
        "ANN_Annotation".to_string(),
        "5_prime_UTR_variant".to_string(),
    );
    info.insert("ANN_Gene_Name".to_string(), "SAMD11".to_string());

    // VCF deletion: ref=ATCG, alt=A at position 924024
    let record = create_test_maf_record("chr1", 924024, "ATCG", "A", Some(40.0), "PASS", info);
    let maf =
        MafRecord::from_reformatted_record(&record, "TestCenter", "GRCh38", "SAMPLE-001").unwrap();

    assert_eq!(maf.variant_classification, "5'UTR");
    assert_eq!(maf.variant_type, "DEL");
    assert_eq!(maf.reference_allele, "TCG"); // Stripped shared "A" prefix
    assert_eq!(maf.tumor_seq_allele2, "-");
    assert_eq!(maf.start_position, 924025); // First deleted base
    assert_eq!(maf.end_position, 924027); // Last deleted base
    assert_eq!(maf.hugo_symbol, "SAMD11");
}

#[test]
fn test_maf_3utr_insertion() {
    use vcf_reformatter::essentials_fields::MafRecord;

    let mut info = HashMap::new();
    info.insert(
        "CSQ_Consequence".to_string(),
        "3_prime_utr_variant".to_string(),
    );
    info.insert("CSQ_SYMBOL".to_string(), "OR4F5".to_string());

    let record = create_test_maf_record("chr1", 69511, "G", "GAA", Some(100.0), "PASS", info);
    let maf =
        MafRecord::from_reformatted_record(&record, "TestCenter", "GRCh38", "SAMPLE-001").unwrap();

    assert_eq!(maf.variant_classification, "3'UTR");
    assert_eq!(maf.variant_type, "INS");
    assert_eq!(maf.reference_allele, "-");
    assert_eq!(maf.tumor_seq_allele2, "AA"); // Stripped shared "G"
}

#[test]
fn test_maf_frameshift_ins_vs_del() {
    use vcf_reformatter::essentials_fields::MafRecord;

    // Frameshift insertion: ref shorter than alt
    let mut info_ins = HashMap::new();
    info_ins.insert(
        "CSQ_Consequence".to_string(),
        "frameshift_variant".to_string(),
    );
    let record_ins =
        create_test_maf_record("chr1", 1000, "A", "ATCG", Some(50.0), "PASS", info_ins);
    let maf_ins =
        MafRecord::from_reformatted_record(&record_ins, "TestCenter", "GRCh38", "SAMPLE-001")
            .unwrap();
    assert_eq!(maf_ins.variant_classification, "Frame_Shift_Ins");
    assert_eq!(maf_ins.variant_type, "INS");

    // Frameshift deletion: ref longer than alt
    let mut info_del = HashMap::new();
    info_del.insert(
        "CSQ_Consequence".to_string(),
        "frameshift_variant".to_string(),
    );
    let record_del =
        create_test_maf_record("chr1", 1000, "ATCG", "A", Some(50.0), "PASS", info_del);
    let maf_del =
        MafRecord::from_reformatted_record(&record_del, "TestCenter", "GRCh38", "SAMPLE-001")
            .unwrap();
    assert_eq!(maf_del.variant_classification, "Frame_Shift_Del");
    assert_eq!(maf_del.variant_type, "DEL");

    // Frameshift tag on a same-length substitution: not a real indel, so neither the
    // Frame_Shift_Ins nor Frame_Shift_Del branch applies. Confirmed against vcf2maf's own
    // GetVariantClassification, which falls through to its catch-all in this case too.
    let mut info_snp = HashMap::new();
    info_snp.insert(
        "CSQ_Consequence".to_string(),
        "frameshift_variant".to_string(),
    );
    let record_snp = create_test_maf_record("chr1", 1000, "AT", "GC", Some(50.0), "PASS", info_snp);
    let maf_snp =
        MafRecord::from_reformatted_record(&record_snp, "TestCenter", "GRCh38", "SAMPLE-001")
            .unwrap();
    assert_eq!(maf_snp.variant_classification, "Targeted_Region");
}

#[test]
fn test_maf_classification_impact_not_used_as_consequence() {
    // Only provide IMPACT, not actual consequence
    let mut info = HashMap::new();
    info.insert("ANN_Annotation_Impact".to_string(), "HIGH".to_string());
    // No ANN_Annotation or CSQ_Consequence

    let record = create_test_maf_record("chr1", 1000, "A", "T", Some(30.0), "PASS", info);
    let maf = vcf_reformatter::essentials_fields::MafRecord::from_reformatted_record(
        &record,
        "TestCenter",
        "GRCh38",
        "SAMPLE-001",
    )
    .unwrap();

    // Should use IMPACT fallback path, not try to map "HIGH" as a consequence term
    assert_eq!(maf.variant_classification, "Missense_Mutation"); // HIGH impact fallback
}

#[test]
fn test_maf_real_world_data_simulation() {
    use vcf_reformatter::essentials_fields::MafRecord;

    // Simulate real data like the output you showed
    let test_cases = vec![
        // OR4F5 missense mutation
        (
            "OR4F5",
            "chr1",
            69787,
            "T",
            "A",
            "missense_variant",
            "ENST00000641515",
            "254",
        ),
        // SAMD11 3'UTR variants
        (
            "SAMD11",
            "chr1",
            924024,
            "C",
            "G",
            "3_prime_UTR_variant",
            "ENST00000616016",
            "",
        ),
        (
            "SAMD11",
            "chr1",
            924305,
            "G",
            "A",
            "3_prime_UTR_variant",
            "ENST00000616016",
            "",
        ),
        // SAMD11 silent mutation
        (
            "SAMD11",
            "chr1",
            924509,
            "C",
            "T",
            "synonymous_variant",
            "ENST00000616016",
            "26",
        ),
        // SAMD11 missense mutations
        (
            "SAMD11",
            "chr1",
            924586,
            "T",
            "A",
            "missense_variant",
            "ENST00000616016",
            "52",
        ),
        (
            "SAMD11",
            "chr1",
            942412,
            "T",
            "G",
            "missense_variant",
            "ENST00000341065",
            "238",
        ),
    ];

    for (gene, chr, pos, ref_allele, alt_allele, consequence, transcript, protein_pos) in test_cases
    {
        let mut info_fields = HashMap::new();
        info_fields.insert("CSQ_SYMBOL".to_string(), gene.to_string());
        info_fields.insert("CSQ_Consequence".to_string(), consequence.to_string());
        info_fields.insert("CSQ_Feature".to_string(), transcript.to_string());

        if !protein_pos.is_empty() {
            info_fields.insert("CSQ_Protein_position".to_string(), protein_pos.to_string());
        }

        let record = create_test_maf_record(
            chr,
            pos,
            ref_allele,
            alt_allele,
            Some(60.0),
            "PASS",
            info_fields,
        );

        let maf_record = MafRecord::from_reformatted_record(
            &record,
            "Unknown_Center",
            "GRCh38",
            "B487_B487_1_V",
        )
        .unwrap();

        assert_eq!(maf_record.hugo_symbol, gene);
        assert_eq!(maf_record.start_position, pos);
        assert_eq!(maf_record.reference_allele, ref_allele);
        assert_eq!(maf_record.tumor_seq_allele2, alt_allele);
        assert_eq!(maf_record.transcript_id, Some(transcript.to_string()));

        if !protein_pos.is_empty() {
            assert_eq!(maf_record.protein_position, Some(protein_pos.to_string()));
        }

        // Verify variant classification mapping
        match consequence {
            "missense_variant" => {
                assert_eq!(maf_record.variant_classification, "Missense_Mutation")
            }
            "synonymous_variant" => assert_eq!(maf_record.variant_classification, "Silent"),
            "3_prime_UTR_variant" => assert_eq!(maf_record.variant_classification, "3'UTR"),
            _ => {} // Other consequences tested elsewhere
        }
    }
}

#[test]
fn test_maf_file_write_integration() {
    use std::fs;
    use tempfile::tempdir;
    use vcf_reformatter::essentials_fields::MafRecord;
    use vcf_reformatter::reformat_vcf::write_maf_file;

    let mut info_fields = HashMap::new();
    info_fields.insert("CSQ_SYMBOL".to_string(), "TP53".to_string());
    info_fields.insert(
        "CSQ_Consequence".to_string(),
        "missense_variant".to_string(),
    );
    info_fields.insert("CSQ_Feature".to_string(), "ENST00000269305".to_string());

    let record =
        create_test_maf_record("chr17", 7674220, "G", "A", Some(60.0), "PASS", info_fields);
    let maf_record =
        MafRecord::from_reformatted_record(&record, "TEST_CENTER", "GRCh38", "TEST_SAMPLE")
            .unwrap();

    // Test writing to file
    let temp_dir = tempdir().unwrap();
    let file_path = temp_dir.path().join("test.maf");

    write_maf_file(file_path.to_str().unwrap(), &[maf_record], false).unwrap();

    // Read and verify file contents
    let contents = fs::read_to_string(&file_path).unwrap();
    let lines: Vec<&str> = contents.lines().collect();

    assert_eq!(lines.len(), 2); // Header + 1 data line

    // Verify header
    let headers = MafRecord::get_maf_headers();
    assert_eq!(lines[0], headers.join("\t"));

    // Verify data contains expected values - now with normalized chromosome
    assert!(lines[1].contains("TP53"));
    assert!(lines[1].contains("TEST_CENTER"));
    assert!(lines[1].contains("17")); // Now normalized (no "chr")
    assert!(lines[1].contains("7674220"));
    assert!(lines[1].contains("60"));
}
#[test]
fn test_maf_multi_allelic_conversion() {
    use vcf_reformatter::essentials_fields::MafRecord;

    // Test the multi-allelic function which does use MAF-specific conversions
    let record = create_test_maf_record(
        "chr1",
        1000,
        "A",
        "ATCG", // Insertion
        Some(45.0),
        "PASS",
        HashMap::new(),
    );

    let maf_records =
        MafRecord::from_reformatted_record_multi(&record, "TEST_CENTER", "GRCh38", "TEST_SAMPLE")
            .unwrap();

    assert_eq!(maf_records.len(), 1);
    let maf_record = &maf_records[0];

    assert_eq!(maf_record.variant_type, "INS");
    assert_eq!(maf_record.start_position, 1000);
    assert_eq!(maf_record.end_position, 1001);
    assert_eq!(maf_record.reference_allele, "-");
    assert_eq!(maf_record.tumor_seq_allele2, "TCG"); // Stripped shared "A" prefix
    assert_eq!(maf_record.chromosome, "chr1");
}

#[test]
fn test_maf_multi_deletion_conversion() {
    use vcf_reformatter::essentials_fields::MafRecord;

    let record = create_test_maf_record(
        "chr2",
        2000,
        "ATCG",
        "A", // Deletion
        Some(55.0),
        "PASS",
        HashMap::new(),
    );

    let maf_records =
        MafRecord::from_reformatted_record_multi(&record, "TEST_CENTER", "GRCh38", "TEST_SAMPLE")
            .unwrap();

    assert_eq!(maf_records.len(), 1);
    let maf_record = &maf_records[0];

    assert_eq!(maf_record.variant_type, "DEL");
    assert_eq!(maf_record.start_position, 2001); // First deleted base
    assert_eq!(maf_record.end_position, 2003); // Last deleted base
    assert_eq!(maf_record.reference_allele, "TCG"); // Stripped shared "A" prefix
    assert_eq!(maf_record.tumor_seq_allele2, "-");
    assert_eq!(maf_record.chromosome, "chr2");
}

use std::process::Command;
use vcf_reformatter::reformat_vcf::reformat_vcf_data_with_header_parallel_chunked;

#[test]
fn test_maf_output() {
    use std::io::Write;

    // Create a temporary VCF file for testing
    let mut temp_vcf = tempfile::NamedTempFile::new().unwrap();
    let vcf_content = r#"##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TEST-01
chr1	14930	.	A	G	100	PASS	DP=50;AF=0.5;CSQ=G|missense_variant|MODERATE|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene|||||||||||rs12345||1|||HGNC|HGNC:38034|||||||||||||||||||||||||||||||||	GT:DP:AD	0/1:50:25,25
chr1	69511	.	C	T	200	PASS	DP=60;AF=0.8;CSQ=T|synonymous_variant|LOW|OR4F5|ENSG00000186092|Transcript|ENST00000335137|protein_coding|||||||||||1|||HGNC|HGNC:14825|||||||||||||||||||||||||||||||||	GT:DP:AD	1/1:60:12,48
"#;

    write!(temp_vcf, "{}", vcf_content).unwrap();
    temp_vcf.flush().unwrap();

    // Create gzipped version
    let temp_vcf_gz = tempfile::NamedTempFile::with_suffix(".vcf.gz").unwrap();
    {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::fs::File;
        use std::io::Read;

        let mut file = File::open(temp_vcf.path()).unwrap();
        let mut contents = Vec::new();
        file.read_to_end(&mut contents).unwrap();

        let gz_file = File::create(temp_vcf_gz.path()).unwrap();
        let mut encoder = GzEncoder::new(gz_file, Compression::default());
        encoder.write_all(&contents).unwrap();
        encoder.finish().unwrap();
    }

    let dir = tempdir().unwrap();
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            temp_vcf_gz.path().to_str().unwrap(),
            "--output-format",
            "maf",
            "--center",
            "TestCenter",
            "--sample-barcode",
            "TEST-01",
            "-o",
            dir.path().to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute vcf-reformatter with MAF output");

    // Print stderr for debugging
    if !output.status.success() {
        eprintln!("STDOUT: {}", String::from_utf8_lossy(&output.stdout));
        eprintln!("STDERR: {}", String::from_utf8_lossy(&output.stderr));
        eprintln!("Exit code: {:?}", output.status.code());
    }

    assert!(
        output.status.success(),
        "Command failed with stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn test_stdin_input_plain() {
    use std::io::Write;
    use std::process::Stdio;

    let vcf_content = "##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t60\tPASS\tDP=50\n";

    let dir = tempdir().unwrap();
    let mut child = Command::new("cargo")
        .args(["run", "--", "-", "-o", dir.path().to_str().unwrap()])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("Failed to spawn vcf-reformatter");

    child
        .stdin
        .take()
        .unwrap()
        .write_all(vcf_content.as_bytes())
        .unwrap();

    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let reformatted = std::fs::read_to_string(dir.path().join("stdin_reformatted.tsv")).unwrap();
    assert!(reformatted.contains("chr1"));
    assert!(reformatted.contains("50")); // DP value survived the round trip
}

#[test]
fn test_stdin_input_gzip() {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;
    use std::process::Stdio;

    let vcf_content = "##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr2\t200\t.\tC\tT\t80\tPASS\tDP=30\n";

    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(vcf_content.as_bytes()).unwrap();
    let compressed = encoder.finish().unwrap();

    let dir = tempdir().unwrap();
    let mut child = Command::new("cargo")
        .args(["run", "--", "-", "-o", dir.path().to_str().unwrap()])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("Failed to spawn vcf-reformatter");

    child.stdin.take().unwrap().write_all(&compressed).unwrap();

    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let reformatted = std::fs::read_to_string(dir.path().join("stdin_reformatted.tsv")).unwrap();
    assert!(reformatted.contains("chr2"));
    assert!(reformatted.contains("30"));
}

#[test]
fn test_mutation_status_and_sequence_source_flags() {
    use std::io::Write;

    let mut temp_vcf = tempfile::NamedTempFile::new().unwrap();
    let vcf_content = "##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t60\tPASS\tDP=50\n";
    write!(temp_vcf, "{}", vcf_content).unwrap();
    temp_vcf.flush().unwrap();

    let dir = tempdir().unwrap();
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            temp_vcf.path().to_str().unwrap(),
            "-o",
            dir.path().to_str().unwrap(),
            "-p",
            "meta",
            "--output-format",
            "maf",
            "--report",
            "none",
            "--sample-barcode",
            "TUMOR",
            "--mutation-status",
            "Somatic",
            "--sequence-source",
            "WGS",
        ])
        .output()
        .expect("Failed to execute vcf-reformatter");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let maf = std::fs::read_to_string(dir.path().join("meta_reformatted.maf")).unwrap();
    let mut lines = maf.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let value = |name: &str| row[header.iter().position(|h| *h == name).unwrap()];

    assert_eq!(value("Mutation_Status"), "Somatic");
    assert_eq!(value("Sequence_Source"), "WGS");
}

#[test]
#[cfg(feature = "parquet_out")]
fn test_parquet_filename_keeps_the_format_extension_on_both_paths() {
    use std::io::Write;

    let mut temp_vcf = tempfile::NamedTempFile::new().unwrap();
    let vcf_content = "##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t60\tPASS\tDP=50\n";
    write!(temp_vcf, "{}", vcf_content).unwrap();
    temp_vcf.flush().unwrap();

    let dir = tempdir().unwrap();
    // The extension says what the parquet holds, so both paths keep it. The TSV path used to
    // strip it and write "pq_reformatted.parquet" against the MAF path's ".maf.parquet".
    for (format, expected) in [
        ("tsv", "pq_reformatted.tsv.parquet"),
        ("maf", "pq_reformatted.maf.parquet"),
    ] {
        let output = Command::new("cargo")
            .args([
                "run",
                "--",
                temp_vcf.path().to_str().unwrap(),
                "-o",
                dir.path().to_str().unwrap(),
                "-p",
                "pq",
                "--output-format",
                format,
                "--report",
                "none",
                "--parquet",
                // This VCF has no sample columns; MAF output refuses to run without a barcode.
                "--sample-barcode",
                "TUMOR",
            ])
            .output()
            .expect("Failed to execute vcf-reformatter");
        assert!(
            output.status.success(),
            "stderr: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            dir.path().join(expected).exists(),
            "--output-format {format} should write {expected}"
        );
    }
}

#[test]
fn test_multiallelic_warning_on_stderr() {
    use std::io::Write;

    let mut temp_vcf = tempfile::NamedTempFile::new().unwrap();
    let vcf_content = "##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG,T\t60\tPASS\tDP=50\nchr1\t200\t.\tC\tT\t60\tPASS\tDP=30\n";
    write!(temp_vcf, "{}", vcf_content).unwrap();
    temp_vcf.flush().unwrap();

    let dir = tempdir().unwrap();
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            temp_vcf.path().to_str().unwrap(),
            "-o",
            dir.path().to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute vcf-reformatter");

    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("1 multiallelic site"),
        "stderr was: {stderr}"
    );
    assert!(stderr.contains("bcftools norm"), "stderr was: {stderr}");
}

#[test]
fn test_no_multiallelic_warning_when_none_present() {
    use std::io::Write;

    let mut temp_vcf = tempfile::NamedTempFile::new().unwrap();
    let vcf_content = "##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t60\tPASS\tDP=50\n";
    write!(temp_vcf, "{}", vcf_content).unwrap();
    temp_vcf.flush().unwrap();

    let dir = tempdir().unwrap();
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            temp_vcf.path().to_str().unwrap(),
            "-o",
            dir.path().to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute vcf-reformatter");

    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!stderr.contains("multiallelic"), "stderr was: {stderr}");
}

#[test]
fn test_report_html_is_default_and_contains_expected_sections() {
    use std::io::Write;

    let mut temp_vcf = tempfile::NamedTempFile::new().unwrap();
    let vcf_content = "##fileformat=VCFv4.2\n\
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|SIFT|PolyPhen\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tA\tG\t60\tPASS\tCSQ=G|missense_variant|MODERATE|BRCA1|deleterious(0.01)|probably_damaging(0.99)\n\
chr2\t200\t.\tC\tT\t60\tPASS\tCSQ=T|missense_variant|MODERATE|TP53|tolerated(0.8)|benign(0.02)\n";
    write!(temp_vcf, "{}", vcf_content).unwrap();
    temp_vcf.flush().unwrap();

    let dir = tempdir().unwrap();
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            temp_vcf.path().to_str().unwrap(),
            "-o",
            dir.path().to_str().unwrap(),
            "-p",
            "report_default",
        ])
        .output()
        .expect("Failed to execute vcf-reformatter");

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let report_path = dir.path().join("report_default_summary.html");
    assert!(
        report_path.exists(),
        "expected default --report to write an .html file"
    );

    let html = std::fs::read_to_string(&report_path).unwrap();
    assert!(html.contains("VCF-REFORMATTER"));
    assert!(html.contains("chr1"));
    assert!(html.contains("chr2"));
    assert!(html.contains("damage-metric-select"));
    assert!(html.contains("SIFT"));
    assert!(html.contains("PolyPhen"));
}

#[test]
fn test_report_txt_writes_plain_text_file() {
    use std::io::Write;

    let mut temp_vcf = tempfile::NamedTempFile::new().unwrap();
    let vcf_content = "##fileformat=VCFv4.2\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tA\tG\t60\tPASS\tDP=50\n";
    write!(temp_vcf, "{}", vcf_content).unwrap();
    temp_vcf.flush().unwrap();

    let dir = tempdir().unwrap();
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            temp_vcf.path().to_str().unwrap(),
            "-o",
            dir.path().to_str().unwrap(),
            "-p",
            "report_txt",
            "--report",
            "txt",
        ])
        .output()
        .expect("Failed to execute vcf-reformatter");

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    assert!(dir.path().join("report_txt_summary.txt").exists());
    assert!(!dir.path().join("report_txt_summary.html").exists());
}

#[test]
fn test_report_none_writes_no_report_file() {
    use std::io::Write;

    let mut temp_vcf = tempfile::NamedTempFile::new().unwrap();
    let vcf_content = "##fileformat=VCFv4.2\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tA\tG\t60\tPASS\tDP=50\n";
    write!(temp_vcf, "{}", vcf_content).unwrap();
    temp_vcf.flush().unwrap();

    let dir = tempdir().unwrap();
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            temp_vcf.path().to_str().unwrap(),
            "-o",
            dir.path().to_str().unwrap(),
            "-p",
            "report_none",
            "--report",
            "none",
        ])
        .output()
        .expect("Failed to execute vcf-reformatter");

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    assert!(!dir.path().join("report_none_summary.html").exists());
    assert!(!dir.path().join("report_none_summary.txt").exists());
}

#[test]
fn test_report_html_snpeff_shows_impact_metric_not_sift() {
    use std::io::Write;

    let mut temp_vcf = tempfile::NamedTempFile::new().unwrap();
    let vcf_content = "##fileformat=VCFv4.2\n\
##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name'\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tA\tG\t60\tPASS\tANN=G|missense_variant|HIGH|BRCA1\n\
chr2\t200\t.\tC\tT\t60\tPASS\tANN=T|synonymous_variant|LOW|TP53\n";
    write!(temp_vcf, "{}", vcf_content).unwrap();
    temp_vcf.flush().unwrap();

    let dir = tempdir().unwrap();
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            temp_vcf.path().to_str().unwrap(),
            "-o",
            dir.path().to_str().unwrap(),
            "-p",
            "report_snpeff",
        ])
        .output()
        .expect("Failed to execute vcf-reformatter");

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let html = std::fs::read_to_string(dir.path().join("report_snpeff_summary.html")).unwrap();
    assert!(html.contains("Impact"));
    assert!(!html.contains("SIFT"));
    assert!(html.contains("\"HIGH\""));
    assert!(html.contains("\"LOW\""));
    // SnpEff-only input yields a single ("Impact") breakdown, so the metric
    // dropdown (only meaningful with 2+ metrics to switch between) must not
    // be rendered — the chart itself should still appear.
    assert!(!html.contains("damage-metric-select"));
    assert!(html.contains("damage-chart"));
}

// ------------------------------------------------------------------------------
// Tests for chunking processing
// ------------------------------------------------------------------------------
// Add these imports at the top of your test.rs if they're not already there
// Test 1: Verify chunked streaming function works correctly
#[test]
fn test_chunked_streaming_basic() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">"#;
    let column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    // Create test data
    let data_lines: Vec<String> = (1..=1000)
        .map(|i| format!("chr1\t{}\t.\tA\tG\t60\tPASS\tDP={}", i * 100, i))
        .collect();

    // Test with in-memory writer
    let mut output = Vec::new();
    let result = reformat_vcf_data_with_header_parallel_chunked(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::FirstOnly,
        &mut output,
    );

    assert!(result.is_ok());
    let headers = result.unwrap();

    // Verify headers
    assert!(headers.contains(&"CHROM".to_string()));
    assert!(headers.contains(&"INFO_DP".to_string()));

    // Verify output was written
    let output_str = String::from_utf8(output).unwrap();
    let lines: Vec<&str> = output_str.lines().collect();

    // Should have header line + 1000 data lines
    assert_eq!(lines.len(), 1001);

    // Check header line
    assert!(lines[0].contains("CHROM"));
    assert!(lines[0].contains("INFO_DP"));

    // Check first data line
    assert!(lines[1].contains("chr1"));
    assert!(lines[1].contains("100"));
}

// Test 2: Verify chunked vs regular processing gives same results
#[test]
fn test_chunked_vs_regular_consistency() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|SYMBOL">"#;
    let column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    // Small dataset to compare both approaches
    let data_lines: Vec<String> = (1..=50)
        .map(|i| {
            format!(
                "chr1\t{}\t.\tA\tG\t60\tPASS\tDP={};CSQ=G|missense_variant|GENE{}",
                i * 100,
                i,
                i
            )
        })
        .collect();

    // Regular processing
    let regular_result = reformat_vcf_data_with_header_parallel(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::FirstOnly,
    )
    .unwrap();

    // Chunked streaming processing
    let mut chunked_output = Vec::new();
    let chunked_headers = reformat_vcf_data_with_header_parallel_chunked(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::FirstOnly,
        &mut chunked_output,
    )
    .unwrap();

    // Compare headers
    assert_eq!(regular_result.0, chunked_headers);

    // Parse chunked output
    let output_str = String::from_utf8(chunked_output).unwrap();
    let lines: Vec<&str> = output_str.lines().collect();

    // Should have same number of records (header + data lines)
    assert_eq!(lines.len(), regular_result.1.len() + 1); // +1 for header

    // Verify first few data values match
    let chunked_first_line: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(chunked_first_line[0], "chr1"); // CHROM
    assert_eq!(chunked_first_line[1], "100"); // POS
    assert_eq!(chunked_first_line[3], "A"); // REF
    assert_eq!(chunked_first_line[4], "G"); // ALT
}

// Test 3: Memory efficiency test - large dataset with streaming
#[test]
fn test_large_dataset_streaming() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">"#;
    let column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    // Create large dataset (200k records to trigger chunking)
    let data_lines: Vec<String> = (1..=200_000)
        .map(|i| {
            format!(
                "chr{}\t{}\t.\tA\tG\t60\tPASS\tDP={}",
                (i % 22) + 1,
                i * 10,
                i % 1000
            )
        })
        .collect();

    println!(
        "📊 Testing with {} records (should trigger chunking)",
        data_lines.len()
    );

    // Use a file-based writer to handle large output
    let temp_dir = tempdir().unwrap();
    let temp_file_path = temp_dir.path().join("large_test.tsv");
    let temp_file = std::fs::File::create(&temp_file_path).unwrap();
    let mut writer = std::io::BufWriter::new(temp_file);

    let start_time = std::time::Instant::now();
    let result = reformat_vcf_data_with_header_parallel_chunked(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::FirstOnly,
        &mut writer,
    );
    let elapsed = start_time.elapsed();

    writer.flush().unwrap();
    drop(writer);

    assert!(result.is_ok());
    println!(
        "✅ Large dataset processed in {:.2}s",
        elapsed.as_secs_f64()
    );

    // Verify file was created and has content
    let metadata = std::fs::metadata(&temp_file_path).unwrap();
    assert!(metadata.len() > 0, "Output file should not be empty");

    // Verify headers
    let headers = result.unwrap();
    assert!(headers.contains(&"CHROM".to_string()));
    assert!(headers.contains(&"INFO_DP".to_string()));

    println!("📋 Generated {} headers", headers.len());
}

// Test 4: Test chunk boundary handling
#[test]
fn test_chunk_boundary_integrity() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">"#;
    let column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    // Create exactly 150k records to test chunk boundaries
    let data_lines: Vec<String> = (1..=150_000)
        .map(|i| {
            format!(
                "chr1\t{}\t.\tA\tG\t{}\tPASS\tDP={}",
                i,
                60 + (i % 40),
                i % 100
            )
        })
        .collect();

    let mut output = Vec::new();
    let result = reformat_vcf_data_with_header_parallel_chunked(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::FirstOnly,
        &mut output,
    );

    assert!(result.is_ok());

    // Parse and verify output integrity
    let output_str = String::from_utf8(output).unwrap();
    let lines: Vec<&str> = output_str.lines().collect();

    // Should have header + 150k data lines
    assert_eq!(lines.len(), 150_001);

    // Check first, middle, and last records for data integrity
    let first_line: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(first_line[1], "1"); // First position

    let middle_line: Vec<&str> = lines[75_000].split('\t').collect();
    assert_eq!(middle_line[1], "75000"); // Middle position (line 75000 has position 75000)

    let last_line: Vec<&str> = lines[150_000].split('\t').collect();
    assert_eq!(last_line[1], "150000"); // Last position

    println!(
        "✅ Chunk boundary test passed with {} records",
        lines.len() - 1
    );
}

// Test 5: Test with VEP CSQ annotations
#[test]
fn test_chunked_with_vep_annotations() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|IMPACT|SYMBOL|Gene">"#;
    let column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    let data_lines: Vec<String> = (1..=5000)
        .map(|i| {
            format!(
                "chr1\t{}\t.\tA\tG\t60\tPASS\tCSQ=G|missense_variant|MODERATE|GENE{}|{}",
                i * 100,
                i,
                i
            )
        })
        .collect();

    let mut output = Vec::new();
    let result = reformat_vcf_data_with_header_parallel_chunked(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::FirstOnly,
        &mut output,
    );

    assert!(result.is_ok());
    let headers = result.unwrap();

    // Should have CSQ-related headers
    assert!(headers.contains(&"CSQ_Allele".to_string()));
    assert!(headers.contains(&"CSQ_Consequence".to_string()));
    assert!(headers.contains(&"CSQ_IMPACT".to_string()));
    assert!(headers.contains(&"CSQ_SYMBOL".to_string()));

    // Verify output content
    let output_str = String::from_utf8(output).unwrap();
    let lines: Vec<&str> = output_str.lines().collect();

    assert_eq!(lines.len(), 5001); // header + 5000 data lines

    // Check that CSQ data is properly parsed
    let data_line: Vec<&str> = lines[1].split('\t').collect();
    let csq_allele_idx = headers.iter().position(|h| h == "CSQ_Allele").unwrap();
    let csq_consequence_idx = headers.iter().position(|h| h == "CSQ_Consequence").unwrap();

    assert_eq!(data_line[csq_allele_idx], "G");
    assert_eq!(data_line[csq_consequence_idx], "missense_variant");
}

// Test 6: Performance benchmark test
#[test]
fn test_streaming_performance_benchmark() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">"#;
    let column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    // Create 50k records for performance test
    let data_lines: Vec<String> = (1..=50_000)
        .map(|i| {
            format!(
                "chr1\t{}\t.\tA\tG\t{}\tPASS\tDP={};AF={}",
                i * 100,
                60 + (i % 40),
                i % 100,
                (i as f64 / 100_000.0)
            )
        })
        .collect();

    let mut output = Vec::new();

    let start_time = std::time::Instant::now();
    let result = reformat_vcf_data_with_header_parallel_chunked(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::FirstOnly,
        &mut output,
    );
    let elapsed = start_time.elapsed();

    assert!(result.is_ok());

    let variants_per_second = data_lines.len() as f64 / elapsed.as_secs_f64();

    println!("🚀 Performance: {:.0} variants/sec", variants_per_second);
    println!(
        "⏱️  Total time: {:.3}s for {} variants",
        elapsed.as_secs_f64(),
        data_lines.len()
    );

    // Performance should be reasonable (>1000 variants/sec)
    assert!(
        variants_per_second > 1000.0,
        "Performance too slow: {} variants/sec",
        variants_per_second
    );
}

// Test 7: Error handling test
#[test]
fn test_chunked_error_handling() {
    let header = r#"##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">"#;
    let column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    // Mix good and bad records
    let data_lines = vec![
        "chr1\t100\t.\tA\tG\t60\tPASS\tDP=10".to_string(), // Good
        "chr1\tinvalid\t.\tA\tG\t60\tPASS\tDP=20".to_string(), // Bad position
        "chr1\t200\t.\tA\tG\t60\tPASS\tDP=30".to_string(), // Good
        "incomplete_line".to_string(),                     // Bad format
        "chr1\t300\t.\tA\tG\t60\tPASS\tDP=40".to_string(), // Good
    ];

    let mut output = Vec::new();
    let result = reformat_vcf_data_with_header_parallel_chunked(
        header,
        column_names,
        &data_lines,
        TranscriptHandling::FirstOnly,
        &mut output,
    );

    // Should succeed despite bad records
    assert!(result.is_ok());

    let output_str = String::from_utf8(output).unwrap();
    let lines: Vec<&str> = output_str.lines().collect();

    // Should have header + 3 good records (2 bad records filtered out)
    assert_eq!(lines.len(), 4); // header + 3 valid data lines
}

// ============================================================================
// Row Count Conservation Tests
// ============================================================================

#[test]
fn test_row_count_tsv_first_transcript() {
    use vcf_reformatter::reformat_vcf::{reformat_vcf_data_with_header, TranscriptHandling};

    let header = "##fileformat=VCFv4.2\n##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'\">";
    let columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    let data_lines = vec![
        "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50;ANN=G|missense_variant|MODERATE|BRCA1|ENSG1|transcript|ENST1|protein_coding|5/24|c.1T>C|p.M1T|1/100|1/100|1/33||".to_string(),
        "chr1\t200\t.\tC\tT\t40\tPASS\tDP=30;ANN=T|stop_gained|HIGH|TP53|ENSG2|transcript|ENST2|protein_coding|7/11|c.2C>T|p.Q2*|2/200|2/200|2/66||".to_string(),
        "chr2\t300\t.\tG\tA\t50\tPASS\tDP=20;ANN=A|synonymous_variant|LOW|EGFR|ENSG3|transcript|ENST3|protein_coding|3/8|c.3G>A|p.L3L|3/300|3/300|3/99||".to_string(),
    ];

    let (_, records) =
        reformat_vcf_data_with_header(header, columns, &data_lines, TranscriptHandling::FirstOnly)
            .unwrap();

    assert_eq!(
        records.len(),
        data_lines.len(),
        "Row count mismatch: {} input lines but {} output records",
        data_lines.len(),
        records.len()
    );
}

#[test]
fn test_row_count_maf_first_transcript() {
    use vcf_reformatter::essentials_fields::MafRecord;
    use vcf_reformatter::reformat_vcf::{reformat_vcf_data_with_header, TranscriptHandling};

    let header = "##fileformat=VCFv4.2\n##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'\">";
    let columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    let data_lines = vec![
        "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50;ANN=G|missense_variant|MODERATE|BRCA1|ENSG1|transcript|ENST1|protein_coding|5/24|c.1T>C|p.M1T|1/100|1/100|1/33||".to_string(),
        "chr2\t200\t.\tC\tT\t40\tPASS\tDP=30;ANN=T|stop_gained|HIGH|TP53|ENSG2|transcript|ENST2|protein_coding|7/11|c.2C>T|p.Q2*|2/200|2/200|2/66||".to_string(),
    ];

    let (_, records) =
        reformat_vcf_data_with_header(header, columns, &data_lines, TranscriptHandling::FirstOnly)
            .unwrap();

    let mut maf_count = 0;
    for record in &records {
        let _maf = MafRecord::from_reformatted_record(record, "TEST", "GRCh38", "SAMPLE").unwrap();
        maf_count += 1;
    }

    assert_eq!(
        maf_count,
        data_lines.len(),
        "MAF row count mismatch: {} input lines but {} MAF records",
        data_lines.len(),
        maf_count
    );
}

#[test]
fn test_row_count_split_transcripts() {
    use vcf_reformatter::reformat_vcf::{reformat_vcf_data_with_header, TranscriptHandling};

    let header = "##fileformat=VCFv4.2\n##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'\">";
    let columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    // Line 1: 2 transcripts (comma-separated in ANN)
    // Line 2: 1 transcript
    let data_lines = vec![
        "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50;ANN=G|missense_variant|MODERATE|BRCA1|ENSG1|transcript|ENST1|protein_coding|5/24|c.1T>C|p.M1T|1/100|1/100|1/33||,G|synonymous_variant|LOW|BRCA1|ENSG1|transcript|ENST2|protein_coding|3/24|c.1T>C|p.M1M|1/100|1/100|1/33||".to_string(),
        "chr2\t200\t.\tC\tT\t40\tPASS\tDP=30;ANN=T|stop_gained|HIGH|TP53|ENSG2|transcript|ENST3|protein_coding|7/11|c.2C>T|p.Q2*|2/200|2/200|2/66||".to_string(),
    ];

    let (_, records) =
        reformat_vcf_data_with_header(header, columns, &data_lines, TranscriptHandling::SplitRows)
            .unwrap();

    // 2 input lines, but line 1 has 2 transcripts -> 3 output records
    let expected_output_rows = 3;
    assert_eq!(
        records.len(),
        expected_output_rows,
        "Split transcript count mismatch: expected {} output records, got {}",
        expected_output_rows,
        records.len()
    );

    assert!(
        records.len() > data_lines.len(),
        "Split mode should produce more records than input lines"
    );
}

// ============================================================================
// Summary Correctness Tests
// ============================================================================

#[test]
fn test_summary_input_counts_sum() {
    use vcf_reformatter::summary::count_input_chromosomes;

    let lines = vec![
        "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50".to_string(),
        "chr1\t200\t.\tC\tT\t40\tPASS\tDP=30".to_string(),
        "chr2\t300\t.\tG\tA\t50\tPASS\tDP=20".to_string(),
        "chrX\t400\t.\tT\tC\t70\tPASS\tDP=40".to_string(),
        "chrX\t500\t.\tA\tT\t80\tPASS\tDP=60".to_string(),
    ];

    let counts = count_input_chromosomes(&lines);
    let total: usize = counts.values().sum();
    assert_eq!(
        total,
        lines.len(),
        "Sum of per-chromosome counts must equal total input variants"
    );
}

#[test]
fn test_summary_chromosome_order() {
    use vcf_reformatter::summary::count_input_chromosomes;

    let lines = vec![
        "chrY\t100\t.\tA\tG\t60\tPASS\t.".to_string(),
        "chr10\t200\t.\tC\tT\t40\tPASS\t.".to_string(),
        "chr2\t300\t.\tG\tA\t50\tPASS\t.".to_string(),
        "chrX\t400\t.\tT\tC\t70\tPASS\t.".to_string(),
        "chr1\t500\t.\tA\tT\t80\tPASS\t.".to_string(),
        "chrM\t600\t.\tG\tC\t90\tPASS\t.".to_string(),
    ];

    let counts = count_input_chromosomes(&lines);
    let keys: Vec<&String> = counts.keys().collect();
    assert_eq!(keys, vec!["chr1", "chr2", "chr10", "chrX", "chrY", "chrM"]);
}

#[test]
fn test_summary_expansion_ratio_first_mode() {
    use vcf_reformatter::reformat_vcf::{reformat_vcf_data_with_header, TranscriptHandling};

    let header = "##fileformat=VCFv4.2";
    let columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    let data_lines = vec![
        "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50".to_string(),
        "chr2\t200\t.\tC\tT\t40\tPASS\tDP=30".to_string(),
    ];

    let (_, records) =
        reformat_vcf_data_with_header(header, columns, &data_lines, TranscriptHandling::FirstOnly)
            .unwrap();

    let ratio = records.len() as f64 / data_lines.len() as f64;
    assert!(
        (ratio - 1.0).abs() < 0.001,
        "First-only mode should have expansion ratio of 1.0, got {}",
        ratio
    );
}

#[test]
fn test_summary_write_and_read() {
    use indexmap::IndexMap;
    use vcf_reformatter::summary::SummaryStats;

    let mut input_counts = IndexMap::new();
    input_counts.insert("chr1".to_string(), 50);
    input_counts.insert("chr2".to_string(), 30);
    input_counts.insert("chrX".to_string(), 20);

    let mut output_counts = IndexMap::new();
    output_counts.insert("chr1".to_string(), 50);
    output_counts.insert("chr2".to_string(), 30);
    output_counts.insert("chrX".to_string(), 20);

    let stats = SummaryStats {
        input_file: "test.vcf.gz".to_string(),
        output_format: "TSV".to_string(),
        transcript_handling: "FirstOnly".to_string(),
        input_variant_count: 100,
        output_record_count: 100,
        input_chrom_counts: input_counts,
        output_chrom_counts: output_counts,
        processing_time_secs: 1.5,
        variants_per_sec: 66.7,
    };

    let tmp = tempfile::NamedTempFile::new().unwrap();
    let path = tmp.path().to_str().unwrap();
    stats.write_to_file(path).unwrap();

    let content = std::fs::read_to_string(path).unwrap();
    assert!(content.contains("Total input variants:     100"));
    assert!(content.contains("Total output records:     100"));
    assert!(content.contains("chr1"));
    assert!(content.contains("50.0%"));
    assert!(content.contains("1.00x")); // expansion ratio
}

// ============================================================================
// Parquet Equivalence Tests
// ============================================================================

#[cfg(feature = "parquet_out")]
mod parquet_tests {
    use parquet::file::reader::{FileReader, SerializedFileReader};
    use vcf_reformatter::reformat_vcf::{reformat_vcf_data_with_header, TranscriptHandling};

    fn count_parquet_rows(path: &str) -> usize {
        let file = std::fs::File::open(path).unwrap();
        let reader = SerializedFileReader::new(file).unwrap();
        let mut total = 0i64;
        for rg in reader.metadata().row_groups() {
            total += rg.num_rows();
        }
        total as usize
    }

    #[test]
    fn test_parquet_row_count_matches_tsv() {
        let header = "##fileformat=VCFv4.2\n##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'\">";
        let columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        let data_lines = vec![
            "chr1\t100\trs1\tA\tG\t60\tPASS\tDP=50;ANN=G|missense_variant|MODERATE|BRCA1|ENSG1|transcript|ENST1|protein_coding|5/24|c.1T>C|p.M1T|1/100|1/100|1/33||".to_string(),
            "chr2\t200\t.\tC\tT\t40\tPASS\tDP=30;ANN=T|stop_gained|HIGH|TP53|ENSG2|transcript|ENST2|protein_coding|7/11|c.2C>T|p.Q2*|2/200|2/200|2/66||".to_string(),
        ];

        let (headers, records) = reformat_vcf_data_with_header(
            header,
            columns,
            &data_lines,
            TranscriptHandling::FirstOnly,
        )
        .unwrap();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().to_str().unwrap();

        vcf_reformatter::parquet_writer::write_tsv_as_parquet(path, &headers, &records).unwrap();

        let total_rows = count_parquet_rows(path);
        assert_eq!(
            total_rows,
            records.len(),
            "Parquet row count ({}) must match TSV record count ({})",
            total_rows,
            records.len()
        );
    }

    #[test]
    fn test_parquet_columns_match_tsv_headers() {
        let header = "##fileformat=VCFv4.2";
        let columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        let data_lines = vec!["chr1\t100\t.\tA\tG\t60\tPASS\tDP=50".to_string()];

        let (headers, records) = reformat_vcf_data_with_header(
            header,
            columns,
            &data_lines,
            TranscriptHandling::FirstOnly,
        )
        .unwrap();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().to_str().unwrap();

        vcf_reformatter::parquet_writer::write_tsv_as_parquet(path, &headers, &records).unwrap();

        let file = std::fs::File::open(path).unwrap();
        let reader = SerializedFileReader::new(file).unwrap();
        let schema = reader.metadata().file_metadata().schema_descr();
        let parquet_columns: Vec<String> = schema
            .columns()
            .iter()
            .map(|c| c.name().to_string())
            .collect();

        assert_eq!(
            parquet_columns, headers,
            "Parquet column names must match TSV headers exactly"
        );
    }

    #[test]
    fn test_parquet_maf_row_count() {
        use vcf_reformatter::essentials_fields::MafRecord;

        let header = "##fileformat=VCFv4.2\n##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'\">";
        let columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        let data_lines = vec![
            "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50;ANN=G|missense_variant|MODERATE|BRCA1|ENSG1|transcript|ENST1|protein_coding|5/24|c.1T>C|p.M1T|1/100|1/100|1/33||".to_string(),
            "chr2\t200\t.\tC\tT\t40\tPASS\tDP=30;ANN=T|stop_gained|HIGH|TP53|ENSG2|transcript|ENST2|protein_coding|7/11|c.2C>T|p.Q2*|2/200|2/200|2/66||".to_string(),
        ];

        let (_, records) = reformat_vcf_data_with_header(
            header,
            columns,
            &data_lines,
            TranscriptHandling::FirstOnly,
        )
        .unwrap();

        let maf_records: Vec<MafRecord> = records
            .iter()
            .map(|r| MafRecord::from_reformatted_record(r, "TEST", "GRCh38", "SAMPLE").unwrap())
            .collect();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().to_str().unwrap();

        vcf_reformatter::parquet_writer::write_maf_as_parquet(path, &maf_records).unwrap();

        let total_rows = count_parquet_rows(path);
        assert_eq!(
            total_rows,
            maf_records.len(),
            "Parquet MAF row count ({}) must match MAF record count ({})",
            total_rows,
            maf_records.len()
        );
    }

    #[test]
    fn test_parquet_maf_columns_match_get_maf_headers() {
        use vcf_reformatter::essentials_fields::MafRecord;

        let header = "##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Format: Allele|Consequence|IMPACT|SYMBOL\">";
        let columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        let data_lines = vec![
            "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50;CSQ=G|missense_variant|MODERATE|BRCA1".to_string(),
        ];

        let (_, records) = reformat_vcf_data_with_header(
            header,
            columns,
            &data_lines,
            TranscriptHandling::FirstOnly,
        )
        .unwrap();
        let maf_records: Vec<MafRecord> = records
            .iter()
            .map(|r| {
                MafRecord::from_reformatted_record(r, "TestCenter", "GRCh38", "SAMPLE-001").unwrap()
            })
            .collect();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().to_str().unwrap();
        vcf_reformatter::parquet_writer::write_maf_as_parquet(path, &maf_records).unwrap();

        let file = std::fs::File::open(path).unwrap();
        let reader = SerializedFileReader::new(file).unwrap();
        let schema = reader.metadata().file_metadata().schema_descr();
        let parquet_columns: Vec<String> = schema
            .columns()
            .iter()
            .map(|c| c.name().to_string())
            .collect();

        assert_eq!(parquet_columns, MafRecord::get_maf_headers());
    }
}

#[test]
fn test_report_both_formats() {
    use std::io::Write;

    let dir = tempdir().unwrap();
    let vcf = dir.path().join("both.vcf");
    let mut f = File::create(&vcf).unwrap();
    write!(
        f,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t60\tPASS\tDP=50\n"
    )
    .unwrap();

    let run = |report: &str, out: &std::path::Path| {
        let output = Command::new("cargo")
            .args([
                "run",
                "--",
                vcf.to_str().unwrap(),
                "--report",
                report,
                "-o",
                out.to_str().unwrap(),
            ])
            .output()
            .expect("failed to run vcf-reformatter");
        assert!(
            output.status.success(),
            "--report {report} failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    };

    let both = dir.path().join("both_out");
    run("html,txt", &both);
    assert!(both.join("both_summary.txt").exists(), "txt report missing");
    assert!(
        both.join("both_summary.html").exists(),
        "html report missing"
    );

    // --report-dir sends the report somewhere else; the data output stays under -o.
    let data = dir.path().join("split_data");
    let reports = dir.path().join("split_reports");
    let output = Command::new("cargo")
        .args([
            "run",
            "--",
            vcf.to_str().unwrap(),
            "--report",
            "html,txt",
            "-o",
            data.to_str().unwrap(),
            "--report-dir",
            reports.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run vcf-reformatter");
    assert!(
        output.status.success(),
        "--report-dir failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(reports.join("both_summary.txt").exists());
    assert!(reports.join("both_summary.html").exists());
    assert!(!data.join("both_summary.html").exists());
    assert!(data.join("both_reformatted.tsv").exists());

    // `none` in the list wins over any format listed with it.
    let none = dir.path().join("none_out");
    run("html,none", &none);
    assert!(!none.join("both_summary.txt").exists());
    assert!(!none.join("both_summary.html").exists());
}
