//! VCF in, TSV records out: reading, INFO/FORMAT parsing, transcript handling, the TSV writer.
mod common;

use std::collections::HashMap;
use std::io::Write;
use tempfile::tempdir;

use vcf_reformatter::extract_sample_info::parse_format_and_samples;
use vcf_reformatter::read_vcf_gz::open_vcf;
use vcf_reformatter::reformat_vcf::{
    parse_info_field, reformat_vcf_data_with_header, reformat_vcf_data_with_header_parallel,
    write_tsv_header, write_tsv_rows, AnnotationFieldType, ReformattedVcfRecord,
    TranscriptHandling,
};

const COLUMNS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
const DP_HEADER: &str =
    "##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">";
const CSQ_HEADER: &str = "##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT\">";

// --- reading ------------------------------------------------------------------------------

#[test]
fn open_vcf_splits_header_from_data_and_skips_blank_lines() {
    let dir = tempdir().unwrap();
    let path = common::write_vcf(
        dir.path(),
        "in.vcf",
        "##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         \n\
         chr1\t100\t.\tA\tG\t60\tPASS\tDP=10\n\
         \x20 \n\
         chr2\t200\trs123\tC\tT\t80\tPASS\tDP=20\n",
    );

    let stream = open_vcf(&path).unwrap();
    assert_eq!(stream.header.lines().count(), 2); // the two ## lines
    assert!(stream.columns_title.starts_with("#CHROM"));
    let data: Vec<String> = stream.lines.map(Result::unwrap).collect();
    assert_eq!(data.len(), 2, "blank and whitespace-only lines are skipped");
    assert!(data[1].starts_with("chr2\t200\trs123"));
}

// --- INFO parsing ---------------------------------------------------------------------------

#[test]
fn parse_info_field_without_annotations() {
    let (records, field_type) =
        parse_info_field("DP=10;AF=0.5", &None, &None, TranscriptHandling::FirstOnly).unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].get("INFO_DP"), Some(&"10".to_string()));
    assert_eq!(records[0].get("INFO_AF"), Some(&"0.5".to_string()));
    assert_eq!(field_type, AnnotationFieldType::None);

    let (records, field_type) =
        parse_info_field("", &None, &None, TranscriptHandling::FirstOnly).unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(field_type, AnnotationFieldType::None);
}

#[test]
fn parse_info_field_with_csq() {
    let csq = Some(vec!["Allele".to_string(), "Consequence".to_string()]);
    let (records, field_type) = parse_info_field(
        "DP=10;CSQ=A|missense_variant,T|synonymous_variant;AF=0.5",
        &csq,
        &None,
        TranscriptHandling::SplitRows,
    )
    .unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].get("CSQ_Allele"), Some(&"A".to_string()));
    assert_eq!(
        records[0].get("CSQ_Consequence"),
        Some(&"missense_variant".to_string())
    );
    assert_eq!(records[0].get("INFO_DP"), Some(&"10".to_string()));
    assert_eq!(field_type, AnnotationFieldType::Csq);
}

#[test]
fn parse_info_field_with_ann() {
    let ann = Some(
        ["Allele", "Annotation", "Annotation_Impact", "Gene_Name"]
            .map(String::from)
            .to_vec(),
    );
    let (records, field_type) = parse_info_field(
        "DP=15;ANN=T|missense_variant|MODERATE|BRCA1,G|synonymous_variant|LOW|TP53;AF=0.3",
        &None,
        &ann,
        TranscriptHandling::SplitRows,
    )
    .unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].get("ANN_Allele"), Some(&"T".to_string()));
    assert_eq!(
        records[0].get("ANN_Annotation"),
        Some(&"missense_variant".to_string())
    );
    assert_eq!(records[0].get("ANN_Gene_Name"), Some(&"BRCA1".to_string()));
    assert_eq!(records[0].get("INFO_DP"), Some(&"15".to_string()));
    assert_eq!(records[0].get("INFO_AF"), Some(&"0.3".to_string()));
    assert_eq!(field_type, AnnotationFieldType::Ann);
}

// --- FORMAT / sample parsing ----------------------------------------------------------------

#[test]
fn missing_trailing_format_values_read_as_dot() {
    let parsed = parse_format_and_samples(
        Some("GT:DP:AD:PL"),
        &["0/1:20:10,10".to_string()], // no PL
        &["SAMPLE1".to_string()],
    )
    .unwrap();
    assert_eq!(parsed.format_keys, vec!["GT", "DP", "AD", "PL"]);
    let s = &parsed.samples[0].format_fields;
    assert_eq!(s.get("AD").unwrap(), "10,10");
    assert_eq!(s.get("PL").unwrap(), ".");
}

#[test]
fn no_format_column_yields_no_samples() {
    let parsed = parse_format_and_samples(
        None,
        &["0/1:20:10,10".to_string()],
        &["SAMPLE1".to_string()],
    )
    .unwrap();
    assert!(parsed.format_keys.is_empty());
    assert!(parsed.samples.is_empty());
}

#[test]
fn sample_names_with_underscores_survive_into_column_headers() {
    // Column names are `<sample>_<key>`; a sample name that itself contains underscores
    // must not be split or mangled on the way through.
    let names = [
        "SAMPLE_WITH_MANY_UNDERSCORES_1",
        "ANOTHER_COMPLEX_SAMPLE_NAME_2",
        "FINAL_TEST_SAMPLE_NAME_3",
    ]
    .map(String::from);
    let fields = ["0/0:100", "0/1:200", "1/1:300"].map(String::from);
    let parsed = parse_format_and_samples(Some("GT:DP"), &fields, &names).unwrap();

    assert_eq!(parsed.samples.len(), 3);
    assert_eq!(parsed.samples[2].sample_name, "FINAL_TEST_SAMPLE_NAME_3");
    assert_eq!(parsed.samples[2].format_fields.get("DP").unwrap(), "300");
    assert_eq!(
        parsed.get_headers_for_samples(),
        vec![
            "SAMPLE_WITH_MANY_UNDERSCORES_1_GT",
            "SAMPLE_WITH_MANY_UNDERSCORES_1_DP",
            "ANOTHER_COMPLEX_SAMPLE_NAME_2_GT",
            "ANOTHER_COMPLEX_SAMPLE_NAME_2_DP",
            "FINAL_TEST_SAMPLE_NAME_3_GT",
            "FINAL_TEST_SAMPLE_NAME_3_DP",
        ]
    );
}

#[test]
fn from_vcf_line_carries_sample_columns() {
    let line = "chr1\t123456\t.\tA\tG\t1000\tPASS\tDP=570;AF=0.5\tGT:AD:AF:DP:F1R2:F2R1:FAD:SB\t0/0:257,4:0.017:261:109,0:87,3:229,3:164,93,3,1\t0/1:303,6:0.020:309:105,3:106,1:245,4:187,116,2,4";
    let columns = [
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
        line,
        &columns,
        &None,
        &None,
        TranscriptHandling::FirstOnly,
    )
    .unwrap();
    assert_eq!(records.len(), 1);
    let r = &records[0];
    assert_eq!((r.chromosome.as_str(), r.position), ("chr1", 123456));

    let samples = &r.format_sample_data.as_ref().unwrap().samples;
    assert_eq!(samples.len(), 2);
    assert_eq!(samples[0].sample_name, "B487_B487_1_cOM");
    assert_eq!(samples[0].format_fields.get("AD").unwrap(), "257,4");
    assert_eq!(samples[1].sample_name, "B487_B487_2_LN");
    assert_eq!(samples[1].format_fields.get("GT").unwrap(), "0/1");
}

// --- whole-file conversion ------------------------------------------------------------------

#[test]
fn reformat_plain_info_fields() {
    let lines = [
        "chr1\t100\t.\tA\tG\t60\tPASS\tDP=10",
        "chr2\t200\trs123\tC\tT\t80\tPASS\tDP=20",
    ]
    .map(String::from);
    let (headers, records) =
        reformat_vcf_data_with_header(DP_HEADER, COLUMNS, &lines, TranscriptHandling::FirstOnly)
            .unwrap();
    assert!(headers.contains(&"CHROM".to_string()));
    assert!(headers.contains(&"INFO_DP".to_string()));
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].chromosome, "chr1");
    assert_eq!(records[0].position, 100);
    assert_eq!(records[1].info_fields.get("INFO_DP").unwrap(), "20");
}

#[test]
fn reformat_csq_first_only_keeps_the_annotators_first_entry() {
    let lines = [
        "chr1\t100\t.\tA\tG\t60\tPASS\tCSQ=G|missense_variant|MODERATE",
        "chr2\t200\trs123\tC\tT\t80\tPASS\tCSQ=T|synonymous_variant|LOW,T|intron_variant|MODIFIER",
    ]
    .map(String::from);
    let (headers, records) =
        reformat_vcf_data_with_header(CSQ_HEADER, COLUMNS, &lines, TranscriptHandling::FirstOnly)
            .unwrap();
    for h in ["CSQ_Allele", "CSQ_Consequence", "CSQ_IMPACT"] {
        assert!(headers.contains(&h.to_string()), "missing {h}");
    }
    assert_eq!(records.len(), 2);
    assert_eq!(
        records[1].info_fields.get("CSQ_Consequence").unwrap(),
        "synonymous_variant"
    );
    assert_eq!(records[1].info_fields.get("CSQ_IMPACT").unwrap(), "LOW");
}

#[test]
fn reformat_csq_split_rows_yields_one_record_per_transcript() {
    let lines = [
        "chr1\t100\t.\tA\tG\t60\tPASS\tCSQ=G|missense_variant|MODERATE,G|intron_variant|MODIFIER",
    ]
    .map(String::from);
    let (_, records) =
        reformat_vcf_data_with_header(CSQ_HEADER, COLUMNS, &lines, TranscriptHandling::SplitRows)
            .unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(
        records[0].info_fields.get("CSQ_Consequence").unwrap(),
        "missense_variant"
    );
    assert_eq!(
        records[1].info_fields.get("CSQ_Consequence").unwrap(),
        "intron_variant"
    );
}

#[test]
fn reformat_csq_most_severe_picks_by_consequence_not_list_order() {
    let lines = ["chr1\t100\t.\tA\tG\t60\tPASS\tCSQ=G|synonymous_variant|LOW,G|stop_gained|HIGH,G|intron_variant|MODIFIER"]
        .map(String::from);
    let (_, records) =
        reformat_vcf_data_with_header(CSQ_HEADER, COLUMNS, &lines, TranscriptHandling::MostSevere)
            .unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(
        records[0].info_fields.get("CSQ_Consequence").unwrap(),
        "stop_gained"
    );
}

#[test]
fn reformat_skips_unparseable_lines_and_keeps_the_rest() {
    let lines = [
        "chr1\t100\t.\tA\tG\t60\tPASS\tDP=10",
        "chr1\tinvalid\t.\tA\tG\t60\tPASS\tDP=20", // bad POS
        "chr1\t200\t.\tA\tG\t60\tPASS\tDP=30",
        "incomplete_line",
        "chr1\t300\t.\tA\tG\t60\tPASS\tDP=40",
    ]
    .map(String::from);
    let (_, records) =
        reformat_vcf_data_with_header(DP_HEADER, COLUMNS, &lines, TranscriptHandling::FirstOnly)
            .unwrap();
    let positions: Vec<u64> = records.iter().map(|r| r.position).collect();
    assert_eq!(positions, vec![100, 200, 300]);
}

#[test]
fn parallel_path_returns_every_record() {
    let lines: Vec<String> = (1..=100)
        .map(|i| format!("chr1\t{}\t.\tA\tG\t60\tPASS\tDP={}", i * 100, i))
        .collect();
    let (_, records) = reformat_vcf_data_with_header_parallel(
        DP_HEADER,
        COLUMNS,
        &lines,
        TranscriptHandling::FirstOnly,
    )
    .unwrap();
    assert_eq!(records.len(), 100);
}

// --- writing --------------------------------------------------------------------------------

#[test]
fn tsv_writer_renders_header_and_rows_in_header_order() {
    let mut r = common::record("chr1", 100, "A", "G", Some(60.0), "PASS", HashMap::new());
    r.id = None;
    r.info_fields = common::info(&[
        ("INFO_DP", "10"),
        ("CSQ_Allele", "G"),
        ("CSQ_Consequence", "missense_variant"),
    ]);
    let headers: Vec<String> = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO_DP",
        "CSQ_Allele",
        "CSQ_Consequence",
    ]
    .map(String::from)
    .to_vec();

    let mut out = Vec::new();
    write_tsv_header(&mut out, &headers).unwrap();
    write_tsv_rows(&mut out, &headers, &[r]).unwrap();
    out.flush().unwrap();

    let text = String::from_utf8(out).unwrap();
    let mut lines = text.lines();
    assert_eq!(lines.next().unwrap(), headers.join("\t"));
    assert_eq!(
        lines.next().unwrap(),
        "chr1\t100\t.\tA\tG\t60\tPASS\t10\tG\tmissense_variant"
    );
    assert!(lines.next().is_none());
}
