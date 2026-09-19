//! Annotation fields: pulling CSQ/ANN out of INFO, reading their layout from the header, and
//! SnpEff-specific parsing.
use vcf_reformatter::extract_csq_and_csq_names::extract_csq_regex;
use vcf_reformatter::get_info_from_header::{
    extract_ann_format_from_header, extract_csq_format_from_header,
};
use vcf_reformatter::reformat_vcf::{
    get_ann_impact_severity, reformat_vcf_data_with_header, sanitize_field_name,
    ReformattedVcfRecord, TranscriptHandling,
};

const SNPEFF_HEADER: &str = "##fileformat=VCFv4.2\n##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'\">\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">";

fn vcf_fields(info: &str) -> Vec<String> {
    ["chr1", "100", ".", "A", "G", "60", "PASS", info]
        .map(String::from)
        .to_vec()
}

fn ann_names() -> Option<Vec<String>> {
    extract_ann_format_from_header(SNPEFF_HEADER)
}

// --- CSQ extraction ---------------------------------------------------------------------------

#[test]
fn csq_is_cut_out_of_info_and_returned() {
    let mut fields = vcf_fields("DP=10;CSQ=A|B|C;AF=0.5");
    assert_eq!(extract_csq_regex(&mut fields).as_deref(), Some("A|B|C"));
    assert_eq!(fields[7], "DP=10;AF=0.5");
}

#[test]
fn csq_absent_or_empty_leaves_info_alone() {
    let mut fields = vcf_fields("DP=10;AF=0.5");
    assert!(extract_csq_regex(&mut fields).is_none());
    assert_eq!(fields[7], "DP=10;AF=0.5");

    let mut fields = vcf_fields("DP=10;CSQ=;AF=0.5");
    assert!(extract_csq_regex(&mut fields).is_none());
    assert_eq!(fields[7], "DP=10;CSQ=;AF=0.5");

    let mut short = vcf_fields("x");
    short.pop();
    assert!(extract_csq_regex(&mut short).is_none());
}

// --- header layout ----------------------------------------------------------------------------

#[test]
fn csq_layout_comes_from_the_format_clause() {
    let header = "##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene\">";
    assert_eq!(
        extract_csq_format_from_header(header).unwrap(),
        vec!["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene"]
    );
    assert!(extract_csq_format_from_header(
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
    )
    .is_none());
    assert!(extract_csq_format_from_header(
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"VEP annotations\">"
    )
    .is_none());
}

#[test]
fn ann_layout_keeps_snpeffs_spaced_names() {
    let fields = ann_names().unwrap();
    assert_eq!(fields.len(), 16);
    assert_eq!(fields[0], "Allele");
    assert_eq!(fields[9], "HGVS.c");
    assert_eq!(fields[11], "cDNA.pos / cDNA.length");
    assert_eq!(fields[15], "ERRORS / WARNINGS / INFO");
}

#[test]
fn field_names_are_sanitised_for_column_headers() {
    assert_eq!(sanitize_field_name("HGVS.c"), "HGVS_c");
    assert_eq!(
        sanitize_field_name("cDNA.pos / cDNA.length"),
        "cDNA_pos___cDNA_length"
    );
    assert_eq!(
        sanitize_field_name("ERRORS / WARNINGS / INFO"),
        "ERRORS___WARNINGS___INFO"
    );
    assert_eq!(sanitize_field_name("Gene_Name"), "Gene_Name");
    assert_eq!(sanitize_field_name(""), "");
}

#[test]
fn ann_impact_ranks_high_to_modifier_case_and_space_insensitively() {
    assert_eq!(get_ann_impact_severity("HIGH"), 4);
    assert_eq!(get_ann_impact_severity("MODERATE"), 3);
    assert_eq!(get_ann_impact_severity("LOW"), 2);
    assert_eq!(get_ann_impact_severity("MODIFIER"), 1);
    assert_eq!(get_ann_impact_severity("UNKNOWN"), 0);
    assert_eq!(get_ann_impact_severity("high"), 4);
    assert_eq!(get_ann_impact_severity(" HIGH "), 4);
}

// --- SnpEff ANN through from_vcf_line ---------------------------------------------------------

#[test]
fn ann_fields_are_parsed_with_sanitised_keys() {
    let line = "chr1\t100\t.\tA\tG\t60\tPASS\tANN=G|missense_variant|MODERATE|BRCA1|ENSG00000012048|transcript|ENST00000357654|protein_coding|5/24|c.181T>C|p.Cys61Arg|181/5592|181/4863|61/1620||\tGT\t0/1";
    let columns = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE1",
    ];
    let records = ReformattedVcfRecord::from_vcf_line(
        line,
        &columns,
        &None,
        &ann_names(),
        TranscriptHandling::FirstOnly,
    )
    .unwrap();
    assert_eq!(records.len(), 1);
    let f = &records[0].info_fields;
    assert_eq!(f.get("ANN_Gene_Name").unwrap(), "BRCA1");
    assert_eq!(f.get("ANN_Feature_ID").unwrap(), "ENST00000357654");
    assert_eq!(f.get("ANN_HGVS_c").unwrap(), "c.181T>C");
    assert_eq!(f.get("ANN_HGVS_p").unwrap(), "p.Cys61Arg");
}

#[test]
fn ann_multiple_transcripts_split_into_rows() {
    let line = "chr1\t100\t.\tA\tG\t60\tPASS\tANN=G|missense_variant|MODERATE|BRCA1|ENSG00000012048|transcript|ENST00000357654|protein_coding|5/24|c.181T>C|p.Cys61Arg|||||,G|synonymous_variant|LOW|BRCA1|ENSG00000012048|transcript|ENST00000123456|protein_coding|6/25|c.200A>G|p.Leu67Leu|||||";
    let columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];
    let records = ReformattedVcfRecord::from_vcf_line(
        line,
        &columns,
        &None,
        &ann_names(),
        TranscriptHandling::SplitRows,
    )
    .unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(
        records[0].info_fields.get("ANN_Feature_ID").unwrap(),
        "ENST00000357654"
    );
    assert_eq!(
        records[1].info_fields.get("ANN_Annotation").unwrap(),
        "synonymous_variant"
    );
}

#[test]
fn ann_most_severe_picks_the_highest_impact() {
    let line = "chr1\t100\t.\tA\tG\t60\tPASS\tANN=G|synonymous_variant|LOW|BRCA1||||||||||||,G|missense_variant|MODERATE|BRCA1||||||||||||,G|stop_gained|HIGH|BRCA1||||||||||||";
    let columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];
    let records = ReformattedVcfRecord::from_vcf_line(
        line,
        &columns,
        &None,
        &ann_names(),
        TranscriptHandling::MostSevere,
    )
    .unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(
        records[0].info_fields.get("ANN_Annotation").unwrap(),
        "stop_gained"
    );
}

#[test]
fn ann_coexists_with_other_info_keys() {
    let line = "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50;AF=0.25;ANN=G|missense_variant|MODERATE|BRCA1||||||||||||;AC=2;AN=4";
    let columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];
    let records = ReformattedVcfRecord::from_vcf_line(
        line,
        &columns,
        &None,
        &ann_names(),
        TranscriptHandling::FirstOnly,
    )
    .unwrap();
    let f = &records[0].info_fields;
    assert_eq!(f.get("INFO_DP").unwrap(), "50");
    assert_eq!(f.get("INFO_AC").unwrap(), "2");
    assert_eq!(f.get("INFO_AN").unwrap(), "4");
    assert_eq!(f.get("ANN_Gene_Name").unwrap(), "BRCA1");
}

#[test]
fn empty_ann_produces_no_ann_columns() {
    let line = "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50;ANN=;AF=0.25";
    let columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];
    let records = ReformattedVcfRecord::from_vcf_line(
        line,
        &columns,
        &None,
        &ann_names(),
        TranscriptHandling::FirstOnly,
    )
    .unwrap();
    let f = &records[0].info_fields;
    assert_eq!(f.get("INFO_DP").unwrap(), "50");
    assert_eq!(f.get("INFO_AF").unwrap(), "0.25");
    assert!(!f.contains_key("ANN_Allele"));
}

#[test]
fn snpeff_file_end_to_end() {
    let columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    let lines = [
        "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50;ANN=G|missense_variant|MODERATE|BRCA1|ENSG00000012048|transcript|ENST00000357654|protein_coding|5/24|c.181T>C|p.Cys61Arg|181/5592|181/4863|61/1620||",
        "chr2\t200\t.\tC\tT\t40\tPASS\tDP=30;ANN=T|stop_gained|HIGH|TP53|ENSG00000141510|transcript|ENST00000269305|protein_coding|7/11|c.916C>T|p.Arg306Ter|916/1182|916/1182|306/393||",
    ]
    .map(String::from);
    let (headers, records) = reformat_vcf_data_with_header(
        SNPEFF_HEADER,
        columns,
        &lines,
        TranscriptHandling::FirstOnly,
    )
    .unwrap();
    for h in ["INFO_DP", "ANN_Allele", "ANN_Gene_Name", "ANN_HGVS_c"] {
        assert!(headers.contains(&h.to_string()), "missing {h}");
    }
    assert_eq!(records.len(), 2);
    assert_eq!(records[1].info_fields.get("ANN_Gene_Name").unwrap(), "TP53");
    assert_eq!(
        records[1].info_fields.get("ANN_Annotation").unwrap(),
        "stop_gained"
    );
}
