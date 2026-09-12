//! # VCF Reformatting Module
//!
//! This module provides functionality to parse and reformat VCF (Variant Call Format) files,
//! with special support for VEP (Variant Effect Predictor) and SnpEff annotations.
//!
//! The main purpose is to convert complex VCF files into flattened tab-separated format
//! that's easier to work with in downstream analysis pipelines.
//!
//! ## Key Features
//!
//! - Parse VEP CSQ and SnpEff ANN annotations
//! - Handle multiple transcripts per variant (first-only, most-severe, or split-rows)
//! - Parallel processing support for large files
//! - Flexible output formatting
//! ```
use crate::essentials_fields::{biotype_priority, MafRecord};
use crate::extract_ann_and_ann_names::extract_ann_regex;
use crate::extract_csq_and_csq_names::extract_csq_regex;
use crate::extract_sample_info::ParsedFormatSample;
use crate::get_info_from_header::{extract_ann_format_from_header, extract_csq_format_from_header};
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::{create_dir_all, File};
use std::io::{BufWriter, Write};
use std::path::Path;

/// Specifies which type of variant annotation to parse from VCF files
///
/// Different annotation tools produce different INFO field formats:
/// - VEP produces CSQ fields
/// - SnpEff produces ANN fields
/// - Auto-detection tries both
#[derive(Debug, Clone, Copy)]
pub enum AnnotationType {
    Vep,
    SnpEff,
    Auto,
}
/// Internal enum to track which annotation field type was found during parsing
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnotationFieldType {
    Csq,
    Ann,
    None,
}

/// Result of parsing annotation fields from a VCF INFO column
/// This internal structure separates annotation data from other INFO fields
/// for more efficient processing.
#[allow(dead_code)]
#[derive(Debug, Clone)]
struct AnnotationParseResult {
    field_type: AnnotationFieldType,
    records: Vec<HashMap<String, String>>,
    remaining_info: String,
}

/// Defines how to handle variants with multiple transcript annotations
///
/// Many variants affect multiple transcripts of the same gene. This enum
/// controls how those multiple annotations are processed.
///
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TranscriptHandling {
    MostSevere,
    FirstOnly,
    SplitRows,
}

/// Defines how to handle variants with multiple transcript annotations
///
/// Many variants affect multiple transcripts of the same gene. This enum
/// controls how those multiple annotations are processed.

#[derive(Debug, Clone)]
pub struct ReformattedVcfRecord {
    pub chromosome: String,
    pub position: u64,
    pub id: Option<String>,
    pub reference: String,
    pub alternate: String,
    pub quality: Option<f64>,
    pub filter: String,
    pub info_fields: HashMap<String, String>,
    pub format_sample_data: Option<ParsedFormatSample>,
    pub annotation_field_type: AnnotationFieldType,
}

impl ReformattedVcfRecord {
    /// Parse a VCF data line into one or more reformatted records
    ///
    /// This is the main entry point for converting a raw VCF line into
    /// structured data. It handles annotation parsing and can generate
    /// multiple output records if multiple transcripts are present.
    ///
    /// # Arguments
    ///
    /// * `line` - A tab-separated VCF data line
    /// * `column_names` - Column headers from the VCF file
    /// * `csq_field_names` - VEP CSQ field names from header, if available
    /// * `ann_field_names` - SnpEff ANN field names from header, if available
    /// * `transcript_handling` - How to handle multiple transcripts
    ///
    /// # Returns
    ///
    /// A vector of reformatted records. Usually contains one record, but may
    /// contain multiple if `TranscriptHandling::SplitRows` is used.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The line has fewer than 8 required VCF columns
    /// - Position or quality fields contain invalid numbers
    /// - Annotation parsing fails
    pub fn from_vcf_line(
        line: &str,
        column_names: &[&str],
        csq_field_names: &Option<Vec<String>>,
        ann_field_names: &Option<Vec<String>>,
        transcript_handling: TranscriptHandling,
    ) -> std::result::Result<Vec<Self>, Box<dyn std::error::Error>> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 8 {
            return Err("Invalid VCF line: insufficient fields".into());
        }

        let chromosome = fields[0].to_string();
        let position: u64 = fields[1]
            .parse()
            .map_err(|e| format!("Invalid position '{}': {}", fields[1], e))?;
        let id = if fields[2] == "." {
            None
        } else {
            Some(fields[2].to_string())
        };
        let reference = fields[3].to_string();
        let alternate = fields[4].to_string();
        let quality = if fields[5] == "." {
            None
        } else {
            Some(
                fields[5]
                    .parse()
                    .map_err(|e| format!("Invalid quality '{}': {}", fields[5], e))?,
            )
        };
        let filter = fields[6].to_string();
        let info = fields[7];

        let format_sample_data = if fields.len() > 8 {
            Some(ParsedFormatSample::from_vcf_fields(
                &fields[8..],
                column_names,
            )?)
        } else {
            None
        };

        let (info_variants, annotation_field_type) =
            parse_info_field(info, csq_field_names, ann_field_names, transcript_handling)?;

        let len = info_variants.len();
        let records: Vec<Self> = if len == 1 {
            // Common case: single record, move values without cloning
            vec![Self {
                chromosome,
                position,
                id,
                reference,
                alternate,
                quality,
                filter,
                info_fields: info_variants.into_iter().next().unwrap(),
                format_sample_data,
                annotation_field_type,
            }]
        } else {
            // Multiple records: clone for all but the last, move for the last
            let mut records = Vec::with_capacity(len);
            let mut variants_iter = info_variants.into_iter().peekable();
            while let Some(info_fields) = variants_iter.next() {
                if variants_iter.peek().is_some() {
                    records.push(Self {
                        chromosome: chromosome.clone(),
                        position,
                        id: id.clone(),
                        reference: reference.clone(),
                        alternate: alternate.clone(),
                        quality,
                        filter: filter.clone(),
                        info_fields,
                        format_sample_data: format_sample_data.clone(),
                        annotation_field_type,
                    });
                } else {
                    records.push(Self {
                        chromosome,
                        position,
                        id,
                        reference,
                        alternate,
                        quality,
                        filter,
                        info_fields,
                        format_sample_data,
                        annotation_field_type,
                    });
                    break;
                }
            }
            records
        };

        Ok(records)
    }
}

/// One flattened record per transcript, plus which annotation field they came from.
type ParsedInfoField = (Vec<HashMap<String, String>>, AnnotationFieldType);

/// Parse the INFO field of a VCF record, extracting and processing annotations
///
/// This function separates annotation data (CSQ/ANN) from standard INFO fields,
/// processes the annotations according to the specified transcript handling,
/// and combines everything into structured output records.
///
/// # Arguments
///
/// * `info` - The INFO field string from a VCF record
/// * `csq_field_names` - VEP CSQ field names, if available
/// * `ann_field_names` - SnpEff ANN field names, if available
/// * `transcript_handling` - How to handle multiple transcripts
///
/// # Returns
///
/// A vector of HashMaps, each representing a flattened record with all
/// INFO and annotation fields as key-value pairs.
///
pub fn parse_info_field(
    info: &str,
    csq_field_names: &Option<Vec<String>>,
    ann_field_names: &Option<Vec<String>>,
    transcript_handling: TranscriptHandling,
) -> Result<ParsedInfoField, Box<dyn std::error::Error>> {
    if info.is_empty() {
        return Ok((vec![HashMap::new()], AnnotationFieldType::None));
    }

    let annotation_result =
        parse_annotation_fields(info, csq_field_names, ann_field_names, transcript_handling)?;

    let field_type = annotation_result.field_type;
    let remaining_info_map = parse_remaining_info_fields(&annotation_result.remaining_info)?;

    if annotation_result.records.is_empty() {
        Ok((vec![remaining_info_map], field_type))
    } else {
        let combined_records =
            combine_annotation_with_info(annotation_result.records, remaining_info_map);
        Ok((combined_records, field_type))
    }
}

/// Parse annotation fields (CSQ or ANN) from an INFO string
///
/// This internal function tries to extract VEP CSQ annotations first,
/// then falls back to SnpEff ANN annotations if CSQ is not found.
fn parse_annotation_fields(
    info: &str,
    csq_field_names: &Option<Vec<String>>,
    ann_field_names: &Option<Vec<String>>,
    transcript_handling: TranscriptHandling,
) -> Result<AnnotationParseResult, Box<dyn std::error::Error>> {
    let mut parsed_lines = create_dummy_vcf_line(info);

    // Try CSQ first (VEP annotations)
    if let Some(csq_value) = extract_csq_regex(&mut parsed_lines) {
        if let Some(field_names) = csq_field_names {
            if !field_names.is_empty() && !csq_value.trim().is_empty() {
                match parse_annotation_field_with_handling(
                    "CSQ",
                    &csq_value,
                    field_names,
                    transcript_handling,
                    find_most_severe_consequence,
                ) {
                    Ok(records) if !records.is_empty() => {
                        return Ok(AnnotationParseResult {
                            field_type: AnnotationFieldType::Csq,
                            records,
                            remaining_info: std::mem::take(&mut parsed_lines[7]),
                        });
                    }
                    Ok(_) => {}
                    Err(e) => {
                        eprintln!("Warning: Failed to parse CSQ field: {e}");
                    }
                }
            }
        }
    }

    // Reset and try ANN (SnpEff annotations)
    parsed_lines[7] = info.to_string();
    if let Some(ann_value) = extract_ann_regex(&mut parsed_lines) {
        if let Some(field_names) = ann_field_names {
            if !field_names.is_empty() && !ann_value.trim().is_empty() {
                match parse_annotation_field_with_handling(
                    "ANN",
                    &ann_value,
                    field_names,
                    transcript_handling,
                    find_most_severe_ann_consequence,
                ) {
                    Ok(records) if !records.is_empty() => {
                        return Ok(AnnotationParseResult {
                            field_type: AnnotationFieldType::Ann,
                            records,
                            remaining_info: std::mem::take(&mut parsed_lines[7]),
                        });
                    }
                    Ok(_) => {}
                    Err(e) => {
                        eprintln!("Warning: Failed to parse ANN field: {e}");
                    }
                }
            }
        }
    }

    // No annotation fields found
    Ok(AnnotationParseResult {
        field_type: AnnotationFieldType::None,
        records: Vec::new(),
        remaining_info: info.to_string(),
    })
}

/// Create a dummy VCF line for use with regex extraction functions
///
/// Some extraction functions expect a full VCF line but we only have
/// the INFO field. This creates a minimal valid VCF line.
fn create_dummy_vcf_line(info: &str) -> Vec<String> {
    vec![
        "chr1".to_string(),
        "100".to_string(),
        ".".to_string(),
        "A".to_string(),
        "G".to_string(),
        "60".to_string(),
        "PASS".to_string(),
        info.to_string(),
    ]
}

/// Split a CSQ/ANN value into individual `|`-delimited transcript annotations
/// and reduce them to output records according to `transcript_handling`.
/// `find_most_severe` implements the annotation-type-specific severity ranking
/// (VEP orders by a fixed consequence-term list; SnpEff ranks by IMPACT).
fn parse_annotation_field_with_handling(
    prefix: &str,
    value: &str,
    field_names: &[String],
    transcript_handling: TranscriptHandling,
    find_most_severe: impl Fn(
        &[&str],
        &[String],
    ) -> Result<HashMap<String, String>, Box<dyn std::error::Error>>,
) -> Result<Vec<HashMap<String, String>>, Box<dyn std::error::Error>> {
    if value.trim().is_empty() {
        return Ok(Vec::new());
    }

    let annotations: Vec<&str> = value.split(',').filter(|s| !s.trim().is_empty()).collect();

    if annotations.is_empty() {
        return Ok(Vec::new());
    }

    match transcript_handling {
        TranscriptHandling::FirstOnly => {
            let first_annotation = annotations
                .first()
                .ok_or("No annotations found after filtering")?;
            let parsed = parse_single_annotation(prefix, first_annotation, field_names)?;
            Ok(vec![parsed])
        }
        TranscriptHandling::MostSevere => {
            let most_severe = find_most_severe(&annotations, field_names)?;
            Ok(vec![most_severe])
        }
        TranscriptHandling::SplitRows => {
            let mut all_annotations = Vec::new();
            for annotation in annotations {
                match parse_single_annotation(prefix, annotation, field_names) {
                    Ok(parsed) => all_annotations.push(parsed),
                    Err(e) => {
                        eprintln!(
                            "Warning: Failed to parse {prefix} annotation '{annotation}': {e}"
                        );
                    }
                }
            }
            Ok(all_annotations)
        }
    }
}

/// Parse a single `|`-delimited CSQ/ANN annotation into a `{prefix}_{field}`
/// map, e.g. `CSQ_Consequence` or `ANN_Gene_Name`. Values beyond the known
/// field names are kept under `{prefix}_EXTRA_N`.
fn parse_single_annotation(
    prefix: &str,
    annotation: &str,
    field_names: &[String],
) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    if annotation.trim().is_empty() {
        return Err(format!("Empty {prefix} annotation").into());
    }

    let values: Vec<&str> = annotation.split('|').collect();
    let mut annotation_map = HashMap::new();

    for (i, field_name) in field_names.iter().enumerate() {
        let value = values.get(i).unwrap_or(&"").trim();
        let formatted_value = if value.is_empty() { "." } else { value };
        annotation_map.insert(
            format!("{prefix}_{}", sanitize_field_name(field_name)),
            formatted_value.to_string(),
        );
    }

    if values.len() > field_names.len() {
        for (i, value) in values.iter().enumerate().skip(field_names.len()) {
            annotation_map.insert(
                format!("{prefix}_EXTRA_{}", i - field_names.len() + 1),
                value.trim().to_string(),
            );
        }
    }

    Ok(annotation_map)
}

fn find_most_severe_ann_consequence(
    annotations: &[&str],
    ann_field_names: &[String],
) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    if annotations.is_empty() {
        return Err("No annotations provided".into());
    }

    let impact_index = ann_field_names
        .iter()
        .position(|name| {
            name.to_lowercase().contains("impact")
                || name.to_lowercase().contains("annotation_impact")
        })
        .unwrap_or(2);

    let mut most_severe: Option<&str> = None;
    let mut highest_severity = 0;

    for annotation in annotations {
        let values: Vec<&str> = annotation.split('|').collect();
        if let Some(impact) = values.get(impact_index) {
            let severity = get_ann_impact_severity(impact);
            if severity > highest_severity {
                highest_severity = severity;
                most_severe = Some(annotation);
            }
        }
    }

    let selected_annotation = most_severe.unwrap_or(annotations[0]);
    parse_single_annotation("ANN", selected_annotation, ann_field_names)
}

/// Convert SnpEff impact level to numeric severity score
///
/// SnpEff uses four impact levels to categorize the severity of variants:
/// - HIGH: Variant has high impact (score: 4)
/// - MODERATE: Variant has moderate impact (score: 3)
/// - LOW: Variant has low impact (score: 2)
/// - MODIFIER: Variant is unlikely to change protein behavior (score: 1)
/// - Unknown impacts get score 0
///
/// # Arguments
///
/// * `impact` - The impact string from SnpEff ANN field
///
/// # Returns
///
/// Numeric severity score (0-4, where 4 is most severe)
pub fn get_ann_impact_severity(impact: &str) -> u8 {
    match impact.trim().to_uppercase().as_str() {
        "HIGH" => 4,
        "MODERATE" => 3,
        "LOW" => 2,
        "MODIFIER" => 1,
        _ => 0,
    }
}
/// Parse standard INFO fields (non-annotation) into a HashMap
fn parse_remaining_info_fields(
    remaining_info: &str,
) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    let mut info_map = HashMap::new();

    if remaining_info.trim().is_empty() {
        return Ok(info_map);
    }

    for pair in remaining_info.split(';') {
        let trimmed_pair = pair.trim();
        if trimmed_pair.is_empty() {
            continue;
        }

        match trimmed_pair.split_once('=') {
            Some((key, value)) => {
                let sanitized_key = sanitize_field_name(key.trim());
                let sanitized_value = value.trim();
                if !sanitized_key.is_empty() {
                    info_map.insert(format!("INFO_{sanitized_key}"), sanitized_value.to_string());
                }
            }
            None => {
                let sanitized_key = sanitize_field_name(trimmed_pair);
                if !sanitized_key.is_empty() {
                    info_map.insert(format!("INFO_{sanitized_key}"), "true".to_string());
                }
            }
        }
    }

    Ok(info_map)
}
/// Combine annotation records with standard INFO fields
fn combine_annotation_with_info(
    annotation_records: Vec<HashMap<String, String>>,
    info_fields: HashMap<String, String>,
) -> Vec<HashMap<String, String>> {
    if annotation_records.is_empty() {
        return vec![info_fields];
    }

    annotation_records
        .into_iter()
        .map(|mut annotation_map| {
            for (key, value) in &info_fields {
                annotation_map.insert(key.clone(), value.clone());
            }
            annotation_map
        })
        .collect()
}

/// Sanitize field names to be safe for use as column headers
///
/// Converts special characters to underscores and removes leading/trailing
/// underscores to create valid, clean column names.
///
/// # Arguments
///
/// * `field_name` - The raw field name to sanitize
///
/// # Returns
///
/// A sanitized field name safe for use as a column header
pub fn sanitize_field_name(field_name: &str) -> String {
    field_name
        .chars()
        .map(|c| {
            if c.is_alphanumeric() || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect::<String>()
        .trim_start_matches('_')
        .trim_end_matches('_')
        .to_string()
}

/// Pick the one CSQ annotation vcf2maf would report, porting `vcf2maf.pl:871-894`:
/// sort by transcript biotype, then consequence severity, then longest transcript; then
/// take the worst-affected *gene* and that gene's canonical isoform. Biotype outranks
/// severity, so a milder consequence on a protein_coding transcript beats a worse one on
/// an lncRNA — and the canonical isoform wins even when a sibling isoform is more severe.
fn find_most_severe_consequence(
    annotations: &[&str],
    csq_field_names: &[String],
) -> std::result::Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    if annotations.is_empty() {
        return Err("No annotations provided".into());
    }

    let idx = |name: &str, default: usize| {
        csq_field_names
            .iter()
            .position(|n| n == name)
            .unwrap_or(default)
    };
    let (i_cons, i_symbol, i_biotype, i_canonical, i_cdna) = (
        idx("Consequence", 1),
        idx("SYMBOL", 3),
        idx("BIOTYPE", 7),
        idx("CANONICAL", 23),
        idx("cDNA_position", 12),
    );

    let field = |a: &str, i: usize| a.split('|').nth(i).unwrap_or("").to_string();

    // vcf2maf.pl:862-863 — Transcript_Length is the denominator of cDNA_position.
    let transcript_length = |a: &str| -> u64 {
        field(a, i_cdna)
            .rsplit_once('/')
            .and_then(|(_, len)| len.parse().ok())
            .unwrap_or(0)
    };
    let severity = |a: &str| {
        MafRecord::effect_priority(&MafRecord::resolve_one_consequence(
            &field(a, i_cons).to_lowercase(),
        ))
    };

    // vcf2maf.pl:871-875. Rust's sort_by is stable, as Perl's sort is, so equal-key
    // annotations keep their input order in both.
    let mut sorted: Vec<&str> = annotations.to_vec();
    sorted.sort_by(|a, b| {
        biotype_priority(&field(a, i_biotype))
            .cmp(&biotype_priority(&field(b, i_biotype)))
            .then_with(|| severity(a).cmp(&severity(b)))
            .then_with(|| transcript_length(b).cmp(&transcript_length(a)))
    });

    let has_symbol = |a: &str| !field(a, i_symbol).is_empty();
    let is_canonical = |a: &str| field(a, i_canonical) == "YES";

    // vcf2maf.pl:878-880 — the worst affected GENE, not the worst effect.
    let maf_gene = sorted
        .iter()
        .find(|a| has_symbol(a))
        .map(|a| field(a, i_symbol));

    // vcf2maf.pl:888, then :891, then :893. The two --custom-enst branches (:883, :886)
    // have no equivalent here — this tool exposes no isoform override.
    let selected = maf_gene
        .as_ref()
        .and_then(|gene| {
            sorted
                .iter()
                .find(|a| &field(a, i_symbol) == gene && is_canonical(a))
        })
        .or_else(|| sorted.iter().find(|a| has_symbol(a) && is_canonical(a)))
        .copied()
        .unwrap_or(sorted[0]);

    parse_single_annotation("CSQ", selected, csq_field_names)
}
/// Reformat VCF data with header information for annotation field extraction
///
/// This is the main processing function that takes raw VCF data and converts
/// it into a flattened, tab-separated format suitable for analysis.
///
/// # Arguments
///
/// * `header` - VCF header string containing metadata and field definitions
/// * `column_names` - Column header line from VCF (starts with #CHROM)
/// * `data_lines` - Vector of VCF data lines to process
/// * `transcript_handling` - How to handle multiple transcripts per variant
///
/// # Returns
///
/// A tuple containing:
/// - Vector of column headers for the output
/// - Vector of reformatted VCF records
pub fn reformat_vcf_data_with_header(
    header: &str,
    column_names: &str,
    data_lines: &[String],
    transcript_handling: TranscriptHandling,
) -> std::result::Result<(Vec<String>, Vec<ReformattedVcfRecord>), Box<dyn std::error::Error>> {
    let csq_field_names = extract_csq_format_from_header(header);
    let ann_field_names = extract_ann_format_from_header(header);

    let column_names_vec: Vec<&str> = column_names.trim_start_matches('#').split('\t').collect();

    let mut all_records = Vec::new();

    for (line_num, line) in data_lines.iter().enumerate() {
        match ReformattedVcfRecord::from_vcf_line(
            line,
            &column_names_vec,
            &csq_field_names,
            &ann_field_names,
            transcript_handling,
        ) {
            Ok(mut records) => {
                all_records.append(&mut records);
            }
            Err(e) => {
                eprintln!(
                    "Warning: Failed to parse line {}: {} ({})",
                    line_num + 1,
                    e,
                    line
                );
            }
        }
    }

    let headers = generate_headers_from_records(&all_records, &column_names_vec, header);

    Ok((headers, all_records))
}

/// The `##INFO=<ID=…` / `##FORMAT=<ID=…` ids a VCF header declares, in declaration order.
fn declared_ids(vcf_header: &str, kind: &str) -> Vec<String> {
    let needle = format!("##{kind}=<ID=");
    vcf_header
        .lines()
        .filter_map(|line| line.strip_prefix(needle.as_str()))
        .map(|rest| {
            rest.split([',', '>'])
                .next()
                .unwrap_or("")
                .trim()
                .to_string()
        })
        .filter(|id| !id.is_empty())
        .collect()
}

/// Generate column headers from the first reformatted record, plus every INFO and FORMAT field
/// the VCF header declares that this record happens not to carry.
///
/// Deriving the column list from record #1 alone silently drops any key absent from it — measured
/// on real files: 73% of rows losing `LOF`/`NMD` (SnpEff), 27% losing `PON`/`STR`/`RU` and 34%
/// losing the `PGT`/`PID`/`PS` phasing trio (Mutect2). Rows are rendered by header name, so the
/// cost of a declared-but-unused field is one column of `.`, never a shifted row.
///
/// CSQ/ANN sub-fields are positional and therefore always complete on every record, so only INFO
/// and FORMAT need the header pass.
/// The column list a VCF declares before any variant is seen: the fixed columns, then
/// the CSQ/ANN sub-fields named in the annotation's own `Format:` string. The INFO and
/// sample blocks are added by the caller, which already reads them from the header.
fn headers_from_declarations_only(vcf_header: &str) -> Vec<String> {
    let mut headers: Vec<String> = FIXED_COLUMNS.iter().map(|c| c.to_string()).collect();
    // Gated on the declaration, because extract_ann_format_from_header falls back to a
    // default SnpEff layout: without this, an unannotated VCF grows 16 phantom columns.
    let declared = declared_ids(vcf_header, "INFO");
    for (prefix, names) in [
        ("CSQ", extract_csq_format_from_header(vcf_header)),
        ("ANN", extract_ann_format_from_header(vcf_header)),
    ] {
        if !declared.iter().any(|id| id == prefix) {
            continue;
        }
        for name in names.unwrap_or_default() {
            headers.push(format!("{prefix}_{}", sanitize_field_name(&name)));
        }
    }
    headers
}

fn generate_headers_from_records(
    records: &[ReformattedVcfRecord],
    column_names_vec: &[&str],
    vcf_header: &str,
) -> Vec<String> {
    let sample_names: Vec<String> = if column_names_vec.len() > 9 {
        column_names_vec[9..]
            .iter()
            .map(|s| s.to_string())
            .collect()
    } else {
        vec![]
    };

    // With no records, every column still comes from the VCF header's own declarations,
    // so a variant-free VCF gets a column header rather than a zero-byte file.
    let base = match records.first() {
        Some(first_record) => generate_headers_from_record(first_record, &sample_names),
        None => headers_from_declarations_only(vcf_header),
    };
    let (fixed, rest) = base.split_at(FIXED_COLUMNS.len());

    // INFO block: the union, still alphabetical, so existing columns keep their position.
    let mut info: Vec<String> = declared_ids(vcf_header, "INFO")
        .into_iter()
        .filter(|id| id != "CSQ" && id != "ANN")
        .map(|id| format!("INFO_{}", sanitize_field_name(&id)))
        .collect();
    info.extend(rest.iter().filter(|h| h.starts_with("INFO_")).cloned());
    info.sort();
    info.dedup();

    let annotation: Vec<String> = rest
        .iter()
        .filter(|h| h.starts_with("CSQ_") || h.starts_with("ANN_"))
        .cloned()
        .collect();

    // Sample block: this record's FORMAT keys first, then the declared ones it lacks.
    let declared_format = declared_ids(vcf_header, "FORMAT");
    let sample_block: Vec<&String> = rest
        .iter()
        .filter(|h| !h.starts_with("INFO_") && !h.starts_with("CSQ_") && !h.starts_with("ANN_"))
        .collect();
    let mut samples = Vec::new();
    for name in &sample_names {
        let prefix = format!("{name}_");
        let existing: Vec<String> = sample_block
            .iter()
            .filter(|h| h.starts_with(&prefix))
            .map(|h| (*h).clone())
            .collect();
        samples.extend(existing.iter().cloned());
        for key in &declared_format {
            let column = format!("{prefix}{key}");
            if !existing.contains(&column) {
                samples.push(column);
            }
        }
    }

    let mut headers = fixed.to_vec();
    headers.extend(info);
    headers.extend(annotation);
    headers.extend(samples);
    headers
}
/// Parallel version of VCF data reformatting for improved performance on large files
///
/// This function works identically to `reformat_vcf_data_with_header` but uses
/// parallel processing via Rayon for better performance on multi-core systems.
///
/// # Arguments
///
/// * `header` - VCF header string containing metadata and field definitions
/// * `column_names` - Column header line from VCF (starts with #CHROM)
/// * `data_lines` - Vector of VCF data lines to process
/// * `transcript_handling` - How to handle multiple transcripts per variant
///
/// # Returns
///
/// A tuple containing:
/// - Vector of column headers for the output
/// - Vector of reformatted VCF records
///
/// # Performance
///
/// Use this function for files with >10,000 variants. The parallel processing
/// overhead isn't worth it for smaller files.
///
pub fn reformat_vcf_data_with_header_parallel(
    header: &str,
    column_names: &str,
    data_lines: &[String],
    transcript_handling: TranscriptHandling,
) -> std::result::Result<(Vec<String>, Vec<ReformattedVcfRecord>), Box<dyn std::error::Error>> {
    let csq_field_names = extract_csq_format_from_header(header);
    let ann_field_names = extract_ann_format_from_header(header);

    let column_names_vec: Vec<&str> = column_names.trim_start_matches('#').split('\t').collect();

    // Fixed: Collect into Vec<Vec<ReformattedVcfRecord>> first, then flatten
    let all_results: Vec<Vec<ReformattedVcfRecord>> = data_lines
        .par_iter()
        .enumerate()
        .map(|(line_num, line)| {
            ReformattedVcfRecord::from_vcf_line(
                line,
                &column_names_vec,
                &csq_field_names,
                &ann_field_names,
                transcript_handling,
            )
            .unwrap_or_else(|e| {
                eprintln!(
                    "Warning: Failed to parse line {}: {} ({})",
                    line_num + 1,
                    e,
                    line
                );
                Vec::new()
            })
        })
        .collect();

    let mut flattened_records = Vec::new();
    for mut records in all_results {
        flattened_records.append(&mut records);
    }

    let headers = generate_headers_from_records(&flattened_records, &column_names_vec, header);

    Ok((headers, flattened_records))
}

// FIXED: Update to handle both CSQ and ANN headers
const FIXED_COLUMNS: [&str; 7] = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"];

fn generate_headers_from_record(
    record: &ReformattedVcfRecord,
    _sample_names: &[String],
) -> Vec<String> {
    let mut headers: Vec<String> = FIXED_COLUMNS.iter().map(|c| c.to_string()).collect();

    let mut info_keys: Vec<String> = record
        .info_fields
        .keys()
        .filter(|k| k.starts_with("INFO_"))
        .cloned()
        .collect();
    info_keys.sort();
    headers.extend(info_keys);

    let mut csq_keys: Vec<String> = record
        .info_fields
        .keys()
        .filter(|k| k.starts_with("CSQ_"))
        .cloned()
        .collect();
    csq_keys.sort();
    headers.extend(csq_keys);

    let mut ann_keys: Vec<String> = record
        .info_fields
        .keys()
        .filter(|k| k.starts_with("ANN_"))
        .cloned()
        .collect();
    ann_keys.sort();
    headers.extend(ann_keys);

    if let Some(ref sample_data) = record.format_sample_data {
        headers.extend(sample_data.get_headers_for_samples());
    }

    headers
}
/// Write reformatted VCF records to a TSV file with optional compression
///
/// This function outputs the reformatted data in tab-separated format,
/// optionally compressing the output with gzip.
///
/// # Arguments
///
/// * `filename` - Output file path
/// * `headers` - Column headers for the output
/// * `records` - Reformatted VCF records to write
/// * `compress` - Whether to compress output with gzip
// Superseded by the streaming path; kept because tests/test.rs still exercises it.
#[allow(dead_code)]
pub fn write_reformatted_vcf(
    filename: &str,
    headers: &[String],
    records: &[ReformattedVcfRecord],
    compress: bool,
) -> std::io::Result<()> {
    if let Some(parent) = Path::new(filename).parent() {
        create_dir_all(parent)?;
    }

    let file = File::create(filename)?;

    if compress {
        let encoder = GzEncoder::new(file, Compression::default());
        let mut writer = BufWriter::new(encoder);
        write_tsv_content(&mut writer, headers, records)?;
        writer.flush()?;
    } else {
        let mut writer = BufWriter::new(file);
        write_tsv_content(&mut writer, headers, records)?;
        writer.flush()?;
    }

    Ok(())
}

#[allow(clippy::collapsible_else_if)]
fn write_tsv_content<W: Write>(
    writer: &mut W,
    headers: &[String],
    records: &[ReformattedVcfRecord],
) -> std::io::Result<()> {
    write_tsv_header(writer, headers)?;
    write_tsv_rows(writer, headers, records)
}

/// The column line. Streaming writers call this once, before the first batch.
pub fn write_tsv_header<W: Write>(writer: &mut W, headers: &[String]) -> std::io::Result<()> {
    writeln!(writer, "{}", headers.join("\t"))
}

/// Rows only, no column line — so a caller can write one batch at a time and never hold the
/// whole file. This is the single rendering both the buffered and the streaming path use.
pub fn write_tsv_rows<W: Write>(
    writer: &mut W,
    headers: &[String],
    records: &[ReformattedVcfRecord],
) -> std::io::Result<()> {
    let dot = ".";

    for record in records {
        let mut row: Vec<Cow<str>> = Vec::with_capacity(headers.len());

        for header in headers {
            let value: Cow<str> = match header.as_str() {
                "CHROM" => Cow::Borrowed(record.chromosome.as_str()),
                "POS" => Cow::Owned(record.position.to_string()),
                "ID" => Cow::Borrowed(record.id.as_deref().unwrap_or(dot)),
                "REF" => Cow::Borrowed(record.reference.as_str()),
                "ALT" => Cow::Borrowed(record.alternate.as_str()),
                "QUAL" => match record.quality {
                    Some(q) => Cow::Owned(q.to_string()),
                    None => Cow::Borrowed(dot),
                },
                "FILTER" => Cow::Borrowed(record.filter.as_str()),
                _ => {
                    if header.starts_with("INFO_")
                        || header.starts_with("CSQ_")
                        || header.starts_with("ANN_")
                    {
                        match record.info_fields.get(header) {
                            Some(v) => Cow::Borrowed(v.as_str()),
                            None => Cow::Borrowed(dot),
                        }
                    } else {
                        if let Some(ref sample_data) = record.format_sample_data {
                            let mut found_value: Option<&str> = None;

                            for sample in &sample_data.samples {
                                for format_key in &sample_data.format_keys {
                                    let expected_header =
                                        format!("{}_{}", sample.sample_name, format_key);
                                    if expected_header == *header {
                                        found_value = sample
                                            .format_fields
                                            .get(format_key)
                                            .map(|s| s.as_str());
                                        break;
                                    }
                                }
                                if found_value.is_some() {
                                    break;
                                }
                            }

                            Cow::Borrowed(found_value.unwrap_or(dot))
                        } else {
                            Cow::Borrowed(dot)
                        }
                    }
                }
            };
            row.push(value);
        }

        writeln!(writer, "{}", row.join("\t"))?;
    }

    Ok(())
}
// Kept for the library API — `tests/test.rs` exercises it. The binary now streams instead.
#[allow(dead_code)]
pub fn write_maf_file(
    filename: &str,
    records: &[MafRecord],
    compress: bool,
) -> std::io::Result<()> {
    if compress {
        let file = std::fs::File::create(filename)?;
        let mut encoder = GzEncoder::new(file, Compression::default());
        write_maf_content(&mut encoder, records)?;
        encoder.finish()?;
    } else {
        let mut file = std::fs::File::create(filename)?;
        write_maf_content(&mut file, records)?;
    }
    Ok(())
}

/// A MAF file being written a batch at a time, so the caller never has to hold every record.
/// The header goes out on `create`; `finish` is required for the gzip trailer.
pub enum MafWriter {
    Plain(std::io::BufWriter<std::fs::File>),
    Gz(Box<GzEncoder<std::io::BufWriter<std::fs::File>>>),
}

impl MafWriter {
    pub fn create(filename: &str, compress: bool) -> std::io::Result<Self> {
        let file = std::io::BufWriter::new(std::fs::File::create(filename)?);
        let mut writer = if compress {
            MafWriter::Gz(Box::new(GzEncoder::new(file, Compression::default())))
        } else {
            MafWriter::Plain(file)
        };
        let headers = MafRecord::get_maf_headers();
        writeln!(writer.inner(), "{}", headers.join("\t"))?;
        Ok(writer)
    }

    fn inner(&mut self) -> &mut dyn Write {
        match self {
            MafWriter::Plain(w) => w,
            MafWriter::Gz(w) => w.as_mut(),
        }
    }

    pub fn write_rows(&mut self, records: &[MafRecord]) -> std::io::Result<()> {
        let writer = self.inner();
        for record in records {
            writeln!(writer, "{}", record.to_tsv_line())?;
        }
        Ok(())
    }

    pub fn finish(self) -> std::io::Result<()> {
        match self {
            MafWriter::Plain(mut w) => w.flush(),
            MafWriter::Gz(w) => w.finish().map(|mut f| f.flush()).and_then(|r| r),
        }
    }
}

// Helper function to write MAF content
#[allow(dead_code)]
fn write_maf_content<W: Write>(writer: &mut W, records: &[MafRecord]) -> std::io::Result<()> {
    // Write MAF header
    let headers = MafRecord::get_maf_headers();
    writeln!(writer, "{}", headers.join("\t"))?;

    // Write MAF records
    for record in records {
        writeln!(writer, "{}", record.to_tsv_line())?;
    }

    Ok(())
}

/// Process VCF data in chunks to avoid memory exhaustion on large files
// Superseded by the streaming path; kept because tests/test.rs still exercises it.
#[allow(dead_code)]
pub fn reformat_vcf_data_with_header_parallel_chunked(
    header: &str,
    column_names: &str,
    data_lines: &[String],
    transcript_handling: TranscriptHandling,
    output_writer: &mut dyn Write,
) -> std::result::Result<Vec<String>, Box<dyn std::error::Error>> {
    let csq_field_names = extract_csq_format_from_header(header);
    let ann_field_names = extract_ann_format_from_header(header);
    let column_names_vec: Vec<&str> = column_names.trim_start_matches('#').split('\t').collect();

    // Calculate chunk size
    let chunk_size = if data_lines.len() > 1_000_000 {
        50_000
    } else if data_lines.len() > 100_000 {
        100_000
    } else {
        data_lines.len()
    };

    let mut headers_generated = false;
    let mut output_headers: Vec<String> = Vec::new();
    let mut total_processed = 0usize;

    println!(
        "🔄 Processing {} lines in chunks of {}",
        data_lines.len(),
        chunk_size
    );

    // Process each chunk and stream output immediately
    for (chunk_idx, chunk) in data_lines.chunks(chunk_size).enumerate() {
        // Process chunk in parallel
        let chunk_results: Vec<Vec<ReformattedVcfRecord>> = chunk
            .par_iter()
            .enumerate()
            .map(|(line_num, line)| {
                ReformattedVcfRecord::from_vcf_line(
                    line,
                    &column_names_vec,
                    &csq_field_names,
                    &ann_field_names,
                    transcript_handling,
                )
                .unwrap_or_else(|e| {
                    let global_line_num = chunk_idx * chunk_size + line_num + 1;
                    eprintln!(
                        "Warning: Failed to parse line {}: {} ({})",
                        global_line_num, e, line
                    );
                    Vec::new()
                })
            })
            .collect();

        // Flatten this chunk's results
        let chunk_records: Vec<ReformattedVcfRecord> =
            chunk_results.into_iter().flatten().collect();

        // Generate headers from first non-empty chunk only
        if !headers_generated && !chunk_records.is_empty() {
            output_headers =
                generate_headers_from_records(&chunk_records, &column_names_vec, header);

            // Write headers to output
            writeln!(output_writer, "{}", output_headers.join("\t"))?;
            headers_generated = true;

            println!("📋 Generated {} column headers", output_headers.len());
        }

        // Stream each record immediately (NO ACCUMULATION!)
        for record in chunk_records {
            let values = extract_values_from_record(&record, &output_headers);
            writeln!(output_writer, "{}", values.join("\t"))?;
        }

        total_processed += chunk.len();

        // Progress logging every 100k lines
        if total_processed.is_multiple_of(100_000) {
            println!("   📊 Streamed {} lines so far...", total_processed);
        }
    }

    println!(
        "✅ Streaming complete! Processed {} total lines",
        total_processed
    );
    Ok(output_headers)
}

/// Extract values from a record in the same order as headers
fn extract_values_from_record<'a>(
    record: &'a ReformattedVcfRecord,
    headers: &[String],
) -> Vec<Cow<'a, str>> {
    let dot = ".";
    headers
        .iter()
        .map(|header| match header.as_str() {
            "CHROM" => Cow::Borrowed(record.chromosome.as_str()),
            "POS" => Cow::Owned(record.position.to_string()),
            "ID" => Cow::Borrowed(record.id.as_deref().unwrap_or(dot)),
            "REF" => Cow::Borrowed(record.reference.as_str()),
            "ALT" => Cow::Borrowed(record.alternate.as_str()),
            "QUAL" => match record.quality {
                Some(q) => Cow::Owned(q.to_string()),
                None => Cow::Borrowed(dot),
            },
            "FILTER" => Cow::Borrowed(record.filter.as_str()),
            _ => {
                if let Some(value) = record.info_fields.get(header) {
                    Cow::Borrowed(value.as_str())
                } else if let Some(sample_data) = &record.format_sample_data {
                    extract_sample_value_for_header_cow(sample_data, header)
                } else {
                    Cow::Borrowed(dot)
                }
            }
        })
        .collect()
}

/// Helper function to extract sample values by header name, returning Cow to avoid cloning
fn extract_sample_value_for_header_cow<'a>(
    sample_data: &'a ParsedFormatSample,
    header: &str,
) -> Cow<'a, str> {
    for sample in &sample_data.samples {
        for format_key in &sample_data.format_keys {
            let expected_header = format!("{}_{}", sample.sample_name, format_key);
            if expected_header == header {
                return match sample.format_fields.get(format_key) {
                    Some(v) => Cow::Borrowed(v.as_str()),
                    None => Cow::Borrowed("."),
                };
            }
        }
    }
    Cow::Borrowed(".")
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A VCF with a header and no variants must still produce a column header, built
    /// from the declarations alone, or the run writes a zero-byte file no reader opens.
    #[test]
    fn test_headers_come_from_declarations_when_there_are_no_records() {
        let header = concat!(
            "##fileformat=VCFv4.2\n",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n",
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Format: Allele|Consequence|SYMBOL\">\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n",
        );
        let columns: Vec<&str> = "CHROM POS ID REF ALT QUAL FILTER INFO FORMAT TUMOR"
            .split(' ')
            .collect();

        let headers = generate_headers_from_records(&[], &columns, header);
        assert_eq!(
            headers,
            vec![
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
                "CSQ_SYMBOL",
                "TUMOR_GT",
            ]
        );
    }

    /// The ANN extractor falls back to a default SnpEff layout when the header declares
    /// none, so the no-record path must gate on the declaration or invent 16 columns.
    #[test]
    fn test_no_records_and_no_annotation_declared_yields_only_fixed_columns() {
        let header = "##fileformat=VCFv4.2\n";
        let columns: Vec<&str> = "CHROM POS ID REF ALT QUAL FILTER INFO".split(' ').collect();

        let headers = generate_headers_from_records(&[], &columns, header);
        assert_eq!(headers, FIXED_COLUMNS.to_vec());
    }

    fn csq_fields() -> Vec<String> {
        [
            "Consequence",
            "SYMBOL",
            "BIOTYPE",
            "CANONICAL",
            "cDNA_position",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect()
    }

    #[test]
    fn most_severe_ranks_terms_the_old_36_term_list_never_had() {
        // splice_donor_5th_base_variant is in EFFECT_PRIORITY but was absent from the
        // hand-written list this function used to carry, so an annotation holding only it
        // scored nothing and could never be selected.
        let annotations = vec![
            "intron_variant|GENEA|protein_coding|YES|100/1000",
            "splice_donor_5th_base_variant|GENEB|protein_coding|YES|200/2000",
        ];
        let picked = find_most_severe_consequence(&annotations, &csq_fields()).unwrap();
        assert_eq!(picked.get("CSQ_SYMBOL").map(String::as_str), Some("GENEB"));
    }

    #[test]
    fn consequence_terms_are_ranked_case_insensitively() {
        // VEP writes NMD_transcript_variant and TFBS_ablation with capitals; the ported
        // table is lowercase, so an unlowered lookup would score them as unknown.
        let annotations = vec![
            "NMD_transcript_variant|GENEA|protein_coding|YES|100/1000",
            "downstream_gene_variant|GENEB|protein_coding|YES|200/2000",
        ];
        let picked = find_most_severe_consequence(&annotations, &csq_fields()).unwrap();
        assert_eq!(picked.get("CSQ_SYMBOL").map(String::as_str), Some("GENEA"));
    }

    #[test]
    fn biotype_outranks_severity_when_choosing_a_transcript() {
        // vcf2maf.pl:871-875 sorts on biotype BEFORE effect. A more severe consequence on a
        // lncRNA loses to a milder one on a protein_coding transcript.
        let annotations = vec![
            "non_coding_transcript_exon_variant|LAMTOR5-AS1|lncRNA|YES|100/1000",
            "intron_variant|LAMTOR5|protein_coding|YES|200/2000",
        ];
        let picked = find_most_severe_consequence(&annotations, &csq_fields()).unwrap();
        assert_eq!(
            picked.get("CSQ_SYMBOL").map(String::as_str),
            Some("LAMTOR5")
        );
    }

    #[test]
    fn canonical_isoform_of_the_worst_gene_wins_over_a_worse_noncanonical_one() {
        // vcf2maf.pl:878-891: pick the worst affected GENE first, then that gene's
        // CANONICAL=YES isoform — even though the non-canonical isoform is more severe.
        let annotations = vec![
            "stop_gained|GENEA|protein_coding||100/1000",
            "missense_variant|GENEA|protein_coding|YES|200/2000",
        ];
        let picked = find_most_severe_consequence(&annotations, &csq_fields()).unwrap();
        assert_eq!(
            picked.get("CSQ_Consequence").map(String::as_str),
            Some("missense_variant")
        );
    }

    #[test]
    fn longest_transcript_breaks_a_biotype_and_severity_tie() {
        // Third sort key, vcf2maf.pl:874 — Transcript_Length is the cDNA_position
        // denominator (vcf2maf.pl:862-863), descending.
        let annotations = vec![
            "missense_variant|GENEA|protein_coding|YES|10/500",
            "missense_variant|GENEB|protein_coding|YES|10/5000",
        ];
        let picked = find_most_severe_consequence(&annotations, &csq_fields()).unwrap();
        assert_eq!(picked.get("CSQ_SYMBOL").map(String::as_str), Some("GENEB"));
    }

    #[test]
    fn falls_back_to_the_worst_effect_when_no_annotation_has_a_symbol() {
        // vcf2maf.pl:893 — $all_effects[0] after the sort, when nothing has a SYMBOL.
        let annotations = vec![
            "intron_variant||lncRNA||100/1000",
            "stop_gained||protein_coding||200/2000",
        ];
        let picked = find_most_severe_consequence(&annotations, &csq_fields()).unwrap();
        assert_eq!(
            picked.get("CSQ_Consequence").map(String::as_str),
            Some("stop_gained")
        );
    }

    #[test]
    fn headers_cover_fields_the_first_record_lacks() {
        // Record 1 carries DP and GT only; the header also declares LOF, PGT and an unused SB.
        let header = "##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n##INFO=<ID=LOF,Number=.,Type=String,Description=\"Loss of function\">\n##INFO=<ID=SB,Number=1,Type=Integer,Description=\"Never used\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=PGT,Number=1,Type=String,Description=\"Phasing\">\n";
        let column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1";
        let lines = vec![
            "chr1\t100\t.\tA\tG\t60\tPASS\tDP=30\tGT\t0/1".to_string(),
            "chr1\t200\t.\tC\tT\t60\tPASS\tDP=40;LOF=(X|X|1|1.00)\tGT:PGT\t0/1:0|1".to_string(),
        ];

        let (headers, _records) = reformat_vcf_data_with_header(
            header,
            column_names,
            &lines,
            TranscriptHandling::FirstOnly,
        )
        .unwrap();

        // Present on record 2 only — dropped entirely before the header pass.
        assert!(headers.contains(&"INFO_LOF".to_string()));
        assert!(headers.contains(&"S1_PGT".to_string()));
        // Declared but never used: one empty column is the price, not a missing one.
        assert!(headers.contains(&"INFO_SB".to_string()));
        // Nothing lost, nothing duplicated.
        assert!(headers.contains(&"INFO_DP".to_string()));
        assert!(headers.contains(&"S1_GT".to_string()));
        let mut sorted = headers.clone();
        sorted.sort();
        sorted.dedup();
        assert_eq!(sorted.len(), headers.len(), "duplicate column names");
    }
}
