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
use crate::essentials_fields::MafRecord;
use crate::extract_ann_and_ann_names::extract_ann_regex;
use crate::extract_csq_and_csq_names::extract_csq_regex;
use crate::extract_sample_info::ParsedFormatSample;
use crate::get_info_from_header::{extract_ann_format_from_header, extract_csq_format_from_header};
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
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

        let info_variants =
            parse_info_field(info, csq_field_names, ann_field_names, transcript_handling)?;

        let records: Vec<Self> = info_variants
            .into_iter()
            .map(|info_fields| Self {
                chromosome: chromosome.clone(),
                position,
                id: id.clone(),
                reference: reference.clone(),
                alternate: alternate.clone(),
                quality,
                filter: filter.clone(),
                info_fields,
                format_sample_data: format_sample_data.clone(),
            })
            .collect();

        Ok(records)
    }
}
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
) -> std::result::Result<Vec<HashMap<String, String>>, Box<dyn std::error::Error>> {
    if info.is_empty() {
        return Ok(vec![HashMap::new()]);
    }

    let annotation_result =
        parse_annotation_fields(info, csq_field_names, ann_field_names, transcript_handling)?;

    let remaining_info_map = parse_remaining_info_fields(&annotation_result.remaining_info)?;
    let combined_records =
        combine_annotation_with_info(annotation_result.records, remaining_info_map.clone());

    if combined_records.is_empty() {
        Ok(vec![remaining_info_map])
    } else {
        Ok(combined_records)
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
                match parse_csq_field_with_handling(&csq_value, field_names, transcript_handling) {
                    Ok(records) if !records.is_empty() => {
                        return Ok(AnnotationParseResult {
                            field_type: AnnotationFieldType::Csq,
                            records,
                            remaining_info: parsed_lines[7].clone(),
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
                match parse_ann_field_with_handling(&ann_value, field_names, transcript_handling) {
                    Ok(records) if !records.is_empty() => {
                        return Ok(AnnotationParseResult {
                            field_type: AnnotationFieldType::Ann,
                            records,
                            remaining_info: parsed_lines[7].clone(),
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

// FIXED: Remove the duplicate function and keep only this safer version
fn parse_csq_field_with_handling(
    csq_value: &str,
    csq_field_names: &[String],
    transcript_handling: TranscriptHandling,
) -> Result<Vec<HashMap<String, String>>, Box<dyn std::error::Error>> {
    if csq_value.trim().is_empty() {
        return Ok(Vec::new());
    }

    let annotations: Vec<&str> = csq_value
        .split(',')
        .filter(|s| !s.trim().is_empty())
        .collect();

    if annotations.is_empty() {
        return Ok(Vec::new());
    }

    match transcript_handling {
        TranscriptHandling::FirstOnly => {
            let first_annotation = annotations
                .first()
                .ok_or("No annotations found after filtering")?;
            let parsed = parse_single_csq_annotation(first_annotation, csq_field_names)?;
            Ok(vec![parsed])
        }
        TranscriptHandling::MostSevere => {
            let most_severe = find_most_severe_consequence(&annotations, csq_field_names)?;
            Ok(vec![most_severe])
        }
        TranscriptHandling::SplitRows => {
            let mut all_annotations = Vec::new();
            for annotation in annotations {
                match parse_single_csq_annotation(annotation, csq_field_names) {
                    Ok(parsed) => all_annotations.push(parsed),
                    Err(e) => {
                        eprintln!("Warning: Failed to parse CSQ annotation '{annotation}': {e}");
                    }
                }
            }
            Ok(all_annotations)
        }
    }
}

fn parse_ann_field_with_handling(
    ann_value: &str,
    ann_field_names: &[String],
    transcript_handling: TranscriptHandling,
) -> Result<Vec<HashMap<String, String>>, Box<dyn std::error::Error>> {
    if ann_value.trim().is_empty() {
        return Ok(Vec::new());
    }

    let annotations: Vec<&str> = ann_value
        .split(',')
        .filter(|s| !s.trim().is_empty())
        .collect();

    if annotations.is_empty() {
        return Ok(Vec::new());
    }

    match transcript_handling {
        TranscriptHandling::FirstOnly => {
            let first_annotation = annotations
                .first()
                .ok_or("No annotations found after filtering")?;
            let parsed = parse_single_ann_annotation(first_annotation, ann_field_names)?;
            Ok(vec![parsed])
        }
        TranscriptHandling::MostSevere => {
            let most_severe = find_most_severe_ann_consequence(&annotations, ann_field_names)?;
            Ok(vec![most_severe])
        }
        TranscriptHandling::SplitRows => {
            let mut all_annotations = Vec::new();
            for annotation in annotations {
                match parse_single_ann_annotation(annotation, ann_field_names) {
                    Ok(parsed) => all_annotations.push(parsed),
                    Err(e) => {
                        eprintln!("Warning: Failed to parse ANN annotation '{annotation}': {e}");
                    }
                }
            }
            Ok(all_annotations)
        }
    }
}

fn parse_single_csq_annotation(
    annotation: &str,
    csq_field_names: &[String],
) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    if annotation.trim().is_empty() {
        return Err("Empty CSQ annotation".into());
    }

    let values: Vec<&str> = annotation.split('|').collect();
    let mut annotation_map = HashMap::new();

    for (i, field_name) in csq_field_names.iter().enumerate() {
        let value = values.get(i).unwrap_or(&"").trim();
        let formatted_value = if value.is_empty() { "." } else { value };
        annotation_map.insert(
            format!("CSQ_{}", sanitize_field_name(field_name)),
            formatted_value.to_string(),
        );
    }

    if values.len() > csq_field_names.len() {
        for (i, value) in values.iter().enumerate().skip(csq_field_names.len()) {
            annotation_map.insert(
                format!("CSQ_EXTRA_{}", i - csq_field_names.len() + 1),
                value.trim().to_string(),
            );
        }
    }

    Ok(annotation_map)
}

fn parse_single_ann_annotation(
    annotation: &str,
    ann_field_names: &[String],
) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    if annotation.trim().is_empty() {
        return Err("Empty ANN annotation".into());
    }

    let values: Vec<&str> = annotation.split('|').collect();
    let mut annotation_map = HashMap::new();

    for (i, field_name) in ann_field_names.iter().enumerate() {
        let value = values.get(i).unwrap_or(&"").trim();
        let formatted_value = if value.is_empty() { "." } else { value };
        annotation_map.insert(
            format!("ANN_{}", sanitize_field_name(field_name)),
            formatted_value.to_string(),
        );
    }

    if values.len() > ann_field_names.len() {
        for (i, value) in values.iter().enumerate().skip(ann_field_names.len()) {
            annotation_map.insert(
                format!("ANN_EXTRA_{}", i - ann_field_names.len() + 1),
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
    parse_single_ann_annotation(selected_annotation, ann_field_names)
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

fn find_most_severe_consequence(
    annotations: &[&str],
    csq_field_names: &[String],
) -> std::result::Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    if annotations.is_empty() {
        return Err("No annotations provided".into());
    }

    let severity_order = vec![
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
        "transcript_amplification",
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant",
        "splice_region_variant",
        "incomplete_terminal_codon_variant",
        "start_retained_variant",
        "stop_retained_variant",
        "synonymous_variant",
        "coding_sequence_variant",
        "mature_miRNA_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "non_coding_transcript_exon_variant",
        "intron_variant",
        "NMD_transcript_variant",
        "non_coding_transcript_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "TFBS_ablation",
        "TFBS_amplification",
        "TF_binding_site_variant",
        "regulatory_region_ablation",
        "regulatory_region_amplification",
        "feature_elongation",
        "regulatory_region_variant",
        "feature_truncation",
        "intergenic_variant",
    ];

    let mut most_severe_annotation = annotations[0];
    let mut best_severity = usize::MAX;

    let consequence_index = csq_field_names
        .iter()
        .position(|name| name == "Consequence")
        .unwrap_or(1);

    for annotation in annotations {
        let values: Vec<&str> = annotation.split('|').collect();
        if let Some(consequence) = values.get(consequence_index) {
            let consequences: Vec<&str> = consequence.split('&').collect();

            for cons in consequences {
                if let Some(severity) = severity_order.iter().position(|&x| x == cons) {
                    if severity < best_severity {
                        best_severity = severity;
                        most_severe_annotation = annotation;
                    }
                }
            }
        }
    }

    parse_single_csq_annotation(most_severe_annotation, csq_field_names)
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

    let headers = generate_headers_from_records(&all_records, &column_names_vec);

    Ok((headers, all_records))
}

/// Generate column headers from the first reformatted record
fn generate_headers_from_records(
    records: &[ReformattedVcfRecord],
    column_names_vec: &[&str],
) -> Vec<String> {
    if let Some(first_record) = records.first() {
        let sample_names: Vec<String> = if column_names_vec.len() > 9 {
            column_names_vec[9..]
                .iter()
                .map(|s| s.to_string())
                .collect()
        } else {
            vec![]
        };
        generate_headers_from_record(first_record, &sample_names)
    } else {
        vec![]
    }
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

    let headers = generate_headers_from_records(&flattened_records, &column_names_vec);

    Ok((headers, flattened_records))
}

// FIXED: Update to handle both CSQ and ANN headers
fn generate_headers_from_record(
    record: &ReformattedVcfRecord,
    _sample_names: &[String],
) -> Vec<String> {
    let mut headers = vec![
        "CHROM".to_string(),
        "POS".to_string(),
        "ID".to_string(),
        "REF".to_string(),
        "ALT".to_string(),
        "QUAL".to_string(),
        "FILTER".to_string(),
    ];

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
    writeln!(writer, "{}", headers.join("\t"))?;

    for record in records {
        let mut row = Vec::new();

        for header in headers {
            let value = match header.as_str() {
                "CHROM" => record.chromosome.clone(),
                "POS" => record.position.to_string(),
                "ID" => record.id.as_ref().unwrap_or(&".".to_string()).clone(),
                "REF" => record.reference.clone(),
                "ALT" => record.alternate.clone(),
                "QUAL" => record
                    .quality
                    .map(|q| q.to_string())
                    .unwrap_or(".".to_string()),
                "FILTER" => record.filter.clone(),
                _ => {
                    if header.starts_with("INFO_")
                        || header.starts_with("CSQ_")
                        || header.starts_with("ANN_")
                    {
                        record
                            .info_fields
                            .get(header)
                            .unwrap_or(&".".to_string())
                            .clone()
                    } else {
                        if let Some(ref sample_data) = record.format_sample_data {
                            let mut found_value = None;

                            for sample in &sample_data.samples {
                                for format_key in &sample_data.format_keys {
                                    let expected_header =
                                        format!("{}_{}", sample.sample_name, format_key);
                                    if expected_header == *header {
                                        found_value = sample.format_fields.get(format_key).cloned();
                                        break;
                                    }
                                }
                                if found_value.is_some() {
                                    break;
                                }
                            }

                            found_value.unwrap_or(".".to_string())
                        } else {
                            ".".to_string()
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
/// Simplified VCF reformatting function with default settings
///
/// This convenience function reformats VCF data using default settings:
/// - No header parsing (empty header)
/// - First transcript only
///
/// Use `reformat_vcf_data_with_header` for more control over the processing.
///
/// # Arguments
///
/// * `column_names` - Column header line from VCF
/// * `data_lines` - VCF data lines to process
///
/// # Returns
///
/// A tuple of (headers, reformatted_records)
pub fn _reformat_vcf_data(
    column_names: &str,
    data_lines: &[String],
) -> std::result::Result<(Vec<String>, Vec<ReformattedVcfRecord>), Box<dyn std::error::Error>> {
    reformat_vcf_data_with_header("", column_names, data_lines, TranscriptHandling::FirstOnly)
}

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

// Helper function to write MAF content
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
            output_headers = generate_headers_from_records(&chunk_records, &column_names_vec);

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
        if total_processed % 100_000 == 0 {
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
fn extract_values_from_record(record: &ReformattedVcfRecord, headers: &[String]) -> Vec<String> {
    headers
        .iter()
        .map(|header| {
            match header.as_str() {
                "CHROM" => record.chromosome.clone(),
                "POS" => record.position.to_string(),
                "ID" => record.id.as_deref().unwrap_or(".").to_string(),
                "REF" => record.reference.clone(),
                "ALT" => record.alternate.clone(),
                "QUAL" => record.quality.map_or(".".to_string(), |q| q.to_string()),
                "FILTER" => record.filter.clone(),
                _ => {
                    // Handle INFO fields, CSQ fields, ANN fields, and sample data
                    if let Some(value) = record.info_fields.get(header) {
                        value.clone()
                    } else if let Some(sample_data) = &record.format_sample_data {
                        // Get sample value for this header
                        extract_sample_value_for_header(sample_data, header)
                    } else {
                        ".".to_string()
                    }
                }
            }
        })
        .collect()
}

/// Helper function to extract sample values by header name
fn extract_sample_value_for_header(sample_data: &ParsedFormatSample, header: &str) -> String {
    // Header format: "SAMPLE_NAME_FORMAT_KEY"
    for sample in &sample_data.samples {
        for format_key in &sample_data.format_keys {
            let expected_header = format!("{}_{}", sample.sample_name, format_key);
            if expected_header == header {
                return sample
                    .format_fields
                    .get(format_key)
                    .cloned()
                    .unwrap_or_else(|| ".".to_string());
            }
        }
    }
    ".".to_string()
}
