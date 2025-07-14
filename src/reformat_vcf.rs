use crate::extract_sample_info::{parse_format_and_samples, ParsedFormatSample};
use crate::get_info_from_header::extract_csq_format_from_header;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::{create_dir_all, File};
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone, Copy)]
pub enum TranscriptHandling {
    MostSevere,
    FirstOnly,
    SplitRows,
}

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
    pub fn from_vcf_line(
        line: &str,
        column_names: &[&str],
        csq_field_names: &Option<Vec<String>>,
        transcript_handling: TranscriptHandling,
    ) -> std::result::Result<Vec<Self>, Box<dyn std::error::Error>> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 8 {
            return Err("VCF line has too few fields".into());
        }

        // Parse basic VCF fields
        let chromosome = fields[0].to_string();
        let position = fields[1].parse::<u64>()?;
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
            Some(fields[5].parse::<f64>()?)
        };
        let filter = fields[6].to_string();
        let info = fields[7].to_string();

        // Parse INFO field and handle CSQ
        let info_variants = parse_info_field(&info, csq_field_names, transcript_handling)?;

        // Parse FORMAT and sample data
        let format_field = if fields.len() > 8 && !fields[8].is_empty() && fields[8] != "." {
            Some(fields[8])
        } else {
            None
        };

        let sample_fields: Vec<String> = if fields.len() > 9 {
            fields[9..].iter().map(|s| s.to_string()).collect()
        } else {
            Vec::new()
        };

        let sample_names: Vec<String> = if column_names.len() > 9 {
            column_names[9..].iter().map(|s| s.to_string()).collect()
        } else {
            Vec::new()
        };

        // Only parse sample data if we have both format field and sample data
        let format_sample_data = if let Some(format_str) = format_field {
            if !sample_fields.is_empty() && !sample_names.is_empty() {
                Some(parse_format_and_samples(
                    Some(format_str),
                    &sample_fields,
                    &sample_names,
                )?)
            } else {
                None
            }
        } else {
            None
        };

        // Create records for each info variant (handling multiple transcripts)
        let mut records = Vec::new();
        for info_fields in info_variants {
            records.push(ReformattedVcfRecord {
                chromosome: chromosome.clone(),
                position,
                id: id.clone(),
                reference: reference.clone(),
                alternate: alternate.clone(),
                quality,
                filter: filter.clone(),
                info_fields,
                format_sample_data: format_sample_data.clone(),
            });
        }

        Ok(records)
    }
}

fn parse_info_field(
    info: &str,
    csq_field_names: &Option<Vec<String>>,
    transcript_handling: TranscriptHandling,
) -> std::result::Result<Vec<HashMap<String, String>>, Box<dyn std::error::Error>> {
    let mut info_map = HashMap::new();
    let mut csq_value = None;

    // Split the INFO field by semicolon and collect all fields
    for pair in info.split(';') {
        if pair.is_empty() {
            continue;
        }

        if let Some(eq_pos) = pair.find('=') {
            let key = pair[..eq_pos].to_string();
            let value = pair[eq_pos + 1..].to_string();

            // Store CSQ separately for special handling
            if key == "CSQ" {
                csq_value = Some(value);
            } else {
                info_map.insert(format!("INFO_{key}"), value);
            }
        } else {
            // Flag without value
            info_map.insert(format!("INFO_{pair}"), "true".to_string());
        }
    }

    // Handle CSQ field if present
    if let Some(csq_val) = csq_value {
        if let Some(csq_names) = csq_field_names {
            let csq_variants =
                parse_csq_field_with_handling(&csq_val, csq_names, transcript_handling)?;

            // Combine INFO fields with CSQ fields for each variant
            let mut results = Vec::new();
            for csq_map in csq_variants {
                let mut combined_map = info_map.clone();
                combined_map.extend(csq_map);
                results.push(combined_map);
            }
            return Ok(results);
        }
    }

    // If no CSQ processing was done, return a single variant
    Ok(vec![info_map])
}

fn parse_csq_field_with_handling(
    csq_value: &str,
    csq_field_names: &[String],
    transcript_handling: TranscriptHandling,
) -> std::result::Result<Vec<HashMap<String, String>>, Box<dyn std::error::Error>> {
    let annotations: Vec<&str> = csq_value.split(',').collect();

    match transcript_handling {
        TranscriptHandling::MostSevere => {
            let most_severe = find_most_severe_consequence(&annotations, csq_field_names)?;
            Ok(vec![most_severe])
        }
        TranscriptHandling::FirstOnly => {
            if let Some(first_annotation) = annotations.first() {
                let parsed = parse_single_csq_annotation(first_annotation, csq_field_names)?;
                Ok(vec![parsed])
            } else {
                Ok(vec![HashMap::new()])
            }
        }
        TranscriptHandling::SplitRows => {
            let mut results = Vec::new();
            for annotation in annotations {
                let parsed = parse_single_csq_annotation(annotation, csq_field_names)?;
                results.push(parsed);
            }
            Ok(results)
        }
    }
}

fn parse_single_csq_annotation(
    annotation: &str,
    csq_field_names: &[String],
) -> std::result::Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    let mut csq_map = HashMap::new();
    let values: Vec<&str> = annotation.split('|').collect();

    for (i, field_name) in csq_field_names.iter().enumerate() {
        let value = values.get(i).unwrap_or(&"").to_string();
        csq_map.insert(format!("CSQ_{field_name}"), value);
    }

    Ok(csq_map)
}

fn find_most_severe_consequence(
    annotations: &[&str],
    csq_field_names: &[String],
) -> std::result::Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    if annotations.is_empty() {
        return Ok(HashMap::new());
    }

    // VEP consequence severity ranking (most severe first)
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

    // Find consequence column index
    let consequence_index = csq_field_names
        .iter()
        .position(|name| name == "Consequence")
        .unwrap_or(0);

    for annotation in annotations {
        let values: Vec<&str> = annotation.split('|').collect();
        if let Some(consequence) = values.get(consequence_index) {
            // Handle multiple consequences separated by &
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

pub fn reformat_vcf_data_with_header(
    header: &str,
    column_names: &str,
    data_lines: &[String],
    transcript_handling: TranscriptHandling,
) -> std::result::Result<(Vec<String>, Vec<ReformattedVcfRecord>), Box<dyn std::error::Error>> {
    // Extract CSQ field names from the header
    let csq_field_names = extract_csq_format_from_header(header);

    // Parse column names
    let columns: Vec<&str> = column_names.split('\t').collect();

    // Extract sample names (columns after FORMAT)
    let sample_names: Vec<String> = if columns.len() > 9 {
        columns[9..].iter().map(|s| s.to_string()).collect()
    } else {
        Vec::new()
    };

    // Process data lines
    let mut records = Vec::new();
    for line in data_lines {
        let line_records = ReformattedVcfRecord::from_vcf_line(
            line,
            &columns,
            &csq_field_names,
            transcript_handling,
        )?;
        records.extend(line_records);
    }

    // Generate headers from the first record
    let headers = if let Some(first_record) = records.first() {
        generate_headers_from_record(first_record, &sample_names)
    } else {
        vec![
            "CHROM".to_string(),
            "POS".to_string(),
            "ID".to_string(),
            "REF".to_string(),
            "ALT".to_string(),
            "QUAL".to_string(),
            "FILTER".to_string(),
        ]
    };

    Ok((headers, records))
}

pub fn reformat_vcf_data_with_header_parallel(
    header: &str,
    column_names: &str,
    data_lines: &[String],
    transcript_handling: TranscriptHandling,
) -> std::result::Result<(Vec<String>, Vec<ReformattedVcfRecord>), Box<dyn std::error::Error>> {
    // Extract CSQ field names from the header
    let csq_field_names = extract_csq_format_from_header(header);

    // Parse column names
    let columns: Vec<&str> = column_names.split('\t').collect();

    // Extract sample names (columns after FORMAT)
    let sample_names: Vec<String> = if columns.len() > 9 {
        columns[9..].iter().map(|s| s.to_string()).collect()
    } else {
        Vec::new()
    };

    // Process data lines in parallel
    let results: std::result::Result<Vec<Vec<ReformattedVcfRecord>>, String> = data_lines
        .par_iter()
        .map(|line| {
            ReformattedVcfRecord::from_vcf_line(
                line,
                &columns,
                &csq_field_names,
                transcript_handling,
            )
            .map_err(|e| format!("Error processing line: {e}"))
        })
        .collect();

    let records: Vec<ReformattedVcfRecord> = results?.into_iter().flatten().collect();

    // Generate headers from the first record
    let headers = if let Some(first_record) = records.first() {
        generate_headers_from_record(first_record, &sample_names)
    } else {
        vec![
            "CHROM".to_string(),
            "POS".to_string(),
            "ID".to_string(),
            "REF".to_string(),
            "ALT".to_string(),
            "QUAL".to_string(),
            "FILTER".to_string(),
        ]
    };

    Ok((headers, records))
}

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

    // Add INFO field headers (sorted for consistency)
    let mut info_keys: Vec<String> = record
        .info_fields
        .keys()
        .filter(|k| k.starts_with("INFO_"))
        .cloned()
        .collect();
    info_keys.sort();
    headers.extend(info_keys);

    // Add CSQ field headers (sorted for consistency)
    let mut csq_keys: Vec<String> = record
        .info_fields
        .keys()
        .filter(|k| k.starts_with("CSQ_"))
        .cloned()
        .collect();
    csq_keys.sort();
    headers.extend(csq_keys);

    // Add sample headers
    if let Some(ref sample_data) = record.format_sample_data {
        headers.extend(sample_data.get_headers_for_samples());
    }

    headers
}

pub fn write_reformatted_vcf(
    filename: &str,
    headers: &[String],
    records: &[ReformattedVcfRecord],
    compress: bool,
) -> std::io::Result<()> {
    // Create the directory structure if it doesn't exist
    if let Some(parent) = Path::new(filename).parent() {
        create_dir_all(parent)?;
    }

    let file = File::create(filename)?;

    if compress {
        // Use gzip compression
        let encoder = GzEncoder::new(file, Compression::default());
        let mut writer = BufWriter::new(encoder);
        write_tsv_content(&mut writer, headers, records)?;
        writer.flush()?;
    } else {
        // Write uncompressed
        let mut writer = BufWriter::new(file);
        write_tsv_content(&mut writer, headers, records)?;
        writer.flush()?;
    }

    Ok(())
}

fn write_tsv_content<W: Write>(
    writer: &mut W,
    headers: &[String],
    records: &[ReformattedVcfRecord],
) -> std::io::Result<()> {
    // Write headers
    writeln!(writer, "{}", headers.join("\t"))?;

    // Write records
    for record in records {
        let mut row = Vec::new();

        // Process each header in order to maintain column alignment
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
                    // Handle INFO, CSQ, and sample fields
                    if header.starts_with("INFO_") || header.starts_with("CSQ_") {
                        record
                            .info_fields
                            .get(header)
                            .unwrap_or(&".".to_string())
                            .clone()
                    } else {
                        // This is likely a sample field (SAMPLE_FORMAT)
                        if let Some(ref sample_data) = record.format_sample_data {
                            // Find the sample and format key by checking all samples
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

pub fn reformat_vcf_data(
    column_names: &str,
    data_lines: &[String],
) -> std::result::Result<(Vec<String>, Vec<ReformattedVcfRecord>), Box<dyn std::error::Error>> {
    // Default to FirstOnly for backward compatibility
    reformat_vcf_data_with_header("", column_names, data_lines, TranscriptHandling::FirstOnly)
}
