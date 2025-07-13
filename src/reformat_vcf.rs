
use std::collections::HashMap;
use regex::Regex;
use crate::extract_sample_info::{parse_format_and_samples, ParsedFormatSample};
use crate::get_info_from_header::extract_csq_format_from_header;

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
    pub fn from_vcf_line(line: &str, column_names: &[&str], csq_field_names: &Option<Vec<String>>) -> Result<Self, Box<dyn std::error::Error>> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 8 {
            return Err("VCF line has too few fields".into());
        }

        // Parse standard VCF fields
        let chromosome = fields[0].to_string();
        let position = fields[1].parse::<u64>()?;
        let id = if fields[2] == "." { None } else { Some(fields[2].to_string()) };
        let reference = fields[3].to_string();
        let alternate = fields[4].to_string();
        let quality = if fields[5] == "." { None } else { Some(fields[5].parse::<f64>()?) };
        let filter = fields[6].to_string();

        // Parse INFO field into individual key-value pairs
        let info_fields = parse_info_field(fields[7], csq_field_names)?;

        // Parse FORMAT and sample fields if they exist
        let format_sample_data = if fields.len() > 8 {
            let format_field = if fields.len() > 8 { Some(fields[8]) } else { None };
            let sample_fields: Vec<String> = if fields.len() > 9 {
                fields[9..].iter().map(|s| s.to_string()).collect()
            } else {
                Vec::new()
            };

            // Extract sample names from column headers
            let sample_names: Vec<String> = if column_names.len() > 9 {
                column_names[9..].iter().map(|s| s.to_string()).collect()
            } else {
                (0..sample_fields.len()).map(|i| format!("SAMPLE_{}", i + 1)).collect()
            };

            if format_field.is_some() && !sample_fields.is_empty() {
                Some(parse_format_and_samples(format_field, &sample_fields, &sample_names)?)
            } else {
                None
            }
        } else {
            None
        };

        Ok(ReformattedVcfRecord {
            chromosome,
            position,
            id,
            reference,
            alternate,
            quality,
            filter,
            info_fields,
            format_sample_data,
        })
    }
}

fn parse_info_field(info: &str, csq_field_names: &Option<Vec<String>>) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    let mut info_map = HashMap::new();

    if info == "." {
        return Ok(info_map);
    }

    // Split by semicolon to get individual INFO entries
    let entries: Vec<&str> = info.split(';').collect();

    for entry in entries {
        if entry.is_empty() {
            continue;
        }

        // Check if entry contains '=' (key=value format)
        if let Some(eq_pos) = entry.find('=') {
            let key = entry[..eq_pos].to_string();
            let value = entry[eq_pos + 1..].to_string();

            // Special handling for CSQ field
            if key == "CSQ" && csq_field_names.is_some() {
                let csq_fields = parse_csq_field(&value, csq_field_names.as_ref().unwrap())?;
                for (csq_key, csq_value) in csq_fields {
                    info_map.insert(format!("CSQ_{}", csq_key), csq_value);
                }
            } else {
                info_map.insert(key, value);
            }
        } else {
            // Flag fields (no value, just presence indicates true)
            info_map.insert(entry.to_string(), "true".to_string());
        }
    }

    Ok(info_map)
}

fn parse_csq_field(csq_value: &str, csq_field_names: &[String]) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    let mut csq_map = HashMap::new();

    // Split CSQ value by comma to get individual annotations (multiple transcripts)
    let annotations: Vec<&str> = csq_value.split(',').collect();

    // Process the first annotation (can be extended to handle multiple)
    if let Some(first_annotation) = annotations.first() {
        let csq_values: Vec<&str> = first_annotation.split('|').collect();

        // Map each CSQ field name to its corresponding value
        for (i, field_name) in csq_field_names.iter().enumerate() {
            let value = csq_values.get(i).unwrap_or(&"").to_string();
            // Use "." for empty values to be consistent with VCF format
            let final_value = if value.is_empty() { ".".to_string() } else { value };
            csq_map.insert(field_name.clone(), final_value);
        }

        // If there are multiple annotations, we could also add a count field
        if annotations.len() > 1 {
            csq_map.insert("ANNOTATION_COUNT".to_string(), annotations.len().to_string());
        }
    }

    Ok(csq_map)
}

pub fn reformat_vcf_data(column_names: &str, data_lines: &[String]) -> Result<(Vec<String>, Vec<ReformattedVcfRecord>), Box<dyn std::error::Error>> {
    let column_names_vec: Vec<&str> = column_names.trim_start_matches('#').split('\t').collect();
    let mut reformatted_records = Vec::new();
    let mut all_info_keys = std::collections::HashSet::new();
    let mut all_format_keys = std::collections::HashSet::new();
    let mut sample_names = std::collections::HashSet::new();

    // For now, we'll need to get header from somewhere - this will need to be passed in
    // This is a temporary solution - header should be passed as parameter
    let csq_field_names: Option<Vec<String>> = None; // Will be populated from header

    // First pass: collect all unique INFO field keys and FORMAT keys
    for line in data_lines {
        let record = ReformattedVcfRecord::from_vcf_line(line, &column_names_vec, &csq_field_names)?;

        // Collect INFO keys
        for key in record.info_fields.keys() {
            all_info_keys.insert(key.clone());
        }

        // Collect FORMAT keys and sample names
        if let Some(ref format_sample) = record.format_sample_data {
            for key in format_sample.get_all_format_keys() {
                all_format_keys.insert(key);
            }
            for name in format_sample.get_sample_names() {
                sample_names.insert(name);
            }
        }

        reformatted_records.push(record);
    }

    // Create new column headers
    let mut new_headers = vec![
        "CHROM".to_string(),
        "POS".to_string(),
        "ID".to_string(),
        "REF".to_string(),
        "ALT".to_string(),
        "QUAL".to_string(),
        "FILTER".to_string(),
    ];

    // Add sorted INFO field keys as separate columns
    let mut sorted_info_keys: Vec<String> = all_info_keys.into_iter().collect();
    sorted_info_keys.sort();

    for key in &sorted_info_keys {
        new_headers.push(format!("INFO_{}", key));
    }

    // Add FORMAT/sample columns if they exist
    if !all_format_keys.is_empty() {
        let mut sorted_format_keys: Vec<String> = all_format_keys.into_iter().collect();
        sorted_format_keys.sort();

        let mut sorted_sample_names: Vec<String> = sample_names.into_iter().collect();
        sorted_sample_names.sort();

        for sample_name in &sorted_sample_names {
            for format_key in &sorted_format_keys {
                new_headers.push(format!("{}_{}", sample_name, format_key));
            }
        }
    }

    Ok((new_headers, reformatted_records))
}

// Updated function signature to accept header
pub fn reformat_vcf_data_with_header(header: &str, column_names: &str, data_lines: &[String]) -> Result<(Vec<String>, Vec<ReformattedVcfRecord>), Box<dyn std::error::Error>> {
    let column_names_vec: Vec<&str> = column_names.trim_start_matches('#').split('\t').collect();
    let mut reformatted_records = Vec::new();
    let mut all_info_keys = std::collections::HashSet::new();
    let mut all_format_keys = std::collections::HashSet::new();
    let mut sample_names = std::collections::HashSet::new();

    // Extract CSQ field names from header
    let csq_field_names = extract_csq_format_from_header(header);

    // First pass: collect all unique INFO field keys and FORMAT keys
    for line in data_lines {
        let record = ReformattedVcfRecord::from_vcf_line(line, &column_names_vec, &csq_field_names)?;

        // Collect INFO keys
        for key in record.info_fields.keys() {
            all_info_keys.insert(key.clone());
        }

        // Collect FORMAT keys and sample names
        if let Some(ref format_sample) = record.format_sample_data {
            for key in format_sample.get_all_format_keys() {
                all_format_keys.insert(key);
            }
            for name in format_sample.get_sample_names() {
                sample_names.insert(name);
            }
        }

        reformatted_records.push(record);
    }

    // Create new column headers
    let mut new_headers = vec![
        "CHROM".to_string(),
        "POS".to_string(),
        "ID".to_string(),
        "REF".to_string(),
        "ALT".to_string(),
        "QUAL".to_string(),
        "FILTER".to_string(),
    ];

    // Add sorted INFO field keys as separate columns
    let mut sorted_info_keys: Vec<String> = all_info_keys.into_iter().collect();
    sorted_info_keys.sort();

    for key in &sorted_info_keys {
        new_headers.push(format!("INFO_{}", key));
    }

    // Add FORMAT/sample columns if they exist
    if !all_format_keys.is_empty() {
        let mut sorted_format_keys: Vec<String> = all_format_keys.into_iter().collect();
        sorted_format_keys.sort();

        let mut sorted_sample_names: Vec<String> = sample_names.into_iter().collect();
        sorted_sample_names.sort();

        for sample_name in &sorted_sample_names {
            for format_key in &sorted_format_keys {
                new_headers.push(format!("{}_{}", sample_name, format_key));
            }
        }
    }

    Ok((new_headers, reformatted_records))
}

pub fn write_reformatted_vcf(filename: &str, headers: &[String], records: &[ReformattedVcfRecord]) -> Result<(), Box<dyn std::error::Error>> {
    use std::io::Write;
    use std::fs::File;

    let mut file = File::create(filename)?;

    // Write headers
    writeln!(file, "{}", headers.join("\t"))?;

    // Get all INFO keys from headers (those that start with "INFO_")
    let info_keys: Vec<String> = headers.iter()
        .filter(|h| h.starts_with("INFO_"))
        .map(|h| h[5..].to_string()) // Remove "INFO_" prefix
        .collect();

    // Get all sample format keys from headers (those that contain '_' but don't start with "INFO_")
    let sample_format_headers: Vec<String> = headers.iter()
        .filter(|h| h.contains('_') && !h.starts_with("INFO_"))
        .map(|h| h.to_string())
        .collect();

    // Write data
    for record in records {
        let mut row = Vec::new();

        // Standard VCF fields
        row.push(record.chromosome.clone());
        row.push(record.position.to_string());
        row.push(record.id.as_ref().unwrap_or(&".".to_string()).clone());
        row.push(record.reference.clone());
        row.push(record.alternate.clone());
        row.push(record.quality.map_or(".".to_string(), |q| q.to_string()));
        row.push(record.filter.clone());

        // INFO fields in the same order as headers
        for key in &info_keys {
            row.push(record.info_fields.get(key).unwrap_or(&".".to_string()).clone());
        }

        // Sample format fields
        for header in &sample_format_headers {
            let mut value = ".".to_string();

            if let Some(ref format_sample) = record.format_sample_data {
                // Parse header to get sample name and format key
                if let Some(underscore_pos) = header.rfind('_') {
                    let sample_name = &header[..underscore_pos];
                    let format_key = &header[underscore_pos + 1..];

                    // Find the sample and get the value
                    for sample in &format_sample.samples {
                        if sample.sample_name == sample_name {
                            value = sample.format_fields.get(format_key).unwrap_or(&".".to_string()).clone();
                            break;
                        }
                    }
                }
            }

            row.push(value);
        }

        writeln!(file, "{}", row.join("\t"))?;
    }

    Ok(())
}