use regex::Regex;
use std::collections::HashMap;

pub fn extract_csq_format_from_header(header: &str) -> Option<Vec<String>> {
    // Multiple regex patterns to handle different VEP CSQ format variations
    let patterns = vec![
        r#"##INFO=<ID=CSQ,.*?Description=".*?Format:\s*([^"]+)""#,
        r#"##INFO=<ID=CSQ,.*?Format:\s*([^"]+)""#,
        r#"##INFO=<ID=CSQ,.*?Format:([^"]+)""#,
    ];

    for pattern in patterns {
        if let Ok(csq_regex) = Regex::new(pattern) {
            if let Some(captures) = csq_regex.captures(header) {
                if let Some(format_str) = captures.get(1) {
                    let field_names: Vec<String> = format_str
                        .as_str()
                        .split('|')
                        .map(|s| s.trim().to_string())
                        .filter(|s| !s.is_empty())
                        .collect();

                    if !field_names.is_empty() {
                        return Some(field_names);
                    }
                }
            }
        }
    }

    None
}

pub fn extract_ann_format_from_header(header: &str) -> Option<Vec<String>> {
    // SnpEff ANN format patterns
    let patterns = vec![
        r#"##INFO=<ID=ANN,.*?Description=".*?Format:\s*([^"]+)""#,
        r#"##INFO=<ID=ANN,.*?Format:\s*([^"]+)""#,
        r#"##INFO=<ID=ANN,.*?Format:([^"]+)""#,
    ];

    for pattern in patterns {
        if let Ok(ann_regex) = Regex::new(pattern) {
            if let Some(captures) = ann_regex.captures(header) {
                if let Some(format_str) = captures.get(1) {
                    let field_names: Vec<String> = format_str
                        .as_str()
                        .split('|')
                        .map(|s| s.trim().to_string())
                        .filter(|s| !s.is_empty())
                        .collect();

                    if !field_names.is_empty() {
                        return Some(field_names);
                    }
                }
            }
        }
    }

    // Default SnpEff ANN field order if not found in header
    Some(vec![
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
        "ERRORS / WARNINGS / INFO".to_string(),
    ])
}

pub fn _extract_all_info_descriptions(header: &str) -> HashMap<String, String> {
    let mut info_descriptions = HashMap::new();

    // More robust regex pattern
    let patterns = vec![
        r#"##INFO=<ID=([^,]+),.*?Description="([^"]+)""#,
        r#"##INFO=<ID=([^,]+),[^>]*Description=([^,>]+)"#,
    ];

    for pattern in patterns {
        if let Ok(info_regex) = Regex::new(pattern) {
            for captures in info_regex.captures_iter(header) {
                if let (Some(id), Some(desc)) = (captures.get(1), captures.get(2)) {
                    info_descriptions.insert(
                        id.as_str().to_string(),
                        desc.as_str().trim_matches('"').to_string(),
                    );
                }
            }
        }
    }

    info_descriptions
}
