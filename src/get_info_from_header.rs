use regex::Regex;

/// Look up the `Format:` field-name list for a given INFO `id` (e.g. "CSQ" or
/// "ANN") in a VCF header, trying a few Description-formatting variations.
fn extract_format_field_names(header: &str, id: &str) -> Option<Vec<String>> {
    let patterns = [
        format!(r#"##INFO=<ID={id},.*?Description=".*?Format:\s*([^"]+)""#),
        format!(r#"##INFO=<ID={id},.*?Format:\s*([^"]+)""#),
        format!(r#"##INFO=<ID={id},.*?Format:([^"]+)""#),
    ];

    for pattern in &patterns {
        if let Ok(regex) = Regex::new(pattern) {
            if let Some(format_str) = regex.captures(header).and_then(|c| c.get(1)) {
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

    None
}

pub fn extract_csq_format_from_header(header: &str) -> Option<Vec<String>> {
    extract_format_field_names(header, "CSQ")
}

pub fn extract_ann_format_from_header(header: &str) -> Option<Vec<String>> {
    extract_format_field_names(header, "ANN").or_else(|| {
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
    })
}
