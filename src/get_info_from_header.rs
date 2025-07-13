
use regex::Regex;
use std::collections::HashMap;

pub fn extract_csq_format_from_header(header: &str) -> Option<Vec<String>> {
    // Look for CSQ INFO field definition in header
    let csq_regex = Regex::new(r#"##INFO=<ID=CSQ,.*?Description=".*?Format:\s*([^"]+)""#).unwrap();

    if let Some(captures) = csq_regex.captures(header) {
        if let Some(format_str) = captures.get(1) {
            // Split the format string by | to get field names
            let field_names: Vec<String> = format_str.as_str()
                .split('|')
                .map(|s| s.trim().to_string())
                .collect();
            return Some(field_names);
        }
    }

    None
}

pub fn extract_all_info_descriptions(header: &str) -> HashMap<String, String> {
    let mut info_descriptions = HashMap::new();
    let info_regex = Regex::new(r#"##INFO=<ID=([^,]+),.*?Description="([^"]+)""#).unwrap();

    for captures in info_regex.captures_iter(header) {
        if let (Some(id), Some(desc)) = (captures.get(1), captures.get(2)) {
            info_descriptions.insert(id.as_str().to_string(), desc.as_str().to_string());
        }
    }

    info_descriptions
}