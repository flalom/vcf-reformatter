use regex::Regex;

pub fn extract_csq_regex(parsed_lines: &mut Vec<String>) -> Option<String> {
    if parsed_lines.len() < 8 {
        return None;
    }

    let re = Regex::new(r"CSQ=([^;]+)").unwrap();
    let info_field = &parsed_lines[7];

    let csq_value = re.find(info_field)
        .map(|m| m.as_str()[4..].to_string()); // Remove "CSQ=" prefix

    // Remove CSQ from INFO field
    parsed_lines[7] = re.replace(info_field, "").replace(";;", ";").trim_end_matches(';').to_string();

    csq_value
}