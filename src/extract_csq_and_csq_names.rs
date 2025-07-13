use regex::Regex;

pub fn extract_csq_regex(parsed_lines: &mut Vec<String>) -> Option<String> {
    if parsed_lines.len() < 8 {
        return None;
    }

    let info_field = &parsed_lines[7];
    let re = Regex::new(r"CSQ=([^;]+)").ok()?;

    if let Some(captures) = re.captures(info_field) {
        if let Some(csq_match) = captures.get(1) {
            let csq_value = csq_match.as_str().to_string();

            // Remove CSQ from the INFO field
            let new_info = re.replace(info_field, "").to_string();

            // Clean up any double semicolons or leading/trailing semicolons
            let cleaned_info = new_info
                .replace(";;", ";")
                .trim_start_matches(';')
                .trim_end_matches(';')
                .to_string();

            parsed_lines[7] = cleaned_info;

            return Some(csq_value);
        }
    }

    None
}

