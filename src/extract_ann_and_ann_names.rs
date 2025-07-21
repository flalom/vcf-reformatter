use regex::Regex;

pub fn extract_ann_regex(parsed_lines: &mut [String]) -> Option<String> {
    if parsed_lines.len() < 8 {
        return None;
    }

    let info_field = &parsed_lines[7];
    let re = Regex::new(r"ANN=([^;]+)").ok()?;

    if let Some(captures) = re.captures(info_field) {
        if let Some(ann_match) = captures.get(1) {
            let ann_value = ann_match.as_str().to_string();

            // Remove ANN from the INFO field
            let new_info = re.replace(info_field, "").to_string();

            // Clean up any double semicolons or leading/trailing semicolons
            let cleaned_info = new_info
                .replace(";;", ";")
                .trim_start_matches(';')
                .trim_end_matches(';')
                .to_string();

            parsed_lines[7] = cleaned_info;

            return Some(ann_value);
        }
    }

    None
}
