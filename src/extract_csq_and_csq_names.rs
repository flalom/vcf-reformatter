use regex::Regex;

/// Extract a `TAG=...` annotation value from the INFO field (parsed_lines[7]),
/// stripping it (and any resulting stray `;`) out of the remaining INFO string.
pub(crate) fn extract_tag_regex(parsed_lines: &mut [String], tag: &str) -> Option<String> {
    if parsed_lines.len() < 8 {
        return None;
    }

    let info_field = &parsed_lines[7];
    let re = Regex::new(&format!(r"{tag}=([^;]+)")).ok()?;

    let captures = re.captures(info_field)?;
    let value = captures.get(1)?.as_str().to_string();

    let new_info = re.replace(info_field, "").to_string();
    let cleaned_info = new_info
        .replace(";;", ";")
        .trim_start_matches(';')
        .trim_end_matches(';')
        .to_string();

    parsed_lines[7] = cleaned_info;

    Some(value)
}

pub fn extract_csq_regex(parsed_lines: &mut [String]) -> Option<String> {
    extract_tag_regex(parsed_lines, "CSQ")
}