use crate::extract_csq_and_csq_names::extract_tag_regex;

pub fn extract_ann_regex(parsed_lines: &mut [String]) -> Option<String> {
    extract_tag_regex(parsed_lines, "ANN")
}
