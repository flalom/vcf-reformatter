use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct ParsedSample {
    pub sample_name: String,
    pub format_fields: HashMap<String, String>,
}

#[derive(Debug, Clone)]
pub struct ParsedFormatSample {
    pub format_keys: Vec<String>,
    pub samples: Vec<ParsedSample>,
}

impl ParsedFormatSample {
    pub fn new() -> Self {
        ParsedFormatSample {
            format_keys: Vec::new(),
            samples: Vec::new(),
        }
    }

    pub fn get_headers_for_samples(&self) -> Vec<String> {
        let mut headers = Vec::new();

        for sample in &self.samples {
            for format_key in &self.format_keys {
                headers.push(format!("{}_{}", sample.sample_name, format_key));
            }
        }

        headers
    }

    pub fn get_values_for_samples(&self) -> Vec<String> {
        let mut values = Vec::new();
        let default_value = ".".to_string();

        for sample in &self.samples {
            for format_key in &self.format_keys {
                let value = sample.format_fields.get(format_key).unwrap_or(&default_value);
                values.push(value.clone());
            }
        }

        values
    }

    pub fn get_all_format_keys(&self) -> Vec<String> {
        self.format_keys.clone()
    }

    pub fn get_sample_names(&self) -> Vec<String> {
        self.samples.iter().map(|s| s.sample_name.clone()).collect()
    }
}

// Alias for backward compatibility with tests
pub type SampleData = ParsedSample;

impl ParsedSample {
    pub fn new(sample_name: String) -> Self {
        ParsedSample {
            sample_name,
            format_fields: HashMap::new(),
        }
    }
}

pub fn parse_format_and_samples(
    format_field: Option<&str>,
    sample_fields: &[String],
    sample_names: &[String],
) -> Result<ParsedFormatSample, Box<dyn std::error::Error>> {
    let mut parsed = ParsedFormatSample::new();

    if let Some(format_str) = format_field {
        // Skip if format string is empty or just whitespace
        if format_str.trim().is_empty() {
            return Ok(parsed);
        }
        
        parsed.format_keys = format_str.split(':').map(|s| s.to_string()).collect();

        for (i, sample_field) in sample_fields.iter().enumerate() {
            let sample_name = sample_names.get(i)
                .map(|s| s.to_string())
                .unwrap_or_else(|| format!("SAMPLE_{}", i + 1));

            let sample_values: Vec<&str> = sample_field.split(':').collect();
            let mut format_fields = HashMap::new();

            for (j, format_key) in parsed.format_keys.iter().enumerate() {
                let value = sample_values.get(j)
                    .unwrap_or(&".")
                    .to_string();
                format_fields.insert(format_key.clone(), value);
            }

            parsed.samples.push(ParsedSample {
                sample_name,
                format_fields,
            });
        }
    }

    Ok(parsed)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_format_and_samples() {
        let format = Some("GT:DP:AD:RO:QR:AO:QA:GL");
        let sample_fields = vec!["1/1:19:0,19:0:0:19:579:-32.2782,-5.71957,0".to_string()];
        let sample_names = vec!["B505_B505_V_1".to_string()];

        let parsed = parse_format_and_samples(format, &sample_fields, &sample_names).unwrap();

        assert_eq!(parsed.format_keys, vec!["GT", "DP", "AD", "RO", "QR", "AO", "QA", "GL"]);
        assert_eq!(parsed.samples.len(), 1);
        assert_eq!(parsed.samples[0].sample_name, "B505_B505_V_1");
        assert_eq!(parsed.samples[0].format_fields.get("GT").unwrap(), "1/1");
        assert_eq!(parsed.samples[0].format_fields.get("DP").unwrap(), "19");
        assert_eq!(parsed.samples[0].format_fields.get("AD").unwrap(), "0,19");
    }

    #[test]
    fn test_get_headers_for_samples() {
        let mut parsed = ParsedFormatSample::new();
        parsed.format_keys = vec!["GT".to_string(), "DP".to_string()];

        let mut sample1 = SampleData::new("SAMPLE1".to_string());
        sample1.format_fields.insert("GT".to_string(), "0/1".to_string());
        sample1.format_fields.insert("DP".to_string(), "20".to_string());

        let mut sample2 = SampleData::new("SAMPLE2".to_string());
        sample2.format_fields.insert("GT".to_string(), "1/1".to_string());
        sample2.format_fields.insert("DP".to_string(), "30".to_string());

        parsed.samples = vec![sample1, sample2];

        let headers = parsed.get_headers_for_samples();
        assert_eq!(headers, vec![
            "SAMPLE1_GT", "SAMPLE1_DP",
            "SAMPLE2_GT", "SAMPLE2_DP"
        ]);
    }
}