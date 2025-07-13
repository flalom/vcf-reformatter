use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct SampleData {
    pub sample_name: String,
    pub format_fields: HashMap<String, String>,
}

impl SampleData {
    pub fn new(sample_name: String) -> Self {
        SampleData {
            sample_name,
            format_fields: HashMap::new(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct ParsedFormatSample {
    pub format_keys: Vec<String>,
    pub samples: Vec<SampleData>,
}

impl ParsedFormatSample {
    pub fn new() -> Self {
        ParsedFormatSample {
            format_keys: Vec::new(),
            samples: Vec::new(),
        }
    }

    pub fn from_vcf_fields(
        format_field: Option<&str>,
        sample_fields: &[String],
        sample_names: &[String],
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let mut parsed = ParsedFormatSample::new();

        // If no FORMAT field, return empty
        let format = match format_field {
            Some(f) => f,
            None => return Ok(parsed),
        };

        // Parse FORMAT field to get keys
        parsed.format_keys = format.split(':').map(|s| s.to_string()).collect();

        // Parse each sample
        for (i, sample_field) in sample_fields.iter().enumerate() {
            let sample_name = sample_names.get(i)
                .unwrap_or(&format!("SAMPLE_{}", i + 1))
                .clone();

            let mut sample_data = SampleData::new(sample_name);

            // Split sample values by ':'
            let sample_values: Vec<&str> = sample_field.split(':').collect();

            // Map FORMAT keys to sample values
            for (j, key) in parsed.format_keys.iter().enumerate() {
                let value = sample_values.get(j).unwrap_or(&".").to_string();
                sample_data.format_fields.insert(key.clone(), value);
            }

            parsed.samples.push(sample_data);
        }

        Ok(parsed)
    }

    pub fn get_all_format_keys(&self) -> Vec<String> {
        self.format_keys.clone()
    }

    pub fn get_sample_names(&self) -> Vec<String> {
        self.samples.iter().map(|s| s.sample_name.clone()).collect()
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
}

pub fn parse_format_and_samples(
    format_field: Option<&str>,
    sample_fields: &[String],
    sample_names: &[String],
) -> Result<ParsedFormatSample, Box<dyn std::error::Error>> {
    ParsedFormatSample::from_vcf_fields(format_field, sample_fields, sample_names)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_format_sample() {
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
}