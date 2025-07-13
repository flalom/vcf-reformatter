#[derive(Debug, Clone)]
pub struct VcfVariant {
    pub chromosome: String,
    pub position: u64,
    pub id: Option<String>,
    pub reference: String,
    pub alternate: String,
    pub quality: Option<f64>,
    pub filter: String,
    pub info: String,
    pub format: Option<String>,
    pub samples: Vec<String>,
}
impl VcfVariant {
    // Constructor from a VCF line and column headers
    pub fn from_line(line: &str, column_names: &[&str]) -> Result<Self, Box<dyn std::error::Error>> {
        //Flexible: Works with Vec<&str>, arrays [&str; N], or any slice
        // No copying: Just borrows references, very efficient
        // No ownership: Function doesn't take ownership of your data
        let fields: Vec<&str> = line.split('\t').collect();

        // Ensure we have enough fields
        if fields.len() < 8 {
            return Err("VCF line has too few fields".into());
        }

        // Parse required fields (first 8 are standard VCF)
        let chromosome = fields[0].to_string();
        let position = fields[1].parse::<u64>()?;
        let id = if fields[2] == "." { None } else { Some(fields[2].to_string()) };
        let reference = fields[3].to_string();
        let alternate = fields[4].to_string();
        let quality = if fields[5] == "." { None } else { Some(fields[5].parse::<f64>()?) };
        let filter = fields[6].to_string();
        let info = fields[7].to_string();

        // Handle optional FORMAT and sample columns
        let format = if fields.len() > 8 { Some(fields[8].to_string()) } else { None };
        let samples = if fields.len() > 9 {
            fields[9..].iter().map(|s| s.to_string()).collect()
        } else {
            Vec::new()
        };

        Ok(VcfVariant {
            chromosome,
            position,
            id,
            reference,
            alternate,
            quality,
            filter,
            info,
            format,
            samples,
        })
    }
}

#[derive(Debug, Clone)]
pub struct MafRecord {
    pub hugo_symbol: String,
    pub entrez_gene_id: Option<u32>,
    pub center: String,
    pub ncbi_build: String,
    pub chromosome: String,
    pub start_position: u64,
    pub end_position: u64,
    pub strand: String,
    pub variant_classification: String,
    pub variant_type: String,
    pub reference_allele: String,
    pub tumor_seq_allele1: String,
    pub tumor_seq_allele2: String,
    pub dbsnp_rs: Option<String>,
    pub dbsnp_val_status: Option<String>,
    pub tumor_sample_barcode: String,
    pub matched_norm_sample_barcode: Option<String>,
    // ... additional MAF fields
}

