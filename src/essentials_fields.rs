use crate::extract_sample_info::ParsedFormatSample;
use crate::reformat_vcf::ReformattedVcfRecord;
use std::collections::HashMap;

#[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn from_line(
        line: &str,
        _column_names: &[&str],
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 8 {
            return Err("VCF line has too few fields".into());
        }

        let chromosome = fields[0].to_string();
        let position = fields[1].parse::<u64>()?;
        let id = if fields[2] == "." {
            None
        } else {
            Some(fields[2].to_string())
        };
        let reference = fields[3].to_string();
        let alternate = fields[4].to_string();
        let quality = if fields[5] == "." {
            None
        } else {
            Some(fields[5].parse::<f64>()?)
        };
        let filter = fields[6].to_string();
        let info = fields[7].to_string();

        let format = if fields.len() > 8 {
            Some(fields[8].to_string())
        } else {
            None
        };
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

    #[allow(dead_code)]
    pub fn get_info_field(&self, field_name: &str) -> Option<String> {
        for pair in self.info.split(';') {
            if let Some((key, value)) = pair.split_once('=') {
                if key == field_name {
                    return Some(value.to_string());
                }
            } else if pair == field_name {
                return Some("true".to_string());
            }
        }
        None
    }
}

#[derive(Debug, Clone, PartialEq)]
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
    pub mutation_status: String,
    pub validation_status: Option<String>,
    pub sequencer: Option<String>,
    pub sequence_source: String,
    pub depth: Option<u32>,
    pub total_depth: Option<u32>,
    pub vaf: Option<f32>,
    pub hgvsp: Option<String>,
    pub hgvsc: Option<String>,
    // NEW FIELDS
    pub qual: Option<f64>,                // VCF QUAL score
    pub filter_status: String,            // VCF FILTER field
    pub transcript_id: Option<String>,    // Transcript ID from annotations
    pub protein_position: Option<String>, // Protein position from annotations
}

impl MafRecord {
    /// Convert from ReformattedVcfRecord to MafRecord
    pub fn from_reformatted_record(
        record: &ReformattedVcfRecord,
        center: &str,
        ncbi_build: &str,
        sample_barcode: &str,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        // Determine variant type and get proper MAF positions/alleles
        let variant_type = Self::determine_variant_type(&record.reference, &record.alternate);
        let (start_pos, end_pos) = Self::calculate_maf_positions(
            record.position,
            &record.reference,
            &record.alternate,
            &variant_type,
        );
        let (ref_allele, tumor_seq_allele1, tumor_seq_allele2) =
            Self::get_maf_alleles(&record.reference, &record.alternate, &variant_type);

        // Extract depth information
        let mut tumor_depth = Self::extract_tumor_depth(&record.info_fields);
        let mut total_depth = Self::extract_total_depth(&record.info_fields);

        // Fallback to sample data if INFO fields don't have depth
        if tumor_depth.is_none() || total_depth.is_none() {
            let (sample_total, sample_alt) =
                Self::extract_depth_from_sample_data(&record.format_sample_data);
            if tumor_depth.is_none() {
                tumor_depth = sample_alt;
            }
            if total_depth.is_none() {
                total_depth = sample_total;
            }
        }

        // Calculate VAF if we have both tumor and total depth
        let vaf = match (tumor_depth, total_depth) {
            (Some(t_depth), Some(total)) if total > 0 => Some(t_depth as f32 / total as f32),
            _ => None,
        };

        Ok(MafRecord {
            hugo_symbol: Self::get_annotation_field(
                &record.info_fields,
                &[
                    "ANN_Gene_Name", // SnpEff primary
                    "CSQ_SYMBOL",
                    "ANN_SYMBOL",
                    "SYMBOL",
                    "Gene_Name",
                    "Gene",
                ],
            )
            .unwrap_or(".".to_string()),
            entrez_gene_id: Self::get_entrez_gene_id(&record.info_fields),
            center: if center.trim().is_empty() {
                "Unknown_Center".to_string()
            } else {
                center.to_string()
            },
            ncbi_build: ncbi_build.to_string(),
            chromosome: Self::normalize_chromosome(&record.chromosome), // Now normalize consistently
            start_position: start_pos,                                  // Use MAF positions
            end_position: end_pos,                                      // Use MAF positions
            strand: Self::get_strand(&record.info_fields),
            variant_classification: Self::get_variant_classification(&record.info_fields),
            variant_type,
            reference_allele: ref_allele, // Use MAF alleles
            tumor_seq_allele1,            // Use MAF alleles
            tumor_seq_allele2,            // Use MAF alleles
            dbsnp_rs: record.id.clone().filter(|id| id != "."),
            dbsnp_val_status: None,
            tumor_sample_barcode: sample_barcode.to_string(),
            matched_norm_sample_barcode: None,
            mutation_status: "Somatic".to_string(),
            validation_status: None,
            sequencer: Self::extract_sequencing_info(&record.info_fields),
            sequence_source: "WXS".to_string(),
            depth: tumor_depth,
            total_depth,
            vaf,
            hgvsp: Self::get_annotation_field(
                &record.info_fields,
                &["CSQ_HGVSp", "ANN_HGVS_p", "HGVSp", "HGVS_p"],
            )
            .filter(|s| s != "."),
            hgvsc: Self::get_annotation_field(
                &record.info_fields,
                &["CSQ_HGVSc", "ANN_HGVS_c", "HGVSc", "HGVS_c"],
            )
            .filter(|s| s != "."),
            qual: record.quality,
            filter_status: record.filter.clone(),
            transcript_id: Self::get_transcript_id(&record.info_fields),
            protein_position: Self::get_protein_position(&record.info_fields),
        })
    }

    fn get_transcript_id(info_fields: &HashMap<String, String>) -> Option<String> {
        Self::get_annotation_field(
            info_fields,
            &[
                "ANN_Feature_ID", // SnpEff primary
                "CSQ_Feature",    // VEP
                "CSQ_Transcript_ID",
                "ANN_Transcript_ID",
            ],
        )
    }

    fn get_protein_position(info_fields: &HashMap<String, String>) -> Option<String> {
        Self::get_annotation_field(
            info_fields,
            &[
                "ANN_AA_pos___AA_length", // SnpEff format from your data
                "ANN_Protein_position",
                "CSQ_Protein_position", // VEP
                "ANN_AA_pos",
                "CSQ_AA_pos",
            ],
        )
    }

    /// Handle multi-allelic variants by creating separate MafRecord for each alternate allele
    pub fn from_reformatted_record_multi(
        record: &ReformattedVcfRecord,
        center: &str,
        ncbi_build: &str,
        sample_barcode: &str,
    ) -> Result<Vec<Self>, Box<dyn std::error::Error>> {
        let alternates: Vec<&str> = record.alternate.split(',').collect();
        let mut maf_records = Vec::new();

        for (alt_index, alternate) in alternates.iter().enumerate() {
            // Create a record for each alternate allele
            let single_alt_record = ReformattedVcfRecord {
                chromosome: record.chromosome.clone(),
                position: record.position,
                id: record.id.clone(),
                reference: record.reference.clone(),
                alternate: alternate.to_string(),
                quality: record.quality,
                filter: record.filter.clone(),
                info_fields: record.info_fields.clone(),
                format_sample_data: record.format_sample_data.clone(),
            };

            // Use the unified conversion logic
            let mut maf_record = Self::from_reformatted_record(
                &single_alt_record,
                center,
                ncbi_build,
                sample_barcode,
            )?;

            // Adjust depth for specific allele if available
            if let Some(allele_depth) =
                Self::extract_tumor_depth_for_allele(&record.info_fields, alt_index)
            {
                maf_record.depth = Some(allele_depth);

                // Recalculate VAF with allele-specific depth
                if let Some(total) = maf_record.total_depth {
                    if total > 0 {
                        maf_record.vaf = Some(allele_depth as f32 / total as f32);
                    }
                }
            }

            maf_records.push(maf_record);
        }

        Ok(maf_records)
    }

    // Extract tumor depth for specific allele index
    fn extract_tumor_depth_for_allele(
        info_fields: &HashMap<String, String>,
        allele_index: usize,
    ) -> Option<u32> {
        if let Some(ao_str) = Self::get_annotation_field(info_fields, &["INFO_AO", "AO"]) {
            let depths: Vec<&str> = ao_str.split(',').collect();
            if let Some(depth_str) = depths.get(allele_index) {
                if let Ok(depth) = depth_str.parse::<u32>() {
                    return Some(depth);
                }
            }
        }

        // Fallback to total tumor depth
        Self::extract_tumor_depth(info_fields)
    }

    fn get_strand(info_fields: &HashMap<String, String>) -> String {
        if let Some(strand) = Self::get_annotation_field(
            info_fields,
            &[
                "CSQ_STRAND", // VEP
                "ANN_Strand", // SnpEff (rare)
                "STRAND",     // Generic
            ],
        ) {
            match strand.as_str() {
                "1" | "+" => "+".to_string(),
                "-1" | "-" => "-".to_string(),
                _ => "+".to_string(),
            }
        } else {
            "*".to_string()
        }
    }

    // Extract Entrez Gene ID from HGNC if available
    fn get_entrez_gene_id(info_fields: &HashMap<String, String>) -> Option<u32> {
        // Try HGNC first (VEP style)
        if let Some(hgnc) = Self::get_annotation_field(info_fields, &["CSQ_HGNC_ID"]) {
            if let Some(id) = hgnc.strip_prefix("HGNC:").and_then(|id| id.parse().ok()) {
                return Some(id);
            }
        }

        Self::get_annotation_field(
            info_fields,
            &[
                "ENTREZ_GENE_ID", // Direct field (rare)
                "ANN_Entrez_ID",  // SnpEff Entrez (rare)
                "CSQ_Gene",       // VEP Gene ID
                "Gene_ID",        // Generic
            ],
        )
        .and_then(|id| id.parse().ok())
    }

    fn extract_depth_from_sample_data(
        format_sample_data: &Option<ParsedFormatSample>,
    ) -> (Option<u32>, Option<u32>) {
        if let Some(sample_data) = format_sample_data {
            let mut total_depth = None;
            let mut alt_depth = None;

            // Check each sample for depth information
            for sample in &sample_data.samples {
                // Total depth (DP field)
                if let Some(dp) = sample.format_fields.get("DP") {
                    if let Ok(depth) = dp.parse::<u32>() {
                        total_depth = Some(depth);
                    }
                }

                // Alternative depth (AD field - usually comma-separated: ref,alt)
                if let Some(ad) = sample.format_fields.get("AD") {
                    let depths: Vec<&str> = ad.split(',').collect();
                    if depths.len() >= 2 {
                        if let Ok(alt_d) = depths[1].parse::<u32>() {
                            alt_depth = Some(alt_d);
                        }
                    }
                }

                // If we found depth in this sample, use it (typically first sample is tumor)
                if total_depth.is_some() {
                    break;
                }
            }

            (total_depth, alt_depth)
        } else {
            (None, None)
        }
    }

    fn extract_total_depth(info_fields: &HashMap<String, String>) -> Option<u32> {
        // Try different field names for total depth
        Self::get_annotation_field(
            info_fields,
            &[
                "INFO_DP",
                "DP",
                "TotalDepth",
                "DEPTH",
                "Total_Depth",
                "INFO_DEPTH",
                // SnpEff specific fields
                "ANN_DP",
                "ANN_TotalDepth",
                // Sample-based depth (might be prefixed)
                "SAMPLE_DP",
                "FORMAT_DP",
            ],
        )
        .and_then(|s| s.parse().ok())
    }

    fn extract_tumor_depth(info_fields: &HashMap<String, String>) -> Option<u32> {
        // Try different field names for tumor/alternative allele depth
        Self::get_annotation_field(
            info_fields,
            &[
                "INFO_AO",
                "AO",
                "AD_ALT",
                "t_alt_count",
                "TumorDepth",
                "ALT_DEPTH",
                // SnpEff specific fields
                "ANN_AO",
                "ANN_AD",
                // Sample-based fields
                "SAMPLE_AO",
                "FORMAT_AO",
                "FORMAT_AD",
                // Alternative depth representations
                "ALT_COUNT",
                "TUMOR_ALT_COUNT",
            ],
        )
        .and_then(|s| {
            // Handle comma-separated values by taking the first one
            let first_value = s.split(',').next().unwrap_or(&s);
            first_value.parse().ok()
        })
    }

    /// Extract sequencing platform or variant calling tool info
    fn extract_sequencing_info(info_fields: &HashMap<String, String>) -> Option<String> {
        // Look for various tool/platform indicators in INFO fields
        for field_name in &[
            "INFO_source",
            "INFO_caller",
            "INFO_platform",
            "source",
            "caller",
        ] {
            if let Some(value) = info_fields.get(*field_name) {
                if !value.is_empty() && value != "." {
                    return Some(value.clone());
                }
            }
        }

        // If no specific field found, return None instead of defaulting
        None
    }

    fn get_annotation_field(
        info_fields: &HashMap<String, String>,
        field_names: &[&str],
    ) -> Option<String> {
        for field_name in field_names {
            if let Some(value) = info_fields.get(*field_name) {
                if !value.is_empty() && value != "." {
                    return Some(value.clone());
                }
            }
        }
        None
    }

    fn calculate_maf_positions(
        position: u64,
        reference: &str,
        _alternate: &str,
        variant_type: &str,
    ) -> (u64, u64) {
        match variant_type {
            "INS" => (position, position + 1),
            "DEL" => {
                let end_pos = position + reference.len() as u64 - 1;
                (position, end_pos.max(position))
            }
            _ => (position, position),
        }
    }

    fn get_maf_alleles(
        reference: &str,
        alternate: &str,
        variant_type: &str,
    ) -> (String, String, String) {
        match variant_type {
            "INS" => ("-".to_string(), "-".to_string(), alternate.to_string()),
            "DEL" => (
                reference.to_string(),
                reference.to_string(),
                "-".to_string(),
            ),
            _ => (
                reference.to_string(),
                reference.to_string(),
                alternate.to_string(),
            ),
        }
    }

    fn determine_variant_type(reference: &str, alternate: &str) -> String {
        if reference.len() == 1 && alternate.len() == 1 {
            "SNP".to_string()
        } else if reference.len() < alternate.len() {
            "INS".to_string()
        } else if reference.len() > alternate.len() {
            "DEL".to_string()
        } else {
            "ONP".to_string()
        }
    }

    fn get_variant_classification(info_fields: &HashMap<String, String>) -> String {
        if let Some(consequence) = Self::get_annotation_field(
            info_fields,
            &[
                "CSQ_Consequence",
                "ANN_Annotation",
                "ANN_Annotation_Impact",
                "Consequence",
                "Annotation",
                "Effect",
                "Variant_Effect",
            ],
        ) {
            Self::map_consequence_to_maf(&consequence, info_fields)
        } else {
            // Fallback: check IMPACT field directly
            if let Some(impact) = Self::get_annotation_field(
                info_fields,
                &["CSQ_IMPACT", "ANN_Annotation_Impact", "IMPACT"],
            ) {
                match impact.to_uppercase().as_str() {
                    "HIGH" => "Missense_Mutation".to_string(),
                    "MODERATE" => "Missense_Mutation".to_string(),
                    "LOW" => "Silent".to_string(),
                    "MODIFIER" => "Silent".to_string(),
                    _ => "Unknown".to_string(),
                }
            } else {
                "Unknown".to_string()
            }
        }
    }

    fn map_consequence_to_maf(consequence: &str, info_fields: &HashMap<String, String>) -> String {
        let consequence_lower = consequence.to_lowercase();

        // Handle multiple consequences separated by &, |, or ,
        let consequences: Vec<&str> = consequence_lower
            .split(&['&', '|', ','][..])
            .map(|s| s.trim())
            .collect();

        // Find the most severe consequence
        for cons in &consequences {
            // High impact
            if cons.contains("stop_gained") || cons.contains("nonsense") {
                return "Nonsense_Mutation".to_string();
            } else if cons.contains("missense") || cons.contains("missense_variant") {
                return "Missense_Mutation".to_string();
            } else if cons.contains("frameshift") {
                return "Frame_Shift_Del".to_string();
            } else if cons.contains("splice")
                && (cons.contains("acceptor") || cons.contains("donor"))
            {
                return "Splice_Site".to_string();
            }
            // Medium impact variants
            else if cons.contains("inframe_insertion") {
                return "In_Frame_Ins".to_string();
            } else if cons.contains("inframe_deletion") {
                return "In_Frame_Del".to_string();
            } else if cons.contains("stop_lost") {
                return "Nonstop_Mutation".to_string();
            } else if cons.contains("start_lost") {
                return "Translation_Start_Site".to_string();
            }
            // Low imp
            else if cons.contains("synonymous") || cons.contains("silent") {
                return "Silent".to_string();
            } else if cons.contains("5_prime_utr")
                || cons.contains("3_prime_utr")
                || cons.contains("utr")
            {
                return "3'UTR".to_string();
            }
            // Splice region variants
            else if cons.contains("splice_region") {
                return "Splice_Region".to_string();
            }
            // Intronic and intergenic
            else if cons.contains("intronic") || cons.contains("intron") {
                return "Intron".to_string();
            } else if cons.contains("intergenic") {
                return "IGR".to_string();
            }
        }

        // If no specific consequence matched, check VEP IMPACT field as fallback
        if let Some(impact) = Self::get_annotation_field(
            info_fields,
            &["CSQ_IMPACT", "ANN_Annotation_Impact", "IMPACT"],
        ) {
            match impact.to_uppercase().as_str() {
                "HIGH" => "Missense_Mutation".to_string(),
                "MODERATE" => "Missense_Mutation".to_string(),
                "LOW" => "Silent".to_string(),
                "MODIFIER" => "Silent".to_string(),
                _ => "Unknown".to_string(),
            }
        } else {
            "Unknown".to_string()
        }
    }

    fn normalize_chromosome(chr: &str) -> String {
        chr.trim_start_matches("chr").to_string()
    }

    pub fn get_maf_headers() -> Vec<String> {
        vec![
            "Hugo_Symbol".to_string(),
            "Entrez_Gene_Id".to_string(),
            "Center".to_string(),
            "NCBI_Build".to_string(),
            "Chromosome".to_string(),
            "Start_Position".to_string(),
            "End_Position".to_string(),
            "Strand".to_string(),
            "Variant_Classification".to_string(),
            "Variant_Type".to_string(),
            "Reference_Allele".to_string(),
            "Tumor_Seq_Allele1".to_string(),
            "Tumor_Seq_Allele2".to_string(),
            "dbSNP_RS".to_string(),
            "dbSNP_Val_Status".to_string(),
            "Tumor_Sample_Barcode".to_string(),
            "Matched_Norm_Sample_Barcode".to_string(),
            "Mutation_Status".to_string(),
            "Validation_Status".to_string(),
            "Sequencer".to_string(),
            "Sequence_Source".to_string(),
            "t_depth".to_string(),
            "total_depth".to_string(),
            "VAF".to_string(),
            "HGVSp".to_string(),
            "HGVSc".to_string(),
            // NEW HEADERS
            "QUAL".to_string(),             // VCF quality score
            "FILTER".to_string(),           // VCF filter status
            "Transcript_ID".to_string(),    // Transcript identifier
            "Protein_Position".to_string(), // Protein position
        ]
    }

    pub fn to_tsv_line(&self) -> String {
        vec![
            self.hugo_symbol.clone(),
            self.entrez_gene_id
                .map_or(".".to_string(), |id| id.to_string()),
            self.center.clone(),
            self.ncbi_build.clone(),
            self.chromosome.clone(),
            self.start_position.to_string(),
            self.end_position.to_string(),
            self.strand.clone(),
            self.variant_classification.clone(),
            self.variant_type.clone(),
            self.reference_allele.clone(),
            self.tumor_seq_allele1.clone(),
            self.tumor_seq_allele2.clone(),
            self.dbsnp_rs.as_ref().unwrap_or(&".".to_string()).clone(),
            self.dbsnp_val_status
                .as_ref()
                .unwrap_or(&".".to_string())
                .clone(),
            self.tumor_sample_barcode.clone(),
            self.matched_norm_sample_barcode
                .as_ref()
                .unwrap_or(&".".to_string())
                .clone(),
            self.mutation_status.clone(),
            self.validation_status
                .as_ref()
                .unwrap_or(&".".to_string())
                .clone(),
            self.sequencer.as_ref().unwrap_or(&".".to_string()).clone(),
            self.sequence_source.clone(),
            self.depth.map_or(".".to_string(), |d| d.to_string()),
            self.total_depth.map_or(".".to_string(), |d| d.to_string()),
            self.vaf.map_or(".".to_string(), |v| format!("{:.4}", v)),
            self.hgvsp.as_ref().unwrap_or(&".".to_string()).clone(),
            self.hgvsc.as_ref().unwrap_or(&".".to_string()).clone(),
            // NEW FIELDS OUTPUT
            self.qual.map_or(".".to_string(), |q| q.to_string()), // QUAL score
            self.filter_status.clone(),                           // FILTER status
            self.transcript_id
                .as_ref()
                .unwrap_or(&".".to_string())
                .clone(), // Transcript ID
            self.protein_position
                .as_ref()
                .unwrap_or(&".".to_string())
                .clone(), // Protein position
        ]
        .join("\t")
    }
}
