use crate::extract_sample_info::ParsedFormatSample;
use crate::reformat_vcf::ReformattedVcfRecord;
use std::collections::HashMap;

type ClassificationRule = (fn(&str) -> bool, &'static str);

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
                &["ANN_Gene_Name", "CSQ_SYMBOL", "ANN_SYMBOL"],
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
            variant_classification: Self::get_variant_classification(&record.info_fields, &variant_type),
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
            hgvsp: Self::get_annotation_field(&record.info_fields, &["CSQ_HGVSp", "ANN_HGVS_p"])
                .filter(|s| s != "."),
            hgvsc: Self::get_annotation_field(&record.info_fields, &["CSQ_HGVSc", "ANN_HGVS_c"])
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
                annotation_field_type: record.annotation_field_type,
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
        if let Some(strand) = Self::get_annotation_field(info_fields, &["CSQ_STRAND", "ANN_Strand"]) {
            match strand.as_str() {
                "1" | "+" => "+".to_string(),
                "-1" | "-" => "-".to_string(),
                _ => "+".to_string(),
            }
        } else {
            "+".to_string() // MAF spec requires + or -, default to +
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

        Self::get_annotation_field(info_fields, &["ANN_Entrez_ID", "CSQ_Gene"])
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
        Self::get_annotation_field(info_fields, &["INFO_DP", "INFO_DEPTH", "ANN_DP", "ANN_TotalDepth"])
            .and_then(|s| s.parse().ok())
    }

    fn extract_tumor_depth(info_fields: &HashMap<String, String>) -> Option<u32> {
        Self::get_annotation_field(info_fields, &["INFO_AO", "ANN_AO", "ANN_AD"]).and_then(|s| {
            // Handle comma-separated values by taking the first one
            let first_value = s.split(',').next().unwrap_or(&s);
            first_value.parse().ok()
        })
    }

    /// Extract sequencing platform or variant calling tool info
    fn extract_sequencing_info(info_fields: &HashMap<String, String>) -> Option<String> {
        for field_name in &["INFO_source", "INFO_caller", "INFO_platform"] {
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
        alternate: &str,
        variant_type: &str,
    ) -> (u64, u64) {
        // VCF indels have a shared anchor base prefix. MAF positions refer to
        // the actual inserted/deleted bases, not the anchor.
        // Note: VCF alleles are ASCII (A/C/G/T), so chars().count() == byte length.
        let shared_prefix_len = reference
            .chars()
            .zip(alternate.chars())
            .take_while(|(r, a)| r == a)
            .count() as u64;

        match variant_type {
            "INS" => {
                // Insertion: start = last shared base, end = start + 1
                // Guard against malformed VCF with no anchor base (shared_prefix_len == 0)
                let start = if shared_prefix_len > 0 {
                    position + shared_prefix_len - 1
                } else {
                    position
                };
                (start, start + 1)
            }
            "DEL" => {
                // Deletion: start = first deleted base, end = last deleted base
                let deleted_len = reference.len() as u64 - shared_prefix_len;
                if deleted_len == 0 {
                    return (position, position);
                }
                let start = position + shared_prefix_len;
                let end = start + deleted_len - 1;
                (start, end.max(start))
            }
            _ => (position, position),
        }
    }

    fn get_maf_alleles(
        reference: &str,
        alternate: &str,
        variant_type: &str,
    ) -> (String, String, String) {
        // VCF indels include a shared anchor base that MAF strips out.
        // E.g., VCF ref=A alt=ATCG → MAF ref="-" alt="TCG"
        //       VCF ref=ATCG alt=A → MAF ref="TCG" alt="-"
        let shared_prefix_len = reference
            .chars()
            .zip(alternate.chars())
            .take_while(|(r, a)| r == a)
            .count();

        match variant_type {
            "INS" => {
                let inserted = &alternate[shared_prefix_len..];
                ("-".to_string(), "-".to_string(), inserted.to_string())
            }
            "DEL" => {
                let deleted = &reference[shared_prefix_len..];
                (deleted.to_string(), deleted.to_string(), "-".to_string())
            }
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

    fn get_variant_classification(info_fields: &HashMap<String, String>, variant_type: &str) -> String {
        match Self::get_annotation_field(info_fields, &["CSQ_Consequence", "ANN_Annotation"]) {
            Some(consequence) => Self::map_consequence_to_maf(&consequence, info_fields, variant_type),
            None => Self::classify_by_impact(info_fields),
        }
    }

    /// Ordered by severity (HIGH > MODERATE > LOW > MODIFIER). `frameshift` is
    /// handled separately since its classification depends on `variant_type`.
    const CLASSIFICATION_RULES: &'static [ClassificationRule] = &[
        (|c| c.contains("stop_gained") || c.contains("nonsense"), "Nonsense_Mutation"),
        (|c| c.contains("splice") && (c.contains("acceptor") || c.contains("donor")), "Splice_Site"),
        (|c| c.contains("start_lost") || c.contains("initiator_codon_variant"), "Translation_Start_Site"),
        (|c| c.contains("stop_lost"), "Nonstop_Mutation"),
        (|c| c.contains("missense") || c.contains("rare_amino_acid_variant"), "Missense_Mutation"),
        (|c| c.contains("inframe_insertion"), "In_Frame_Ins"),
        (|c| c.contains("inframe_deletion"), "In_Frame_Del"),
        (|c| c.contains("synonymous") || c.contains("silent") || c.contains("stop_retained_variant"), "Silent"),
        (|c| c.contains("splice_region"), "Splice_Region"),
        // Specific UTR types are checked before anything more generic.
        (|c| c.contains("5_prime_utr"), "5'UTR"),
        (|c| c.contains("3_prime_utr"), "3'UTR"),
        (|c| c.contains("upstream_gene_variant"), "5'Flank"),
        (|c| c.contains("downstream_gene_variant"), "3'Flank"),
        (|c| c.contains("non_coding_transcript"), "RNA"),
        (|c| c.contains("regulatory_region") || c.contains("tf_binding_site"), "Targeted_Region"),
        (|c| c.contains("intronic") || c.contains("intron"), "Intron"),
        (|c| c.contains("intergenic"), "IGR"),
    ];

    fn map_consequence_to_maf(consequence: &str, info_fields: &HashMap<String, String>, variant_type: &str) -> String {
        let consequence_lower = consequence.to_lowercase();

        // A variant can carry multiple consequences separated by &, |, or ,;
        // the first one that matches a rule determines severity.
        for cons in consequence_lower.split(&['&', '|', ','][..]).map(|s| s.trim()) {
            if cons.contains("frameshift") {
                let classification = if variant_type == "INS" { "Frame_Shift_Ins" } else { "Frame_Shift_Del" };
                return classification.to_string();
            }
            if let Some((_, classification)) =
                Self::CLASSIFICATION_RULES.iter().find(|(predicate, _)| predicate(cons))
            {
                return classification.to_string();
            }
        }

        Self::classify_by_impact(info_fields)
    }

    fn classify_by_impact(info_fields: &HashMap<String, String>) -> String {
        let impact = Self::get_annotation_field(info_fields, &["CSQ_IMPACT", "ANN_Annotation_Impact"]);
        match impact.as_deref().map(|s| s.to_uppercase()).as_deref() {
            Some("HIGH") | Some("MODERATE") => "Missense_Mutation".to_string(),
            Some("LOW") | Some("MODIFIER") => "Silent".to_string(),
            _ => "Unknown".to_string(),
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
        let dot = ".";
        let entrez = self.entrez_gene_id.map(|id| id.to_string());
        let start = self.start_position.to_string();
        let end = self.end_position.to_string();
        let depth = self.depth.map(|d| d.to_string());
        let total_depth = self.total_depth.map(|d| d.to_string());
        let vaf = self.vaf.map(|v| format!("{:.4}", v));
        let qual = self.qual.map(|q| q.to_string());

        [
            self.hugo_symbol.as_str(),
            entrez.as_deref().unwrap_or(dot),
            self.center.as_str(),
            self.ncbi_build.as_str(),
            self.chromosome.as_str(),
            start.as_str(),
            end.as_str(),
            self.strand.as_str(),
            self.variant_classification.as_str(),
            self.variant_type.as_str(),
            self.reference_allele.as_str(),
            self.tumor_seq_allele1.as_str(),
            self.tumor_seq_allele2.as_str(),
            self.dbsnp_rs.as_deref().unwrap_or(dot),
            self.dbsnp_val_status.as_deref().unwrap_or(dot),
            self.tumor_sample_barcode.as_str(),
            self.matched_norm_sample_barcode.as_deref().unwrap_or(dot),
            self.mutation_status.as_str(),
            self.validation_status.as_deref().unwrap_or(dot),
            self.sequencer.as_deref().unwrap_or(dot),
            self.sequence_source.as_str(),
            depth.as_deref().unwrap_or(dot),
            total_depth.as_deref().unwrap_or(dot),
            vaf.as_deref().unwrap_or(dot),
            self.hgvsp.as_deref().unwrap_or(dot),
            self.hgvsc.as_deref().unwrap_or(dot),
            qual.as_deref().unwrap_or(dot),
            self.filter_status.as_str(),
            self.transcript_id.as_deref().unwrap_or(dot),
            self.protein_position.as_deref().unwrap_or(dot),
        ]
        .join("\t")
    }
}
