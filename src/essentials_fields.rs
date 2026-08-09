use crate::extract_sample_info::ParsedFormatSample;
use crate::reformat_vcf::ReformattedVcfRecord;
use std::collections::HashMap;

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
        let inframe = record.reference.len().abs_diff(record.alternate.len()) % 3 == 0;
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
        let mut tumor_ref_depth = Self::extract_ref_depth(&record.info_fields);

        // Fallback to sample data if INFO fields don't have depth
        if tumor_depth.is_none() || total_depth.is_none() || tumor_ref_depth.is_none() {
            let (sample_total, sample_ref, sample_alt) =
                Self::extract_depth_from_sample_data(&record.format_sample_data);
            if tumor_depth.is_none() {
                tumor_depth = sample_alt;
            }
            if total_depth.is_none() {
                total_depth = sample_total;
            }
            if tumor_ref_depth.is_none() {
                tumor_ref_depth = sample_ref;
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
            variant_classification: Self::get_variant_classification(
                &record.info_fields,
                &variant_type,
                inframe,
            ),
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

    /// VEP's `EXON` and SnpEff's `Rank` are both formatted `<rank>/<total>` and are
    /// vcf2maf's source for the MAF `Exon_Number` column.
    fn get_exon_number(info_fields: &HashMap<String, String>) -> Option<String> {
        Self::get_annotation_field(info_fields, &["CSQ_EXON", "ANN_Rank"])
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
    ) -> (Option<u32>, Option<u32>, Option<u32>) {
        if let Some(sample_data) = format_sample_data {
            let mut total_depth = None;
            let mut ref_depth = None;
            let mut alt_depth = None;

            // Check each sample for depth information
            for sample in &sample_data.samples {
                // Total depth (DP field)
                if let Some(dp) = sample.format_fields.get("DP") {
                    if let Ok(depth) = dp.parse::<u32>() {
                        total_depth = Some(depth);
                    }
                }

                // Reference/alternative depth (AD field - usually comma-separated: ref,alt)
                if let Some(ad) = sample.format_fields.get("AD") {
                    let depths: Vec<&str> = ad.split(',').collect();
                    if let Some(r) = depths.first().and_then(|d| d.parse::<u32>().ok()) {
                        ref_depth = Some(r);
                    }
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

            (total_depth, ref_depth, alt_depth)
        } else {
            (None, None, None)
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

    /// freebayes-style INFO `RO` (Reference Observation count) is this tool's only
    /// direct source for t_ref_count; sample-level `AD[0]` is the fallback (see
    /// `extract_depth_from_sample_data`), mirroring how `extract_tumor_depth` already
    /// falls back from INFO `AO` to sample `AD[1]`.
    fn extract_ref_depth(info_fields: &HashMap<String, String>) -> Option<u32> {
        Self::get_annotation_field(info_fields, &["INFO_RO"]).and_then(|s| s.parse().ok())
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
        if reference.len() < alternate.len() {
            "INS".to_string()
        } else if reference.len() > alternate.len() {
            "DEL".to_string()
        } else {
            match reference.len() {
                1 => "SNP",
                2 => "DNP",
                3 => "TNP",
                _ => "ONP",
            }
            .to_string()
        }
    }

    fn get_variant_classification(
        info_fields: &HashMap<String, String>,
        variant_type: &str,
        inframe: bool,
    ) -> String {
        match Self::get_annotation_field(info_fields, &["CSQ_Consequence", "ANN_Annotation"]) {
            Some(consequence) => Self::map_consequence_to_maf(&consequence, variant_type, inframe),
            None => Self::classify_by_impact(info_fields),
        }
    }

    /// Sequence Ontology term severity ranks, ported from vcf2maf's `%effectPriority`
    /// (mskcc/vcf2maf, GetEffectPriority) so consequence resolution matches the reference
    /// tool exactly rather than picking whichever `&`-joined term happens to appear first.
    /// Lower number = more severe. Unrecognized terms default to 20, same as vcf2maf.
    const EFFECT_PRIORITY: &'static [(&'static str, u8)] = &[
        ("transcript_ablation", 1),
        ("exon_loss_variant", 1),
        ("sequence_feature + exon_loss_variant", 1),
        ("feature_ablation", 1),
        ("chromosome_number_variation", 1),
        ("bidirectional_gene_fusion", 2),
        ("duplication", 2),
        ("gene_fusion", 2),
        ("inversion", 2),
        ("splice_donor_variant", 2),
        ("splice_acceptor_variant", 2),
        ("stop_gained", 3),
        ("frameshift_variant", 3),
        ("stop_lost", 3),
        ("initiator_codon_variant+non_canonical_start_codon", 4),
        ("rearranged_at_dna_level", 4),
        ("start_lost", 4),
        ("initiator_codon_variant", 4),
        ("transcript_amplification", 4),
        ("feature_elongation", 4),
        ("feature_truncation", 4),
        ("disruptive_inframe_insertion", 5),
        ("disruptive_inframe_deletion", 5),
        ("conservative_inframe_insertion", 5),
        ("conservative_inframe_deletion", 5),
        ("inframe_insertion", 5),
        ("inframe_deletion", 5),
        ("protein_altering_variant", 5),
        ("missense_variant", 6),
        ("conservative_missense_variant", 6),
        ("rare_amino_acid_variant", 6),
        ("5_prime_utr_truncation + exon_loss_variant", 8),
        ("protein_protein_contact", 8),
        ("3_prime_utr_truncation + exon_loss", 8),
        ("structural_interaction_variant", 8),
        ("splice_branch_variant", 8),
        ("splice_region_variant", 8),
        ("splice_donor_5th_base_variant", 8),
        ("splice_donor_region_variant", 8),
        ("splice_polypyrimidine_tract_variant", 8),
        ("start_retained_variant", 9),
        ("stop_retained_variant", 9),
        ("synonymous_variant", 9),
        ("start_retained", 9),
        ("incomplete_terminal_codon_variant", 10),
        ("coding_sequence_variant", 11),
        ("mature_mirna_variant", 11),
        ("exon_variant", 11),
        ("transcript_variant", 11),
        ("5_prime_utr_variant", 12),
        ("5_prime_utr_premature_start_codon_gain_variant", 12),
        ("3_prime_utr_variant", 12),
        ("non_coding_exon_variant", 13),
        ("non_coding_transcript_exon_variant", 13),
        ("non_coding_transcript_variant", 14),
        ("nc_transcript_variant", 14),
        ("intron_variant", 14),
        ("intragenic_variant", 14),
        ("intragenic", 14),
        ("nmd_transcript_variant", 15),
        ("coding_transcript_variant", 15),
        ("upstream_gene_variant", 16),
        ("downstream_gene_variant", 16),
        ("tfbs_ablation", 17),
        ("tfbs_amplification", 17),
        ("tf_binding_site_variant", 17),
        ("regulatory_region_ablation", 17),
        ("regulatory_region_amplification", 17),
        ("regulatory_region_variant", 17),
        ("regulatory_region", 17),
        ("mirna", 17),
        ("intergenic_variant", 19),
        ("intergenic_region", 19),
        ("sequence_feature", 19),
        ("conserved_intron_variant", 19),
        ("gene_variant", 19),
        ("conserved_intergenic_variant", 20),
        ("sequence_variant", 20),
        ("custom", 20),
    ];

    fn effect_priority(term: &str) -> u8 {
        Self::EFFECT_PRIORITY
            .iter()
            .find(|(name, _)| *name == term)
            .map(|(_, priority)| *priority)
            .unwrap_or(20)
    }

    /// Resolve a (possibly `&`/`,`/`|`-joined) multi-term consequence string down to the
    /// single most severe term, mirroring vcf2maf's sort-by-priority-then-take-first.
    fn resolve_one_consequence(consequence: &str) -> String {
        consequence
            .split(&['&', '|', ','][..])
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .min_by_key(|term| Self::effect_priority(term))
            .unwrap_or("intergenic_variant")
            .to_string()
    }

    /// Ported from vcf2maf's `GetVariantClassification`: classifies the single most-severe
    /// consequence term, using variant_type/inframe only to disambiguate the terms whose MAF
    /// classification depends on them (frameshift_variant, protein_altering_variant).
    fn map_consequence_to_maf(consequence: &str, variant_type: &str, inframe: bool) -> String {
        let term = Self::resolve_one_consequence(&consequence.to_lowercase());

        if matches!(term.as_str(), "splice_acceptor_variant" | "splice_donor_variant") {
            return "Splice_Site".to_string();
        }
        if term == "stop_gained" {
            return "Nonsense_Mutation".to_string();
        }
        let is_frameshift_like = term == "frameshift_variant"
            || (term == "protein_altering_variant" && !inframe);
        if is_frameshift_like && variant_type == "DEL" {
            return "Frame_Shift_Del".to_string();
        }
        if is_frameshift_like && variant_type == "INS" {
            return "Frame_Shift_Ins".to_string();
        }
        if term == "stop_lost" {
            return "Nonstop_Mutation".to_string();
        }
        if matches!(term.as_str(), "initiator_codon_variant" | "start_lost") {
            return "Translation_Start_Site".to_string();
        }
        let is_inframe_ins = term.ends_with("inframe_insertion")
            || (term == "protein_altering_variant" && inframe && variant_type == "INS");
        if is_inframe_ins {
            return "In_Frame_Ins".to_string();
        }
        let is_inframe_del = term.ends_with("inframe_deletion")
            || (term == "protein_altering_variant" && inframe && variant_type == "DEL");
        if is_inframe_del {
            return "In_Frame_Del".to_string();
        }
        if matches!(
            term.as_str(),
            "missense_variant" | "coding_sequence_variant" | "conservative_missense_variant" | "rare_amino_acid_variant"
        ) {
            return "Missense_Mutation".to_string();
        }
        if matches!(
            term.as_str(),
            "transcript_amplification" | "intron_variant" | "intragenic" | "intragenic_variant"
        ) {
            return "Intron".to_string();
        }
        if matches!(
            term.as_str(),
            "splice_region_variant"
                | "splice_donor_5th_base_variant"
                | "splice_donor_region_variant"
                | "splice_polypyrimidine_tract_variant"
        ) {
            return "Splice_Region".to_string();
        }
        if matches!(
            term.as_str(),
            "incomplete_terminal_codon_variant"
                | "synonymous_variant"
                | "stop_retained_variant"
                | "start_retained_variant"
                | "nmd_transcript_variant"
        ) {
            return "Silent".to_string();
        }
        if matches!(
            term.as_str(),
            "mature_mirna_variant"
                | "exon_variant"
                | "non_coding_exon_variant"
                | "non_coding_transcript_exon_variant"
                | "non_coding_transcript_variant"
                | "nc_transcript_variant"
        ) {
            return "RNA".to_string();
        }
        if matches!(
            term.as_str(),
            "5_prime_utr_variant" | "5_prime_utr_premature_start_codon_gain_variant"
        ) {
            return "5'UTR".to_string();
        }
        if term == "3_prime_utr_variant" {
            return "3'UTR".to_string();
        }
        if matches!(
            term.as_str(),
            "tf_binding_site_variant"
                | "regulatory_region_variant"
                | "regulatory_region"
                | "intergenic_variant"
                | "intergenic_region"
        ) {
            return "IGR".to_string();
        }
        if term == "upstream_gene_variant" {
            return "5'Flank".to_string();
        }
        if term == "downstream_gene_variant" {
            return "3'Flank".to_string();
        }

        // Everything else (TFBS/regulatory ablation/amplification, feature_elongation/
        // truncation, coding_transcript_variant, sequence_variant, ...): vcf2maf's own
        // catch-all.
        "Targeted_Region".to_string()
    }

    /// 3-letter → 1-letter amino acid code table, verbatim from vcf2maf's `%aa3to1`
    /// (mskcc/vcf2maf, vcf2maf.pl) so HGVSp_Short matches the reference tool exactly.
    const AA_3_TO_1: &'static [(&'static str, &'static str)] = &[
        ("Ala", "A"), ("Arg", "R"), ("Asn", "N"), ("Asp", "D"), ("Asx", "B"),
        ("Cys", "C"), ("Glu", "E"), ("Gln", "Q"), ("Glx", "Z"), ("Gly", "G"),
        ("His", "H"), ("Ile", "I"), ("Leu", "L"), ("Lys", "K"), ("Met", "M"),
        ("Phe", "F"), ("Pro", "P"), ("Ser", "S"), ("Thr", "T"), ("Trp", "W"),
        ("Tyr", "Y"), ("Val", "V"), ("Xxx", "X"), ("Ter", "*"),
    ];

    /// Convert an HGVSp protein-change string to its short form, e.g.
    /// "p.Val600Glu" -> "p.V600E". Ported from vcf2maf's `%aa3to1` substitution.
    fn hgvsp_to_short(hgvsp: &str) -> String {
        let mut short = hgvsp.to_string();
        for (three, one) in Self::AA_3_TO_1 {
            short = short.replace(three, one);
        }
        short
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hgvsp_to_short_converts_three_letter_codes() {
        assert_eq!(MafRecord::hgvsp_to_short("p.Val600Glu"), "p.V600E");
        assert_eq!(MafRecord::hgvsp_to_short("p.Trp24Ter"), "p.W24*");
        assert_eq!(MafRecord::hgvsp_to_short("p.Gly12Asp"), "p.G12D");
    }

    #[test]
    fn test_hgvsp_to_short_passes_through_non_matching_input() {
        assert_eq!(MafRecord::hgvsp_to_short(""), "");
        assert_eq!(MafRecord::hgvsp_to_short("."), ".");
    }

    #[test]
    fn test_get_exon_number_prefers_vep_exon_field() {
        let mut info = HashMap::new();
        info.insert("CSQ_EXON".to_string(), "3/10".to_string());
        info.insert("ANN_Rank".to_string(), "4/12".to_string());
        assert_eq!(MafRecord::get_exon_number(&info), Some("3/10".to_string()));
    }

    #[test]
    fn test_get_exon_number_falls_back_to_snpeff_rank() {
        let mut info = HashMap::new();
        info.insert("ANN_Rank".to_string(), "4/12".to_string());
        assert_eq!(MafRecord::get_exon_number(&info), Some("4/12".to_string()));
    }

    #[test]
    fn test_get_exon_number_none_when_absent() {
        let info = HashMap::new();
        assert_eq!(MafRecord::get_exon_number(&info), None);
    }

    #[test]
    fn test_extract_ref_depth_from_info_ro() {
        let mut info = HashMap::new();
        info.insert("INFO_RO".to_string(), "42".to_string());
        assert_eq!(MafRecord::extract_ref_depth(&info), Some(42));
    }

    #[test]
    fn test_extract_ref_depth_none_when_absent() {
        let info = HashMap::new();
        assert_eq!(MafRecord::extract_ref_depth(&info), None);
    }

    #[test]
    fn test_extract_depth_from_sample_data_returns_ref_and_alt() {
        use crate::extract_sample_info::{ParsedFormatSample, ParsedSample};
        let mut format_fields = HashMap::new();
        format_fields.insert("DP".to_string(), "50".to_string());
        format_fields.insert("AD".to_string(), "30,20".to_string());
        let sample_data = Some(ParsedFormatSample {
            format_keys: vec!["DP".to_string(), "AD".to_string()],
            samples: vec![ParsedSample {
                sample_name: "SAMPLE-001".to_string(),
                format_fields,
            }],
        });
        let (total, refc, alt) = MafRecord::extract_depth_from_sample_data(&sample_data);
        assert_eq!(total, Some(50));
        assert_eq!(refc, Some(30));
        assert_eq!(alt, Some(20));
    }
}
