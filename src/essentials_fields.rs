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
    pub validation_status: Option<String>,
    pub mutation_status: String,
    pub sequence_source: String,
    pub sequencer: Option<String>,
    pub hgvsc: Option<String>,
    pub hgvsp: Option<String>,
    pub hgvsp_short: Option<String>,
    pub transcript_id: Option<String>,
    pub exon_number: Option<String>,
    pub t_depth: Option<u32>,      // was `total_depth` — INFO DP, unchanged extraction
    pub t_ref_count: Option<u32>,  // new — INFO RO / sample AD[0]
    pub t_alt_count: Option<u32>,  // was `depth` — INFO AO / sample AD[1], unchanged extraction
    pub n_depth: Option<u32>,      // matched normal, from --normal-id's FORMAT/DP
    pub n_ref_count: Option<u32>,  // matched normal, from --normal-id's AD[0]
    pub n_alt_count: Option<u32>,  // matched normal, from --normal-id's AD[1]
    // Trailing custom columns (this tool's own additions, not part of vcf2maf's core 46):
    pub filter_status: String,
    pub qual: Option<f64>,
    pub vaf: Option<f32>,
    pub protein_position: Option<String>,
}

/// vcf2maf's `%biotype_priority` (vcf2maf.pl:144-227), transcribed verbatim. Lower is
/// better. An unrecognized biotype gets 10, as vcf2maf.pl:228-232 does (it also warns;
/// we stay silent — a per-variant warning on a 92k-variant file is unusable).
const BIOTYPE_PRIORITY: &[(&str, u8)] = &[
    ("protein_coding", 1),
    ("LRG_gene", 2),
    ("IG_C_gene", 2),
    ("IG_D_gene", 2),
    ("IG_J_gene", 2),
    ("IG_LV_gene", 2),
    ("IG_V_gene", 2),
    ("TR_C_gene", 2),
    ("TR_D_gene", 2),
    ("TR_J_gene", 2),
    ("TR_V_gene", 2),
    ("miRNA", 3),
    ("snRNA", 3),
    ("snoRNA", 3),
    ("ribozyme", 3),
    ("tRNA", 3),
    ("sRNA", 3),
    ("scaRNA", 3),
    ("rRNA", 3),
    ("scRNA", 3),
    ("lincRNA", 3),
    ("lncRNA", 3),
    ("bidirectional_promoter_lncrna", 3),
    ("bidirectional_promoter_lncRNA", 3),
    ("known_ncrna", 4),
    ("vaultRNA", 4),
    ("vault_RNA", 4),
    ("macro_lncRNA", 4),
    ("Mt_tRNA", 4),
    ("Mt_rRNA", 4),
    ("antisense", 5),
    ("antisense_RNA", 5),
    ("sense_intronic", 5),
    ("sense_overlapping", 5),
    ("3prime_overlapping_ncrna", 5),
    ("3prime_overlapping_ncRNA", 5),
    ("misc_RNA", 5),
    ("non_coding", 5),
    ("regulatory_region", 6),
    ("disrupted_domain", 6),
    ("processed_transcript", 6),
    ("protein_coding_CDS_not_defined", 6),
    ("TEC", 6),
    ("TF_binding_site", 7),
    ("CTCF_binding_site", 7),
    ("promoter_flanking_region", 7),
    ("enhancer", 7),
    ("promoter", 7),
    ("open_chromatin_region", 7),
    ("retained_intron", 7),
    ("nonsense_mediated_decay", 7),
    ("non_stop_decay", 7),
    ("ambiguous_orf", 7),
    ("pseudogene", 8),
    ("processed_pseudogene", 8),
    ("polymorphic_pseudogene", 8),
    ("protein_coding_LoF", 8),
    ("retrotransposed", 8),
    ("translated_processed_pseudogene", 8),
    ("translated_unprocessed_pseudogene", 8),
    ("transcribed_processed_pseudogene", 8),
    ("transcribed_unprocessed_pseudogene", 8),
    ("transcribed_unitary_pseudogene", 8),
    ("unitary_pseudogene", 8),
    ("unprocessed_pseudogene", 8),
    ("Mt_tRNA_pseudogene", 8),
    ("tRNA_pseudogene", 8),
    ("snoRNA_pseudogene", 8),
    ("snRNA_pseudogene", 8),
    ("scRNA_pseudogene", 8),
    ("rRNA_pseudogene", 8),
    ("misc_RNA_pseudogene", 8),
    ("miRNA_pseudogene", 8),
    ("IG_pseudogene", 8),
    ("IG_C_pseudogene", 8),
    ("IG_D_pseudogene", 8),
    ("IG_J_pseudogene", 8),
    ("IG_V_pseudogene", 8),
    ("TR_J_pseudogene", 8),
    ("TR_V_pseudogene", 8),
    ("artifact", 9),
    ("", 10),
];

/// Transcript biotype rank. Lower is better; anything unlisted is 10.
pub fn biotype_priority(biotype: &str) -> u8 {
    BIOTYPE_PRIORITY
        .iter()
        .find(|(name, _)| *name == biotype)
        .map(|(_, priority)| *priority)
        .unwrap_or(10)
}

impl MafRecord {
    /// Convert from ReformattedVcfRecord to MafRecord
    // ponytail: unused by the binary since sample selection landed; the lib's tests are the callers.
    #[allow(dead_code)]
    pub fn from_reformatted_record(
        record: &ReformattedVcfRecord,
        center: &str,
        ncbi_build: &str,
        sample_barcode: &str,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        Self::from_reformatted_record_for_samples(record, center, ncbi_build, sample_barcode, None, None)
    }

    /// As `from_reformatted_record`, but reading the depth columns from explicitly named
    /// samples. `tumor` is the sample the t_* columns describe; `normal` populates the
    /// matched-normal columns. Passing `None` for either keeps the historical behaviour of
    /// using the first sample that declares `DP`.
    pub fn from_reformatted_record_for_samples(
        record: &ReformattedVcfRecord,
        center: &str,
        ncbi_build: &str,
        sample_barcode: &str,
        tumor: Option<&str>,
        normal: Option<&str>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        // Everything positional is decided on the trimmed alleles, exactly as vcf2maf does.
        let (pos, vcf_ref, vcf_alt) =
            Self::trim_shared_prefix(record.position, &record.reference, &record.alternate);
        let variant_type = Self::determine_variant_type(&vcf_ref, &vcf_alt);
        // Kept character-for-character as vcf2maf.pl:769 writes it:
        //   $inframe = ( abs( $ref_length - $var_length ) % 3 == 0 ? 1 : 0 );
        // clippy would rewrite this as .is_multiple_of(3); that is semantically identical but
        // breaks the line-for-line correspondence with the reference implementation.
        #[allow(clippy::manual_is_multiple_of)]
        let inframe = vcf_ref.len().abs_diff(vcf_alt.len()) % 3 == 0;
        let (start_pos, end_pos) = Self::calculate_maf_positions(pos, &vcf_ref);
        let (ref_allele, tumor_seq_allele1, tumor_seq_allele2) = Self::get_maf_alleles(
            &vcf_ref,
            &vcf_alt,
            Self::tumor_genotype(&record.format_sample_data, tumor).as_deref(),
        );

        // All three tumor counts must describe the SAME reads. They previously did not:
        // t_depth came from the sample's FORMAT/DP (one sample) while t_ref_count/t_alt_count
        // came from INFO RO/AO, which sum every sample in the VCF. On a tumor/normal file that
        // credited the normal's reads to the tumor and produced t_ref + t_alt > t_depth on
        // 98.4% of rows. Sample values now win for all three; INFO is only a fallback, and
        // only when no sample was named.
        let (sample_total, sample_ref, sample_alt) =
            Self::extract_depth_for_sample(&record.format_sample_data, tumor);

        // With a tumor named, its columns are the only source, as in vcf2maf. INFO is pooled
        // over every sample, so falling back to it credits the normal's reads to a tumor that
        // is absent from the VCF — or present with no call at all (".:.:.:."), which is how
        // 75 rows of the B487 tumor/normal file reported the normal's depth as the tumor's.
        let (t_depth, t_ref_count, t_alt_count) = if tumor.is_some() {
            (sample_total, sample_ref, sample_alt)
        } else {
            (
                // vcf2maf.pl:936 reads the sample's FORMAT/DP, not INFO/DP, which counts every
                // read the caller saw at the locus and runs higher.
                sample_total.or_else(|| Self::extract_total_depth(&record.info_fields)),
                sample_ref.or_else(|| Self::extract_ref_depth(&record.info_fields)),
                sample_alt.or_else(|| Self::extract_tumor_depth(&record.info_fields)),
            )
        };

        let (n_depth, n_ref_count, n_alt_count) = match normal {
            Some(name) => Self::extract_depth_for_sample(&record.format_sample_data, Some(name)),
            None => (None, None, None),
        };

        let vaf = match (t_alt_count, t_depth) {
            (Some(alt), Some(total)) if total > 0 => Some(alt as f32 / total as f32),
            _ => None,
        };

        let hgvsc = Self::get_annotation_field(&record.info_fields, &["CSQ_HGVSc", "ANN_HGVS_c"])
            .filter(|s| s != ".")
            .map(|s| Self::strip_accession(&s));
        let hgvsp = Self::get_annotation_field(&record.info_fields, &["CSQ_HGVSp", "ANN_HGVS_p"])
            .filter(|s| s != ".")
            .map(|s| Self::strip_accession(&s));
        // The splice rule runs after the 3->1 conversion in vcf2maf.pl and assigns
        // unconditionally, so it wins over a real HGVSp on the rare rows that have both.
        let hgvsp_short = Self::splice_hgvsp_short(&record.info_fields, hgvsc.as_deref())
            .or_else(|| hgvsp.as_deref().map(Self::hgvsp_to_short));
        let transcript_id = Self::get_transcript_id(&record.info_fields);

        Ok(MafRecord {
            hugo_symbol: Self::get_hugo_symbol(&record.info_fields, &transcript_id),
            entrez_gene_id: Self::get_entrez_gene_id(&record.info_fields),
            center: if center.trim().is_empty() {
                "Unknown_Center".to_string()
            } else {
                center.to_string()
            },
            ncbi_build: ncbi_build.to_string(),
            // The VCF's own naming is passed through unchanged (user's call
            // 2026-08-30, reversing the earlier strip-to-bare-name behaviour). vcf2maf
            // does the same, so this also removes a whole diff class against it.
            chromosome: record.chromosome.clone(),
            start_position: start_pos,
            end_position: end_pos,
            // vcf2maf.pl:907 — per the MAF definition, only "+" is an accepted value here.
            // VEP's transcript strand stays available in the TSV output's CSQ_STRAND column.
            strand: "+".to_string(),
            variant_classification: Self::get_variant_classification(
                &record.info_fields,
                &variant_type,
                inframe,
            ),
            variant_type,
            reference_allele: ref_allele,
            tumor_seq_allele1,
            tumor_seq_allele2,
            dbsnp_rs: Self::get_dbsnp_rs(&record.info_fields, record.id.as_deref()),
            dbsnp_val_status: None,
            tumor_sample_barcode: sample_barcode.to_string(),
            matched_norm_sample_barcode: normal.map(str::to_string),
            validation_status: None,
            // Neither is derivable from a VCF; main.rs overwrites them when the user passes
            // --mutation-status / --sequence-source.
            mutation_status: String::new(),
            sequence_source: String::new(),
            sequencer: Self::extract_sequencing_info(&record.info_fields),
            hgvsc,
            hgvsp,
            hgvsp_short,
            transcript_id: transcript_id.clone(),
            exon_number: Self::get_exon_number(&record.info_fields),
            t_depth,
            t_ref_count,
            t_alt_count,
            n_depth,
            n_ref_count,
            n_alt_count,
            filter_status: record.filter.clone(),
            qual: record.quality,
            vaf,
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
    /// vcf2maf.pl:854-860 — dbSNP_RS comes from VEP's `Existing_variation`, keeping only real
    /// rs IDs; a variant known only to COSMIC et al. leaves the column blank, and "novel" means
    /// VEP looked it up and found nothing. SnpEff's ANN carries no equivalent field, so on that
    /// path the VCF ID column is all there is.
    fn get_dbsnp_rs(info_fields: &HashMap<String, String>, id: Option<&str>) -> Option<String> {
        let Some(existing) = info_fields.get("CSQ_Existing_variation") else {
            return id.filter(|id| *id != ".").map(str::to_string);
        };
        if existing.is_empty() || existing == "." {
            return Some("novel".to_string());
        }
        // VEP joins co-located variants with "&"; vcf2maf.pl:784 rewrites those to "," first.
        let rs_ids: Vec<&str> = existing
            .split(['&', ','])
            .filter(|t| {
                t.strip_prefix("rs")
                    .is_some_and(|d| !d.is_empty() && d.bytes().all(|b| b.is_ascii_digit()))
            })
            .collect();
        (!rs_ids.is_empty()).then(|| rs_ids.join(","))
    }

    fn get_exon_number(info_fields: &HashMap<String, String>) -> Option<String> {
        Self::get_annotation_field(info_fields, &["CSQ_EXON", "ANN_Rank"])
    }

    /// vcf2maf never leaves Hugo_Symbol blank: upstream/downstream/intronic-lncRNA/IGR variants
    /// often carry a Gene ID with no HGNC symbol mapping, so VEP's SYMBOL field is empty. vcf2maf
    /// falls back to the transcript ID when one is associated, or the literal "Unknown" when
    /// there's no transcript either (pure IGR) — matched here for exact parity.
    fn get_hugo_symbol(
        info_fields: &HashMap<String, String>,
        transcript_id: &Option<String>,
    ) -> String {
        Self::get_annotation_field(info_fields, &["ANN_Gene_Name", "CSQ_SYMBOL", "ANN_SYMBOL"])
            .or_else(|| transcript_id.clone())
            .unwrap_or_else(|| "Unknown".to_string())
    }

    /// Handle multi-allelic variants by creating separate MafRecord for each alternate allele
    #[allow(dead_code)]
    pub fn from_reformatted_record_multi(
        record: &ReformattedVcfRecord,
        center: &str,
        ncbi_build: &str,
        sample_barcode: &str,
    ) -> Result<Vec<Self>, Box<dyn std::error::Error>> {
        Self::from_reformatted_record_multi_for_samples(
            record, center, ncbi_build, sample_barcode, None, None,
        )
    }

    /// As `from_reformatted_record_multi`, with the depth columns read from explicitly
    /// named tumor and normal samples.
    pub fn from_reformatted_record_multi_for_samples(
        record: &ReformattedVcfRecord,
        center: &str,
        ncbi_build: &str,
        sample_barcode: &str,
        tumor: Option<&str>,
        normal: Option<&str>,
    ) -> Result<Vec<Self>, Box<dyn std::error::Error>> {
        let alternates: Vec<&str> = record.alternate.split(',').collect();
        let mut maf_records = Vec::new();

        for (alt_index, alternate) in alternates.iter().enumerate() {
            // The annotation we kept belongs to one specific allele. On a multiallelic line,
            // carrying it onto the other ALTs would report a gene and consequence for an allele
            // the annotator never described, so those rows go out unannotated instead.
            let mut info_fields = record.info_fields.clone();
            if alternates.len() > 1
                && !Self::annotation_describes_allele(
                    &info_fields,
                    record.position,
                    &record.reference,
                    alternate,
                    alt_index,
                )
            {
                info_fields.retain(|k, _| !(k.starts_with("CSQ_") || k.starts_with("ANN_")));
            }

            // Create a record for each alternate allele
            let single_alt_record = ReformattedVcfRecord {
                chromosome: record.chromosome.clone(),
                position: record.position,
                id: record.id.clone(),
                reference: record.reference.clone(),
                alternate: alternate.to_string(),
                quality: record.quality,
                filter: record.filter.clone(),
                info_fields,
                format_sample_data: record.format_sample_data.clone(),
                annotation_field_type: record.annotation_field_type,
            };

            // Use the unified conversion logic
            let mut maf_record = Self::from_reformatted_record_for_samples(
                &single_alt_record,
                center,
                ncbi_build,
                sample_barcode,
                tumor,
                normal,
            )?;

            // Tumor_Seq_Allele1 needs the sibling ALTs a single-ALT record no longer carries:
            // a 1/2 genotype reports the *other* ALT (vcf2maf.pl:918-921).
            if alternates.len() > 1 {
                if let Some(genotype) = Self::tumor_genotype(&record.format_sample_data, tumor) {
                    if let Some(allele1) = Self::genotype_allele1(
                        &record.reference,
                        &alternates,
                        alt_index,
                        &genotype,
                    ) {
                        maf_record.tumor_seq_allele1 = allele1;
                    }
                }
            }

            // Per-allele alt count. With a tumor named it comes from that sample's AD only;
            // INFO AO pools every sample, so on a tumor/normal VCF it credits the normal's
            // reads to the tumor.
            let allele_alt_count = if tumor.is_some() {
                Self::sample_alt_depth_for_allele(&record.format_sample_data, tumor, alt_index)
            } else {
                Self::extract_tumor_depth_for_allele(&record.info_fields, alt_index)
            };

            if let Some(allele_depth) = allele_alt_count {
                maf_record.t_alt_count = Some(allele_depth);

                // Recalculate VAF with allele-specific depth
                if let Some(total) = maf_record.t_depth {
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
    /// The alt-observation count for one specific ALT, from a sample's `AD` — which carries
    /// one entry per allele (`ref,alt1,alt2,...`), so ALT number `n` is `AD[n + 1]`.
    /// This is the per-sample counterpart of `extract_tumor_depth_for_allele`, which reads
    /// INFO `AO` and therefore sums every sample in a multi-sample VCF.
    fn sample_alt_depth_for_allele(
        format_sample_data: &Option<ParsedFormatSample>,
        sample: Option<&str>,
        allele_index: usize,
    ) -> Option<u32> {
        format_sample_data
            .as_ref()?
            .samples
            .iter()
            .filter(|s| sample.is_none_or(|name| s.sample_name == name))
            .find_map(|s| s.format_fields.get("AD"))?
            .split(',')
            .nth(allele_index + 1)?
            .parse()
            .ok()
    }

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

    /// Depth, ref count and alt count from one sample. With `sample` set, only that sample is
    /// considered (and an absent name yields nothing); with `None`, the historical behaviour
    /// of taking the first sample that declares `DP` is kept.
    fn extract_depth_for_sample(
        format_sample_data: &Option<ParsedFormatSample>,
        sample: Option<&str>,
    ) -> (Option<u32>, Option<u32>, Option<u32>) {
        if let Some(sample_data) = format_sample_data {
            let mut total_depth = None;
            let mut ref_depth = None;
            let mut alt_depth = None;

            // Check each candidate sample for depth information
            let candidates = sample_data
                .samples
                .iter()
                .filter(|s| sample.is_none_or(|name| s.sample_name == name));
            for sample in candidates {
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

    /// Does the annotation we kept actually describe this ALT?
    ///
    /// VEP's `ALLELE_NUM` settles it outright when present — that is the key vcf2maf.pl:867
    /// uses — but VEP only emits it under `--allele_number`, which plenty of real files were
    /// not annotated with. Without it, fall back to the allele name in the first CSQ/ANN
    /// subfield, which VEP writes in any of three forms: the ALT verbatim, the ALT with its
    /// anchor base removed, or the fully prefix-trimmed allele ("-" when that leaves nothing).
    /// No allele named at all — unannotated input — counts as a match, so such records are
    /// left exactly as they were.
    fn annotation_describes_allele(
        info_fields: &HashMap<String, String>,
        position: u64,
        reference: &str,
        alternate: &str,
        alt_index: usize,
    ) -> bool {
        if let Some(allele_num) = Self::get_annotation_field(info_fields, &["CSQ_ALLELE_NUM"]) {
            return allele_num.parse::<usize>() == Ok(alt_index + 1);
        }
        let Some(annotated) = Self::get_annotation_field(info_fields, &["CSQ_Allele", "ANN_Allele"])
        else {
            return true;
        };
        let (_, _, trimmed) = Self::trim_shared_prefix(position, reference, alternate);
        let dash = |s: &str| if s.is_empty() { "-" } else { s }.to_string();
        [
            alternate.to_string(),
            dash(&trimmed),
            dash(&alternate[1.min(alternate.len())..]),
        ]
        .contains(&annotated)
    }

    /// Strip the bases REF and ALT share at the front, moving the position along with them.
    /// vcf2maf.pl:749-752 does this for *every* variant type, not just indels — so an
    /// un-normalized `CCCCA>CCCCC` is reported as the A>C SNP it actually is.
    /// Note: VCF alleles are ASCII (A/C/G/T), so char count == byte length.
    fn shared_prefix_len(reference: &str, alternate: &str) -> usize {
        // vcf2maf's loop is guarded on `$ref ne $var`, so identical alleles are left alone.
        if reference == alternate {
            return 0;
        }
        reference
            .chars()
            .zip(alternate.chars())
            .take_while(|(r, a)| r == a)
            .count()
    }

    fn trim_shared_prefix(position: u64, reference: &str, alternate: &str) -> (u64, String, String) {
        let shared = Self::shared_prefix_len(reference, alternate);
        (
            position + shared as u64,
            reference[shared..].to_string(),
            alternate[shared..].to_string(),
        )
    }

    /// Takes the *trimmed* REF. An empty REF is an insertion, which MAF anchors between the
    /// two flanking bases; everything else spans the reference bases it replaces or deletes.
    fn calculate_maf_positions(position: u64, reference: &str) -> (u64, u64) {
        if reference.is_empty() {
            (position.saturating_sub(1), position)
        } else {
            (position, position + reference.len() as u64 - 1)
        }
    }

    /// Takes the *trimmed* alleles. MAF writes a bare "-" where a trimmed allele is empty.
    ///
    /// vcf2maf.pl:913-921 picks Tumor_Seq_Allele1 as the first genotype allele that isn't the
    /// variant, so a homozygous-alt call reports the ALT twice; with no usable GT it assumes a
    /// ref/var heterozygote. A `1/2` call needs the sibling ALTs, which only
    /// `from_reformatted_record_multi` still has — it patches the result afterwards.
    fn get_maf_alleles(
        reference: &str,
        alternate: &str,
        tumor_genotype: Option<&str>,
    ) -> (String, String, String) {
        let dash = |s: &str| if s.is_empty() { "-" } else { s }.to_string();
        let hom_alt = tumor_genotype.is_some_and(|gt| {
            let mut indices = gt.split(['/', '|']).peekable();
            indices.peek().is_some()
                && indices.all(|i| matches!(i.parse::<u32>(), Ok(index) if index > 0))
        });
        let allele1 = if hom_alt { alternate } else { reference };
        (dash(reference), dash(allele1), dash(alternate))
    }

    /// Tumor_Seq_Allele1 for one ALT of a multiallelic line, ported from vcf2maf.pl:918-921:
    /// the first GT allele that isn't this row's variant, so `1/2` reports the *sibling* ALT.
    /// Alleles are trimmed by the prefix this row's REF/ALT share (vcf2maf.pl:749-752 trims the
    /// whole allele list together), and a sibling shorter than that trim collapses to "-", the
    /// same as Perl's `substr` past the end of the string.
    fn genotype_allele1(
        reference: &str,
        alternates: &[&str],
        alt_index: usize,
        genotype: &str,
    ) -> Option<String> {
        let variant = alternates.get(alt_index)?;
        let shared = Self::shared_prefix_len(reference, variant);
        let trim = |allele: &str| match allele.get(shared..) {
            Some(trimmed) if !trimmed.is_empty() => trimmed.to_string(),
            _ => "-".to_string(),
        };
        let allele_at = |index: usize| match index.checked_sub(1) {
            None => Some(trim(reference)),
            Some(alt) => alternates.get(alt).map(|a| trim(a)),
        };

        let mut indices = genotype.split(['/', '|']).map(|i| i.parse::<usize>().ok());
        let first = indices.next().flatten()?;
        // "If GT was monoploid, then $idx2 will be undefined, and we should set it equal to $idx1"
        let second = indices.next().flatten().unwrap_or(first);

        let allele1 = allele_at(first)?;
        if allele1 != trim(variant) {
            Some(allele1)
        } else {
            allele_at(second)
        }
    }

    /// The tumor sample's GT, taken from the first sample that declares one — the same
    /// "first sample is the tumor" assumption `extract_depth_from_sample_data` makes.
    fn tumor_genotype(
        format_sample_data: &Option<ParsedFormatSample>,
        sample: Option<&str>,
    ) -> Option<String> {
        format_sample_data
            .as_ref()?
            .samples
            .iter()
            .filter(|s| sample.is_none_or(|name| s.sample_name == name))
            .find_map(|s| s.format_fields.get("GT").cloned())
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

    pub fn effect_priority(term: &str) -> u8 {
        Self::EFFECT_PRIORITY
            .iter()
            .find(|(name, _)| *name == term)
            .map(|(_, priority)| *priority)
            .unwrap_or(20)
    }

    /// Resolve a (possibly `&`/`,`/`|`-joined) multi-term consequence string down to the
    /// single most severe term, mirroring vcf2maf's sort-by-priority-then-take-first.
    pub fn resolve_one_consequence(consequence: &str) -> String {
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

    /// Drop the reference-sequence accession from an HGVS string, as vcf2maf.pl:786-787 does
    /// with `s/^.*://`. VEP writes full HGVS ("ENST00000641515.2:c.760T>A"); the MAF column
    /// carries only the change, the accession being already present in Transcript_ID. SnpEff
    /// writes the change bare, with no colon, so it passes through untouched.
    fn strip_accession(hgvs: &str) -> String {
        match hgvs.rfind(':') {
            Some(colon) => hgvs[colon + 1..].to_string(),
            None => hgvs.to_string(),
        }
    }

    /// The synthetic `p.X{codon}_splice` protein change vcf2maf.pl:824-834 builds for splice
    /// acceptor/donor variants. They sit in an intron, so VEP reports no protein change at all;
    /// vcf2maf derives a codon number from the cDNA position so the MAF still carries a protein
    /// coordinate. Returns None for every other consequence, and for an HGVSc whose position is
    /// not a plain number (a 5' UTR `c.-14+1G>T` or a 3' UTR `c.*91G>T` never matches
    /// vcf2maf's `m/^c.(\d+)/`, so it leaves HGVSp_Short alone).
    fn splice_hgvsp_short(
        info_fields: &HashMap<String, String>,
        hgvsc: Option<&str>,
    ) -> Option<String> {
        let consequence =
            Self::get_annotation_field(info_fields, &["CSQ_Consequence", "ANN_Annotation"])?;
        // vcf2maf gates on One_Consequence — the most severe term, not the first one listed.
        let term = Self::resolve_one_consequence(&consequence.to_lowercase());
        if !matches!(
            term.as_str(),
            "splice_acceptor_variant" | "splice_donor_variant"
        ) {
            return None;
        }

        let digits: String = hgvsc?
            .strip_prefix("c.")?
            .chars()
            .take_while(char::is_ascii_digit)
            .collect();
        // vcf2maf.pl:828 guards against cDNA positions below 1 before dividing.
        let c_pos = digits.parse::<u64>().ok()?.max(1);
        // vcf2maf.pl:829 — sprintf( "%.0f", ( $c_pos + $c_pos % 3 ) / 3 ). Perl divides in
        // floating point and rounds; integer division would truncate instead. The quotient is
        // always whole, x.333 or x.667, never a .5 tie, so the rounding mode does not matter.
        let p_pos = ((c_pos + c_pos % 3) as f64 / 3.0).round() as u64;
        Some(format!("p.X{p_pos}_splice"))
    }

    fn classify_by_impact(info_fields: &HashMap<String, String>) -> String {
        let impact = Self::get_annotation_field(info_fields, &["CSQ_IMPACT", "ANN_Annotation_Impact"]);
        match impact.as_deref().map(|s| s.to_uppercase()).as_deref() {
            Some("HIGH") | Some("MODERATE") => "Missense_Mutation".to_string(),
            Some("LOW") | Some("MODIFIER") => "Silent".to_string(),
            // vcf2maf.pl:1050 — "Targeted_Region" when there is no effect to go on at all.
            _ => "Targeted_Region".to_string(),
        }
    }

    pub fn get_maf_headers() -> Vec<String> {
        [
            "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
            "Start_Position", "End_Position", "Strand", "Variant_Classification",
            "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode",
            "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
            "Tumor_Validation_Allele1", "Tumor_Validation_Allele2",
            "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2",
            "Verification_Status", "Validation_Status", "Mutation_Status",
            "Sequencing_Phase", "Sequence_Source", "Validation_Method", "Score",
            "BAM_File", "Sequencer", "Tumor_Sample_UUID", "Matched_Norm_Sample_UUID",
            "HGVSc", "HGVSp", "HGVSp_Short", "Transcript_ID", "Exon_Number",
            "t_depth", "t_ref_count", "t_alt_count", "n_depth", "n_ref_count",
            "n_alt_count", "all_effects",
            "FILTER", "QUAL", "VAF", "Protein_Position",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect()
    }

    pub fn to_tsv_line(&self) -> String {
        // vcf2maf leaves a cell it has no data for empty. A "." there was our own invention
        // and the single largest diff class against it. Cells whose "." came from the VCF
        // itself (FILTER, ID) are untouched — that dot is content, not absence.
        let empty = "";
        let entrez = self.entrez_gene_id.map(|id| id.to_string());
        let start = self.start_position.to_string();
        let end = self.end_position.to_string();
        let t_depth = self.t_depth.map(|d| d.to_string());
        let t_ref_count = self.t_ref_count.map(|d| d.to_string());
        let t_alt_count = self.t_alt_count.map(|d| d.to_string());
        let n_depth = self.n_depth.map(|d| d.to_string());
        let n_ref_count = self.n_ref_count.map(|d| d.to_string());
        let n_alt_count = self.n_alt_count.map(|d| d.to_string());
        let vaf = self.vaf.map(|v| format!("{:.4}", v));
        let qual = self.qual.map(|q| q.to_string());

        [
            self.hugo_symbol.as_str(),
            entrez.as_deref().unwrap_or(empty),
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
            self.dbsnp_rs.as_deref().unwrap_or(empty),
            self.dbsnp_val_status.as_deref().unwrap_or(empty),
            self.tumor_sample_barcode.as_str(),
            self.matched_norm_sample_barcode.as_deref().unwrap_or(empty),
            empty, // 18 Match_Norm_Seq_Allele1 — the normal's GT is not read, only its depths
            empty, // 19 Match_Norm_Seq_Allele2
            empty, // 20 Tumor_Validation_Allele1
            empty, // 21 Tumor_Validation_Allele2
            empty, // 22 Match_Norm_Validation_Allele1
            empty, // 23 Match_Norm_Validation_Allele2
            empty, // 24 Verification_Status
            self.validation_status.as_deref().unwrap_or(empty),
            self.mutation_status.as_str(),
            empty, // 27 Sequencing_Phase
            self.sequence_source.as_str(),
            empty, // 29 Validation_Method
            empty, // 30 Score
            empty, // 31 BAM_File
            self.sequencer.as_deref().unwrap_or(empty),
            empty, // 33 Tumor_Sample_UUID
            empty, // 34 Matched_Norm_Sample_UUID
            self.hgvsc.as_deref().unwrap_or(empty),
            self.hgvsp.as_deref().unwrap_or(empty),
            self.hgvsp_short.as_deref().unwrap_or(empty),
            self.transcript_id.as_deref().unwrap_or(empty),
            self.exon_number.as_deref().unwrap_or(empty),
            t_depth.as_deref().unwrap_or(empty),
            t_ref_count.as_deref().unwrap_or(empty),
            t_alt_count.as_deref().unwrap_or(empty),
            n_depth.as_deref().unwrap_or(empty),
            n_ref_count.as_deref().unwrap_or(empty),
            n_alt_count.as_deref().unwrap_or(empty),
            empty, // 46 all_effects — see Global Constraints: deferred, needs full transcript list
            self.filter_status.as_str(),
            qual.as_deref().unwrap_or(empty),
            vaf.as_deref().unwrap_or(empty),
            self.protein_position.as_deref().unwrap_or(empty),
        ]
        .join("\t")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn maf_from(position: u64, reference: &str, alternate: &str) -> MafRecord {
        let record = ReformattedVcfRecord {
            chromosome: "chr1".to_string(),
            position,
            id: None,
            reference: reference.to_string(),
            alternate: alternate.to_string(),
            quality: Some(60.0),
            filter: "PASS".to_string(),
            info_fields: HashMap::new(),
            format_sample_data: None,
            annotation_field_type: crate::reformat_vcf::AnnotationFieldType::None,
        };
        MafRecord::from_reformatted_record(&record, "test", "GRCh38", "sample").unwrap()
    }

    #[test]
    fn test_dnp_end_position_spans_both_bases() {
        // vcf2maf.pl:755 — ( $start, $stop ) = ( $pos, $pos + $var_length - 1 )
        let maf = maf_from(100, "AC", "GT");
        assert_eq!(maf.variant_type, "DNP");
        assert_eq!((maf.start_position, maf.end_position), (100, 101));
    }

    #[test]
    fn test_onp_end_position_spans_all_bases() {
        let maf = maf_from(100, "ACGT", "TGCA");
        assert_eq!(maf.variant_type, "ONP");
        assert_eq!((maf.start_position, maf.end_position), (100, 103));
    }

    #[test]
    fn test_untrimmed_equal_length_alleles_reduce_to_snp() {
        // Real un-normalized freebayes site: chr1:8324505 CCCCA>CCCCC is an A>C SNP at 8324509.
        // vcf2maf.pl:749-752 strips shared leading bases for every variant type, not just indels.
        let maf = maf_from(8324505, "CCCCA", "CCCCC");
        assert_eq!(maf.variant_type, "SNP");
        assert_eq!((maf.start_position, maf.end_position), (8324509, 8324509));
        assert_eq!(maf.reference_allele, "A");
        assert_eq!(maf.tumor_seq_allele2, "C");
    }

    #[test]
    fn test_insertion_positions_and_alleles_unchanged() {
        let maf = maf_from(100, "A", "ATCG");
        assert_eq!(maf.variant_type, "INS");
        assert_eq!((maf.start_position, maf.end_position), (100, 101));
        assert_eq!(maf.reference_allele, "-");
        assert_eq!(maf.tumor_seq_allele2, "TCG");
    }

    #[test]
    fn test_deletion_positions_and_alleles_unchanged() {
        let maf = maf_from(100, "ATCG", "A");
        assert_eq!(maf.variant_type, "DEL");
        assert_eq!((maf.start_position, maf.end_position), (101, 103));
        assert_eq!(maf.reference_allele, "TCG");
        assert_eq!(maf.tumor_seq_allele2, "-");
    }

    fn annotated_record(
        position: u64,
        reference: &str,
        alternate: &str,
        csq_allele: &str,
    ) -> ReformattedVcfRecord {
        let mut info_fields = HashMap::new();
        info_fields.insert("CSQ_Allele".to_string(), csq_allele.to_string());
        info_fields.insert("CSQ_SYMBOL".to_string(), "SLC45A1".to_string());
        info_fields.insert(
            "CSQ_Consequence".to_string(),
            "missense_variant".to_string(),
        );
        ReformattedVcfRecord {
            chromosome: "chr1".to_string(),
            position,
            id: None,
            reference: reference.to_string(),
            alternate: alternate.to_string(),
            quality: Some(60.0),
            filter: "PASS".to_string(),
            info_fields,
            format_sample_data: None,
            annotation_field_type: crate::reformat_vcf::AnnotationFieldType::Csq,
        }
    }

    #[test]
    fn test_strand_is_always_plus() {
        // vcf2maf.pl:907 — "Per MAF definition, only the positive strand is an accepted value".
        // The MAF Strand column is genomic; VEP's transcript strand does not belong in it.
        let mut record = annotated_record(100, "A", "G", "G");
        record
            .info_fields
            .insert("CSQ_STRAND".to_string(), "-1".to_string());
        let maf = MafRecord::from_reformatted_record(&record, "c", "GRCh38", "s").unwrap();

        assert_eq!(maf.strand, "+");
    }

    #[test]
    fn test_unasserted_metadata_columns_default_to_empty() {
        // Nothing in a VCF says the calls are somatic or that the library was an exome.
        let maf = maf_from(100, "A", "G");
        assert_eq!(maf.mutation_status, "");
        assert_eq!(maf.sequence_source, "");
    }

    #[test]
    fn test_unannotated_record_uses_vcf2maf_fallback_classification() {
        // vcf2maf.pl:1050 — return "Targeted_Region" if( not defined $effect or not $effect );
        // "Unknown" is not a value the MAF spec allows in this column.
        let maf = maf_from(100, "A", "G");
        assert_eq!(maf.variant_classification, "Targeted_Region");
    }

    #[test]
    fn test_multiallelic_annotation_stays_on_its_own_allele() {
        // VEP annotated GCCCC only; CCCCC must not inherit its gene and consequence.
        let record = annotated_record(8324505, "CCCCA", "GCCCC,CCCCC", "GCCCC");
        let rows =
            MafRecord::from_reformatted_record_multi(&record, "c", "GRCh38", "s").unwrap();

        assert_eq!(rows.len(), 2, "one row per ALT");
        assert_eq!(rows[0].hugo_symbol, "SLC45A1");
        assert_eq!(rows[0].variant_classification, "Missense_Mutation");
        assert_eq!(rows[1].hugo_symbol, "Unknown");
        assert_ne!(rows[1].variant_classification, "Missense_Mutation");
    }

    #[test]
    fn test_multiallelic_annotation_kept_when_it_names_the_second_allele() {
        let record = annotated_record(8324505, "CCCCA", "GCCCC,CCCCC", "CCCCC");
        let rows =
            MafRecord::from_reformatted_record_multi(&record, "c", "GRCh38", "s").unwrap();

        assert_eq!(rows[0].hugo_symbol, "Unknown");
        assert_eq!(rows[1].hugo_symbol, "SLC45A1");
    }

    #[test]
    fn test_multiallelic_matches_vep_minimal_indel_allele() {
        // VEP reports deletions as "-": REF=AT ALT=A is a deletion of T.
        let record = annotated_record(100, "AT", "A,ATT", "-");
        let rows =
            MafRecord::from_reformatted_record_multi(&record, "c", "GRCh38", "s").unwrap();

        assert_eq!(rows[0].hugo_symbol, "SLC45A1", "deletion allele is the annotated one");
        assert_eq!(rows[1].hugo_symbol, "Unknown");
    }

    #[test]
    fn test_multiallelic_matches_vep_anchor_stripped_insertion_allele() {
        // Real site: REF=TGGAGGA ALT=T,TGGAGGAGGA — VEP names the insertion allele
        // "GGAGGAGGA", i.e. the ALT with only its anchor base removed.
        let record = annotated_record(73385903, "TGGAGGA", "T,TGGAGGAGGA", "GGAGGAGGA");
        let rows =
            MafRecord::from_reformatted_record_multi(&record, "c", "GRCh38", "s").unwrap();

        assert_eq!(rows[0].hugo_symbol, "Unknown", "deletion allele was not annotated");
        assert_eq!(rows[1].hugo_symbol, "SLC45A1");
    }

    #[test]
    fn test_allele_num_decides_when_vep_provides_it() {
        // vcf2maf.pl:867 skips effects whose ALLELE_NUM is not this ALT's 1-based index.
        // Present only when VEP ran with --allele_number, and authoritative when it is.
        let mut record = annotated_record(100, "A", "G,T", "does_not_match");
        record
            .info_fields
            .insert("CSQ_ALLELE_NUM".to_string(), "2".to_string());
        let rows =
            MafRecord::from_reformatted_record_multi(&record, "c", "GRCh38", "s").unwrap();

        assert_eq!(rows[0].hugo_symbol, "Unknown");
        assert_eq!(rows[1].hugo_symbol, "SLC45A1");
    }

    #[test]
    fn test_single_allele_annotation_is_never_stripped() {
        // Guard: allele filtering must not touch the ordinary one-ALT case.
        let record = annotated_record(100, "A", "G", "does_not_match");
        let rows =
            MafRecord::from_reformatted_record_multi(&record, "c", "GRCh38", "s").unwrap();

        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].hugo_symbol, "SLC45A1");
    }

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

    /// Build a MafRecord from just the annotation fields the HGVS columns are derived from.
    fn maf_with_hgvs(consequence: &str, hgvsc: Option<&str>, hgvsp: Option<&str>) -> MafRecord {
        let mut info_fields = HashMap::new();
        info_fields.insert("CSQ_Consequence".to_string(), consequence.to_string());
        if let Some(c) = hgvsc {
            info_fields.insert("CSQ_HGVSc".to_string(), c.to_string());
        }
        if let Some(p) = hgvsp {
            info_fields.insert("CSQ_HGVSp".to_string(), p.to_string());
        }
        let record =
            create_test_maf_record("chr1", 100, "A", "G", Some(60.0), "PASS", info_fields);
        MafRecord::from_reformatted_record(&record, "test", "GRCh38", "sample").unwrap()
    }

    #[test]
    fn test_hgvs_strips_the_reference_sequence_accession() {
        // vcf2maf.pl:786-787 — s/^.*:// on both. VEP writes full HGVS ("ENST...:c.760T>A");
        // the MAF column carries only the change, since Transcript_ID sits right beside it.
        let maf = maf_with_hgvs(
            "missense_variant",
            Some("ENST00000641515.2:c.760T>A"),
            Some("ENSP00000493376.2:p.Leu254Met"),
        );
        assert_eq!(maf.hgvsc.as_deref(), Some("c.760T>A"));
        assert_eq!(maf.hgvsp.as_deref(), Some("p.Leu254Met"));
    }

    #[test]
    fn test_hgvsp_short_is_built_from_the_stripped_hgvsp() {
        // The 3->1 conversion must run on the change alone, or the accession is carried along.
        let maf = maf_with_hgvs(
            "missense_variant",
            Some("ENST00000288602.11:c.1799T>A"),
            Some("ENSP00000288602.6:p.Val600Glu"),
        );
        assert_eq!(maf.hgvsp_short.as_deref(), Some("p.V600E"));
    }

    #[test]
    fn test_hgvs_without_an_accession_is_left_alone() {
        // SnpEff's ANN writes the change bare, with no accession and no colon.
        let maf = maf_with_hgvs("synonymous_variant", Some("c.*91G>T"), Some("p.Pro34Pro"));
        assert_eq!(maf.hgvsc.as_deref(), Some("c.*91G>T"));
        assert_eq!(maf.hgvsp.as_deref(), Some("p.Pro34Pro"));
        assert_eq!(maf.hgvsp_short.as_deref(), Some("p.P34P"));
    }

    #[test]
    fn test_splice_site_gets_a_synthetic_hgvsp_short() {
        // vcf2maf.pl:824-834 — splice variants are intronic, so VEP reports no protein change.
        // vcf2maf synthesizes one from the cDNA position: p.X{codon}_splice.
        // c.756+1G>T -> c_pos 756, 756 % 3 == 0 -> 756 / 3 = 252.
        let maf = maf_with_hgvs(
            "splice_donor_variant",
            Some("ENST00000380152.8:c.756+1G>T"),
            None,
        );
        assert_eq!(maf.hgvsp_short.as_deref(), Some("p.X252_splice"));
        assert_eq!(maf.variant_classification, "Splice_Site");
    }

    #[test]
    fn test_synthetic_splice_position_rounds_the_codon_up() {
        // ( c_pos + c_pos % 3 ) / 3 rounded: 757 -> 758/3 = 252.67 -> 253; 758 -> 760/3 -> 253.
        let a = maf_with_hgvs("splice_acceptor_variant", Some("c.757-2A>G"), None);
        assert_eq!(a.hgvsp_short.as_deref(), Some("p.X253_splice"));
        let b = maf_with_hgvs("splice_acceptor_variant", Some("c.758-1A>G"), None);
        assert_eq!(b.hgvsp_short.as_deref(), Some("p.X253_splice"));
    }

    #[test]
    fn test_splice_rule_keys_on_the_most_severe_consequence() {
        // vcf2maf gates on One_Consequence, i.e. after sorting by severity — not on whichever
        // term VEP happened to list first. splice_donor (2) outranks intron_variant (14).
        let maf = maf_with_hgvs("intron_variant&splice_donor_variant", Some("c.300+1G>A"), None);
        assert_eq!(maf.hgvsp_short.as_deref(), Some("p.X100_splice"));
    }

    #[test]
    fn test_splice_rule_needs_a_numeric_cdna_position() {
        // vcf2maf.pl:826 matches /^c.(\d+)/, which a 5' UTR position like c.-14+1 never
        // satisfies; the rewrite is skipped and HGVSp_Short stays empty.
        let maf = maf_with_hgvs("splice_donor_variant", Some("c.-14+1G>T"), None);
        assert_eq!(maf.hgvsp_short, None);
    }

    #[test]
    fn test_non_splice_consequences_keep_their_real_hgvsp_short() {
        let maf = maf_with_hgvs(
            "missense_variant",
            Some("ENST00000288602.11:c.1799T>A"),
            Some("ENSP00000288602.6:p.Val600Glu"),
        );
        assert_eq!(maf.hgvsp_short.as_deref(), Some("p.V600E"));
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
    fn test_get_hugo_symbol_prefers_annotated_symbol() {
        let mut info = HashMap::new();
        info.insert("CSQ_SYMBOL".to_string(), "BRCA1".to_string());
        assert_eq!(
            MafRecord::get_hugo_symbol(&info, &Some("ENST00000123456".to_string())),
            "BRCA1"
        );
    }

    #[test]
    fn test_get_hugo_symbol_falls_back_to_transcript_id_when_symbol_missing() {
        let info = HashMap::new();
        assert_eq!(
            MafRecord::get_hugo_symbol(&info, &Some("ENST00000620188".to_string())),
            "ENST00000620188"
        );
    }

    #[test]
    fn test_get_hugo_symbol_falls_back_to_unknown_when_no_symbol_or_transcript() {
        let info = HashMap::new();
        assert_eq!(MafRecord::get_hugo_symbol(&info, &None), "Unknown");
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

    fn record_with_genotype(reference: &str, alternate: &str, gt: &str) -> ReformattedVcfRecord {
        use crate::extract_sample_info::{ParsedFormatSample, ParsedSample};
        let mut format_fields = HashMap::new();
        format_fields.insert("GT".to_string(), gt.to_string());
        let mut record = annotated_record(100, reference, alternate, alternate);
        record.format_sample_data = Some(ParsedFormatSample {
            format_keys: vec!["GT".to_string()],
            samples: vec![ParsedSample {
                sample_name: "TUMOR".to_string(),
                format_fields,
            }],
        });
        record
    }

    #[test]
    fn test_tumor_seq_allele1_is_the_alt_when_genotype_is_hom_alt() {
        // vcf2maf.pl:913-921 — Tumor_Seq_Allele1 is the first GT allele that isn't the variant,
        // so a 1/1 call reports the ALT twice rather than pretending the site is heterozygous.
        let record = record_with_genotype("T", "C", "1/1");
        let maf = MafRecord::from_reformatted_record(&record, "c", "GRCh38", "s").unwrap();

        assert_eq!(maf.reference_allele, "T");
        assert_eq!(maf.tumor_seq_allele1, "C");
        assert_eq!(maf.tumor_seq_allele2, "C");
    }

    #[test]
    fn test_tumor_seq_allele1_hom_alt_deletion_uses_the_dash_form() {
        let record = record_with_genotype("ATCG", "A", "1|1");
        let maf = MafRecord::from_reformatted_record(&record, "c", "GRCh38", "s").unwrap();

        assert_eq!(maf.reference_allele, "TCG");
        assert_eq!(maf.tumor_seq_allele1, "-");
    }

    #[test]
    fn test_tumor_seq_allele1_stays_reference_for_het_and_missing_genotypes() {
        // Guard: only hom-alt changes. vcf2maf assumes ref/var het when GT is absent or "./.".
        for gt in ["0/1", "0|1", "./.", "."] {
            let record = record_with_genotype("T", "C", gt);
            let maf = MafRecord::from_reformatted_record(&record, "c", "GRCh38", "s").unwrap();
            assert_eq!(maf.tumor_seq_allele1, "T", "GT was {gt}");
        }
        let maf = maf_from(100, "T", "C");
        assert_eq!(maf.tumor_seq_allele1, "T", "no sample columns at all");
    }

    #[test]
    fn test_multiallelic_genotype_reports_the_sibling_allele() {
        // Real site: chr1:240207640 REF=CT ALT=TC,CC GT=1/2. vcf2maf.pl:921 takes the first GT
        // allele that isn't this row's variant, so the TC row reports CC and vice versa.
        let record = record_with_genotype("CT", "TC,CC", "1/2");
        let rows =
            MafRecord::from_reformatted_record_multi(&record, "c", "GRCh38", "s").unwrap();

        assert_eq!(rows[0].tumor_seq_allele2, "TC");
        assert_eq!(rows[0].tumor_seq_allele1, "CC");
        // The CC row trims the C it shares with REF=CT, and its siblings trim with it.
        assert_eq!(rows[1].tumor_seq_allele2, "C");
        assert_eq!(rows[1].tumor_seq_allele1, "C");
    }

    #[test]
    fn test_multiallelic_genotype_with_reference_allele_reports_reference() {
        let record = record_with_genotype("CT", "TC,CC", "0/2");
        let rows =
            MafRecord::from_reformatted_record_multi(&record, "c", "GRCh38", "s").unwrap();

        assert_eq!(rows[1].tumor_seq_allele2, "C");
        assert_eq!(rows[1].tumor_seq_allele1, "T", "GT names the reference allele");
    }

    #[test]
    fn test_multiallelic_sibling_shorter_than_the_trim_becomes_dash() {
        // REF=TGGAGGA ALT=T,TGGAGGAGGA GT=1/2: for the insertion row the trim eats 7 bases, more
        // than the sibling deletion allele has, and vcf2maf's substr loop leaves it as "-".
        let record = record_with_genotype("TGGAGGA", "T,TGGAGGAGGA", "1/2");
        let rows =
            MafRecord::from_reformatted_record_multi(&record, "c", "GRCh38", "s").unwrap();

        assert_eq!(rows[1].tumor_seq_allele2, "GGA", "insertion row");
        assert_eq!(rows[1].tumor_seq_allele1, "-");
    }

    #[test]
    fn test_t_depth_prefers_the_sample_over_info_dp() {
        // vcf2maf.pl:936 takes t_depth from the tumor sample's FORMAT/DP. INFO/DP counts reads
        // the caller saw at the locus, which on mutect2 output is consistently higher.
        use crate::extract_sample_info::{ParsedFormatSample, ParsedSample};
        let mut format_fields = HashMap::new();
        format_fields.insert("DP".to_string(), "90".to_string());
        format_fields.insert("AD".to_string(), "60,30".to_string());

        let mut record = annotated_record(100, "A", "G", "G");
        record
            .info_fields
            .insert("INFO_DP".to_string(), "100".to_string());
        record.format_sample_data = Some(ParsedFormatSample {
            format_keys: vec!["DP".to_string(), "AD".to_string()],
            samples: vec![ParsedSample {
                sample_name: "TUMOR".to_string(),
                format_fields,
            }],
        });

        let maf = MafRecord::from_reformatted_record(&record, "c", "GRCh38", "s").unwrap();
        assert_eq!(maf.t_depth, Some(90));
        assert_eq!(maf.t_ref_count, Some(60));
        assert_eq!(maf.t_alt_count, Some(30));
    }

    #[test]
    fn test_t_depth_falls_back_to_info_dp_without_sample_columns() {
        let mut record = annotated_record(100, "A", "G", "G");
        record
            .info_fields
            .insert("INFO_DP".to_string(), "100".to_string());

        let maf = MafRecord::from_reformatted_record(&record, "c", "GRCh38", "s").unwrap();
        assert_eq!(maf.t_depth, Some(100));
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
        let (total, refc, alt) = MafRecord::extract_depth_for_sample(&sample_data, None);
        assert_eq!(total, Some(50));
        assert_eq!(refc, Some(30));
        assert_eq!(alt, Some(20));
    }

    /// Two samples with the pooled INFO counts a real tumor/normal freebayes VCF carries.
    /// Modelled on the first variant of B487_1_V_vs_B487_1_cOM: INFO DP=7 RO=5 AO=2 is the
    /// SUM over both samples, while the tumor itself has DP=1, AD=1,0.
    fn tumor_normal_record() -> ReformattedVcfRecord {
        use crate::extract_sample_info::{ParsedFormatSample, ParsedSample};
        let sample = |name: &str, dp: &str, ad: &str| ParsedSample {
            sample_name: name.to_string(),
            format_fields: HashMap::from([
                ("DP".to_string(), dp.to_string()),
                ("AD".to_string(), ad.to_string()),
            ]),
        };
        let mut record = create_test_maf_record(
            "chr1", 69787, "T", "A", Some(50.0), "PASS",
            HashMap::from([
                ("INFO_DP".to_string(), "7".to_string()),
                ("INFO_RO".to_string(), "5".to_string()),
                ("INFO_AO".to_string(), "2".to_string()),
            ]),
        );
        record.format_sample_data = Some(ParsedFormatSample {
            format_keys: vec!["DP".to_string(), "AD".to_string()],
            samples: vec![
                sample("B487_1_V", "1", "1,0"),
                sample("B487_1_cOM", "6", "4,2"),
            ],
        });
        record
    }

    #[test]
    fn per_allele_alt_count_comes_from_the_tumor_sample_ad_not_pooled_info_ao() {
        // The multi path overrode t_alt_count from INFO/AO after conversion, which on a
        // tumor/normal VCF pools every sample. The tumor's own AD carries one entry per
        // allele (ref,alt1,alt2), and that is what each row must report.
        use crate::extract_sample_info::{ParsedFormatSample, ParsedSample};
        let mut record = create_test_maf_record(
            "chr1", 1000, "A", "G,T", Some(50.0), "PASS",
            HashMap::from([("INFO_AO".to_string(), "30,40".to_string())]),
        );
        record.format_sample_data = Some(ParsedFormatSample {
            format_keys: vec!["DP".to_string(), "AD".to_string()],
            samples: vec![
                ParsedSample {
                    sample_name: "TUMOR".to_string(),
                    format_fields: HashMap::from([
                        ("DP".to_string(), "20".to_string()),
                        ("AD".to_string(), "10,3,7".to_string()),
                    ]),
                },
                ParsedSample {
                    sample_name: "NORMAL".to_string(),
                    format_fields: HashMap::from([
                        ("DP".to_string(), "60".to_string()),
                        ("AD".to_string(), "5,27,33".to_string()),
                    ]),
                },
            ],
        });

        let mafs = MafRecord::from_reformatted_record_multi_for_samples(
            &record, "c", "GRCh38", "s", Some("TUMOR"), None,
        )
        .unwrap();

        assert_eq!(mafs.len(), 2);
        assert_eq!(mafs[0].t_alt_count, Some(3), "first ALT takes the tumor's AD[1]");
        assert_eq!(mafs[1].t_alt_count, Some(7), "second ALT takes the tumor's AD[2]");
        for m in &mafs {
            let (d, r, a) = (m.t_depth.unwrap(), m.t_ref_count.unwrap(), m.t_alt_count.unwrap());
            assert!(r + a <= d, "t_ref {r} + t_alt {a} exceeds t_depth {d}");
        }
    }

    #[test]
    fn a_tumor_with_no_call_reports_no_depths_rather_than_the_pooled_info() {
        // freebayes writes ".:.:.:." for a sample it made no call in. INFO DP/RO/AO still
        // carry the other sample's reads, and falling back to them reported the normal's
        // depth as the tumor's on 75 rows of B487_1_V_vs_B487_1_cOM. vcf2maf leaves the
        // whole depth family empty here, and so must we.
        use crate::extract_sample_info::{ParsedFormatSample, ParsedSample};
        let mut record = create_test_maf_record(
            "chr1", 21899250, "G", "C", Some(50.0), "PASS",
            HashMap::from([
                ("INFO_DP".to_string(), "2".to_string()),
                ("INFO_RO".to_string(), "0".to_string()),
                ("INFO_AO".to_string(), "2".to_string()),
            ]),
        );
        record.format_sample_data = Some(ParsedFormatSample {
            format_keys: vec!["DP".to_string(), "AD".to_string()],
            samples: vec![
                ParsedSample {
                    sample_name: "B487_1_V".to_string(),
                    format_fields: HashMap::from([
                        ("DP".to_string(), ".".to_string()),
                        ("AD".to_string(), ".".to_string()),
                    ]),
                },
                ParsedSample {
                    sample_name: "B487_1_cOM".to_string(),
                    format_fields: HashMap::from([
                        ("DP".to_string(), "2".to_string()),
                        ("AD".to_string(), "0,2".to_string()),
                    ]),
                },
            ],
        });

        let maf = MafRecord::from_reformatted_record_for_samples(
            &record, "c", "GRCh38", "s", Some("B487_1_V"), None,
        )
        .unwrap();
        assert_eq!(maf.t_depth, None);
        assert_eq!(maf.t_ref_count, None);
        assert_eq!(maf.t_alt_count, None);
    }

    #[test]
    fn tumor_counts_all_come_from_the_tumor_sample_not_pooled_info() {
        // The bug: t_depth was read from the sample while t_ref_count/t_alt_count came from
        // INFO RO/AO, which sum every sample. That produced t_ref + t_alt = 7 against a
        // t_depth of 1, and credited the normal's 2 alt reads to the tumor.
        let record = tumor_normal_record();
        let maf = MafRecord::from_reformatted_record_for_samples(
            &record, "c", "GRCh38", "s", Some("B487_1_V"), None,
        )
        .unwrap();

        assert_eq!(maf.t_depth, Some(1));
        assert_eq!(maf.t_ref_count, Some(1));
        assert_eq!(maf.t_alt_count, Some(0));
    }

    #[test]
    fn tumor_ref_and_alt_counts_never_exceed_tumor_depth() {
        let record = tumor_normal_record();
        let maf = MafRecord::from_reformatted_record_for_samples(
            &record, "c", "GRCh38", "s", Some("B487_1_V"), None,
        )
        .unwrap();
        let (d, r, a) = (
            maf.t_depth.unwrap(),
            maf.t_ref_count.unwrap(),
            maf.t_alt_count.unwrap(),
        );
        assert!(r + a <= d, "t_ref_count {r} + t_alt_count {a} exceeds t_depth {d}");
    }

    #[test]
    fn naming_the_normal_sample_populates_the_matched_normal_columns() {
        let record = tumor_normal_record();
        let maf = MafRecord::from_reformatted_record_for_samples(
            &record, "c", "GRCh38", "s", Some("B487_1_V"), Some("B487_1_cOM"),
        )
        .unwrap();

        assert_eq!(maf.n_depth, Some(6));
        assert_eq!(maf.n_ref_count, Some(4));
        assert_eq!(maf.n_alt_count, Some(2));
        assert_eq!(maf.matched_norm_sample_barcode.as_deref(), Some("B487_1_cOM"));
    }

    #[test]
    fn matched_normal_columns_stay_empty_when_no_normal_is_named() {
        let record = tumor_normal_record();
        let maf = MafRecord::from_reformatted_record_for_samples(
            &record, "c", "GRCh38", "s", Some("B487_1_V"), None,
        )
        .unwrap();

        assert_eq!(maf.n_depth, None);
        assert_eq!(maf.n_ref_count, None);
        assert_eq!(maf.n_alt_count, None);
    }

    #[test]
    fn an_unknown_sample_name_yields_nothing_rather_than_the_first_sample() {
        // Falling back to sample 1 would reproduce the guess this parameter exists to remove.
        let record = tumor_normal_record();
        let maf = MafRecord::from_reformatted_record_for_samples(
            &record, "c", "GRCh38", "s", Some("NOT_IN_THIS_VCF"), None,
        )
        .unwrap();

        assert_eq!(maf.t_ref_count, None);
        assert_eq!(maf.t_alt_count, None);
    }

    #[test]
    fn without_a_named_tumor_the_first_sample_is_still_used() {
        // Existing single-sample behaviour must not change; only naming a sample changes it.
        let record = tumor_normal_record();
        let maf =
            MafRecord::from_reformatted_record(&record, "c", "GRCh38", "s").unwrap();
        assert_eq!(maf.t_depth, Some(1));
    }

    // Local equivalent of tests/test.rs's `create_test_maf_record` helper — that helper lives
    // in the integration-test crate and isn't reachable from this unit-test module.
    fn create_test_maf_record(
        chromosome: &str,
        position: u64,
        reference: &str,
        alternate: &str,
        quality: Option<f64>,
        filter: &str,
        info_fields: HashMap<String, String>,
    ) -> ReformattedVcfRecord {
        ReformattedVcfRecord {
            chromosome: chromosome.to_string(),
            position,
            id: Some("rs123456".to_string()),
            reference: reference.to_string(),
            alternate: alternate.to_string(),
            quality,
            filter: filter.to_string(),
            info_fields,
            format_sample_data: None,
            annotation_field_type: crate::reformat_vcf::AnnotationFieldType::None,
        }
    }

    #[test]
    fn test_maf_headers_match_vcf2maf_core_46_plus_custom() {
        let headers = MafRecord::get_maf_headers();
        assert_eq!(headers.len(), 50);
        let expected = [
            "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
            "Start_Position", "End_Position", "Strand", "Variant_Classification",
            "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode",
            "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
            "Tumor_Validation_Allele1", "Tumor_Validation_Allele2",
            "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2",
            "Verification_Status", "Validation_Status", "Mutation_Status",
            "Sequencing_Phase", "Sequence_Source", "Validation_Method", "Score",
            "BAM_File", "Sequencer", "Tumor_Sample_UUID", "Matched_Norm_Sample_UUID",
            "HGVSc", "HGVSp", "HGVSp_Short", "Transcript_ID", "Exon_Number",
            "t_depth", "t_ref_count", "t_alt_count", "n_depth", "n_ref_count",
            "n_alt_count", "all_effects",
            "FILTER", "QUAL", "VAF", "Protein_Position",
        ];
        assert_eq!(headers, expected.to_vec());
    }

    #[test]
    fn test_to_tsv_line_leaves_unsupported_columns_empty() {
        let record =
            create_test_maf_record("chr1", 100, "A", "G", Some(60.0), "PASS", HashMap::new());
        let maf =
            MafRecord::from_reformatted_record(&record, "TestCenter", "GRCh38", "SAMPLE-001")
                .unwrap();
        let tsv = maf.to_tsv_line();
        let fields: Vec<&str> = tsv.split('\t').collect();
        assert_eq!(fields.len(), 50);
        // Match_Norm_Seq_Allele1 (idx 17), Verification_Status (idx 23), Score (idx 29),
        // n_depth (idx 42), all_effects (idx 45) — no data source, so empty as vcf2maf
        // writes them. A "." here was our own invention and the largest diff class.
        for idx in [17, 23, 29, 42, 45] {
            assert_eq!(fields[idx], "", "column {idx} should be empty");
        }
    }

    fn dbsnp_from(id: Option<&str>, existing_variation: Option<&str>) -> Option<String> {
        let mut info_fields = HashMap::new();
        if let Some(ev) = existing_variation {
            info_fields.insert("CSQ_Existing_variation".to_string(), ev.to_string());
        }
        let record = ReformattedVcfRecord {
            chromosome: "chr1".to_string(),
            position: 100,
            id: id.map(|s| s.to_string()),
            reference: "A".to_string(),
            alternate: "G".to_string(),
            quality: Some(60.0),
            filter: "PASS".to_string(),
            info_fields,
            format_sample_data: None,
            annotation_field_type: crate::reformat_vcf::AnnotationFieldType::None,
        };
        MafRecord::from_reformatted_record(&record, "test", "GRCh38", "sample")
            .unwrap()
            .dbsnp_rs
    }

    #[test]
    fn test_dbsnp_rs_keeps_only_rs_ids_from_existing_variation() {
        // Real CSQ value from B505_1_V.mutect2.filtered_VEP.ann.vcf.gz. VEP joins co-located
        // variants with "&"; vcf2maf.pl:784 rewrites those to "," before filtering on /^rs\d+$/.
        assert_eq!(
            dbsnp_from(None, Some("rs992327&COSV57258734")),
            Some("rs992327".to_string())
        );
    }

    #[test]
    fn test_dbsnp_rs_joins_multiple_rs_ids_with_commas() {
        // vcf2maf.pl:856 — join( ",", grep{m/^rs\d+$/} ... )
        assert_eq!(
            dbsnp_from(None, Some("rs123&COSV1&rs456")),
            Some("rs123,rs456".to_string())
        );
    }

    #[test]
    fn test_dbsnp_rs_blank_when_variant_known_only_outside_dbsnp() {
        // vcf2maf.pl:855 — "If seen in a DB other than dbSNP, this field will remain blank"
        assert_eq!(dbsnp_from(None, Some("COSV57258734")), None);
    }

    #[test]
    fn test_dbsnp_rs_is_novel_when_vep_found_no_existing_variation() {
        // vcf2maf.pl:858-860 — VEP looked and came back empty, which is a result, not missing data.
        assert_eq!(dbsnp_from(None, Some(".")), Some("novel".to_string()));
    }

    #[test]
    fn test_dbsnp_rs_falls_back_to_vcf_id_when_csq_absent() {
        // SnpEff's ANN has no Existing_variation equivalent, so the ID column is all we have.
        assert_eq!(dbsnp_from(Some("rs123"), None), Some("rs123".to_string()));
    }

    #[test]
    fn test_dbsnp_rs_none_without_annotation_or_id() {
        assert_eq!(dbsnp_from(None, None), None);
    }
}
