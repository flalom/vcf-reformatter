use crate::reformat_vcf::ReformattedVcfRecord;
use indexmap::IndexMap;
use std::fmt::Write as FmtWrite;
use std::io::Write;

/// Holds all statistics for a processing run
pub struct SummaryStats {
    pub input_file: String,
    pub output_format: String,
    pub transcript_handling: String,
    pub input_variant_count: usize,
    pub output_record_count: usize,
    pub input_chrom_counts: IndexMap<String, usize>,
    pub output_chrom_counts: IndexMap<String, usize>,
    pub processing_time_secs: f64,
    pub variants_per_sec: f64,
}

impl SummaryStats {
    pub fn write_to_file(&self, path: &str) -> std::io::Result<()> {
        let mut file = std::fs::File::create(path)?;
        let timestamp = chrono::Local::now().format("%Y-%m-%d %H:%M:%S");

        writeln!(file, "VCF REFORMATTER - PROCESSING SUMMARY")?;
        writeln!(file, "=====================================")?;
        writeln!(file, "Input file:          {}", self.input_file)?;
        writeln!(file, "Output format:       {}", self.output_format)?;
        writeln!(file, "Transcript handling: {}", self.transcript_handling)?;
        writeln!(file, "Date:                {}", timestamp)?;
        writeln!(file)?;

        writeln!(file, "INPUT STATISTICS")?;
        writeln!(file, "----------------")?;
        writeln!(
            file,
            "Total input variants:     {}",
            self.input_variant_count
        )?;
        writeln!(file, "Variants per chromosome:")?;
        write!(
            file,
            "{}",
            format_chrom_table(&self.input_chrom_counts, self.input_variant_count)
        )?;
        writeln!(file)?;

        writeln!(file, "OUTPUT STATISTICS")?;
        writeln!(file, "-----------------")?;
        writeln!(
            file,
            "Total output records:     {}",
            self.output_record_count
        )?;
        writeln!(file, "Records per chromosome:")?;
        write!(
            file,
            "{}",
            format_chrom_table(&self.output_chrom_counts, self.output_record_count)
        )?;
        writeln!(file)?;

        let expansion = if self.input_variant_count > 0 {
            self.output_record_count as f64 / self.input_variant_count as f64
        } else {
            0.0
        };

        writeln!(file, "PROCESSING")?;
        writeln!(file, "----------")?;
        writeln!(file, "Expansion ratio:          {:.2}x", expansion)?;
        writeln!(
            file,
            "Processing time:          {:.2}s",
            self.processing_time_secs
        )?;
        writeln!(
            file,
            "Processing rate:          {:.0} variants/sec",
            self.variants_per_sec
        )?;

        Ok(())
    }
}

/// Count VCF data lines whose ALT field (5th tab-separated column) lists
/// more than one comma-separated allele.
pub fn count_multiallelic_sites(data_lines: &[String]) -> usize {
    data_lines
        .iter()
        .filter(|line| line.split('\t').nth(4).is_some_and(|alt| alt.contains(',')))
        .count()
}

/// Count annotation entries whose `|`-separated field count falls short of the count the
/// header's `Format:` string declares — a truncated CSQ/ANN entry. The variant still converts
/// (absent fields read as empty, so the row degrades to `Hugo_Symbol=Unknown` rather than
/// disappearing), which is why this is a note on stderr and not a failure: one malformed
/// annotation must not cost the other 92,000 variants on the file.
pub fn count_malformed_annotation_entries(
    data_lines: &[String],
    key: &str,
    expected_fields: usize,
) -> usize {
    if expected_fields == 0 {
        return 0;
    }
    data_lines
        .iter()
        .filter_map(|line| line.split('\t').nth(7))
        .flat_map(|info| info.split(';'))
        .filter_map(|entry| entry.strip_prefix(key).and_then(|v| v.strip_prefix('=')))
        .flat_map(|value| value.split(','))
        .filter(|annotation| annotation.split('|').count() < expected_fields)
        .count()
}

/// Count VCF data lines whose CSQ (VEP) or ANN (SnpEff) annotation field lists more than one
/// comma-separated transcript entry. These are the only sites where the transcript-handling
/// mode changes what gets reported: `first` takes the annotator's own first entry without
/// re-ranking it, so the consequence it reports need not be the most damaging one present.
/// VEP guarantees a single entry per variant only when run with `--pick`.
pub fn count_multi_transcript_sites(data_lines: &[String]) -> usize {
    data_lines
        .iter()
        .filter(|line| {
            line.split('\t')
                .nth(7)
                .and_then(|info| {
                    info.split(';').find_map(|field| {
                        field
                            .strip_prefix("CSQ=")
                            .or_else(|| field.strip_prefix("ANN="))
                    })
                })
                .is_some_and(|annotation| annotation.contains(','))
        })
        .count()
}

/// Count variants per chromosome from raw VCF data lines.
/// Each line starts with the chromosome name followed by a tab.
// Superseded by the streaming path; kept because tests/test.rs still exercises it.
#[allow(dead_code)]
pub fn count_input_chromosomes(data_lines: &[String]) -> IndexMap<String, usize> {
    let mut counts = IndexMap::new();
    for line in data_lines {
        if let Some(chrom) = line.split('\t').next() {
            *counts.entry(chrom.to_string()).or_insert(0) += 1;
        }
    }
    sort_chromosomes(counts)
}

/// Sort chromosome keys in natural order: 1-22, X, Y, M/MT, then others alphabetically.
pub fn sort_chromosomes(counts: IndexMap<String, usize>) -> IndexMap<String, usize> {
    sort_by_chromosome(counts.into_iter().collect())
        .into_iter()
        .collect()
}

/// Sort any chromosome-keyed entries in natural order (see `chrom_sort_key`).
/// Shared by `sort_chromosomes` and `compute_damage_breakdowns`.
fn sort_by_chromosome<V>(mut entries: Vec<(String, V)>) -> Vec<(String, V)> {
    entries.sort_by_key(|(a, _)| chrom_sort_key(a));
    entries
}

/// Generate a sort key for a chromosome name.
/// Numeric chromosomes sort first (by number), then X, Y, M/MT, then everything else.
fn chrom_sort_key(chrom: &str) -> (u8, u32, String) {
    let name = chrom.strip_prefix("chr").unwrap_or(chrom);
    if let Ok(num) = name.parse::<u32>() {
        (0, num, String::new()) // Numeric: sort group 0, by number
    } else {
        match name.to_uppercase().as_str() {
            "X" => (1, 0, String::new()),
            "Y" => (1, 1, String::new()),
            "M" | "MT" => (1, 2, String::new()),
            _ => (2, 0, name.to_string()), // Non-standard: sort group 2, alphabetically
        }
    }
}

/// Format a chromosome count table as a string for the summary report.
pub fn format_chrom_table(counts: &IndexMap<String, usize>, total: usize) -> String {
    let mut output = String::new();
    for (chrom, count) in counts {
        let pct = if total > 0 {
            *count as f64 / total as f64 * 100.0
        } else {
            0.0
        };
        writeln!(output, "  {:<10} {:>8}  ({:.1}%)", chrom, count, pct).unwrap();
    }
    output
}

/// Per-chromosome breakdown of a single annotation-severity metric (SIFT,
/// PolyPhen, or SnpEff's Impact), used to drive the HTML report's
/// stacked-bar chart.
pub struct DamageBreakdown {
    pub metric_name: String,
    /// Ordered severity-first (most severe category first).
    pub categories: Vec<String>,
    pub per_chrom_counts: IndexMap<String, IndexMap<String, usize>>,
}

const SIFT_CATEGORIES: [&str; 4] = [
    "deleterious",
    "deleterious_low_confidence",
    "tolerated_low_confidence",
    "tolerated",
];
const POLYPHEN_CATEGORIES: [&str; 4] = [
    "probably_damaging",
    "possibly_damaging",
    "benign",
    "unknown",
];
const IMPACT_CATEGORIES: [&str; 4] = ["HIGH", "MODERATE", "LOW", "MODIFIER"];

/// Parse VEP's `"category(score)"` format (e.g. `"deleterious(0.02)"`) into
/// just the category. Returns `None` for missing/empty/`"."` values.
fn parse_bracketed_category(raw: &str) -> Option<String> {
    let raw = raw.trim();
    if raw.is_empty() || raw == "." {
        return None;
    }
    let category = raw.split('(').next().unwrap_or(raw).trim();
    if category.is_empty() {
        None
    } else {
        Some(category.to_string())
    }
}

fn parse_impact_category(raw: &str) -> Option<String> {
    let raw = raw.trim();
    if raw.is_empty() || raw == "." {
        None
    } else {
        Some(raw.to_uppercase())
    }
}

fn build_metric_breakdown(
    records: &[ReformattedVcfRecord],
    info_key: &str,
    metric_name: &str,
    known_categories: &[&str],
    parse: fn(&str) -> Option<String>,
) -> Option<DamageBreakdown> {
    let mut per_chrom_counts: IndexMap<String, IndexMap<String, usize>> = IndexMap::new();
    let mut seen_categories: Vec<String> = Vec::new();

    for record in records {
        let Some(raw) = record.info_fields.get(info_key) else {
            continue;
        };
        let Some(category) = parse(raw) else {
            continue;
        };
        if !seen_categories.contains(&category) {
            seen_categories.push(category.clone());
        }
        *per_chrom_counts
            .entry(record.chromosome.clone())
            .or_default()
            .entry(category)
            .or_insert(0) += 1;
    }

    if seen_categories.is_empty() {
        return None;
    }

    let mut categories: Vec<String> = known_categories
        .iter()
        .map(|s| s.to_string())
        .filter(|c| seen_categories.contains(c))
        .collect();
    let mut extra: Vec<String> = seen_categories
        .into_iter()
        .filter(|c| !categories.contains(c))
        .collect();
    extra.sort();
    categories.extend(extra);

    let per_chrom_counts = sort_by_chromosome(per_chrom_counts.into_iter().collect())
        .into_iter()
        .collect();

    Some(DamageBreakdown {
        metric_name: metric_name.to_string(),
        categories,
        per_chrom_counts,
    })
}

/// Compute per-chromosome annotation-severity breakdowns from whichever
/// metrics are actually present in the data. VEP input yields SIFT and/or
/// PolyPhen entries; SnpEff input yields an Impact entry. Neither present
/// yields an empty `Vec`, in which case the HTML report omits the chart
/// section entirely.
/// Fold one chunk's breakdowns into an accumulator, so a file can be summarised without ever
/// holding all of its records at once. Summing counts per chromosome per category is the whole
/// operation — the result must equal `compute_damage_breakdowns` over the concatenated input.
pub fn merge_damage_breakdowns(acc: &mut Vec<DamageBreakdown>, next: Vec<DamageBreakdown>) {
    for incoming in next {
        match acc
            .iter_mut()
            .find(|b| b.metric_name == incoming.metric_name)
        {
            Some(existing) => {
                for cat in incoming.categories {
                    if !existing.categories.contains(&cat) {
                        existing.categories.push(cat);
                    }
                }
                for (chrom, counts) in incoming.per_chrom_counts {
                    let entry = existing.per_chrom_counts.entry(chrom).or_default();
                    for (cat, n) in counts {
                        *entry.entry(cat).or_insert(0) += n;
                    }
                }
            }
            None => acc.push(incoming),
        }
    }
}

pub fn compute_damage_breakdowns(records: &[ReformattedVcfRecord]) -> Vec<DamageBreakdown> {
    let mut result = Vec::new();

    if let Some(b) = build_metric_breakdown(
        records,
        "CSQ_SIFT",
        "SIFT",
        &SIFT_CATEGORIES,
        parse_bracketed_category,
    ) {
        result.push(b);
    }
    if let Some(b) = build_metric_breakdown(
        records,
        "CSQ_PolyPhen",
        "PolyPhen",
        &POLYPHEN_CATEGORIES,
        parse_bracketed_category,
    ) {
        result.push(b);
    }
    if let Some(b) = build_metric_breakdown(
        records,
        "ANN_Annotation_Impact",
        "Impact",
        &IMPACT_CATEGORIES,
        parse_impact_category,
    ) {
        result.push(b);
    }

    result
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_count_malformed_annotation_entries() {
        let expected = 8; // Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE
        let lines = vec![
            // full entry: 8 fields, fine
            "chr1\t1\t.\tA\tG\t.\tPASS\tDP=3;CSQ=G|missense_variant|MODERATE|G1|E1|Transcript|T1|protein_coding".to_string(),
            // three fields: truncated
            "chr1\t2\t.\tC\tT\t.\tPASS\tDP=3;CSQ=T|synonymous_variant|LOW".to_string(),
            // empty value counts as one field, so also short
            "chr1\t3\t.\tG\tA\t.\tPASS\tDP=3;CSQ=".to_string(),
            // two comma-separated entries, only the second short
            "chr1\t4\t.\tT\tC\t.\tPASS\tCSQ=C|a|b|c|d|e|f|g,C|short".to_string(),
        ];
        assert_eq!(
            count_malformed_annotation_entries(&lines, "CSQ", expected),
            3
        );
    }

    #[test]
    fn test_count_malformed_annotation_entries_ignores_other_info_keys() {
        // A key that merely ends in CSQ, and a pipe inside an unrelated INFO value, must not
        // be mistaken for the annotation field.
        let lines = vec!["chr1\t1\t.\tA\tG\t.\tPASS\tMY_CSQ=x|y;OTHER=a|b".to_string()];
        assert_eq!(count_malformed_annotation_entries(&lines, "CSQ", 8), 0);
    }

    #[test]
    fn test_count_malformed_annotation_entries_no_format_declared() {
        let lines = vec!["chr1\t1\t.\tA\tG\t.\tPASS\tCSQ=G|x".to_string()];
        assert_eq!(count_malformed_annotation_entries(&lines, "CSQ", 0), 0);
    }

    use super::*;

    #[test]
    fn test_count_multiallelic_sites_none() {
        let lines = vec![
            "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50".to_string(),
            "chr1\t200\t.\tC\tT\t40\tPASS\tDP=30".to_string(),
        ];
        assert_eq!(count_multiallelic_sites(&lines), 0);
    }

    #[test]
    fn test_count_multiallelic_sites_some() {
        let lines = vec![
            "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50".to_string(),
            "chr1\t200\t.\tC\tT,A\t40\tPASS\tDP=30".to_string(),
            "chr2\t300\t.\tG\tA,T,C\t50\tPASS\tDP=20".to_string(),
        ];
        assert_eq!(count_multiallelic_sites(&lines), 2);
    }

    #[test]
    fn test_count_multiallelic_sites_ignores_commas_outside_alt() {
        // A comma inside INFO (e.g. an AF list) must not be mistaken for a
        // multiallelic ALT field.
        let lines = vec!["chr1\t100\t.\tA\tG\t60\tPASS\tAF=0.1,0.2".to_string()];
        assert_eq!(count_multiallelic_sites(&lines), 0);
    }

    #[test]
    fn test_count_multi_transcript_sites() {
        let lines = vec![
            // single CSQ entry
            "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50;CSQ=G|missense_variant|MODERATE|BRAF".to_string(),
            // two CSQ entries
            "chr1\t200\t.\tC\tT\t40\tPASS\tCSQ=T|intron_variant||X,T|missense_variant||X"
                .to_string(),
            // two SnpEff ANN entries
            "chr2\t300\t.\tG\tA\t50\tPASS\tANN=A|synonymous_variant||Y,A|stop_gained||Y"
                .to_string(),
            // comma in a different INFO field must not count
            "chr3\t400\t.\tT\tC\t50\tPASS\tAF=0.1,0.2;CSQ=C|intron_variant||Z".to_string(),
            // no annotation at all
            "chr4\t500\t.\tT\tC\t50\tPASS\tDP=10".to_string(),
        ];
        assert_eq!(count_multi_transcript_sites(&lines), 2);
    }

    #[test]
    fn test_count_input_chromosomes() {
        let lines = vec![
            "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50".to_string(),
            "chr1\t200\t.\tC\tT\t40\tPASS\tDP=30".to_string(),
            "chr2\t300\t.\tG\tA\t50\tPASS\tDP=20".to_string(),
        ];
        let counts = count_input_chromosomes(&lines);
        assert_eq!(counts.get("chr1"), Some(&2));
        assert_eq!(counts.get("chr2"), Some(&1));
        assert_eq!(counts.len(), 2);
    }

    #[test]
    fn test_chromosome_sort_order() {
        let mut counts = IndexMap::new();
        counts.insert("chrX".to_string(), 10);
        counts.insert("chr2".to_string(), 20);
        counts.insert("chr10".to_string(), 5);
        counts.insert("chr1".to_string(), 30);
        counts.insert("chrY".to_string(), 2);
        counts.insert("chrM".to_string(), 1);

        let sorted = sort_chromosomes(counts);
        let keys: Vec<&String> = sorted.keys().collect();
        assert_eq!(keys, vec!["chr1", "chr2", "chr10", "chrX", "chrY", "chrM"]);
    }

    #[test]
    fn test_chromosome_sort_without_chr_prefix() {
        let mut counts = IndexMap::new();
        counts.insert("X".to_string(), 10);
        counts.insert("2".to_string(), 20);
        counts.insert("10".to_string(), 5);
        counts.insert("1".to_string(), 30);
        counts.insert("Y".to_string(), 2);

        let sorted = sort_chromosomes(counts);
        let keys: Vec<&String> = sorted.keys().collect();
        assert_eq!(keys, vec!["1", "2", "10", "X", "Y"]);
    }

    #[test]
    fn test_summary_counts_sum_to_total() {
        let lines = vec![
            "chr1\t100\t.\tA\tG\t60\tPASS\tDP=50".to_string(),
            "chr1\t200\t.\tC\tT\t40\tPASS\tDP=30".to_string(),
            "chr2\t300\t.\tG\tA\t50\tPASS\tDP=20".to_string(),
            "chr3\t400\t.\tT\tC\t70\tPASS\tDP=40".to_string(),
        ];
        let counts = count_input_chromosomes(&lines);
        let total: usize = counts.values().sum();
        assert_eq!(total, lines.len());
    }

    #[test]
    fn test_format_chrom_table() {
        let mut counts = IndexMap::new();
        counts.insert("chr1".to_string(), 75);
        counts.insert("chr2".to_string(), 25);
        let total = 100;
        let table = format_chrom_table(&counts, total);
        assert!(table.contains("chr1"));
        assert!(table.contains("75.0%"));
        assert!(table.contains("chr2"));
        assert!(table.contains("25.0%"));
    }

    fn make_record(chrom: &str, info: &[(&str, &str)]) -> ReformattedVcfRecord {
        use crate::reformat_vcf::AnnotationFieldType;
        use std::collections::HashMap;

        let mut info_fields = HashMap::new();
        for (k, v) in info {
            info_fields.insert(k.to_string(), v.to_string());
        }
        ReformattedVcfRecord {
            chromosome: chrom.to_string(),
            position: 100,
            id: None,
            reference: "A".to_string(),
            alternate: "G".to_string(),
            quality: None,
            filter: "PASS".to_string(),
            info_fields,
            format_sample_data: None,
            annotation_field_type: AnnotationFieldType::None,
        }
    }

    #[test]
    fn test_parse_bracketed_category() {
        assert_eq!(
            parse_bracketed_category("deleterious(0.02)"),
            Some("deleterious".to_string())
        );
        assert_eq!(
            parse_bracketed_category("probably_damaging(0.967)"),
            Some("probably_damaging".to_string())
        );
        assert_eq!(parse_bracketed_category("."), None);
        assert_eq!(parse_bracketed_category(""), None);
    }

    #[test]
    fn test_compute_damage_breakdowns_sift_and_polyphen() {
        let records = vec![
            make_record(
                "chr1",
                &[
                    ("CSQ_SIFT", "deleterious(0.01)"),
                    ("CSQ_PolyPhen", "probably_damaging(0.99)"),
                ],
            ),
            make_record(
                "chr1",
                &[
                    ("CSQ_SIFT", "tolerated(0.8)"),
                    ("CSQ_PolyPhen", "benign(0.05)"),
                ],
            ),
            // No PolyPhen value on this one — must not count toward PolyPhen totals.
            make_record("chr2", &[("CSQ_SIFT", "deleterious(0.02)")]),
        ];

        let breakdowns = compute_damage_breakdowns(&records);
        assert_eq!(breakdowns.len(), 2);

        let sift = breakdowns.iter().find(|b| b.metric_name == "SIFT").unwrap();
        assert_eq!(sift.categories, vec!["deleterious", "tolerated"]);
        assert_eq!(sift.per_chrom_counts["chr1"]["deleterious"], 1);
        assert_eq!(sift.per_chrom_counts["chr1"]["tolerated"], 1);
        assert_eq!(sift.per_chrom_counts["chr2"]["deleterious"], 1);

        let polyphen = breakdowns
            .iter()
            .find(|b| b.metric_name == "PolyPhen")
            .unwrap();
        assert_eq!(polyphen.categories, vec!["probably_damaging", "benign"]);
        assert!(!polyphen.per_chrom_counts.contains_key("chr2"));
    }

    #[test]
    fn merging_two_chunks_equals_computing_over_the_whole() {
        // The chunked MAF path never holds every record, so the merged counts must match what
        // a single pass over all of them would have produced.
        let all = vec![
            make_record("chr1", &[("CSQ_SIFT", "deleterious(0.01)")]),
            make_record("chr1", &[("CSQ_SIFT", "tolerated(0.4)")]),
            make_record("chr2", &[("CSQ_SIFT", "deleterious(0.02)")]),
            make_record("chr2", &[("CSQ_SIFT", "deleterious(0.03)")]),
        ];
        let whole = compute_damage_breakdowns(&all);

        let mut merged = Vec::new();
        for chunk in all.chunks(2) {
            merge_damage_breakdowns(&mut merged, compute_damage_breakdowns(chunk));
        }

        assert_eq!(merged.len(), whole.len());
        for (m, w) in merged.iter().zip(whole.iter()) {
            assert_eq!(m.metric_name, w.metric_name);
            assert_eq!(m.per_chrom_counts, w.per_chrom_counts);
        }
    }

    #[test]
    fn test_compute_damage_breakdowns_impact_only_for_snpeff() {
        let records = vec![
            make_record("chr1", &[("ANN_Annotation_Impact", "HIGH")]),
            make_record("chr1", &[("ANN_Annotation_Impact", "LOW")]),
            make_record("chr2", &[("ANN_Annotation_Impact", "MODERATE")]),
        ];

        let breakdowns = compute_damage_breakdowns(&records);
        assert_eq!(breakdowns.len(), 1);
        let impact = &breakdowns[0];
        assert_eq!(impact.metric_name, "Impact");
        assert_eq!(impact.categories, vec!["HIGH", "MODERATE", "LOW"]);
        assert_eq!(impact.per_chrom_counts["chr1"]["HIGH"], 1);
        assert_eq!(impact.per_chrom_counts["chr1"]["LOW"], 1);
        assert_eq!(impact.per_chrom_counts["chr2"]["MODERATE"], 1);
    }

    #[test]
    fn test_compute_damage_breakdowns_empty_when_no_metrics_present() {
        let records = vec![make_record(
            "chr1",
            &[("CSQ_Consequence", "missense_variant")],
        )];
        assert!(compute_damage_breakdowns(&records).is_empty());
    }
}
