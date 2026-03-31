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
        writeln!(file, "Total input variants:     {}", self.input_variant_count)?;
        writeln!(file, "Variants per chromosome:")?;
        write!(file, "{}", format_chrom_table(&self.input_chrom_counts, self.input_variant_count))?;
        writeln!(file)?;

        writeln!(file, "OUTPUT STATISTICS")?;
        writeln!(file, "-----------------")?;
        writeln!(file, "Total output records:     {}", self.output_record_count)?;
        writeln!(file, "Records per chromosome:")?;
        write!(file, "{}", format_chrom_table(&self.output_chrom_counts, self.output_record_count))?;
        writeln!(file)?;

        let expansion = if self.input_variant_count > 0 {
            self.output_record_count as f64 / self.input_variant_count as f64
        } else {
            0.0
        };

        writeln!(file, "PROCESSING")?;
        writeln!(file, "----------")?;
        writeln!(file, "Expansion ratio:          {:.2}x", expansion)?;
        writeln!(file, "Processing time:          {:.2}s", self.processing_time_secs)?;
        writeln!(file, "Processing rate:          {:.0} variants/sec", self.variants_per_sec)?;

        Ok(())
    }
}

/// Count variants per chromosome from raw VCF data lines.
/// Each line starts with the chromosome name followed by a tab.
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
    let mut entries: Vec<(String, usize)> = counts.into_iter().collect();
    entries.sort_by(|(a, _), (b, _)| chrom_sort_key(a).cmp(&chrom_sort_key(b)));
    entries.into_iter().collect()
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

#[cfg(test)]
mod tests {
    use super::*;

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
}