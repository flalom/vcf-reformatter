//! `--report`: which files appear, what the HTML page charts, and the txt summary layout.
mod common;

use indexmap::IndexMap;
use tempfile::tempdir;

use common::{bin, write_vcf, TINY_VCF};
use vcf_reformatter::summary::SummaryStats;

#[test]
fn html_is_the_default_and_charts_sift_and_polyphen_for_vep() {
    let dir = tempdir().unwrap();
    let vcf = write_vcf(
        dir.path(),
        "in.vcf",
        "##fileformat=VCFv4.2\n\
         ##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|SIFT|PolyPhen\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         chr1\t100\t.\tA\tG\t60\tPASS\tCSQ=G|missense_variant|MODERATE|BRCA1|deleterious(0.01)|probably_damaging(0.99)\n\
         chr2\t200\t.\tC\tT\t60\tPASS\tCSQ=T|missense_variant|MODERATE|TP53|tolerated(0.8)|benign(0.02)\n",
    );
    let o = bin()
        .args([&vcf, "-o", dir.path().to_str().unwrap(), "-p", "r"])
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&o.stderr)
    );

    let html = std::fs::read_to_string(dir.path().join("r_summary.html")).unwrap();
    assert!(html.contains("VCF-REFORMATTER"));
    assert!(html.contains("chr1") && html.contains("chr2"));
    assert!(html.contains("damage-metric-select"));
    assert!(html.contains("SIFT") && html.contains("PolyPhen"));
}

#[test]
fn html_for_snpeff_charts_impact_and_hides_the_metric_selector() {
    let dir = tempdir().unwrap();
    let vcf = write_vcf(
        dir.path(),
        "in.vcf",
        "##fileformat=VCFv4.2\n\
         ##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name'\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         chr1\t100\t.\tA\tG\t60\tPASS\tANN=G|missense_variant|HIGH|BRCA1\n\
         chr2\t200\t.\tC\tT\t60\tPASS\tANN=T|synonymous_variant|LOW|TP53\n",
    );
    let o = bin()
        .args([&vcf, "-o", dir.path().to_str().unwrap(), "-p", "r"])
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&o.stderr)
    );

    let html = std::fs::read_to_string(dir.path().join("r_summary.html")).unwrap();
    assert!(html.contains("Impact"));
    assert!(!html.contains("SIFT"));
    assert!(html.contains("\"HIGH\"") && html.contains("\"LOW\""));
    // One breakdown only, so there is nothing to switch between; the chart itself stays.
    assert!(!html.contains("damage-metric-select"));
    assert!(html.contains("damage-chart"));
}

#[test]
fn report_formats_txt_none_both_and_report_dir() {
    let dir = tempdir().unwrap();
    let vcf = write_vcf(dir.path(), "both.vcf", TINY_VCF);
    let run = |extra: &[&str]| {
        let o = bin()
            .args([&vcf, "-p", "both"])
            .args(extra)
            .output()
            .unwrap();
        assert!(
            o.status.success(),
            "{extra:?}: {}",
            String::from_utf8_lossy(&o.stderr)
        );
    };
    let exists = |d: &std::path::Path, f: &str| d.join(f).exists();

    let txt = dir.path().join("txt");
    run(&["--report", "txt", "-o", txt.to_str().unwrap()]);
    assert!(exists(&txt, "both_summary.txt") && !exists(&txt, "both_summary.html"));

    let none = dir.path().join("none");
    run(&["--report", "none", "-o", none.to_str().unwrap()]);
    assert!(!exists(&none, "both_summary.txt") && !exists(&none, "both_summary.html"));

    let both = dir.path().join("both");
    run(&["--report", "html,txt", "-o", both.to_str().unwrap()]);
    assert!(exists(&both, "both_summary.txt") && exists(&both, "both_summary.html"));

    // `none` anywhere in the list wins.
    let mixed = dir.path().join("mixed");
    run(&["--report", "html,none", "-o", mixed.to_str().unwrap()]);
    assert!(!exists(&mixed, "both_summary.txt") && !exists(&mixed, "both_summary.html"));

    // --report-dir moves the reports; the data output stays under -o.
    let data = dir.path().join("data");
    let reports = dir.path().join("reports");
    run(&[
        "--report",
        "html,txt",
        "-o",
        data.to_str().unwrap(),
        "--report-dir",
        reports.to_str().unwrap(),
    ]);
    assert!(exists(&reports, "both_summary.txt") && exists(&reports, "both_summary.html"));
    assert!(!exists(&data, "both_summary.html"));
    assert!(exists(&data, "both_reformatted.tsv"));
}

#[test]
fn txt_summary_lists_totals_and_per_chromosome_shares() {
    let counts: IndexMap<String, usize> = [("chr1", 50), ("chr2", 30), ("chrX", 20)]
        .map(|(k, v)| (k.to_string(), v))
        .into_iter()
        .collect();
    let stats = SummaryStats {
        input_file: "test.vcf.gz".to_string(),
        output_format: "TSV".to_string(),
        transcript_handling: "FirstOnly".to_string(),
        input_variant_count: 100,
        output_record_count: 100,
        input_chrom_counts: counts.clone(),
        output_chrom_counts: counts,
        processing_time_secs: 1.5,
        variants_per_sec: 66.7,
    };
    let tmp = tempfile::NamedTempFile::new().unwrap();
    stats.write_to_file(tmp.path().to_str().unwrap()).unwrap();

    let text = std::fs::read_to_string(tmp.path()).unwrap();
    assert!(text.contains("Total input variants:     100"));
    assert!(text.contains("Total output records:     100"));
    assert!(text.contains("chr1"));
    assert!(text.contains("50.0%"));
    assert!(text.contains("1.00x")); // expansion ratio
}
