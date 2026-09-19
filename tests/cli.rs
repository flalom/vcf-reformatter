//! The binary as a user runs it: input routes, MAF metadata flags, stderr notes.
mod common;

use std::io::Write;
use std::process::Stdio;
use tempfile::tempdir;

use common::{bin, write_vcf, TINY_VCF};

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    enc.write_all(bytes).unwrap();
    enc.finish().unwrap()
}

#[test]
fn maf_output_from_a_gzipped_vep_file() {
    let vcf = "##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTEST-01
chr1\t14930\t.\tA\tG\t100\tPASS\tDP=50;CSQ=G|missense_variant|MODERATE|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene\tGT:DP:AD\t0/1:50:25,25
chr1\t69511\t.\tC\tT\t200\tPASS\tDP=60;CSQ=T|synonymous_variant|LOW|OR4F5|ENSG00000186092|Transcript|ENST00000335137|protein_coding\tGT:DP:AD\t1/1:60:12,48
";
    let dir = tempdir().unwrap();
    let gz = dir.path().join("in.vcf.gz");
    std::fs::write(&gz, gzip(vcf.as_bytes())).unwrap();

    let o = bin()
        .args([
            gz.to_str().unwrap(),
            "--output-format",
            "maf",
            "--center",
            "TestCenter",
        ])
        .args([
            "--sample-barcode",
            "TEST-01",
            "-o",
            dir.path().to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&o.stderr)
    );

    let maf = std::fs::read_to_string(dir.path().join("in_reformatted.maf")).unwrap();
    let lines: Vec<&str> = maf.lines().collect();
    assert_eq!(lines.len(), 3, "header + 2 rows");
    assert!(lines[1].starts_with("WASH7P\t"));
    assert!(lines[2].starts_with("OR4F5\t"));
}

#[test]
fn stdin_plain_and_gzipped_are_read_identically() {
    let vcf = "##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr2\t200\t.\tC\tT\t80\tPASS\tDP=30\n";
    let mut outputs = Vec::new();
    for payload in [vcf.as_bytes().to_vec(), gzip(vcf.as_bytes())] {
        let dir = tempdir().unwrap();
        let mut child = bin()
            .args(["-", "-o", dir.path().to_str().unwrap(), "--report", "none"])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        child.stdin.take().unwrap().write_all(&payload).unwrap();
        let o = child.wait_with_output().unwrap();
        assert!(
            o.status.success(),
            "stderr: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        outputs.push(std::fs::read_to_string(dir.path().join("stdin_reformatted.tsv")).unwrap());
    }
    assert_eq!(outputs[0], outputs[1]);
    assert!(outputs[0].contains("chr2\t200"));
    assert!(outputs[0].contains("\t30")); // DP survived the round trip
}

#[test]
fn mutation_status_and_sequence_source_flags_land_in_the_maf() {
    let dir = tempdir().unwrap();
    let vcf = write_vcf(dir.path(), "in.vcf", TINY_VCF);
    let o = bin()
        .args([
            &vcf,
            "-o",
            dir.path().to_str().unwrap(),
            "-p",
            "meta",
            "--output-format",
            "maf",
        ])
        .args(["--report", "none", "--sample-barcode", "TUMOR"])
        .args(["--mutation-status", "Somatic", "--sequence-source", "WGS"])
        .output()
        .unwrap();
    assert!(
        o.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&o.stderr)
    );

    let maf = std::fs::read_to_string(dir.path().join("meta_reformatted.maf")).unwrap();
    let mut lines = maf.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let at = |name: &str| row[header.iter().position(|h| *h == name).unwrap()];
    assert_eq!(at("Mutation_Status"), "Somatic");
    assert_eq!(at("Sequence_Source"), "WGS");
}

#[test]
fn multiallelic_note_fires_only_when_a_site_has_several_alts() {
    let dir = tempdir().unwrap();
    let multi = write_vcf(
        dir.path(),
        "multi.vcf",
        "##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG,T\t60\tPASS\tDP=50\nchr1\t200\t.\tC\tT\t60\tPASS\tDP=30\n",
    );
    let o = bin()
        .args([&multi, "-o", dir.path().to_str().unwrap()])
        .output()
        .unwrap();
    assert!(o.status.success());
    let stderr = String::from_utf8_lossy(&o.stderr);
    assert!(
        stderr.contains("1 multiallelic site"),
        "stderr was: {stderr}"
    );
    assert!(stderr.contains("bcftools norm"), "stderr was: {stderr}");

    let single = write_vcf(dir.path(), "single.vcf", TINY_VCF);
    let o = bin()
        .args([&single, "-o", dir.path().to_str().unwrap()])
        .output()
        .unwrap();
    assert!(o.status.success());
    let stderr = String::from_utf8_lossy(&o.stderr);
    assert!(!stderr.contains("multiallelic"), "stderr was: {stderr}");
}
