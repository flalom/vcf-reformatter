//! Shared helpers for the integration tests. Each `tests/*.rs` file is its own crate, so
//! anything two of them need lives here.
#![allow(dead_code)] // not every test crate uses every helper

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;
use std::process::Command;

use vcf_reformatter::essentials_fields::MafRecord;
use vcf_reformatter::reformat_vcf::{AnnotationFieldType, ReformattedVcfRecord};

/// A minimal VCF with one variant and no sample columns.
pub const TINY_VCF: &str = "##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t60\tPASS\tDP=50\n";

/// The compiled binary under test, so CLI tests spawn it directly instead of going through
/// `cargo run` (which re-resolves the build on every call).
pub fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_vcf-reformatter"))
}

/// Write `content` to `dir/name` and return the path as a `String` for argv.
pub fn write_vcf(dir: &Path, name: &str, content: &str) -> String {
    let path = dir.join(name);
    std::fs::File::create(&path)
        .unwrap()
        .write_all(content.as_bytes())
        .unwrap();
    path.to_str().unwrap().to_string()
}

/// A sample-free `ReformattedVcfRecord` with the given coordinates and INFO/CSQ/ANN fields.
pub fn record(
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
        annotation_field_type: AnnotationFieldType::None,
    }
}

/// `info_fields` from `(key, value)` pairs.
pub fn info(pairs: &[(&str, &str)]) -> HashMap<String, String> {
    pairs
        .iter()
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect()
}

/// Convert a single-ALT record with no tumour/normal sample named.
pub fn maf(r: &ReformattedVcfRecord) -> MafRecord {
    MafRecord::from_reformatted_record_for_samples(
        r,
        "TEST_CENTER",
        "GRCh38",
        "TEST_SAMPLE",
        None,
        None,
    )
    .unwrap()
}

/// Convert a possibly multi-ALT record, one `MafRecord` per allele.
pub fn maf_multi(r: &ReformattedVcfRecord) -> Vec<MafRecord> {
    MafRecord::from_reformatted_record_multi_for_samples(
        r,
        "TEST_CENTER",
        "GRCh38",
        "TEST_SAMPLE",
        None,
        None,
    )
    .unwrap()
}
