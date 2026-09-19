//! MafRecord: VCF record to MAF row. Column-level behaviour (HGVS, depths, dbSNP, alleles) is
//! unit-tested next to the code in `src/essentials_fields.rs`; these cover the conversion as a
//! whole and the writer.
mod common;

use std::collections::HashMap;

use common::{info, maf, maf_multi, record};
use vcf_reformatter::essentials_fields::MafRecord;
use vcf_reformatter::reformat_vcf::{AnnotationFieldType, MafWriter};

#[test]
fn snp_with_vep_annotation_fills_the_core_columns() {
    let r = record(
        "chr17",
        7674220,
        "G",
        "A",
        Some(60.0),
        "PASS",
        info(&[
            ("CSQ_SYMBOL", "TP53"),
            ("CSQ_Consequence", "missense_variant"),
            ("CSQ_Feature", "ENST00000269305"),
            ("CSQ_Protein_position", "175"),
            ("CSQ_HGVSp", "ENSP00000269305.4:p.Arg175His"),
            ("CSQ_HGVSc", "ENST00000269305.8:c.524G>A"),
            ("INFO_DP", "100"),
        ]),
    );
    let m = maf(&r);
    assert_eq!(m.hugo_symbol, "TP53");
    assert_eq!(m.chromosome, "chr17"); // VCF naming passes through unchanged
    assert_eq!((m.start_position, m.end_position), (7674220, 7674220));
    assert_eq!(m.variant_type, "SNP");
    assert_eq!(m.variant_classification, "Missense_Mutation");
    assert_eq!(m.reference_allele, "G");
    assert_eq!(m.tumor_seq_allele2, "A");
    assert_eq!(m.center, "TEST_CENTER");
    assert_eq!(m.ncbi_build, "GRCh38");
    assert_eq!(m.tumor_sample_barcode, "TEST_SAMPLE");
    assert_eq!(m.qual, Some(60.0));
    assert_eq!(m.filter_status, "PASS");
    assert_eq!(m.transcript_id.as_deref(), Some("ENST00000269305"));
    assert_eq!(m.protein_position.as_deref(), Some("175"));
    // The accession is already in Transcript_ID, so it is stripped from the HGVS strings.
    assert_eq!(m.hgvsp.as_deref(), Some("p.Arg175His"));
    assert_eq!(m.hgvsc.as_deref(), Some("c.524G>A"));
    assert_eq!(m.hgvsp_short.as_deref(), Some("p.R175H"));
}

#[test]
fn insertion_and_deletion_positions_and_alleles() {
    // The shared anchor base is stripped; the MAF spells the missing side as "-".
    let ins = maf(&record(
        "chr1",
        1000,
        "A",
        "ATCG",
        Some(45.0),
        "PASS",
        HashMap::new(),
    ));
    assert_eq!(ins.variant_type, "INS");
    assert_eq!((ins.start_position, ins.end_position), (1000, 1001));
    assert_eq!(
        (
            ins.reference_allele.as_str(),
            ins.tumor_seq_allele2.as_str()
        ),
        ("-", "TCG")
    );

    let del = maf(&record(
        "chr2",
        2000,
        "ATCG",
        "A",
        Some(55.0),
        "PASS",
        HashMap::new(),
    ));
    assert_eq!(del.variant_type, "DEL");
    assert_eq!((del.start_position, del.end_position), (2001, 2003));
    assert_eq!(
        (
            del.reference_allele.as_str(),
            del.tumor_seq_allele2.as_str()
        ),
        ("TCG", "-")
    );

    // Long indels: End_Position spans the whole deleted run.
    let long_del = "A".to_string() + &"G".repeat(50);
    let m = maf(&record(
        "chr1",
        2000,
        &long_del,
        "A",
        Some(40.0),
        "PASS",
        HashMap::new(),
    ));
    assert_eq!(m.reference_allele, "G".repeat(50));
    assert_eq!((m.start_position, m.end_position), (2001, 2050));
    let long_ins = "A".to_string() + &"T".repeat(100);
    let m = maf(&record(
        "chr1",
        1000,
        "A",
        &long_ins,
        Some(30.0),
        "PASS",
        HashMap::new(),
    ));
    assert_eq!(m.tumor_seq_allele2, "T".repeat(100));
}

#[test]
fn vep_consequence_terms_map_to_maf_classes() {
    for (term, expected) in [
        ("missense_variant", "Missense_Mutation"),
        ("synonymous_variant", "Silent"),
        ("stop_gained", "Nonsense_Mutation"),
        ("intron_variant", "Intron"),
        ("5_prime_UTR_variant", "5'UTR"),
        ("3_prime_UTR_variant", "3'UTR"),
    ] {
        let m = maf(&record(
            "chr1",
            1000,
            "G",
            "A",
            Some(60.0),
            "PASS",
            info(&[("CSQ_Consequence", term)]),
        ));
        assert_eq!(m.variant_classification, expected, "for {term}");
    }
}

#[test]
fn snpeff_consequence_terms_map_to_maf_classes() {
    for (term, expected) in [
        ("disruptive_inframe_insertion", "In_Frame_Ins"),
        ("conservative_inframe_insertion", "In_Frame_Ins"),
        ("disruptive_inframe_deletion", "In_Frame_Del"),
        ("conservative_inframe_deletion", "In_Frame_Del"),
        ("initiator_codon_variant", "Translation_Start_Site"),
        ("rare_amino_acid_variant", "Missense_Mutation"),
        ("stop_retained_variant", "Silent"),
        ("upstream_gene_variant", "5'Flank"),
        ("downstream_gene_variant", "3'Flank"),
        ("non_coding_transcript_exon_variant", "RNA"),
        ("non_coding_transcript_variant", "RNA"),
        // case-insensitive
        ("STOP_GAINED", "Nonsense_Mutation"),
        ("Missense_Variant", "Missense_Mutation"),
        ("SPLICE_DONOR_VARIANT", "Splice_Site"),
    ] {
        let m = maf(&record(
            "chr1",
            1000,
            "A",
            "T",
            Some(30.0),
            "PASS",
            info(&[("ANN_Annotation", term)]),
        ));
        assert_eq!(m.variant_classification, expected, "for {term}");
    }
}

#[test]
fn frameshift_class_follows_the_indel_direction() {
    let fs = info(&[("CSQ_Consequence", "frameshift_variant")]);
    let ins = maf(&record(
        "chr1",
        1000,
        "A",
        "ATCG",
        Some(50.0),
        "PASS",
        fs.clone(),
    ));
    assert_eq!(
        (
            ins.variant_classification.as_str(),
            ins.variant_type.as_str()
        ),
        ("Frame_Shift_Ins", "INS")
    );
    let del = maf(&record(
        "chr1",
        1000,
        "ATCG",
        "A",
        Some(50.0),
        "PASS",
        fs.clone(),
    ));
    assert_eq!(
        (
            del.variant_classification.as_str(),
            del.variant_type.as_str()
        ),
        ("Frame_Shift_Del", "DEL")
    );
    // A frameshift tag on a same-length substitution is not an indel; vcf2maf's
    // GetVariantClassification falls through to its catch-all here too.
    let snp = maf(&record("chr1", 1000, "AT", "GC", Some(50.0), "PASS", fs));
    assert_eq!(snp.variant_classification, "Targeted_Region");
}

#[test]
fn impact_alone_is_a_fallback_not_a_consequence_term() {
    let m = maf(&record(
        "chr1",
        1000,
        "A",
        "T",
        Some(30.0),
        "PASS",
        info(&[("ANN_Annotation_Impact", "HIGH")]),
    ));
    assert_eq!(m.variant_classification, "Missense_Mutation");
}

#[test]
fn no_annotation_at_all_falls_back_like_vcf2maf() {
    let m = maf(&record("chr22", 42000, "C", "G", None, ".", HashMap::new()));
    assert_eq!(m.hugo_symbol, "Unknown");
    assert_eq!(m.qual, None);
    assert_eq!(m.filter_status, ".");
    assert_eq!(m.transcript_id, None);
    assert_eq!(m.protein_position, None);
    assert_eq!(m.hgvsp, None);
    assert_eq!(m.hgvsc, None);
}

#[test]
fn snpeff_ann_fields_feed_the_same_columns_as_vep_csq() {
    let mut r = record(
        "chr17",
        7675088,
        "C",
        "T",
        Some(100.0),
        "PASS",
        info(&[
            ("ANN_Gene_Name", "TP53"),
            ("ANN_Annotation", "missense_variant"),
            ("ANN_Feature_ID", "NM_000546.6"),
            ("ANN_HGVS_p", "p.Arg175His"),
            ("ANN_HGVS_c", "c.524G>A"),
            ("ANN_Protein_position", "175"),
        ]),
    );
    r.annotation_field_type = AnnotationFieldType::Ann;
    let m = maf(&r);
    assert_eq!(m.hugo_symbol, "TP53");
    assert_eq!(m.variant_classification, "Missense_Mutation");
    assert_eq!(m.transcript_id.as_deref(), Some("NM_000546.6"));
    assert_eq!(m.hgvsp.as_deref(), Some("p.Arg175His"));
    assert_eq!(m.hgvsc.as_deref(), Some("c.524G>A"));
    assert_eq!(m.protein_position.as_deref(), Some("175"));
}

#[test]
fn multiallelic_record_yields_one_row_per_alt() {
    let r = record(
        "chr17",
        43094692,
        "G",
        "A,T",
        Some(80.0),
        "PASS",
        info(&[
            ("CSQ_SYMBOL", "BRCA1"),
            ("CSQ_Consequence", "missense_variant"),
        ]),
    );
    let rows = maf_multi(&r);
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0].tumor_seq_allele2, "A");
    assert_eq!(rows[1].tumor_seq_allele2, "T");
    assert_eq!(rows[0].start_position, rows[1].start_position);
    assert_eq!(rows[0].qual, rows[1].qual);
}

#[test]
fn tsv_line_has_one_field_per_header_in_order() {
    let r = record(
        "chr7",
        55191822,
        "T",
        "G",
        Some(99.9),
        "PASS",
        info(&[
            ("CSQ_SYMBOL", "EGFR"),
            ("CSQ_Feature", "ENST00000275493"),
            ("CSQ_Protein_position", "858"),
        ]),
    );
    let m = MafRecord::from_reformatted_record_for_samples(
        &r,
        "TCGA",
        "GRCh38",
        "TCGA-01-0001",
        None,
        None,
    )
    .unwrap();
    let line = m.to_tsv_line();
    let fields: Vec<&str> = line.split('\t').collect();
    let headers = MafRecord::get_maf_headers();
    assert_eq!(fields.len(), headers.len());
    let at = |name: &str| fields[headers.iter().position(|h| h == name).unwrap()];
    assert_eq!(at("Hugo_Symbol"), "EGFR");
    assert_eq!(at("Center"), "TCGA");
    assert_eq!(at("Chromosome"), "chr7");
    assert_eq!(at("Start_Position"), "55191822");
    assert_eq!(at("Transcript_ID"), "ENST00000275493");
    assert_eq!(at("FILTER"), "PASS");
    assert_eq!(at("QUAL"), "99.9");
    assert_eq!(at("Protein_Position"), "858");
}

#[test]
fn maf_writer_emits_header_then_rows() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("out.maf");
    let m = maf(&record(
        "chr17",
        7674220,
        "G",
        "A",
        Some(60.0),
        "PASS",
        info(&[
            ("CSQ_SYMBOL", "TP53"),
            ("CSQ_Consequence", "missense_variant"),
        ]),
    ));

    let mut w = MafWriter::create(path.to_str().unwrap(), false).unwrap();
    w.write_rows(&[m]).unwrap();
    w.finish().unwrap();

    let text = std::fs::read_to_string(&path).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 2);
    assert_eq!(lines[0], MafRecord::get_maf_headers().join("\t"));
    assert!(lines[1].starts_with("TP53\t"));
    assert!(lines[1].contains("\tchr17\t7674220\t"));
}
