#![cfg(feature = "parquet_out")]

use crate::essentials_fields::MafRecord;
use crate::reformat_vcf::ReformattedVcfRecord;
use arrow::array::{ArrayRef, Float64Builder, StringArray, StringBuilder, UInt64Builder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use std::borrow::Cow;
use std::sync::Arc;

/// Write ReformattedVcfRecords to a parquet file.
/// Numeric columns (POS, QUAL) use native types; all others are strings.
pub fn write_tsv_as_parquet(
    path: &str,
    headers: &[String],
    records: &[ReformattedVcfRecord],
) -> Result<(), Box<dyn std::error::Error>> {
    let schema = build_tsv_schema(headers);
    let batch = build_tsv_batch(&schema, headers, records)?;

    let file = std::fs::File::create(path)?;
    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .build();
    let mut writer = ArrowWriter::try_new(file, Arc::new(schema), Some(props))?;
    writer.write(&batch)?;
    writer.close()?;

    Ok(())
}

fn build_tsv_schema(headers: &[String]) -> Schema {
    let fields: Vec<Field> = headers
        .iter()
        .map(|h| match h.as_str() {
            "POS" => Field::new(h, DataType::UInt64, false),
            "QUAL" => Field::new(h, DataType::Float64, true),
            _ => Field::new(h, DataType::Utf8, true),
        })
        .collect();
    Schema::new(fields)
}

fn build_tsv_batch(
    schema: &Schema,
    headers: &[String],
    records: &[ReformattedVcfRecord],
) -> Result<RecordBatch, Box<dyn std::error::Error>> {
    let mut columns: Vec<ArrayRef> = Vec::with_capacity(headers.len());

    for header in headers {
        match header.as_str() {
            "POS" => {
                let mut builder = UInt64Builder::with_capacity(records.len());
                for r in records {
                    builder.append_value(r.position);
                }
                columns.push(Arc::new(builder.finish()));
            }
            "QUAL" => {
                let mut builder = Float64Builder::with_capacity(records.len());
                for r in records {
                    match r.quality {
                        Some(q) => builder.append_value(q),
                        None => builder.append_null(),
                    }
                }
                columns.push(Arc::new(builder.finish()));
            }
            _ => {
                let mut builder = StringBuilder::with_capacity(records.len(), records.len() * 8);
                for r in records {
                    let value = get_tsv_string_value(r, header);
                    builder.append_value(value.as_ref());
                }
                columns.push(Arc::new(builder.finish()));
            }
        }
    }

    Ok(RecordBatch::try_new(Arc::new(schema.clone()), columns)?)
}

fn get_tsv_string_value<'a>(record: &'a ReformattedVcfRecord, header: &str) -> Cow<'a, str> {
    match header {
        "CHROM" => Cow::Borrowed(record.chromosome.as_str()),
        "ID" => Cow::Borrowed(record.id.as_deref().unwrap_or(".")),
        "REF" => Cow::Borrowed(record.reference.as_str()),
        "ALT" => Cow::Borrowed(record.alternate.as_str()),
        "FILTER" => Cow::Borrowed(record.filter.as_str()),
        _ => match record.info_fields.get(header) {
            Some(v) => Cow::Borrowed(v.as_str()),
            None => {
                if let Some(ref sd) = record.format_sample_data {
                    for sample in &sd.samples {
                        for fk in &sd.format_keys {
                            if format!("{}_{}", sample.sample_name, fk) == header {
                                return match sample.format_fields.get(fk) {
                                    Some(v) => Cow::Borrowed(v.as_str()),
                                    None => Cow::Borrowed("."),
                                };
                            }
                        }
                    }
                }
                Cow::Borrowed(".")
            }
        },
    }
}

/// Write MafRecords to a parquet file.
/// Numeric columns use native types; all others are strings.
pub fn write_maf_as_parquet(
    path: &str,
    records: &[MafRecord],
) -> Result<(), Box<dyn std::error::Error>> {
    let schema = build_maf_schema();
    let batch = build_maf_batch(&schema, records)?;

    let file = std::fs::File::create(path)?;
    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .build();
    let mut writer = ArrowWriter::try_new(file, Arc::new(schema), Some(props))?;
    writer.write(&batch)?;
    writer.close()?;

    Ok(())
}

fn build_maf_schema() -> Schema {
    Schema::new(vec![
        Field::new("Hugo_Symbol", DataType::Utf8, false),
        Field::new("Entrez_Gene_Id", DataType::Utf8, true),
        Field::new("Center", DataType::Utf8, false),
        Field::new("NCBI_Build", DataType::Utf8, false),
        Field::new("Chromosome", DataType::Utf8, false),
        Field::new("Start_Position", DataType::UInt64, false),
        Field::new("End_Position", DataType::UInt64, false),
        Field::new("Strand", DataType::Utf8, false),
        Field::new("Variant_Classification", DataType::Utf8, false),
        Field::new("Variant_Type", DataType::Utf8, false),
        Field::new("Reference_Allele", DataType::Utf8, false),
        Field::new("Tumor_Seq_Allele1", DataType::Utf8, false),
        Field::new("Tumor_Seq_Allele2", DataType::Utf8, false),
        Field::new("dbSNP_RS", DataType::Utf8, true),
        Field::new("dbSNP_Val_Status", DataType::Utf8, true),
        Field::new("Tumor_Sample_Barcode", DataType::Utf8, false),
        Field::new("Matched_Norm_Sample_Barcode", DataType::Utf8, true),
        Field::new("Mutation_Status", DataType::Utf8, false),
        Field::new("Validation_Status", DataType::Utf8, true),
        Field::new("Sequencer", DataType::Utf8, true),
        Field::new("Sequence_Source", DataType::Utf8, false),
        Field::new("t_depth", DataType::UInt64, true),
        Field::new("total_depth", DataType::UInt64, true),
        Field::new("VAF", DataType::Float64, true),
        Field::new("HGVSp", DataType::Utf8, true),
        Field::new("HGVSc", DataType::Utf8, true),
        Field::new("QUAL", DataType::Float64, true),
        Field::new("FILTER", DataType::Utf8, false),
        Field::new("Transcript_ID", DataType::Utf8, true),
        Field::new("Protein_Position", DataType::Utf8, true),
    ])
}

fn build_maf_batch(
    schema: &Schema,
    records: &[MafRecord],
) -> Result<RecordBatch, Box<dyn std::error::Error>> {
    let len = records.len();
    let dot = ".";

    // String columns
    let hugo: StringArray = records.iter().map(|r| Some(r.hugo_symbol.as_str())).collect();
    let entrez_strings: Vec<String> = records
        .iter()
        .map(|r| {
            r.entrez_gene_id
                .map_or(".".to_string(), |id| id.to_string())
        })
        .collect();
    let entrez: StringArray = entrez_strings.iter().map(|s| Some(s.as_str())).collect();

    let center: StringArray = records.iter().map(|r| Some(r.center.as_str())).collect();
    let ncbi: StringArray = records.iter().map(|r| Some(r.ncbi_build.as_str())).collect();
    let chrom: StringArray = records.iter().map(|r| Some(r.chromosome.as_str())).collect();
    let strand: StringArray = records.iter().map(|r| Some(r.strand.as_str())).collect();
    let var_class: StringArray = records
        .iter()
        .map(|r| Some(r.variant_classification.as_str()))
        .collect();
    let var_type: StringArray = records
        .iter()
        .map(|r| Some(r.variant_type.as_str()))
        .collect();
    let ref_allele: StringArray = records
        .iter()
        .map(|r| Some(r.reference_allele.as_str()))
        .collect();
    let tsa1: StringArray = records
        .iter()
        .map(|r| Some(r.tumor_seq_allele1.as_str()))
        .collect();
    let tsa2: StringArray = records
        .iter()
        .map(|r| Some(r.tumor_seq_allele2.as_str()))
        .collect();
    let dbsnp: StringArray = records
        .iter()
        .map(|r| r.dbsnp_rs.as_deref().or(Some(dot)))
        .collect();
    let dbsnp_val: StringArray = records
        .iter()
        .map(|r| r.dbsnp_val_status.as_deref().or(Some(dot)))
        .collect();
    let barcode: StringArray = records
        .iter()
        .map(|r| Some(r.tumor_sample_barcode.as_str()))
        .collect();
    let matched_norm: StringArray = records
        .iter()
        .map(|r| r.matched_norm_sample_barcode.as_deref().or(Some(dot)))
        .collect();
    let mutation_status: StringArray = records
        .iter()
        .map(|r| Some(r.mutation_status.as_str()))
        .collect();
    let validation: StringArray = records
        .iter()
        .map(|r| r.validation_status.as_deref().or(Some(dot)))
        .collect();
    let sequencer: StringArray = records
        .iter()
        .map(|r| r.sequencer.as_deref().or(Some(dot)))
        .collect();
    let seq_source: StringArray = records
        .iter()
        .map(|r| Some(r.sequence_source.as_str()))
        .collect();
    let hgvsp: StringArray = records
        .iter()
        .map(|r| r.hgvsp.as_deref().or(Some(dot)))
        .collect();
    let hgvsc: StringArray = records
        .iter()
        .map(|r| r.hgvsc.as_deref().or(Some(dot)))
        .collect();
    let filter: StringArray = records
        .iter()
        .map(|r| Some(r.filter_status.as_str()))
        .collect();
    let transcript_id: StringArray = records
        .iter()
        .map(|r| r.transcript_id.as_deref().or(Some(dot)))
        .collect();
    let protein_pos: StringArray = records
        .iter()
        .map(|r| r.protein_position.as_deref().or(Some(dot)))
        .collect();

    // Numeric columns
    let mut start_b = UInt64Builder::with_capacity(len);
    let mut end_b = UInt64Builder::with_capacity(len);
    let mut depth_b = UInt64Builder::with_capacity(len);
    let mut total_depth_b = UInt64Builder::with_capacity(len);
    let mut vaf_b = Float64Builder::with_capacity(len);
    let mut qual_b = Float64Builder::with_capacity(len);

    for r in records {
        start_b.append_value(r.start_position);
        end_b.append_value(r.end_position);
        match r.depth {
            Some(d) => depth_b.append_value(d as u64),
            None => depth_b.append_null(),
        }
        match r.total_depth {
            Some(d) => total_depth_b.append_value(d as u64),
            None => total_depth_b.append_null(),
        }
        match r.vaf {
            Some(v) => vaf_b.append_value(v as f64),
            None => vaf_b.append_null(),
        }
        match r.qual {
            Some(q) => qual_b.append_value(q),
            None => qual_b.append_null(),
        }
    }

    let columns: Vec<ArrayRef> = vec![
        Arc::new(hugo),
        Arc::new(entrez),
        Arc::new(center),
        Arc::new(ncbi),
        Arc::new(chrom),
        Arc::new(start_b.finish()),
        Arc::new(end_b.finish()),
        Arc::new(strand),
        Arc::new(var_class),
        Arc::new(var_type),
        Arc::new(ref_allele),
        Arc::new(tsa1),
        Arc::new(tsa2),
        Arc::new(dbsnp),
        Arc::new(dbsnp_val),
        Arc::new(barcode),
        Arc::new(matched_norm),
        Arc::new(mutation_status),
        Arc::new(validation),
        Arc::new(sequencer),
        Arc::new(seq_source),
        Arc::new(depth_b.finish()),
        Arc::new(total_depth_b.finish()),
        Arc::new(vaf_b.finish()),
        Arc::new(hgvsp),
        Arc::new(hgvsc),
        Arc::new(qual_b.finish()),
        Arc::new(filter),
        Arc::new(transcript_id),
        Arc::new(protein_pos),
    ];

    Ok(RecordBatch::try_new(Arc::new(schema.clone()), columns)?)
}