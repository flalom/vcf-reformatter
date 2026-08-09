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
    let fields: Vec<Field> = MafRecord::get_maf_headers()
        .iter()
        .map(|h| match h.as_str() {
            "Start_Position" | "End_Position" | "t_depth" | "t_ref_count" | "t_alt_count" => {
                Field::new(h, DataType::UInt64, true)
            }
            "QUAL" | "VAF" => Field::new(h, DataType::Float64, true),
            _ => Field::new(h, DataType::Utf8, true),
        })
        .collect();
    Schema::new(fields)
}

fn build_maf_batch(
    schema: &Schema,
    records: &[MafRecord],
) -> Result<RecordBatch, Box<dyn std::error::Error>> {
    let headers = MafRecord::get_maf_headers();
    let lines: Vec<String> = records.iter().map(|r| r.to_tsv_line()).collect();
    let rows: Vec<Vec<&str>> = lines.iter().map(|l| l.split('\t').collect()).collect();

    let columns: Vec<ArrayRef> = headers
        .iter()
        .enumerate()
        .map(|(col_idx, header)| -> ArrayRef {
            match header.as_str() {
                "Start_Position" => Arc::new(
                    records.iter().map(|r| r.start_position).collect::<arrow::array::UInt64Array>(),
                ),
                "End_Position" => Arc::new(
                    records.iter().map(|r| r.end_position).collect::<arrow::array::UInt64Array>(),
                ),
                "t_depth" => Arc::new(
                    records.iter().map(|r| r.t_depth.map(|v| v as u64)).collect::<arrow::array::UInt64Array>(),
                ),
                "t_ref_count" => Arc::new(
                    records.iter().map(|r| r.t_ref_count.map(|v| v as u64)).collect::<arrow::array::UInt64Array>(),
                ),
                "t_alt_count" => Arc::new(
                    records.iter().map(|r| r.t_alt_count.map(|v| v as u64)).collect::<arrow::array::UInt64Array>(),
                ),
                "QUAL" => Arc::new(records.iter().map(|r| r.qual).collect::<arrow::array::Float64Array>()),
                "VAF" => Arc::new(
                    records.iter().map(|r| r.vaf.map(|v| v as f64)).collect::<arrow::array::Float64Array>(),
                ),
                _ => Arc::new(
                    rows.iter()
                        .map(|row| Some(row[col_idx]))
                        .collect::<StringArray>(),
                ),
            }
        })
        .collect();

    Ok(RecordBatch::try_new(Arc::new(schema.clone()), columns)?)
}