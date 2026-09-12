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

type BoxErr = Box<dyn std::error::Error>;

/// Rows per row group. The writer buffers a whole row group before flushing it, so this
/// is what bounds parquet's memory while the caller streams chunks at it.
const ROW_GROUP_ROWS: usize = 65_536;

/// An open parquet file that takes records a chunk at a time, so nothing upstream has to
/// hold the whole file. Built either for the TSV or the MAF layout; write with the
/// matching method and finish with `close`.
pub struct ParquetSink {
    writer: ArrowWriter<std::fs::File>,
    schema: Arc<Schema>,
    /// TSV column names; empty on the MAF path, whose columns are fixed.
    headers: Vec<String>,
}

impl ParquetSink {
    pub fn create_tsv(path: &str, headers: &[String]) -> Result<Self, BoxErr> {
        let schema = Arc::new(build_tsv_schema(headers));
        Ok(Self {
            writer: new_writer(path, schema.clone())?,
            schema,
            headers: headers.to_vec(),
        })
    }

    pub fn create_maf(path: &str) -> Result<Self, BoxErr> {
        let schema = Arc::new(build_maf_schema());
        Ok(Self {
            writer: new_writer(path, schema.clone())?,
            schema,
            headers: Vec::new(),
        })
    }

    pub fn write_tsv(&mut self, records: &[ReformattedVcfRecord]) -> Result<(), BoxErr> {
        if records.is_empty() {
            return Ok(());
        }
        let batch = build_tsv_batch(&self.schema, &self.headers, records)?;
        self.writer.write(&batch)?;
        Ok(())
    }

    pub fn write_maf(&mut self, records: &[MafRecord]) -> Result<(), BoxErr> {
        if records.is_empty() {
            return Ok(());
        }
        let batch = build_maf_batch(&self.schema, records)?;
        self.writer.write(&batch)?;
        Ok(())
    }

    pub fn close(self) -> Result<(), BoxErr> {
        self.writer.close()?;
        Ok(())
    }
}

fn new_writer(path: &str, schema: Arc<Schema>) -> Result<ArrowWriter<std::fs::File>, BoxErr> {
    let file = std::fs::File::create(path)?;
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .set_max_row_group_size(ROW_GROUP_ROWS)
        .build();
    Ok(ArrowWriter::try_new(file, schema, Some(props))?)
}

// One-shot wrappers: unused by the binary, which streams through ParquetSink, but they
// are the library API the integration tests use.
#[allow(dead_code)]
/// Write ReformattedVcfRecords to a parquet file in one go.
/// Numeric columns (POS, QUAL) use native types; all others are strings.
pub fn write_tsv_as_parquet(
    path: &str,
    headers: &[String],
    records: &[ReformattedVcfRecord],
) -> Result<(), BoxErr> {
    let mut sink = ParquetSink::create_tsv(path, headers)?;
    sink.write_tsv(records)?;
    sink.close()
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
    schema: &Arc<Schema>,
    headers: &[String],
    records: &[ReformattedVcfRecord],
) -> Result<RecordBatch, BoxErr> {
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

    Ok(RecordBatch::try_new(schema.clone(), columns)?)
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

#[allow(dead_code)]
/// Write MafRecords to a parquet file in one go.
/// Numeric columns use native types; all others are strings.
pub fn write_maf_as_parquet(path: &str, records: &[MafRecord]) -> Result<(), BoxErr> {
    let mut sink = ParquetSink::create_maf(path)?;
    sink.write_maf(records)?;
    sink.close()
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

fn build_maf_batch(schema: &Arc<Schema>, records: &[MafRecord]) -> Result<RecordBatch, BoxErr> {
    let headers = MafRecord::get_maf_headers();
    let lines: Vec<String> = records.iter().map(|r| r.to_tsv_line()).collect();
    let rows: Vec<Vec<&str>> = lines.iter().map(|l| l.split('\t').collect()).collect();

    let columns: Vec<ArrayRef> = headers
        .iter()
        .enumerate()
        .map(|(col_idx, header)| -> ArrayRef {
            match header.as_str() {
                "Start_Position" => Arc::new(
                    records
                        .iter()
                        .map(|r| r.start_position)
                        .collect::<arrow::array::UInt64Array>(),
                ),
                "End_Position" => Arc::new(
                    records
                        .iter()
                        .map(|r| r.end_position)
                        .collect::<arrow::array::UInt64Array>(),
                ),
                "t_depth" => Arc::new(
                    records
                        .iter()
                        .map(|r| r.t_depth.map(|v| v as u64))
                        .collect::<arrow::array::UInt64Array>(),
                ),
                "t_ref_count" => Arc::new(
                    records
                        .iter()
                        .map(|r| r.t_ref_count.map(|v| v as u64))
                        .collect::<arrow::array::UInt64Array>(),
                ),
                "t_alt_count" => Arc::new(
                    records
                        .iter()
                        .map(|r| r.t_alt_count.map(|v| v as u64))
                        .collect::<arrow::array::UInt64Array>(),
                ),
                "QUAL" => Arc::new(
                    records
                        .iter()
                        .map(|r| r.qual)
                        .collect::<arrow::array::Float64Array>(),
                ),
                "VAF" => Arc::new(
                    records
                        .iter()
                        .map(|r| r.vaf.map(|v| v as f64))
                        .collect::<arrow::array::Float64Array>(),
                ),
                // An unpopulated MAF cell is an empty string in the text output; parquet is
                // typed, so it gets a real NULL rather than an empty string.
                _ => Arc::new(
                    rows.iter()
                        .map(|row| Some(row[col_idx]).filter(|v| !v.is_empty()))
                        .collect::<StringArray>(),
                ),
            }
        })
        .collect();

    Ok(RecordBatch::try_new(schema.clone(), columns)?)
}
