use crate::summary::{DamageBreakdown, SummaryStats};
use std::fmt::Write as FmtWrite;

const CARD_COLORS: [&str; 4] = ["#2E7D46", "#2563A6", "#C9820A", "#B5342E"];
const SEVERITY_COLORS: [&str; 4] = ["#B5342E", "#C9820A", "#2563A6", "#2E7D46"];
const FALLBACK_COLOR: &str = "#6B7280";

const CSS: &str = r#"
:root {
  color-scheme: light dark;
  --bg: #F7F8FA;
  --ink: #1B2430;
  --card-bg: #FFFFFF;
  --border: #E2E5EA;
  --mono: ui-monospace, "SF Mono", "Cascadia Code", Menlo, monospace;
  --sans: ui-sans-serif, system-ui, -apple-system, sans-serif;
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #12151A;
    --ink: #E8EAED;
    --card-bg: #1B1F26;
    --border: #2A2F38;
  }
}
* { box-sizing: border-box; }
body {
  margin: 0;
  padding: 2.5rem 1.5rem;
  background: var(--bg);
  color: var(--ink);
  font-family: var(--sans);
}
.report { max-width: 760px; margin: 0 auto; }
.eyebrow {
  font-family: var(--mono);
  font-size: 0.75rem;
  letter-spacing: 0.08em;
  opacity: 0.6;
  text-transform: uppercase;
}
h1 { font-family: var(--mono); font-size: 1.5rem; margin: 0.25rem 0; }
.meta { font-family: var(--mono); font-size: 0.8rem; opacity: 0.7; margin-bottom: 1.5rem; }
.cards {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
  gap: 1rem;
  margin-bottom: 2rem;
}
.card {
  background: var(--card-bg);
  border: 1px solid var(--border);
  border-top: 3px solid var(--accent);
  border-radius: 6px;
  padding: 1rem;
}
.card .value {
  font-family: var(--mono);
  font-size: 1.5rem;
  font-weight: 600;
  color: var(--accent);
  font-variant-numeric: tabular-nums;
}
.card .label { font-size: 0.75rem; opacity: 0.7; margin-top: 0.25rem; }
h2 { font-size: 1rem; font-family: var(--mono); margin: 2rem 0 0.75rem; }
table.chrom-table { width: 100%; border-collapse: collapse; font-size: 0.85rem; }
table.chrom-table th, table.chrom-table td {
  text-align: right;
  padding: 0.35rem 0.5rem;
  border-bottom: 1px solid var(--border);
  font-variant-numeric: tabular-nums;
}
table.chrom-table th:first-child, table.chrom-table td:first-child { text-align: left; }
td.chrom { font-family: var(--mono); }
.damage-controls { margin-bottom: 0.75rem; font-size: 0.85rem; }
.damage-controls select { font-family: var(--sans); margin-left: 0.5rem; }
#damage-chart svg { width: 100%; height: auto; }
"#;

const DAMAGE_CHART_JS: &str = r#"
    function escapeHtml(s) {
      return String(s)
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;')
        .replace(/'/g, '&#39;');
    }
    function renderDamageChart(metric) {
      const data = DAMAGE_DATA[metric];
      const container = document.getElementById('damage-chart');
      const width = 480;
      const barHeight = 22;
      const gap = 8;
      let svg = '<svg viewBox="0 0 ' + width + ' ' + ((barHeight + gap) * data.chroms.length) + '" xmlns="http://www.w3.org/2000/svg">';
      data.chroms.forEach((row, i) => {
        const total = row.counts.reduce((a, b) => a + b, 0) || 1;
        let x = 90;
        const y = i * (barHeight + gap);
        svg += '<text x="0" y="' + (y + barHeight / 2 + 4) + '" font-size="12" font-family="ui-monospace, monospace">' + escapeHtml(row.chrom) + '</text>';
        row.counts.forEach((count, ci) => {
          const segWidth = (count / total) * (width - 90);
          if (segWidth > 0) {
            svg += '<rect x="' + x + '" y="' + y + '" width="' + segWidth + '" height="' + barHeight + '" fill="' + data.colors[ci] + '"><title>' + escapeHtml(data.categories[ci]) + ': ' + count + '</title></rect>';
            x += segWidth;
          }
        });
      });
      svg += '</svg>';
      container.innerHTML = svg;
    }
"#;

fn escape_html(input: &str) -> String {
    input
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}

fn escape_js_string(input: &str) -> String {
    // Escape `<` as < so a value containing "</script" can't break out of
    // the surrounding <script> block, even though it decodes back to `<` at
    // JS-string-literal runtime. Standard mitigation for JSON-in-<script>
    // (same approach Rails/Django use).
    input
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('<', "\\u003C")
}

fn render_chrom_rows(stats: &SummaryStats) -> String {
    let mut rows = String::new();
    let mut seen = std::collections::HashSet::new();

    for (chrom, input_count) in &stats.input_chrom_counts {
        let output_count = stats.output_chrom_counts.get(chrom).copied().unwrap_or(0);
        writeln!(
            rows,
            "      <tr><td class=\"chrom\">{}</td><td>{}</td><td>{}</td></tr>",
            escape_html(chrom),
            input_count,
            output_count
        )
        .unwrap();
        seen.insert(chrom.clone());
    }
    for (chrom, output_count) in &stats.output_chrom_counts {
        if !seen.contains(chrom) {
            writeln!(
                rows,
                "      <tr><td class=\"chrom\">{}</td><td>0</td><td>{}</td></tr>",
                escape_html(chrom),
                output_count
            )
            .unwrap();
        }
    }
    rows
}

fn render_damage_json(breakdowns: &[DamageBreakdown]) -> String {
    let mut json = String::from("{");
    for (bi, b) in breakdowns.iter().enumerate() {
        if bi > 0 {
            json.push(',');
        }
        write!(json, "\"{}\":{{", escape_js_string(&b.metric_name)).unwrap();

        json.push_str("\"categories\":[");
        for (i, cat) in b.categories.iter().enumerate() {
            if i > 0 {
                json.push(',');
            }
            write!(json, "\"{}\"", escape_js_string(cat)).unwrap();
        }
        json.push(']');

        json.push_str(",\"colors\":[");
        for i in 0..b.categories.len() {
            if i > 0 {
                json.push(',');
            }
            write!(
                json,
                "\"{}\"",
                SEVERITY_COLORS.get(i).copied().unwrap_or(FALLBACK_COLOR)
            )
            .unwrap();
        }
        json.push(']');

        json.push_str(",\"chroms\":[");
        for (ci, (chrom, cat_counts)) in b.per_chrom_counts.iter().enumerate() {
            if ci > 0 {
                json.push(',');
            }
            write!(json, "{{\"chrom\":\"{}\",\"counts\":[", escape_js_string(chrom)).unwrap();
            for (i, cat) in b.categories.iter().enumerate() {
                if i > 0 {
                    json.push(',');
                }
                write!(json, "{}", cat_counts.get(cat).copied().unwrap_or(0)).unwrap();
            }
            json.push_str("]}");
        }
        json.push(']');
        json.push('}');
    }
    json.push('}');
    json
}

fn render_damage_section(breakdowns: &[DamageBreakdown]) -> String {
    if breakdowns.is_empty() {
        return String::new();
    }

    let mut section = String::new();
    section.push_str("\n  <h2>Annotation severity by chromosome</h2>\n");
    if breakdowns.len() > 1 {
        section.push_str("  <div class=\"damage-controls\">\n");
        section.push_str("    <label for=\"damage-metric-select\">Metric:</label>\n");
        section.push_str(
        "    <select id=\"damage-metric-select\" onchange=\"renderDamageChart(this.value)\">\n",
        );
        for (i, b) in breakdowns.iter().enumerate() {
            writeln!(
                section,
                "      <option value=\"{name}\"{selected}>{name}</option>",
                name = escape_html(&b.metric_name),
                selected = if i == 0 { " selected" } else { "" }
            )
            .unwrap();
        }
        section.push_str("    </select>\n");
        section.push_str("  </div>\n");
    }
    section.push_str("  <div id=\"damage-chart\"></div>\n");
    section.push_str("  <script>\n");
    writeln!(
        section,
        "    const DAMAGE_DATA = {};",
        render_damage_json(breakdowns)
    )
    .unwrap();
    section.push_str(DAMAGE_CHART_JS);
    writeln!(
        section,
        "    renderDamageChart(\"{}\");",
        escape_js_string(&breakdowns[0].metric_name)
    )
    .unwrap();
    section.push_str("  </script>\n");
    section
}

/// Render a self-contained HTML processing report. Pure function — no I/O —
/// so it can be tested directly against the returned string.
pub fn render(stats: &SummaryStats, breakdowns: &[DamageBreakdown], timestamp: &str) -> String {
    let expansion = if stats.input_variant_count > 0 {
        stats.output_record_count as f64 / stats.input_variant_count as f64
    } else {
        0.0
    };

    let mut html = String::new();
    write!(
        html,
        r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>VCF Reformatter Report</title>
<style>
{css}
</style>
</head>
<body>
<div class="report">
  <div class="eyebrow">VCF-REFORMATTER &middot; PROCESSING REPORT</div>
  <h1>{input_file} &rarr; {output_format}</h1>
  <div class="meta">{timestamp} &middot; {transcript_handling} transcripts</div>

  <div class="cards">
    <div class="card" style="--accent: {c0}">
      <div class="value">{input_count}</div>
      <div class="label">Input variants</div>
    </div>
    <div class="card" style="--accent: {c1}">
      <div class="value">{output_count}</div>
      <div class="label">Output records</div>
    </div>
    <div class="card" style="--accent: {c2}">
      <div class="value">{expansion:.2}x</div>
      <div class="label">Expansion ratio</div>
    </div>
    <div class="card" style="--accent: {c3}">
      <div class="value">{rate:.0}/s</div>
      <div class="label">Processing rate</div>
    </div>
  </div>

  <h2>Chromosome breakdown</h2>
  <table class="chrom-table">
    <thead><tr><th>Chromosome</th><th>Input</th><th>Output</th></tr></thead>
    <tbody>
{chrom_rows}    </tbody>
  </table>
{damage_section}</div>
</body>
</html>
"#,
        css = CSS,
        input_file = escape_html(&stats.input_file),
        output_format = escape_html(&stats.output_format),
        timestamp = escape_html(timestamp),
        transcript_handling = escape_html(&stats.transcript_handling),
        c0 = CARD_COLORS[0],
        c1 = CARD_COLORS[1],
        c2 = CARD_COLORS[2],
        c3 = CARD_COLORS[3],
        input_count = stats.input_variant_count,
        output_count = stats.output_record_count,
        expansion = expansion,
        rate = stats.variants_per_sec,
        chrom_rows = render_chrom_rows(stats),
        damage_section = render_damage_section(breakdowns),
    )
    .unwrap();

    html
}

/// Write the rendered HTML report to `path`, using the current local time
/// as the report's timestamp.
pub fn write_html_report(
    stats: &SummaryStats,
    breakdowns: &[DamageBreakdown],
    path: &str,
) -> std::io::Result<()> {
    let timestamp = chrono::Local::now().format("%Y-%m-%d %H:%M:%S").to_string();
    let html = render(stats, breakdowns, &timestamp);
    std::fs::write(path, html)
}

#[cfg(test)]
mod tests {
    use super::*;
    use indexmap::IndexMap;

    fn sample_stats() -> SummaryStats {
        let mut input = IndexMap::new();
        input.insert("chr1".to_string(), 75usize);
        input.insert("chr2".to_string(), 25usize);
        let output = input.clone();

        SummaryStats {
            input_file: "sample.vcf.gz".to_string(),
            output_format: "MAF".to_string(),
            transcript_handling: "first".to_string(),
            input_variant_count: 100,
            output_record_count: 100,
            input_chrom_counts: input,
            output_chrom_counts: output,
            processing_time_secs: 2.5,
            variants_per_sec: 40.0,
        }
    }

    #[test]
    fn test_render_contains_stat_values() {
        let stats = sample_stats();
        let html = render(&stats, &[], "2026-07-21 14:30:00");

        assert!(html.contains("sample.vcf.gz"));
        assert!(html.contains("MAF"));
        assert!(html.contains("chr1"));
        assert!(html.contains("chr2"));
        assert!(html.contains("2026-07-21 14:30:00"));
    }

    #[test]
    fn test_render_omits_damage_section_when_empty() {
        let stats = sample_stats();
        let html = render(&stats, &[], "2026-07-21 14:30:00");
        assert!(!html.contains("damage-metric-select"));
    }

    #[test]
    fn test_render_includes_damage_section_with_metric_names() {
        let stats = sample_stats();
        let mut chr1 = IndexMap::new();
        chr1.insert("deleterious".to_string(), 3usize);
        chr1.insert("tolerated".to_string(), 7usize);
        let mut per_chrom = IndexMap::new();
        per_chrom.insert("chr1".to_string(), chr1);

        let sift = DamageBreakdown {
            metric_name: "SIFT".to_string(),
            categories: vec!["deleterious".to_string(), "tolerated".to_string()],
            per_chrom_counts: per_chrom,
        };
        // Two breakdowns: the metric <select> dropdown is only meaningful
        // (and only rendered) when there's more than one metric to switch
        // between.
        let polyphen = DamageBreakdown {
            metric_name: "PolyPhen".to_string(),
            categories: vec!["probably_damaging".to_string()],
            per_chrom_counts: IndexMap::new(),
        };

        let html = render(&stats, &[sift, polyphen], "2026-07-21 14:30:00");
        assert!(html.contains("damage-metric-select"));
        assert!(html.contains("SIFT"));
        assert!(html.contains("PolyPhen"));
        assert!(html.contains("\"deleterious\""));
        assert!(html.contains("\"tolerated\""));
    }

    #[test]
    fn test_render_omits_metric_select_with_single_breakdown() {
        // Single-metric input (e.g. SnpEff-only data, which only has an
        // "Impact" breakdown) must still render the chart, just without a
        // pointless one-option dropdown to switch metrics that don't exist.
        let stats = sample_stats();
        let mut chr1 = IndexMap::new();
        chr1.insert("HIGH".to_string(), 3usize);
        chr1.insert("LOW".to_string(), 7usize);
        let mut per_chrom = IndexMap::new();
        per_chrom.insert("chr1".to_string(), chr1);

        let breakdown = DamageBreakdown {
            metric_name: "Impact".to_string(),
            categories: vec!["HIGH".to_string(), "LOW".to_string()],
            per_chrom_counts: per_chrom,
        };

        let html = render(&stats, &[breakdown], "2026-07-21 14:30:00");
        assert!(!html.contains("damage-metric-select"));
        assert!(html.contains("Annotation severity by chromosome"));
        assert!(html.contains("damage-chart"));
        assert!(html.contains("DAMAGE_DATA"));
        assert!(html.contains("renderDamageChart(\"Impact\")"));
        assert!(html.contains("\"HIGH\""));
        assert!(html.contains("\"LOW\""));
    }

    #[test]
    fn test_escape_html_escapes_special_characters() {
        assert_eq!(escape_html("<a>&\"'"), "&lt;a&gt;&amp;&quot;&#39;");
    }

    #[test]
    fn test_escape_js_string_neutralizes_script_close_tag() {
        let escaped = escape_js_string("</script><script>alert(1)</script>");
        assert!(!escaped.contains("</script"));
        assert!(escaped.contains("\\u003C"));
    }

    #[test]
    fn test_render_neutralizes_script_injection_in_damage_data() {
        let stats = sample_stats();
        let mut chr1 = IndexMap::new();
        chr1.insert("</script><script>alert(1)</script>".to_string(), 3usize);
        let mut per_chrom = IndexMap::new();
        per_chrom.insert(
            "</script><script>alert(2)</script>".to_string(),
            chr1,
        );

        let breakdown = DamageBreakdown {
            metric_name: "</script><script>alert(3)</script>".to_string(),
            categories: vec!["</script><script>alert(1)</script>".to_string()],
            per_chrom_counts: per_chrom,
        };

        let html = render(&stats, &[breakdown], "2026-07-21 14:30:00");

        // The report legitimately contains its own <script> / </script> tags;
        // what must never appear is an *injected* close tag coming from
        // attacker-controlled data breaking out of the JSON payload.
        let script_close_count = html.matches("</script>").count();
        assert_eq!(
            script_close_count, 1,
            "expected exactly one (the report's own) </script> tag, found {script_close_count}"
        );
        assert!(html.contains("\\u003C"));
    }

    #[test]
    fn test_damage_chart_js_escapes_chrom_and_category_before_svg_concat() {
        let stats = sample_stats();
        let mut chr1 = IndexMap::new();
        chr1.insert("SNP".to_string(), 3usize);
        let mut per_chrom = IndexMap::new();
        per_chrom.insert(
            "</svg><img src=x onerror=alert(1)>".to_string(),
            chr1,
        );

        let breakdown = DamageBreakdown {
            metric_name: "Impact".to_string(),
            categories: vec!["</svg><img src=x onerror=alert(2)>".to_string()],
            per_chrom_counts: per_chrom,
        };

        let html = render(&stats, &[breakdown], "2026-07-21 14:30:00");

        // The client-side escaping helper must be present and must be the
        // thing wrapping row.chrom / data.categories[ci] before they're
        // concatenated into the SVG markup string.
        assert!(
            html.contains("function escapeHtml(s)"),
            "expected the JS escapeHtml helper to be embedded in the report"
        );
        assert!(
            html.contains("escapeHtml(row.chrom)"),
            "row.chrom must be passed through escapeHtml before SVG concatenation"
        );
        assert!(
            html.contains("escapeHtml(data.categories[ci])"),
            "data.categories[ci] must be passed through escapeHtml before SVG concatenation"
        );
        assert!(
            !html.contains("+ row.chrom +"),
            "row.chrom must not be concatenated unescaped"
        );
        assert!(
            !html.contains("+ data.categories[ci] +"),
            "data.categories[ci] must not be concatenated unescaped"
        );
    }
}