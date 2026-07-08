//! Minimal VCF v4.2 writer so called variants interoperate with standard tools
//! (IGV, bcftools, …). Positions are stored 0-based and emitted 1-based.

use crate::analysis::variant::Variant;

#[derive(Debug, Clone, PartialEq)]
pub struct VcfRecord {
    pub chrom: String,
    /// 0-based position; written as `pos + 1` (VCF is 1-based).
    pub pos: usize,
    pub id: String,
    pub reference: String,
    pub alternate: String,
    pub qual: Option<f64>,
    pub filter: String,
    /// Pre-formatted INFO column, e.g. `DP=30;AF=0.500`.
    pub info: String,
}

fn or_dot(s: &str) -> &str {
    if s.is_empty() {
        "."
    } else {
        s
    }
}

/// Render records as a VCF document (header + one line per record).
pub fn write_vcf(records: &[VcfRecord]) -> String {
    let mut out = String::new();
    out.push_str("##fileformat=VCFv4.2\n");
    out.push_str("##source=bio_gen\n");
    out.push_str("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">\n");
    out.push_str("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">\n");
    out.push_str("##INFO=<ID=SB,Number=2,Type=Integer,Description=\"Alt support forward,reverse\">\n");
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for r in records {
        let qual = r.qual.map_or_else(|| ".".to_string(), |q| format!("{q:.1}"));
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.chrom,
            r.pos + 1,
            or_dot(&r.id),
            r.reference,
            r.alternate,
            qual,
            or_dot(&r.filter),
            or_dot(&r.info),
        ));
    }
    out
}

/// Build VCF records from alignment-free substitution calls.
pub fn from_substitutions(chrom: &str, variants: &[Variant]) -> Vec<VcfRecord> {
    variants
        .iter()
        .map(|v| VcfRecord {
            chrom: chrom.to_string(),
            pos: v.position,
            id: String::new(),
            reference: (v.reference as char).to_string(),
            alternate: (v.alternate as char).to_string(),
            qual: None,
            filter: "PASS".to_string(),
            info: String::new(),
        })
        .collect()
}
