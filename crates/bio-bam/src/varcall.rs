use bio_core::sequence::Sequence;

use crate::pileup::PileupColumn;

/// A single-nucleotide variant called from a pileup against a reference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PileupVariant {
    pub ref_pos: i32,
    pub reference: char,
    pub alternate: char,
    pub depth: u32,
    pub alt_count: u32,
    /// `alt_count / depth`.
    pub allele_freq: f64,
}

fn normalize(base: u8) -> u8 {
    if base == b'U' {
        b'T'
    } else {
        base.to_ascii_uppercase()
    }
}

/// Call SNVs: for each column at or above `min_depth`, pick the most common
/// non-reference base and emit it when its frequency reaches `min_freq`.
///
/// `ref_offset` is the reference position of `reference[0]` (e.g. the region
/// start). Positions where the reference base is not `A/C/G/T` are skipped.
pub fn call_variants(
    columns: &[PileupColumn],
    reference: &Sequence,
    ref_offset: i32,
    min_depth: u32,
    min_freq: f64,
) -> Vec<PileupVariant> {
    let refb = reference.as_bytes();
    let mut out = Vec::new();

    for col in columns {
        if col.depth < min_depth {
            continue;
        }
        let idx = col.ref_pos - ref_offset;
        if idx < 0 || idx as usize >= refb.len() {
            continue;
        }
        let ref_base = refb[idx as usize];
        let ref_norm = normalize(ref_base);
        if !matches!(ref_norm, b'A' | b'C' | b'G' | b'T') {
            continue;
        }

        let best = [(b'A', col.a), (b'C', col.c), (b'G', col.g), (b'T', col.t)]
            .into_iter()
            .filter(|&(base, count)| base != ref_norm && count > 0)
            .max_by_key(|&(_, count)| count);

        if let Some((alt, alt_count)) = best {
            let freq = alt_count as f64 / col.depth as f64;
            if freq >= min_freq {
                out.push(PileupVariant {
                    ref_pos: col.ref_pos,
                    reference: ref_base as char,
                    alternate: alt as char,
                    depth: col.depth,
                    alt_count,
                    allele_freq: freq,
                });
            }
        }
    }
    out
}
