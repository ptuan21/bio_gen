use std::collections::HashMap;

use crate::bam::BamRecord;

/// Base counts on one strand.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct StrandCount {
    pub a: u32,
    pub c: u32,
    pub g: u32,
    pub t: u32,
    pub n: u32,
}

/// Per-reference-position coverage: total depth, base breakdown, and the same
/// broken down by read strand (for strand-bias assessment).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PileupColumn {
    pub ref_pos: i32,
    pub depth: u32,
    pub a: u32,
    pub c: u32,
    pub g: u32,
    pub t: u32,
    pub n: u32,
    /// Reads carrying a deletion (`D`) at this position.
    pub del: u32,
    pub fwd: StrandCount,
    pub rev: StrandCount,
}

impl PileupColumn {
    fn add_base(&mut self, base: u8, reverse: bool) {
        self.depth += 1;
        match base.to_ascii_uppercase() {
            b'A' => { self.a += 1; if reverse { self.rev.a += 1 } else { self.fwd.a += 1 } }
            b'C' => { self.c += 1; if reverse { self.rev.c += 1 } else { self.fwd.c += 1 } }
            b'G' => { self.g += 1; if reverse { self.rev.g += 1 } else { self.fwd.g += 1 } }
            b'T' | b'U' => { self.t += 1; if reverse { self.rev.t += 1 } else { self.fwd.t += 1 } }
            _ => { self.n += 1; if reverse { self.rev.n += 1 } else { self.fwd.n += 1 } }
        }
    }

    /// Forward and reverse read support for `base` at this column.
    pub fn strand_counts(&self, base: u8) -> (u32, u32) {
        match base.to_ascii_uppercase() {
            b'A' => (self.fwd.a, self.rev.a),
            b'C' => (self.fwd.c, self.rev.c),
            b'G' => (self.fwd.g, self.rev.g),
            b'T' | b'U' => (self.fwd.t, self.rev.t),
            _ => (self.fwd.n, self.rev.n),
        }
    }

    /// The most frequent base and its count (ignoring deletions); `None` when
    /// no base is present.
    pub fn consensus(&self) -> Option<(char, u32)> {
        [('A', self.a), ('C', self.c), ('G', self.g), ('T', self.t), ('N', self.n)]
            .into_iter()
            .filter(|&(_, n)| n > 0)
            .max_by_key(|&(_, n)| n)
    }
}

fn cigar_ops(cigar: &str) -> Vec<(u32, char)> {
    if cigar == "*" {
        return Vec::new();
    }
    let mut ops = Vec::new();
    let mut len = 0u32;
    for ch in cigar.chars() {
        match ch.to_digit(10) {
            Some(d) => len = len * 10 + d,
            None => {
                ops.push((len, ch));
                len = 0;
            }
        }
    }
    ops
}

/// Build a coverage pileup over `[beg, end)` from aligned records, walking each
/// CIGAR to project read bases onto reference positions. Only covered positions
/// are returned, sorted by position. Insertions/soft-clips are skipped; `N`
/// (reference skip) does not add depth, while `D` (deletion) does.
///
/// Bases whose Phred quality is below `min_qual` are not counted (reads with no
/// quality string, `0xFF`, are always counted). Deletions carry no base quality
/// and are unaffected.
pub fn pileup(records: &[BamRecord], beg: i32, end: i32, min_qual: u8) -> Vec<PileupColumn> {
    let mut columns: HashMap<i32, PileupColumn> = HashMap::new();

    for record in records {
        if record.pos < 0 {
            continue;
        }
        let seq = record.seq.as_bytes();
        let qual = &record.qual;
        let reverse = record.flag & 0x10 != 0;
        let mut ref_pos = record.pos;
        let mut query_pos = 0usize;

        for (len, op) in cigar_ops(&record.cigar) {
            let len = len as usize;
            match op {
                'M' | '=' | 'X' => {
                    for i in 0..len {
                        let rp = ref_pos + i as i32;
                        let qpos = query_pos + i;
                        if rp >= beg
                            && rp < end
                            && qual.get(qpos).is_none_or(|&q| q >= min_qual)
                        {
                            let base = seq.get(qpos).copied().unwrap_or(b'N');
                            columns
                                .entry(rp)
                                .or_insert(PileupColumn { ref_pos: rp, ..Default::default() })
                                .add_base(base, reverse);
                        }
                    }
                    ref_pos += len as i32;
                    query_pos += len;
                }
                'I' | 'S' => query_pos += len,
                'D' => {
                    for i in 0..len {
                        let rp = ref_pos + i as i32;
                        if rp >= beg && rp < end {
                            let col = columns
                                .entry(rp)
                                .or_insert(PileupColumn { ref_pos: rp, ..Default::default() });
                            col.depth += 1;
                            col.del += 1;
                        }
                    }
                    ref_pos += len as i32;
                }
                'N' => ref_pos += len as i32,
                _ => {} // H, P consume neither reference nor query
            }
        }
    }

    let mut out: Vec<PileupColumn> = columns.into_values().collect();
    out.sort_by_key(|c| c.ref_pos);
    out
}
