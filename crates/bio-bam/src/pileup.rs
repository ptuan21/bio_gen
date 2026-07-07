use std::collections::HashMap;

use crate::bam::BamRecord;

/// Per-reference-position coverage: total depth and a breakdown by base.
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
}

impl PileupColumn {
    fn add_base(&mut self, base: u8) {
        self.depth += 1;
        match base.to_ascii_uppercase() {
            b'A' => self.a += 1,
            b'C' => self.c += 1,
            b'G' => self.g += 1,
            b'T' | b'U' => self.t += 1,
            _ => self.n += 1,
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
pub fn pileup(records: &[BamRecord], beg: i32, end: i32) -> Vec<PileupColumn> {
    let mut columns: HashMap<i32, PileupColumn> = HashMap::new();

    for record in records {
        if record.pos < 0 {
            continue;
        }
        let seq = record.seq.as_bytes();
        let mut ref_pos = record.pos;
        let mut query_pos = 0usize;

        for (len, op) in cigar_ops(&record.cigar) {
            let len = len as usize;
            match op {
                'M' | '=' | 'X' => {
                    for i in 0..len {
                        let rp = ref_pos + i as i32;
                        if rp >= beg && rp < end {
                            let base = seq.get(query_pos + i).copied().unwrap_or(b'N');
                            columns
                                .entry(rp)
                                .or_insert(PileupColumn { ref_pos: rp, ..Default::default() })
                                .add_base(base);
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
