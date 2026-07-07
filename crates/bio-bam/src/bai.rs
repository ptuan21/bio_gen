use std::collections::HashMap;

use crate::error::{BamError, Result};
use crate::reader::Cursor;

/// A `[beg, end)` pair of BGZF virtual offsets pointing at BAM record bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Chunk {
    pub beg: u64,
    pub end: u64,
}

#[derive(Debug, Clone, Default)]
struct RefIndex {
    bins: HashMap<u32, Vec<Chunk>>,
    /// Linear index: smallest virtual offset touching each 16 kbp window.
    intervals: Vec<u64>,
}

/// A parsed BAI index: one entry per reference in the companion BAM.
#[derive(Debug, Clone, Default)]
pub struct Bai {
    references: Vec<RefIndex>,
}

/// The pseudo-bin BAI uses to store per-reference metadata; not a real region.
const META_BIN: u32 = 37450;
const LINEAR_SHIFT: i32 = 14;

/// The bins that can contain features overlapping `[beg, end)` (SAM spec).
fn reg2bins(beg: i32, end: i32) -> Vec<u32> {
    let mut bins = vec![0u32];
    if end <= beg {
        return bins;
    }
    let b = beg.max(0) as u32;
    let e = (end - 1) as u32;
    for (offset, shift) in [(1u32, 26u32), (9, 23), (73, 20), (585, 17), (4681, 14)] {
        for k in (offset + (b >> shift))..=(offset + (e >> shift)) {
            bins.push(k);
        }
    }
    bins
}

pub fn parse_bai(data: &[u8]) -> Result<Bai> {
    let mut c = Cursor::new(data);
    if c.take(4)? != b"BAI\x01" {
        return Err(BamError::BadBaiMagic);
    }

    let n_ref = c.u32()? as usize;
    let mut references = Vec::with_capacity(n_ref);
    for _ in 0..n_ref {
        let mut bins: HashMap<u32, Vec<Chunk>> = HashMap::new();
        let n_bin = c.u32()?;
        for _ in 0..n_bin {
            let bin = c.u32()?;
            let n_chunk = c.u32()? as usize;
            let mut chunks = Vec::with_capacity(n_chunk);
            for _ in 0..n_chunk {
                chunks.push(Chunk {
                    beg: c.u64()?,
                    end: c.u64()?,
                });
            }
            if bin != META_BIN {
                bins.insert(bin, chunks);
            }
        }
        let n_intv = c.u32()? as usize;
        let mut intervals = Vec::with_capacity(n_intv);
        for _ in 0..n_intv {
            intervals.push(c.u64()?);
        }
        references.push(RefIndex { bins, intervals });
    }
    Ok(Bai { references })
}

impl Bai {
    /// Candidate chunks for `ref_id:[beg, end)`, pruned by the linear index and
    /// sorted by start offset. Empty if the reference is absent.
    pub fn query_chunks(&self, ref_id: usize, beg: i32, end: i32) -> Vec<Chunk> {
        let Some(reference) = self.references.get(ref_id) else {
            return Vec::new();
        };
        let beg = beg.max(0);
        let min_offset = reference
            .intervals
            .get((beg >> LINEAR_SHIFT) as usize)
            .copied()
            .unwrap_or(0);

        let mut chunks: Vec<Chunk> = reg2bins(beg, end)
            .iter()
            .filter_map(|bin| reference.bins.get(bin))
            .flatten()
            .copied()
            .filter(|chunk| chunk.end > min_offset)
            .collect();
        chunks.sort_by_key(|c| c.beg);
        chunks
    }
}
