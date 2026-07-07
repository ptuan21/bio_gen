use crate::sequence::Sequence;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BaseCounts {
    pub a: usize,
    pub c: usize,
    pub g: usize,
    /// `T` for DNA, `U` for RNA.
    pub t: usize,
    pub n: usize,
    pub other: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeqStats {
    pub length: usize,
    pub counts: BaseCounts,
    pub gc_content: f64,
}

/// Base composition and GC content. GC is computed over `A C G T/U` only,
/// ignoring `N`, ambiguity codes and gaps.
pub fn stats(seq: &Sequence) -> SeqStats {
    let mut counts = BaseCounts::default();
    for &b in seq.as_bytes() {
        match b {
            b'A' => counts.a += 1,
            b'C' => counts.c += 1,
            b'G' => counts.g += 1,
            b'T' | b'U' => counts.t += 1,
            b'N' => counts.n += 1,
            _ => counts.other += 1,
        }
    }
    let known = counts.a + counts.c + counts.g + counts.t;
    let gc_content = if known == 0 {
        0.0
    } else {
        (counts.g + counts.c) as f64 / known as f64
    };
    SeqStats {
        length: seq.len(),
        counts,
        gc_content,
    }
}

/// GC skew `(G - C) / (G + C)` over sliding windows of `window` bases advanced
/// by `step`. Windows with no G or C yield `0.0`. Empty when parameters don't
/// fit the sequence.
pub fn gc_skew(seq: &Sequence, window: usize, step: usize) -> Vec<f64> {
    let bytes = seq.as_bytes();
    if window == 0 || step == 0 || window > bytes.len() {
        return Vec::new();
    }
    let mut out = Vec::new();
    let mut start = 0;
    while start + window <= bytes.len() {
        let (mut g, mut c) = (0i64, 0i64);
        for &b in &bytes[start..start + window] {
            match b {
                b'G' => g += 1,
                b'C' => c += 1,
                _ => {}
            }
        }
        let denom = g + c;
        out.push(if denom == 0 {
            0.0
        } else {
            (g - c) as f64 / denom as f64
        });
        start += step;
    }
    out
}
