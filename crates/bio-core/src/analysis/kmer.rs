use std::collections::HashMap;

use crate::sequence::Sequence;

/// Count overlapping k-mers, sorted by descending frequency then lexically.
/// Returns an empty vector when `k` is 0 or larger than the sequence.
pub fn kmer_counts(seq: &Sequence, k: usize) -> Vec<(String, usize)> {
    let bytes = seq.as_bytes();
    if k == 0 || k > bytes.len() {
        return Vec::new();
    }
    let mut counts: HashMap<&[u8], usize> = HashMap::new();
    for window in bytes.windows(k) {
        *counts.entry(window).or_insert(0) += 1;
    }
    let mut out: Vec<(String, usize)> = counts
        .into_iter()
        .map(|(kmer, n)| (String::from_utf8_lossy(kmer).into_owned(), n))
        .collect();
    out.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    out
}
