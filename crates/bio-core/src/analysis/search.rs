use crate::sequence::{iupac_matches, Sequence};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Match {
    /// 0-based start on the forward strand.
    pub start: usize,
    /// Exclusive end on the forward strand.
    pub end: usize,
    pub strand: Strand,
}

fn find_all(haystack: &[u8], pattern: &[u8], strand: Strand, out: &mut Vec<Match>) {
    if pattern.is_empty() || pattern.len() > haystack.len() {
        return;
    }
    let last = haystack.len() - pattern.len();
    for start in 0..=last {
        let hit = pattern
            .iter()
            .zip(&haystack[start..])
            .all(|(&p, &b)| iupac_matches(p, b));
        if hit {
            out.push(Match {
                start,
                end: start + pattern.len(),
                strand,
            });
        }
    }
}

/// Search a motif (IUPAC ambiguity codes allowed) against a sequence.
/// With `both_strands`, reverse-strand hits are reported in forward
/// coordinates. Results are sorted by start position.
pub fn search(seq: &Sequence, pattern: &str, both_strands: bool) -> Vec<Match> {
    let haystack = seq.as_bytes();
    let fwd = pattern.trim().to_ascii_uppercase().into_bytes();
    let mut matches = Vec::new();
    find_all(haystack, &fwd, Strand::Forward, &mut matches);

    if both_strands && !fwd.is_empty() {
        let rc: Vec<u8> = fwd
            .iter()
            .rev()
            .map(|&b| crate::sequence::complement(seq.kind(), b))
            .collect();
        if rc != fwd {
            find_all(haystack, &rc, Strand::Reverse, &mut matches);
        }
    }

    matches.sort_by_key(|m| (m.start, m.strand == Strand::Reverse));
    matches
}
