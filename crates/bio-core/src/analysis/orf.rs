use crate::analysis::search::Strand;
use crate::analysis::translate::codon_to_aa;
use crate::sequence::Sequence;

/// An open reading frame: `ATG` start to the next in-frame stop codon.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Orf {
    /// 0-based start on the forward strand (inclusive).
    pub start: usize,
    /// End on the forward strand (exclusive, includes the stop codon).
    pub end: usize,
    pub strand: Strand,
    /// +1..+3 (forward) or -1..-3 (reverse).
    pub frame: i8,
    /// Translated protein, starting with `M`, excluding the stop.
    pub protein: String,
    /// Coding nucleotide sequence (5'->3'), including the stop codon.
    pub dna: String,
}

fn map_coords(strand: Strand, start: usize, end: usize, n: usize) -> (usize, usize) {
    match strand {
        Strand::Forward => (start, end),
        Strand::Reverse => (n - end, n - start),
    }
}

fn scan(work: &[u8], strand: Strand, min_aa: usize, out: &mut Vec<Orf>) {
    let n = work.len();
    for f in 0..3 {
        let mut protein = String::new();
        let mut in_orf = false;
        let mut orf_start = 0;
        let mut c = f;
        while c + 3 <= n {
            let codon = &work[c..c + 3];
            if !in_orf {
                if codon == b"ATG" || codon == b"AUG" {
                    in_orf = true;
                    orf_start = c;
                    protein.clear();
                    protein.push('M');
                }
            } else {
                let aa = codon_to_aa(codon) as char;
                if aa == '*' {
                    if protein.len() >= min_aa {
                        let (start, end) = map_coords(strand, orf_start, c + 3, n);
                        let dna = String::from_utf8_lossy(&work[orf_start..c + 3]).into_owned();
                        out.push(Orf {
                            start,
                            end,
                            strand,
                            frame: strand_frame(strand, f),
                            protein: std::mem::take(&mut protein),
                            dna,
                        });
                    }
                    in_orf = false;
                } else {
                    protein.push(aa);
                }
            }
            c += 3;
        }
    }
}

fn strand_frame(strand: Strand, f: usize) -> i8 {
    let base = (f as i8) + 1;
    match strand {
        Strand::Forward => base,
        Strand::Reverse => -base,
    }
}

/// Find complete ORFs on both strands, keeping those with at least `min_aa`
/// residues. Results are sorted by forward-strand start then end.
pub fn find_orfs(seq: &Sequence, min_aa: usize) -> Vec<Orf> {
    let mut out = Vec::new();
    scan(seq.as_bytes(), Strand::Forward, min_aa, &mut out);
    let rc = seq.reverse_complement();
    scan(rc.as_bytes(), Strand::Reverse, min_aa, &mut out);
    out.sort_by_key(|o| (o.start, o.end));
    out
}
