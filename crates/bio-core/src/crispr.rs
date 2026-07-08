//! CRISPR guide RNA design: PAM scanning on both strands, on-target scoring,
//! in-sequence off-target counting, and HDR donor templates.
//!
//! Scoring is a transparent heuristic (GC balance, homopolymer / poly-T
//! penalties), not a trained model. Off-target counting is limited to the
//! loaded sequence — genome-wide search needs an external index.

use crate::analysis::search::Strand;
use crate::error::{BioError, Result};
use crate::sequence::{complement, iupac_matches, SeqKind, Sequence};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Enzyme {
    pub name: &'static str,
    /// Recognition PAM (IUPAC codes allowed).
    pub pam: &'static str,
    /// `true` when the PAM sits 3' of the protospacer (Cas9); `false` = 5' (Cas12a).
    pub pam_3prime: bool,
    pub spacer_len: usize,
    /// Cut position measured into the protospacer from its PAM-proximal end.
    pub cut_into_spacer: usize,
}

pub const ENZYMES: &[Enzyme] = &[
    Enzyme { name: "SpCas9", pam: "NGG", pam_3prime: true, spacer_len: 20, cut_into_spacer: 3 },
    Enzyme { name: "SaCas9", pam: "NNGRRT", pam_3prime: true, spacer_len: 21, cut_into_spacer: 3 },
    Enzyme { name: "Cas12a", pam: "TTTV", pam_3prime: false, spacer_len: 23, cut_into_spacer: 18 },
];

pub fn enzyme_by_name(name: &str) -> Option<&'static Enzyme> {
    ENZYMES.iter().find(|e| e.name.eq_ignore_ascii_case(name))
}

#[derive(Debug, Clone, PartialEq)]
pub struct Guide {
    pub strand: Strand,
    /// Protospacer start on the forward strand (0-based, inclusive).
    pub start: usize,
    pub end: usize,
    /// Guide (spacer) sequence, 5'->3'.
    pub spacer: String,
    pub pam: String,
    /// PAM start on the forward strand (inclusive).
    pub pam_start: usize,
    /// Predicted double-strand break position, forward coordinate.
    pub cut_site: usize,
    pub gc: f64,
    /// Heuristic on-target quality, 0..=100.
    pub score: u32,
    /// Near-identical PAM-adjacent sites elsewhere in the sequence.
    pub off_targets: u32,
}

fn pam_matches(window: &[u8], pam: &[u8]) -> bool {
    window.len() == pam.len() && window.iter().zip(pam).all(|(&b, &p)| iupac_matches(p, b))
}

/// Raw protospacer hits within `work`: (start, end, cut, pam_start, spacer, pam).
fn scan(work: &[u8], enzyme: &Enzyme) -> Vec<(usize, usize, usize, usize, String, String)> {
    let pam = enzyme.pam.as_bytes();
    let (plen, slen) = (pam.len(), enzyme.spacer_len);
    let n = work.len();
    let mut out = Vec::new();
    if n < plen {
        return out;
    }
    for j in 0..=(n - plen) {
        if !pam_matches(&work[j..j + plen], pam) {
            continue;
        }
        let (pstart, pend, cut) = if enzyme.pam_3prime {
            if j < slen {
                continue;
            }
            (j - slen, j, j - enzyme.cut_into_spacer)
        } else {
            let pstart = j + plen;
            if pstart + slen > n {
                continue;
            }
            (pstart, pstart + slen, pstart + enzyme.cut_into_spacer)
        };
        let spacer = String::from_utf8_lossy(&work[pstart..pend]).into_owned();
        let pam_str = String::from_utf8_lossy(&work[j..j + plen]).into_owned();
        out.push((pstart, pend, cut, j, spacer, pam_str));
    }
    out
}

fn max_homopolymer(s: &[u8]) -> usize {
    let (mut best, mut run) = (0, 0);
    let mut prev = 0u8;
    for &b in s {
        run = if b == prev { run + 1 } else { 1 };
        prev = b;
        best = best.max(run);
    }
    best
}

fn has_poly_t(s: &[u8]) -> bool {
    s.windows(4).any(|w| w == b"TTTT")
}

/// Heuristic on-target score in 0..=100 plus the spacer GC fraction.
fn score_spacer(spacer: &[u8]) -> (u32, f64) {
    let n = spacer.len().max(1);
    let gc = spacer.iter().filter(|&&b| b == b'G' || b == b'C').count() as f64 / n as f64;
    let mut score = 100.0_f64;
    score -= (gc - 0.5).abs() * 120.0; // ideal GC ~40-60%
    if max_homopolymer(spacer) >= 4 {
        score -= 25.0;
    }
    if has_poly_t(spacer) {
        score -= 15.0; // TTTT terminates Pol III transcription
    }
    (score.clamp(0.0, 100.0) as u32, gc)
}

fn hamming(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).filter(|(x, y)| x != y).count()
}

#[allow(clippy::too_many_arguments)]
fn make_guide(strand: Strand, start: usize, end: usize, cut: usize, pam_start: usize, spacer: String, pam: String) -> Guide {
    Guide { strand, start, end, spacer, pam, pam_start, cut_site: cut, gc: 0.0, score: 0, off_targets: 0 }
}

/// Find all guides for `enzyme` on both strands, scored and annotated with the
/// count of near-identical off-target sites (Hamming distance `<= max_off_mismatch`).
/// Sorted by descending score.
pub fn find_guides(seq: &Sequence, enzyme: &Enzyme, max_off_mismatch: usize) -> Vec<Guide> {
    let fwd = seq.as_bytes();
    let n = fwd.len();
    let rc = seq.reverse_complement();

    let plen = enzyme.pam.len();
    let mut guides = Vec::new();
    for (ps, pe, cut, pam_j, spacer, pam) in scan(fwd, enzyme) {
        guides.push(make_guide(Strand::Forward, ps, pe, cut, pam_j, spacer, pam));
    }
    for (ps, pe, cut, pam_j, spacer, pam) in scan(rc.as_bytes(), enzyme) {
        guides.push(make_guide(Strand::Reverse, n - pe, n - ps, n - cut, n - (pam_j + plen), spacer, pam));
    }

    for g in &mut guides {
        let (score, gc) = score_spacer(g.spacer.as_bytes());
        g.score = score;
        g.gc = gc;
    }

    let spacers: Vec<Vec<u8>> = guides.iter().map(|g| g.spacer.clone().into_bytes()).collect();
    for (i, g) in guides.iter_mut().enumerate() {
        let q = &spacers[i];
        let matches = spacers
            .iter()
            .filter(|s| s.len() == q.len() && hamming(q, s) <= max_off_mismatch)
            .count() as u32;
        g.off_targets = matches.saturating_sub(1);
    }

    guides.sort_by(|a, b| b.score.cmp(&a.score).then(a.start.cmp(&b.start)));
    guides
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HdrTemplate {
    pub cut_site: usize,
    pub arm_len: usize,
    pub left_arm: String,
    pub insert: String,
    pub right_arm: String,
    /// left_arm + insert + right_arm.
    pub donor: String,
}

/// Design an HDR donor: homology arms of `arm_len` flanking `cut_site`, with
/// `insert` (a knock-in cassette, or empty) placed at the cut.
pub fn design_hdr(reference: &Sequence, cut_site: usize, insert: &str, arm_len: usize) -> Result<HdrTemplate> {
    let bytes = reference.as_bytes();
    if cut_site > bytes.len() {
        return Err(BioError::OutOfBounds { position: cut_site, length: bytes.len() });
    }
    let left_start = cut_site.saturating_sub(arm_len);
    let right_end = (cut_site + arm_len).min(bytes.len());
    let left_arm = String::from_utf8_lossy(&bytes[left_start..cut_site]).into_owned();
    let right_arm = String::from_utf8_lossy(&bytes[cut_site..right_end]).into_owned();
    let insert = insert.to_ascii_uppercase();
    let donor = format!("{left_arm}{insert}{right_arm}");
    Ok(HdrTemplate { cut_site, arm_len, left_arm, insert, right_arm, donor })
}

/// A single-base substitution introduced into a donor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Edit {
    /// Forward reference position of the changed base.
    pub pos: usize,
    pub from: char,
    pub to: char,
}

#[derive(Debug, Clone, PartialEq)]
pub struct KnockinDesign {
    pub template: HdrTemplate,
    /// The PAM-disrupting edit applied to the donor, if requested and possible.
    pub pam_edit: Option<Edit>,
}

fn revcomp_bytes(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| complement(SeqKind::Dna, b)).collect()
}

/// Find the smallest single-base change to the (forward) PAM bases that stops
/// the enzyme recognising the PAM on `strand`. Returns `(index_in_pam, new_base)`.
fn find_pam_disruption(pam_ref: &[u8], enzyme: &Enzyme, strand: Strand) -> Option<(usize, u8)> {
    let pat = enzyme.pam.as_bytes();
    for i in 0..pam_ref.len() {
        for &nb in b"ACGT" {
            if nb == pam_ref[i] {
                continue;
            }
            let mut trial = pam_ref.to_vec();
            trial[i] = nb;
            let still_pam = match strand {
                Strand::Forward => pam_matches(&trial, pat),
                Strand::Reverse => pam_matches(&revcomp_bytes(&trial), pat),
            };
            if !still_pam {
                return Some((i, nb));
            }
        }
    }
    None
}

/// Design a knock-in donor around a guide's cut site. With `disrupt_pam`, a
/// single base in the PAM is changed within the donor so the repaired allele is
/// no longer a Cas substrate (not guaranteed synonymous).
pub fn design_knockin(
    reference: &Sequence,
    enzyme: &Enzyme,
    guide: &Guide,
    insert: &str,
    arm_len: usize,
    disrupt_pam: bool,
) -> Result<KnockinDesign> {
    let mut template = design_hdr(reference, guide.cut_site, insert, arm_len)?;
    let mut pam_edit = None;

    if disrupt_pam {
        let refb = reference.as_bytes();
        let plen = guide.pam.len();
        if guide.pam_start + plen <= refb.len() {
            let pam_ref = &refb[guide.pam_start..guide.pam_start + plen];
            if let Some((i, nb)) = find_pam_disruption(pam_ref, enzyme, guide.strand) {
                let p = guide.pam_start + i;
                let cut = guide.cut_site;
                let left_start = cut.saturating_sub(arm_len);
                let edit = |arm: &mut String, idx: usize| -> char {
                    let mut bytes = std::mem::take(arm).into_bytes();
                    let from = bytes[idx];
                    bytes[idx] = nb;
                    *arm = String::from_utf8(bytes).unwrap();
                    from as char
                };
                if p >= left_start && p < cut {
                    let from = edit(&mut template.left_arm, p - left_start);
                    pam_edit = Some(Edit { pos: p, from, to: nb as char });
                } else if p >= cut && p < (cut + arm_len).min(refb.len()) {
                    let from = edit(&mut template.right_arm, p - cut);
                    pam_edit = Some(Edit { pos: p, from, to: nb as char });
                }
                if pam_edit.is_some() {
                    template.donor =
                        format!("{}{}{}", template.left_arm, template.insert, template.right_arm);
                }
            }
        }
    }

    Ok(KnockinDesign { template, pam_edit })
}
