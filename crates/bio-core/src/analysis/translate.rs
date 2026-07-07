use crate::error::{BioError, Result};
use crate::sequence::Sequence;

/// Translate a single codon (standard genetic code) to a one-letter amino acid.
/// `U` is normalised to `T`; anything ambiguous yields `X`, stop codons `*`.
pub fn codon_to_aa(codon: &[u8]) -> u8 {
    if codon.len() != 3 {
        return b'X';
    }
    let n = |b: u8| if b == b'U' { b'T' } else { b.to_ascii_uppercase() };
    let c = [n(codon[0]), n(codon[1]), n(codon[2])];
    match &c {
        b"TTT" | b"TTC" => b'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
        b"ATT" | b"ATC" | b"ATA" => b'I',
        b"ATG" => b'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => b'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
        b"TAT" | b"TAC" => b'Y',
        b"TAA" | b"TAG" | b"TGA" => b'*',
        b"CAT" | b"CAC" => b'H',
        b"CAA" | b"CAG" => b'Q',
        b"AAT" | b"AAC" => b'N',
        b"AAA" | b"AAG" => b'K',
        b"GAT" | b"GAC" => b'D',
        b"GAA" | b"GAG" => b'E',
        b"TGT" | b"TGC" => b'C',
        b"TGG" => b'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => b'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => b'G',
        _ => b'X',
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeneticCode {
    Standard,
    /// Vertebrate mitochondrial: AGA/AGG = stop, ATA = M, TGA = W.
    VertebrateMito,
}

/// Translate a codon under the given genetic code.
pub fn codon_to_aa_with(codon: &[u8], code: GeneticCode) -> u8 {
    let aa = codon_to_aa(codon);
    match code {
        GeneticCode::Standard => aa,
        GeneticCode::VertebrateMito if codon.len() == 3 => {
            let n = |b: u8| if b == b'U' { b'T' } else { b.to_ascii_uppercase() };
            match &[n(codon[0]), n(codon[1]), n(codon[2])] {
                b"AGA" | b"AGG" => b'*',
                b"ATA" => b'M',
                b"TGA" => b'W',
                _ => aa,
            }
        }
        GeneticCode::VertebrateMito => aa,
    }
}

fn translate_bytes(bytes: &[u8]) -> String {
    bytes
        .chunks_exact(3)
        .map(|c| codon_to_aa(c) as char)
        .collect()
}

/// Translate a sequence in frame 0. A trailing partial codon is dropped.
pub fn translate(seq: &Sequence) -> String {
    translate_bytes(seq.as_bytes())
}

/// Translate in frame 0 under a specific genetic code.
pub fn translate_with(seq: &Sequence, code: GeneticCode) -> String {
    seq.as_bytes()
        .chunks_exact(3)
        .map(|c| codon_to_aa_with(c, code) as char)
        .collect()
}

/// Translate starting at `offset` (0, 1 or 2). Other offsets wrap into 0..3.
pub fn translate_frame(seq: &Sequence, offset: usize) -> String {
    let bytes = seq.as_bytes();
    let start = (offset % 3).min(bytes.len());
    translate_bytes(&bytes[start..])
}

/// One reading frame; `frame` is +1..+3 (forward) or -1..-3 (reverse strand).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FrameTranslation {
    pub frame: i8,
    pub protein: String,
}

/// All six reading frames (three forward, three on the reverse complement).
pub fn six_frames(seq: &Sequence) -> Vec<FrameTranslation> {
    let rc = seq.reverse_complement();
    let mut out = Vec::with_capacity(6);
    for f in 0..3 {
        out.push(FrameTranslation {
            frame: (f as i8) + 1,
            protein: translate_frame(seq, f),
        });
    }
    for f in 0..3 {
        out.push(FrameTranslation {
            frame: -((f as i8) + 1),
            protein: translate_frame(&rc, f),
        });
    }
    out
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MutationEffect {
    Silent { residue: char },
    Missense { from: char, to: char },
    Nonsense { from: char },
    StopLost { to: char },
}

/// Effect of a single-base substitution at `pos` (0-based, frame 0) to `alt`.
pub fn point_mutation_effect(seq: &Sequence, pos: usize, alt: u8) -> Result<MutationEffect> {
    let bytes = seq.as_bytes();
    if pos >= bytes.len() {
        return Err(BioError::OutOfBounds {
            position: pos,
            length: bytes.len(),
        });
    }
    let codon_start = (pos / 3) * 3;
    if codon_start + 3 > bytes.len() {
        return Err(BioError::OutOfBounds {
            position: codon_start + 2,
            length: bytes.len(),
        });
    }

    let mut mutated = [bytes[codon_start], bytes[codon_start + 1], bytes[codon_start + 2]];
    mutated[pos - codon_start] = alt.to_ascii_uppercase();

    let from = codon_to_aa(&bytes[codon_start..codon_start + 3]) as char;
    let to = codon_to_aa(&mutated) as char;

    Ok(match (from, to) {
        (f, t) if f == t => MutationEffect::Silent { residue: f },
        (f, '*') if f != '*' => MutationEffect::Nonsense { from: f },
        ('*', t) => MutationEffect::StopLost { to: t },
        (f, t) => MutationEffect::Missense { from: f, to: t },
    })
}
