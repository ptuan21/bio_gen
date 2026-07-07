//! Alphabet helpers shared by DNA and RNA, including IUPAC ambiguity codes.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqKind {
    Dna,
    Rna,
}

impl SeqKind {
    pub fn thymine(self) -> u8 {
        match self {
            SeqKind::Dna => b'T',
            SeqKind::Rna => b'U',
        }
    }
}

/// 4-bit mask over the canonical bases `A C G T` (`U` is treated as `T`).
/// Ambiguity codes expand to the union of the bases they represent, so two
/// symbols match iff their masks share a bit.
pub fn base_mask(byte: u8) -> u8 {
    match byte.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        b'R' => 0b0101,
        b'Y' => 0b1010,
        b'S' => 0b0110,
        b'W' => 0b1001,
        b'K' => 0b1100,
        b'M' => 0b0011,
        b'B' => 0b1110,
        b'D' => 0b1101,
        b'H' => 0b1011,
        b'V' => 0b0111,
        b'N' => 0b1111,
        _ => 0b0000,
    }
}

/// True when a pattern symbol can match a target base under IUPAC rules.
pub fn iupac_matches(pattern: u8, base: u8) -> bool {
    let p = base_mask(pattern);
    p != 0 && (p & base_mask(base)) != 0
}

/// A valid symbol for `kind`: canonical base, ambiguity code, `N`, or gap.
pub fn is_valid(kind: SeqKind, byte: u8) -> bool {
    match byte.to_ascii_uppercase() {
        b'U' => kind == SeqKind::Rna,
        b'T' => kind == SeqKind::Dna,
        b'A' | b'C' | b'G' | b'N' | b'-' => true,
        b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' => true,
        _ => false,
    }
}

pub fn complement(kind: SeqKind, byte: u8) -> u8 {
    let t = kind.thymine();
    match byte.to_ascii_uppercase() {
        b'A' => t,
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other,
    }
}
