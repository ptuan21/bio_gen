mod nucleotide;
mod record;

pub use nucleotide::{base_mask, complement, is_valid, iupac_matches, SeqKind};
pub use record::{split_header, SeqRecord};

use crate::error::{BioError, Result};

/// A validated nucleotide sequence stored as upper-case ASCII symbols.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    kind: SeqKind,
    data: Vec<u8>,
}

impl Sequence {
    /// Validate and normalise `data` (upper-cased) as a sequence of `kind`.
    pub fn new(kind: SeqKind, mut data: Vec<u8>) -> Result<Self> {
        for (i, byte) in data.iter_mut().enumerate() {
            byte.make_ascii_uppercase();
            if !is_valid(kind, *byte) {
                return Err(BioError::InvalidNucleotide {
                    symbol: *byte as char,
                    position: i,
                });
            }
        }
        Ok(Self { kind, data })
    }

    pub fn dna(data: impl Into<Vec<u8>>) -> Result<Self> {
        Self::new(SeqKind::Dna, data.into())
    }

    pub fn rna(data: impl Into<Vec<u8>>) -> Result<Self> {
        Self::new(SeqKind::Rna, data.into())
    }

    pub fn kind(&self) -> SeqKind {
        self.kind
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.data
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn complement(&self) -> Sequence {
        let data = self
            .data
            .iter()
            .map(|&b| complement(self.kind, b))
            .collect();
        Sequence { kind: self.kind, data }
    }

    pub fn reverse_complement(&self) -> Sequence {
        let data = self
            .data
            .iter()
            .rev()
            .map(|&b| complement(self.kind, b))
            .collect();
        Sequence { kind: self.kind, data }
    }

    /// DNA → RNA (`T` → `U`). RNA sequences are returned unchanged.
    pub fn transcribe(&self) -> Sequence {
        match self.kind {
            SeqKind::Rna => self.clone(),
            SeqKind::Dna => {
                let data = self
                    .data
                    .iter()
                    .map(|&b| if b == b'T' { b'U' } else { b })
                    .collect();
                Sequence { kind: SeqKind::Rna, data }
            }
        }
    }
}

impl std::fmt::Display for Sequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // `data` is validated ASCII, so this is always valid UTF-8.
        f.write_str(std::str::from_utf8(&self.data).unwrap())
    }
}
