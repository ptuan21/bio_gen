use crate::error::{BioError, Result};
use crate::sequence::Sequence;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantKind {
    Substitution,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Variant {
    /// 0-based position on the reference.
    pub position: usize,
    pub reference: u8,
    pub alternate: u8,
    pub kind: VariantKind,
}

/// Call substitutions between an equal-length reference and sample.
///
/// This is alignment-free; indels require a proper aligner and are out of
/// scope for the MVP, so differing lengths are rejected.
pub fn call_substitutions(reference: &Sequence, sample: &Sequence) -> Result<Vec<Variant>> {
    if reference.len() != sample.len() {
        return Err(BioError::LengthMismatch {
            expected: reference.len(),
            found: sample.len(),
        });
    }
    let variants = reference
        .as_bytes()
        .iter()
        .zip(sample.as_bytes())
        .enumerate()
        .filter(|(_, (r, s))| r != s)
        .map(|(position, (&reference, &alternate))| Variant {
            position,
            reference,
            alternate,
            kind: VariantKind::Substitution,
        })
        .collect();
    Ok(variants)
}
