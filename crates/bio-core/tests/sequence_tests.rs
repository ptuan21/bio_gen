use bio_core::sequence::{SeqKind, Sequence};

#[test]
fn reverse_complement_dna() {
    let seq = Sequence::dna("ATGC").unwrap();
    assert_eq!(seq.reverse_complement().to_string(), "GCAT");
}

#[test]
fn complement_preserves_length_and_order() {
    let seq = Sequence::dna("AACG").unwrap();
    assert_eq!(seq.complement().to_string(), "TTGC");
}

#[test]
fn transcribe_dna_to_rna() {
    let seq = Sequence::dna("ATGC").unwrap();
    let rna = seq.transcribe();
    assert_eq!(rna.kind(), SeqKind::Rna);
    assert_eq!(rna.to_string(), "AUGC");
}

#[test]
fn transcribe_rna_is_identity() {
    let seq = Sequence::rna("AUGC").unwrap();
    assert_eq!(seq.transcribe().to_string(), "AUGC");
}

#[test]
fn ambiguity_codes_are_accepted_and_complemented() {
    let seq = Sequence::dna("RYSWN").unwrap();
    assert_eq!(seq.reverse_complement().to_string(), "NWSRY");
}

#[test]
fn rna_rejects_thymine() {
    assert!(Sequence::rna("ACGT").is_err());
}

#[test]
fn empty_sequence_is_valid() {
    let seq = Sequence::dna("").unwrap();
    assert!(seq.is_empty());
    assert_eq!(seq.len(), 0);
}
