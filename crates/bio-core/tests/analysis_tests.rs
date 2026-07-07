use bio_core::analysis::search::{search, Strand};
use bio_core::analysis::stats::stats;
use bio_core::analysis::translate::{point_mutation_effect, translate, MutationEffect};
use bio_core::analysis::variant::{call_substitutions, VariantKind};
use bio_core::error::BioError;
use bio_core::sequence::Sequence;

#[test]
fn gc_content_and_counts() {
    let seq = Sequence::dna("GGCCATN").unwrap();
    let s = stats(&seq);
    assert_eq!(s.length, 7);
    assert_eq!(s.counts.g, 2);
    assert_eq!(s.counts.c, 2);
    assert_eq!(s.counts.n, 1);
    // GC over ACGT only: 4 / 6
    assert!((s.gc_content - 4.0 / 6.0).abs() < 1e-9);
}

#[test]
fn gc_content_of_all_n_is_zero() {
    let seq = Sequence::dna("NNNN").unwrap();
    assert_eq!(stats(&seq).gc_content, 0.0);
}

#[test]
fn exact_motif_search_forward() {
    let seq = Sequence::dna("ACGTACGTAC").unwrap();
    let hits = search(&seq, "ACGT", false);
    assert_eq!(hits.len(), 2);
    assert_eq!(hits[0].start, 0);
    assert_eq!(hits[1].start, 4);
}

#[test]
fn iupac_motif_matches_any_base() {
    let seq = Sequence::dna("AAGTAACT").unwrap();
    // "AAN" matches AAG (0) and AAC (4)
    let hits = search(&seq, "AAN", false);
    let starts: Vec<usize> = hits.iter().map(|m| m.start).collect();
    assert_eq!(starts, vec![0, 4]);
}

#[test]
fn both_strand_search_finds_reverse() {
    // reverse complement of "AATT" is "AATT" (palindrome) -> use non-palindrome
    let seq = Sequence::dna("GGGGATCCCC").unwrap();
    // "GGGG" forward at 0; its revcomp "CCCC" at position 6 counts as reverse
    let hits = search(&seq, "GGGG", true);
    assert!(hits.iter().any(|m| m.strand == Strand::Forward && m.start == 0));
    assert!(hits.iter().any(|m| m.strand == Strand::Reverse && m.start == 6));
}

#[test]
fn empty_pattern_yields_no_matches() {
    let seq = Sequence::dna("ACGT").unwrap();
    assert!(search(&seq, "", true).is_empty());
}

#[test]
fn substitution_calling() {
    let reference = Sequence::dna("ACGTACGT").unwrap();
    let sample = Sequence::dna("ACCTACGA").unwrap();
    let variants = call_substitutions(&reference, &sample).unwrap();
    assert_eq!(variants.len(), 2);
    assert_eq!(variants[0].position, 2);
    assert_eq!(variants[0].reference, b'G');
    assert_eq!(variants[0].alternate, b'C');
    assert_eq!(variants[0].kind, VariantKind::Substitution);
    assert_eq!(variants[1].position, 7);
}

#[test]
fn substitution_length_mismatch_rejected() {
    let reference = Sequence::dna("ACGT").unwrap();
    let sample = Sequence::dna("ACG").unwrap();
    assert!(matches!(
        call_substitutions(&reference, &sample),
        Err(BioError::LengthMismatch { .. })
    ));
}

#[test]
fn translation_of_start_and_stop() {
    // ATG=M, GCC=A, TAA=stop
    let seq = Sequence::dna("ATGGCCTAA").unwrap();
    assert_eq!(translate(&seq), "MA*");
}

#[test]
fn translation_drops_partial_codon() {
    let seq = Sequence::dna("ATGGC").unwrap();
    assert_eq!(translate(&seq), "M");
}

#[test]
fn point_mutation_silent() {
    // GCC (A) -> GCA (A) at position 2
    let seq = Sequence::dna("ATGGCC").unwrap();
    let effect = point_mutation_effect(&seq, 5, b'A').unwrap();
    assert_eq!(effect, MutationEffect::Silent { residue: 'A' });
}

#[test]
fn point_mutation_missense() {
    // ATG (M) -> ACG (T) at position 1
    let seq = Sequence::dna("ATGGCC").unwrap();
    let effect = point_mutation_effect(&seq, 1, b'C').unwrap();
    assert_eq!(effect, MutationEffect::Missense { from: 'M', to: 'T' });
}

#[test]
fn point_mutation_nonsense() {
    // TAC (Y) -> TAA (stop) at position 2
    let seq = Sequence::dna("TAC").unwrap();
    let effect = point_mutation_effect(&seq, 2, b'A').unwrap();
    assert_eq!(effect, MutationEffect::Nonsense { from: 'Y' });
}

#[test]
fn point_mutation_out_of_bounds() {
    let seq = Sequence::dna("ATG").unwrap();
    assert!(matches!(
        point_mutation_effect(&seq, 5, b'A'),
        Err(BioError::OutOfBounds { .. })
    ));
}

#[test]
fn point_mutation_incomplete_codon() {
    // position 4 sits in a partial codon (len 5)
    let seq = Sequence::dna("ATGGC").unwrap();
    assert!(matches!(
        point_mutation_effect(&seq, 4, b'A'),
        Err(BioError::OutOfBounds { .. })
    ));
}
