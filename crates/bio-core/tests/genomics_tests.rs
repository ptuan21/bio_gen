use bio_core::analysis::restriction::{digest, find_by_name, find_sites, ENZYMES};
use bio_core::analysis::stats::gc_skew;
use bio_core::analysis::translate::{translate_with, GeneticCode};
use bio_core::sequence::Sequence;

#[test]
fn mito_code_reassigns_codons() {
    // TGA: stop in standard, W in vertebrate mito.
    let seq = Sequence::dna("TGA").unwrap();
    assert_eq!(translate_with(&seq, GeneticCode::Standard), "*");
    assert_eq!(translate_with(&seq, GeneticCode::VertebrateMito), "W");

    // ATA: I standard, M mito. AGA: R standard, stop mito.
    let seq = Sequence::dna("ATAAGA").unwrap();
    assert_eq!(translate_with(&seq, GeneticCode::Standard), "IR");
    assert_eq!(translate_with(&seq, GeneticCode::VertebrateMito), "M*");
}

#[test]
fn gc_skew_windows() {
    // GGGG -> (4-0)/4 = 1.0 ; CCCC -> -1.0
    let seq = Sequence::dna("GGGGCCCC").unwrap();
    let skew = gc_skew(&seq, 4, 4);
    assert_eq!(skew, vec![1.0, -1.0]);
}

#[test]
fn gc_skew_neutral_window_is_zero() {
    let seq = Sequence::dna("ATAT").unwrap();
    assert_eq!(gc_skew(&seq, 4, 4), vec![0.0]);
}

#[test]
fn gc_skew_invalid_params_empty() {
    let seq = Sequence::dna("ACGT").unwrap();
    assert!(gc_skew(&seq, 0, 1).is_empty());
    assert!(gc_skew(&seq, 5, 1).is_empty());
    assert!(gc_skew(&seq, 2, 0).is_empty());
}

#[test]
fn restriction_finds_ecori_cut() {
    let enzyme = find_by_name("ecori").unwrap();
    let seq = Sequence::dna("AAAGAATTCAAA").unwrap();
    let hits = find_sites(&seq, enzyme);
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].start, 3);
    assert_eq!(hits[0].cut, 4); // G^AATTC
}

#[test]
fn digest_sorted_by_cut() {
    // BamHI GGATCC at 0, EcoRI GAATTC at 6
    let seq = Sequence::dna("GGATCCGAATTC").unwrap();
    let hits = digest(&seq, ENZYMES);
    assert_eq!(hits.len(), 2);
    assert_eq!(hits[0].enzyme, "BamHI");
    assert_eq!(hits[1].enzyme, "EcoRI");
    assert!(hits[0].cut < hits[1].cut);
}

#[test]
fn unknown_enzyme_is_none() {
    assert!(find_by_name("NoSuchEnzyme").is_none());
}
