use bio_core::analysis::search::Strand;
use bio_core::crispr::{design_hdr, design_knockin, enzyme_by_name, find_guides, Edit};
use bio_core::error::BioError;
use bio_core::sequence::Sequence;

fn spcas9() -> &'static bio_core::crispr::Enzyme {
    enzyme_by_name("SpCas9").unwrap()
}

#[test]
fn finds_forward_guide_with_pam_and_cut() {
    // 20 bp protospacer immediately 5' of an NGG PAM (CGG).
    let seq = Sequence::dna(format!("{}CGG", "A".repeat(20))).unwrap();
    let guides = find_guides(&seq, spcas9(), 3);
    assert_eq!(guides.len(), 1);
    let g = &guides[0];
    assert_eq!(g.strand, Strand::Forward);
    assert_eq!(g.start, 0);
    assert_eq!(g.end, 20);
    assert_eq!(g.spacer, "A".repeat(20));
    assert_eq!(g.pam, "CGG");
    assert_eq!(g.cut_site, 17); // 3 bp 5' of the PAM
    assert_eq!(g.pam_start, 20);
    assert_eq!(g.off_targets, 0);
}

#[test]
fn knockin_disrupts_forward_pam() {
    // 20 A + CGG: PAM CGG at [20,23). Disrupting NGG changes the first G (pos 21).
    let seq = Sequence::dna(format!("{}CGG", "A".repeat(20))).unwrap();
    let enzyme = spcas9();
    let guide = find_guides(&seq, enzyme, 3).into_iter().next().unwrap();
    let d = design_knockin(&seq, enzyme, &guide, "", 6, true).unwrap();
    assert_eq!(d.pam_edit, Some(Edit { pos: 21, from: 'G', to: 'A' }));
    assert!(d.template.donor.contains("CAG")); // PAM CGG -> CAG in the donor
    assert!(!d.template.donor.contains("CGG"));
}

#[test]
fn knockin_disrupts_reverse_pam() {
    // Reverse-strand guide; forward PAM footprint is CCG at [0,3).
    let seq = Sequence::dna(format!("CCG{}", "T".repeat(20))).unwrap();
    let enzyme = spcas9();
    let guide = find_guides(&seq, enzyme, 3).into_iter().next().unwrap();
    assert_eq!(guide.strand, Strand::Reverse);
    assert_eq!(guide.pam_start, 0);
    let d = design_knockin(&seq, enzyme, &guide, "", 6, true).unwrap();
    let edit = d.pam_edit.unwrap();
    assert_eq!(edit.pos, 0); // a forward C of the CC that forms the reverse GG
    assert_ne!(edit.to, edit.from);
}

#[test]
fn knockin_without_disruption_is_plain_hdr() {
    let seq = Sequence::dna(format!("{}CGG", "A".repeat(20))).unwrap();
    let enzyme = spcas9();
    let guide = find_guides(&seq, enzyme, 3).into_iter().next().unwrap();
    let d = design_knockin(&seq, enzyme, &guide, "GG", 6, false).unwrap();
    assert!(d.pam_edit.is_none());
    let hdr = design_hdr(&seq, guide.cut_site, "GG", 6).unwrap();
    assert_eq!(d.template.donor, hdr.donor);
}

#[test]
fn finds_reverse_strand_guide_mapped_to_forward_coords() {
    // Reverse complement of the forward case: "CCG" + 20 T. The guide lives on
    // the reverse strand and its footprint maps to forward positions [3, 23).
    let seq = Sequence::dna(format!("CCG{}", "T".repeat(20))).unwrap();
    let guides = find_guides(&seq, spcas9(), 3);
    assert_eq!(guides.len(), 1);
    let g = &guides[0];
    assert_eq!(g.strand, Strand::Reverse);
    assert_eq!((g.start, g.end), (3, 23));
    assert_eq!(g.spacer, "A".repeat(20));
}

#[test]
fn scoring_penalises_poly_a_rewards_balanced() {
    let poly = Sequence::dna(format!("{}CGG", "A".repeat(20))).unwrap();
    let low = find_guides(&poly, spcas9(), 3)[0].score;

    // Balanced 50% GC spacer, no homopolymer runs.
    let balanced = Sequence::dna(format!("{}TGG", "ACGT".repeat(5))).unwrap();
    let high = find_guides(&balanced, spcas9(), 3)[0].score;

    assert!(high > low, "balanced {high} should beat poly-A {low}");
    assert_eq!(high, 100);
}

#[test]
fn off_target_counts_duplicate_sites() {
    // Two identical protospacer+PAM units -> each has one off-target (the other).
    let unit = format!("{}TGG", "ACGT".repeat(5));
    let seq = Sequence::dna(format!("{unit}CATCAT{unit}")).unwrap();
    let guides = find_guides(&seq, spcas9(), 3);
    let block = "ACGT".repeat(5);
    let dupes: Vec<_> = guides
        .iter()
        .filter(|g| g.strand == Strand::Forward && g.spacer == block)
        .collect();
    assert_eq!(dupes.len(), 2);
    assert!(dupes.iter().all(|g| g.off_targets == 1));
}

#[test]
fn cas12a_uses_5prime_pam() {
    let enzyme = enzyme_by_name("Cas12a").unwrap();
    // TTTV PAM (TTTA) followed by a 23 nt protospacer.
    let seq = Sequence::dna(format!("TTTA{}", "ACGTACGTACGTACGTACGTACG")).unwrap();
    let guides = find_guides(&seq, enzyme, 3);
    let fwd: Vec<_> = guides.iter().filter(|g| g.strand == Strand::Forward).collect();
    assert_eq!(fwd.len(), 1);
    assert_eq!(fwd[0].pam, "TTTA");
    assert_eq!(fwd[0].start, 4);
    assert_eq!(fwd[0].spacer.len(), 23);
}

#[test]
fn hdr_template_builds_donor() {
    let reference = Sequence::dna("AAAACCCCGGGGTTTT").unwrap();
    let hdr = design_hdr(&reference, 8, "gaga", 4).unwrap();
    assert_eq!(hdr.left_arm, "CCCC");
    assert_eq!(hdr.right_arm, "GGGG");
    assert_eq!(hdr.insert, "GAGA");
    assert_eq!(hdr.donor, "CCCCGAGAGGGG");
}

#[test]
fn hdr_out_of_bounds_errors() {
    let reference = Sequence::dna("ACGT").unwrap();
    assert!(matches!(
        design_hdr(&reference, 99, "A", 4),
        Err(BioError::OutOfBounds { .. })
    ));
}

#[test]
fn unknown_enzyme_is_none() {
    assert!(enzyme_by_name("NotAnEnzyme").is_none());
}
