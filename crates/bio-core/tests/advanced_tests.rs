use bio_core::analysis::kmer::kmer_counts;
use bio_core::analysis::orf::find_orfs;
use bio_core::analysis::search::Strand;
use bio_core::analysis::translate::six_frames;
use bio_core::sequence::Sequence;

#[test]
fn six_frames_has_six_entries_with_signed_labels() {
    let seq = Sequence::dna("ATGGCCTAA").unwrap();
    let frames = six_frames(&seq);
    let labels: Vec<i8> = frames.iter().map(|f| f.frame).collect();
    assert_eq!(labels, vec![1, 2, 3, -1, -2, -3]);
    assert_eq!(frames[0].protein, "MA*");
}

#[test]
fn six_frames_handles_short_sequence() {
    let seq = Sequence::dna("AT").unwrap();
    let frames = six_frames(&seq);
    assert_eq!(frames.len(), 6);
    assert!(frames.iter().all(|f| f.protein.is_empty()));
}

#[test]
fn find_forward_orf() {
    let seq = Sequence::dna("ATGGCCTAA").unwrap();
    let orfs = find_orfs(&seq, 2);
    assert_eq!(orfs.len(), 1);
    let orf = &orfs[0];
    assert_eq!(orf.start, 0);
    assert_eq!(orf.end, 9);
    assert_eq!(orf.strand, Strand::Forward);
    assert_eq!(orf.frame, 1);
    assert_eq!(orf.protein, "MA");
    assert_eq!(orf.dna, "ATGGCCTAA");
}

#[test]
fn find_reverse_orf_mapped_to_forward_coords() {
    // reverse complement of "ATGAAATAA" (M K *) is "TTATTTCAT"
    let seq = Sequence::dna("TTATTTCAT").unwrap();
    let orfs = find_orfs(&seq, 2);
    assert_eq!(orfs.len(), 1);
    let orf = &orfs[0];
    assert_eq!(orf.strand, Strand::Reverse);
    assert_eq!(orf.frame, -1);
    assert_eq!((orf.start, orf.end), (0, 9));
    assert_eq!(orf.protein, "MK");
}

#[test]
fn orf_min_length_filter() {
    let seq = Sequence::dna("ATGTAA").unwrap();
    assert!(find_orfs(&seq, 2).is_empty());
    let orfs = find_orfs(&seq, 1);
    assert_eq!(orfs.len(), 1);
    assert_eq!(orfs[0].protein, "M");
}

#[test]
fn orf_without_stop_is_ignored() {
    let seq = Sequence::dna("ATGGCCGCC").unwrap();
    assert!(find_orfs(&seq, 1).is_empty());
}

#[test]
fn kmer_counts_sorted_by_frequency() {
    let seq = Sequence::dna("ATATAT").unwrap();
    let counts = kmer_counts(&seq, 2);
    assert_eq!(counts, vec![("AT".to_string(), 3), ("TA".to_string(), 2)]);
}

#[test]
fn kmer_larger_than_sequence_is_empty() {
    let seq = Sequence::dna("AT").unwrap();
    assert!(kmer_counts(&seq, 5).is_empty());
    assert!(kmer_counts(&seq, 0).is_empty());
}
