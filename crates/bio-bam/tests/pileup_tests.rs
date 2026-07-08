use bio_bam::bam::BamRecord;
use bio_bam::pileup::pileup;

fn rec(pos: i32, cigar: &str, seq: &str) -> BamRecord {
    BamRecord {
        name: "r".into(),
        flag: 0,
        ref_name: Some("ref".into()),
        pos,
        ref_span: 0, // unused by pileup, which walks the CIGAR
        mapq: 60,
        cigar: cigar.into(),
        seq: seq.into(),
        qual: Vec::new(),
    }
}

#[test]
fn simple_match_pileup() {
    let cols = pileup(&[rec(5, "4M", "ACGT")], 0, 100, 0);
    assert_eq!(cols.len(), 4);
    assert_eq!((cols[0].ref_pos, cols[0].a, cols[0].depth), (5, 1, 1));
    assert_eq!((cols[1].ref_pos, cols[1].c), (6, 1));
    assert_eq!((cols[2].ref_pos, cols[2].g), (7, 1));
    assert_eq!((cols[3].ref_pos, cols[3].t), (8, 1));
}

#[test]
fn overlapping_reads_increase_depth() {
    let cols = pileup(&[rec(0, "3M", "ACG"), rec(1, "3M", "CGT")], 0, 100, 0);
    // positions: 0(A) 1(C,C) 2(G,G) 3(T)
    let depth_at = |p: i32| cols.iter().find(|c| c.ref_pos == p).unwrap().depth;
    assert_eq!(depth_at(0), 1);
    assert_eq!(depth_at(1), 2);
    assert_eq!(depth_at(2), 2);
    assert_eq!(depth_at(3), 1);
    let col1 = cols.iter().find(|c| c.ref_pos == 1).unwrap();
    assert_eq!(col1.c, 2);
    assert_eq!(col1.consensus(), Some(('C', 2)));
}

#[test]
fn deletion_marks_column() {
    // 2M1D2M over "ACGT": ref 0=A 1=C 2=del 3=G 4=T
    let cols = pileup(&[rec(0, "2M1D2M", "ACGT")], 0, 100, 0);
    let col2 = cols.iter().find(|c| c.ref_pos == 2).unwrap();
    assert_eq!(col2.del, 1);
    assert_eq!(col2.depth, 1);
    assert_eq!(col2.a + col2.c + col2.g + col2.t, 0);
    assert_eq!(cols.iter().find(|c| c.ref_pos == 3).unwrap().g, 1);
}

#[test]
fn insertion_does_not_add_columns() {
    // 2M2I2M over "ACGTAA": M=AC, I=GT (skipped), M=AA at ref 2,3
    let cols = pileup(&[rec(0, "2M2I2M", "ACGTAA")], 0, 100, 0);
    assert_eq!(cols.len(), 4);
    assert_eq!(cols.iter().map(|c| c.ref_pos).collect::<Vec<_>>(), vec![0, 1, 2, 3]);
    assert_eq!(cols[2].a, 1);
    assert_eq!(cols[3].a, 1);
}

#[test]
fn ref_skip_does_not_cover() {
    // 2M3N2M: N (intron) leaves ref 2,3,4 uncovered
    let cols = pileup(&[rec(0, "2M3N2M", "ACGT")], 0, 100, 0);
    let positions: Vec<i32> = cols.iter().map(|c| c.ref_pos).collect();
    assert_eq!(positions, vec![0, 1, 5, 6]);
}

#[test]
fn region_clips_columns() {
    let cols = pileup(&[rec(5, "4M", "ACGT")], 6, 8, 0);
    assert_eq!(cols.iter().map(|c| c.ref_pos).collect::<Vec<_>>(), vec![6, 7]);
}

#[test]
fn unmapped_and_star_cigar_ignored() {
    let cols = pileup(&[rec(-1, "*", "ACGT"), rec(3, "*", "ACGT")], 0, 100, 0);
    assert!(cols.is_empty());
}

fn rec_q(pos: i32, seq: &str, quals: Vec<u8>) -> BamRecord {
    BamRecord {
        name: "r".into(),
        flag: 0,
        ref_name: Some("ref".into()),
        pos,
        ref_span: 0,
        mapq: 60,
        cigar: format!("{}M", seq.len()),
        seq: seq.into(),
        qual: quals,
    }
}

#[test]
fn low_quality_bases_are_filtered() {
    // Phred: pos0=40 (keep), pos1=10 (drop at min_qual 20), pos2=30, pos3=5 (drop)
    let record = rec_q(0, "ACGT", vec![40, 10, 30, 5]);
    let cols = pileup(&[record], 0, 100, 20);
    // Only positions 0 and 2 survive.
    assert_eq!(cols.iter().map(|c| c.ref_pos).collect::<Vec<_>>(), vec![0, 2]);
    assert_eq!(cols[0].a, 1);
    assert_eq!(cols[1].g, 1);
}

#[test]
fn min_qual_zero_keeps_everything() {
    let record = rec_q(0, "ACGT", vec![0, 1, 2, 3]);
    let cols = pileup(&[record], 0, 100, 0);
    assert_eq!(cols.len(), 4);
}

#[test]
fn missing_quality_string_always_counts() {
    // Empty qual (0xFF fill) -> no filtering even at a high threshold.
    let record = rec_q(0, "ACGT", Vec::new());
    let cols = pileup(&[record], 0, 100, 40);
    assert_eq!(cols.len(), 4);
}

#[test]
fn tracks_strand_from_flag() {
    let mut fwd = rec(0, "2M", "AC");
    let mut rev = rec(0, "2M", "AC");
    rev.flag = 0x10; // reverse strand
    fwd.flag = 0;
    let cols = pileup(&[fwd, rev], 0, 100, 0);
    let col0 = &cols[0];
    assert_eq!(col0.a, 2);
    assert_eq!(col0.strand_counts(b'A'), (1, 1)); // one forward, one reverse
}
