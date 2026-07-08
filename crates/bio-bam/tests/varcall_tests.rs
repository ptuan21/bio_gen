use bio_bam::pileup::{PileupColumn, StrandCount};
use bio_bam::varcall::call_variants;
use bio_core::sequence::Sequence;

fn col(pos: i32, a: u32, c: u32, g: u32, t: u32) -> PileupColumn {
    // All support placed on the forward strand (strand bias unused in most tests).
    PileupColumn {
        ref_pos: pos,
        depth: a + c + g + t,
        a,
        c,
        g,
        t,
        n: 0,
        del: 0,
        fwd: StrandCount { a, c, g, t, n: 0 },
        rev: StrandCount::default(),
    }
}

/// Column at `pos` with reference C and alt G, split across strands.
fn col_g(pos: i32, c_fwd: u32, c_rev: u32, g_fwd: u32, g_rev: u32) -> PileupColumn {
    let (c, g) = (c_fwd + c_rev, g_fwd + g_rev);
    PileupColumn {
        ref_pos: pos,
        depth: c + g,
        a: 0,
        c,
        g,
        t: 0,
        n: 0,
        del: 0,
        fwd: StrandCount { a: 0, c: c_fwd, g: g_fwd, t: 0, n: 0 },
        rev: StrandCount { a: 0, c: c_rev, g: g_rev, t: 0, n: 0 },
    }
}

#[test]
fn strand_bias_filters_one_sided_alt() {
    let reference = Sequence::dna("ACGT").unwrap(); // pos 1 = C

    // Balanced alt G (4 fwd / 4 rev) passes a 0.2 strand threshold.
    let balanced = vec![col_g(1, 1, 1, 4, 4)];
    let v = call_variants(&balanced, &reference, 0, 3, 0.2, 0.2);
    assert_eq!(v.len(), 1);
    assert_eq!((v[0].alt_fwd, v[0].alt_rev), (4, 4));

    // One-strand alt G (8 fwd / 0 rev) is filtered as strand-biased...
    let biased = vec![col_g(1, 1, 1, 8, 0)];
    assert!(call_variants(&biased, &reference, 0, 3, 0.2, 0.2).is_empty());

    // ...but is still called with the filter disabled.
    assert_eq!(call_variants(&biased, &reference, 0, 3, 0.2, 0.0).len(), 1);
}

#[test]
fn pileup_variants_render_to_vcf() {
    use bio_bam::varcall::pileup_variants_to_vcf;
    let reference = Sequence::dna("ACGT").unwrap();
    let cols = vec![col_g(1, 1, 1, 4, 4)]; // ref C, alt G, DP 10, 4 fwd / 4 rev
    let variants = call_variants(&cols, &reference, 0, 3, 0.2, 0.0);
    let vcf = pileup_variants_to_vcf("chr1", &variants);
    assert!(vcf.starts_with("##fileformat=VCFv4.2"));
    // pos 1 -> POS 2, INFO carries depth/AF/strand support
    assert!(vcf.contains("chr1\t2\t.\tC\tG\t.\tPASS\tDP=10;AF=0.800;SB=4,4\n"), "{vcf}");
}

#[test]
fn calls_snv_above_thresholds() {
    let reference = Sequence::dna("ACGT").unwrap();
    // pos 1 is reference C; reads show G=8, C=2 -> C>G at freq 0.8
    let cols = vec![col(1, 0, 2, 8, 0)];
    let variants = call_variants(&cols, &reference, 0, 3, 0.2, 0.0);
    assert_eq!(variants.len(), 1);
    let v = variants[0];
    assert_eq!((v.ref_pos, v.reference, v.alternate), (1, 'C', 'G'));
    assert_eq!((v.depth, v.alt_count), (10, 8));
    assert!((v.allele_freq - 0.8).abs() < 1e-9);
}

#[test]
fn no_variant_when_matches_reference() {
    let reference = Sequence::dna("ACGT").unwrap();
    let cols = vec![col(0, 10, 0, 0, 0)]; // all A, reference A
    assert!(call_variants(&cols, &reference, 0, 3, 0.2, 0.0).is_empty());
}

#[test]
fn skips_low_depth() {
    let reference = Sequence::dna("ACGT").unwrap();
    let cols = vec![col(2, 0, 0, 0, 2)]; // depth 2 < min_depth 3
    assert!(call_variants(&cols, &reference, 0, 3, 0.2, 0.0).is_empty());
}

#[test]
fn skips_low_frequency() {
    let reference = Sequence::dna("ACGT").unwrap();
    // reference T at pos 3; A=1, T=9 -> alt A freq 0.1 < 0.2
    let cols = vec![col(3, 1, 0, 0, 9)];
    assert!(call_variants(&cols, &reference, 0, 3, 0.2, 0.0).is_empty());
}

#[test]
fn skips_ambiguous_reference() {
    let reference = Sequence::dna("ANGT").unwrap();
    let cols = vec![col(1, 9, 0, 1, 0)]; // reference N at pos 1
    assert!(call_variants(&cols, &reference, 0, 3, 0.2, 0.0).is_empty());
}

#[test]
fn honours_ref_offset() {
    // reference covers positions 100.. ; column at 101 is reference C
    let reference = Sequence::dna("ACGT").unwrap();
    let cols = vec![col(101, 0, 2, 8, 0)];
    let variants = call_variants(&cols, &reference, 100, 3, 0.2, 0.0);
    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0].ref_pos, 101);
}

#[test]
fn rna_reference_uracil_matches_t_reads() {
    // reference U should equal T reads (no false variant)
    let reference = Sequence::rna("ACGU").unwrap();
    let cols = vec![col(3, 0, 0, 0, 10)]; // all T at RNA-U position
    assert!(call_variants(&cols, &reference, 0, 3, 0.2, 0.0).is_empty());
}
