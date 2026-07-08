use bio_core::analysis::variant::call_substitutions;
use bio_core::sequence::Sequence;
use bio_core::vcf::{from_substitutions, write_vcf, VcfRecord};

#[test]
fn writes_header_and_one_based_position() {
    let rec = VcfRecord {
        chrom: "chr1".into(),
        pos: 99, // 0-based -> POS 100
        id: String::new(),
        reference: "A".into(),
        alternate: "G".into(),
        qual: Some(30.0),
        filter: "PASS".into(),
        info: "DP=12;AF=0.750".into(),
    };
    let vcf = write_vcf(&[rec]);
    assert!(vcf.starts_with("##fileformat=VCFv4.2\n"));
    assert!(vcf.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
    assert!(vcf.contains("chr1\t100\t.\tA\tG\t30.0\tPASS\tDP=12;AF=0.750\n"));
}

#[test]
fn empty_fields_render_as_dot() {
    let rec = VcfRecord {
        chrom: "c".into(),
        pos: 0,
        id: String::new(),
        reference: "T".into(),
        alternate: "C".into(),
        qual: None,
        filter: String::new(),
        info: String::new(),
    };
    let vcf = write_vcf(&[rec]);
    assert!(vcf.contains("c\t1\t.\tT\tC\t.\t.\t.\n"));
}

#[test]
fn substitutions_round_trip_to_vcf() {
    let reference = Sequence::dna("ACGTACGT").unwrap();
    let sample = Sequence::dna("ACCTACGA").unwrap(); // diffs at pos 2 (G>C), pos 7 (T>A)
    let variants = call_substitutions(&reference, &sample).unwrap();
    let vcf = write_vcf(&from_substitutions("gene", &variants));
    assert!(vcf.contains("gene\t3\t.\tG\tC\t.\tPASS\t.\n")); // pos 2 -> POS 3
    assert!(vcf.contains("gene\t8\t.\tT\tA\t.\tPASS\t.\n")); // pos 7 -> POS 8
}
