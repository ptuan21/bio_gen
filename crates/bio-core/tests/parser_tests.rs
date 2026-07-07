use bio_core::error::BioError;
use bio_core::parser::{parse_fasta_str, parse_fastq_str, FastaReader};
use bio_core::sequence::SeqKind;

#[test]
fn parses_multiple_records_with_wrapped_lines() {
    let input = ">seq1 first sequence\nACGT\nACGT\n>seq2\nTTTT\n";
    let records = parse_fasta_str(input, SeqKind::Dna).unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].id, "seq1");
    assert_eq!(records[0].description, "first sequence");
    assert_eq!(records[0].sequence.to_string(), "ACGTACGT");
    assert_eq!(records[1].id, "seq2");
    assert_eq!(records[1].sequence.len(), 4);
}

#[test]
fn skips_blank_lines_and_uppercases() {
    let input = "\n\n>lower\nacgt\n\n";
    let records = parse_fasta_str(input, SeqKind::Dna).unwrap();
    assert_eq!(records[0].sequence.to_string(), "ACGT");
}

#[test]
fn empty_input_is_rejected() {
    assert_eq!(parse_fasta_str("   \n", SeqKind::Dna), Err(BioError::EmptyInput));
}

#[test]
fn missing_header_is_malformed() {
    let input = "ACGT\n>seq\nACGT\n";
    let err = parse_fasta_str(input, SeqKind::Dna).unwrap_err();
    assert!(matches!(err, BioError::MalformedFasta { line: 1, .. }));
}

#[test]
fn invalid_nucleotide_reports_position() {
    let input = ">bad\nACGZ\n";
    let err = parse_fasta_str(input, SeqKind::Dna).unwrap_err();
    assert_eq!(err, BioError::InvalidNucleotide { symbol: 'Z', position: 3 });
}

#[test]
fn dna_reader_rejects_uracil() {
    let input = ">rna_in_dna\nACGU\n";
    let err = parse_fasta_str(input, SeqKind::Dna).unwrap_err();
    assert!(matches!(err, BioError::InvalidNucleotide { symbol: 'U', .. }));
}

#[test]
fn streaming_reader_yields_lazily() {
    let input = ">a\nAAAA\n>b\nCCCC\n>c\nGGGG\n";
    let ids: Vec<String> = FastaReader::new(input.as_bytes(), SeqKind::Dna)
        .map(|r| r.unwrap().id)
        .collect();
    assert_eq!(ids, vec!["a", "b", "c"]);
}

#[test]
fn parses_valid_fastq_with_quality() {
    let input = "@read1 sample\nACGT\n+\nIIII\n";
    let records = parse_fastq_str(input, SeqKind::Dna).unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].sequence.to_string(), "ACGT");
    // 'I' == 73, Phred offset 33 -> score 40
    assert_eq!(records[0].quality, Some(vec![40, 40, 40, 40]));
}

#[test]
fn fastq_quality_length_mismatch_is_rejected() {
    let input = "@r\nACGT\n+\nII\n";
    let err = parse_fastq_str(input, SeqKind::Dna).unwrap_err();
    assert_eq!(err, BioError::LengthMismatch { expected: 4, found: 2 });
}

#[test]
fn fastq_missing_separator_is_rejected() {
    let input = "@r\nACGT\nX\nIIII\n";
    let err = parse_fastq_str(input, SeqKind::Dna).unwrap_err();
    assert!(matches!(err, BioError::MalformedFastq { .. }));
}
