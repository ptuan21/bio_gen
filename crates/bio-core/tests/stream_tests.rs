use bio_core::parser::{FastaStreamer, RecordSummary};

fn run(data: &[u8], chunk: usize) -> Vec<RecordSummary> {
    let mut s = FastaStreamer::new();
    let mut out = Vec::new();
    for c in data.chunks(chunk.max(1)) {
        out.extend(s.push(c));
    }
    out.extend(s.finish());
    out
}

const FASTA: &[u8] = b">a first record\nACGT\nACGT\n>b\nGGCC\n>c empty tail\n";

#[test]
fn summaries_match_reference_values() {
    let recs = run(FASTA, FASTA.len());
    assert_eq!(recs.len(), 3);

    assert_eq!(recs[0].id, "a");
    assert_eq!(recs[0].description, "first record");
    assert_eq!(recs[0].length, 8);
    assert_eq!(recs[0].counts.a, 2);
    assert_eq!(recs[0].counts.g, 2);
    assert!((recs[0].gc_content - 0.5).abs() < 1e-12);

    assert_eq!(recs[1].id, "b");
    assert_eq!(recs[1].length, 4);
    assert!((recs[1].gc_content - 1.0).abs() < 1e-12);

    assert_eq!(recs[2].id, "c");
    assert_eq!(recs[2].length, 0);
    assert_eq!(recs[2].gc_content, 0.0);
}

#[test]
fn result_is_independent_of_chunk_size() {
    let reference = run(FASTA, FASTA.len());
    // Every chunk size, including 1 byte (splits mid-header and mid-sequence),
    // must yield identical summaries.
    for chunk in [1, 2, 3, 5, 7, 11, 13, FASTA.len()] {
        assert_eq!(run(FASTA, chunk), reference, "chunk size {chunk}");
    }
}

#[test]
fn handles_crlf_and_lowercase() {
    let data = b">x\r\nacgt\r\nAC\r\n";
    let recs = run(data, 4);
    assert_eq!(recs.len(), 1);
    assert_eq!(recs[0].length, 6);
    assert_eq!(recs[0].counts.a, 2);
    assert_eq!(recs[0].counts.c, 2);
    assert_eq!(recs[0].counts.t, 1);
}

#[test]
fn empty_input_yields_nothing() {
    assert!(run(b"", 4).is_empty());
    assert!(run(b"\n\n", 1).is_empty());
}

#[test]
fn sequence_before_header_is_ignored() {
    // Leading bases with no header belong to no record and are dropped.
    let recs = run(b"ACGT\n>a\nGG\n", 3);
    assert_eq!(recs.len(), 1);
    assert_eq!(recs[0].id, "a");
    assert_eq!(recs[0].length, 2);
}

#[test]
fn counts_unknown_bases_as_other() {
    let recs = run(b">x\nACGTNXYZ\n", 2);
    assert_eq!(recs[0].counts.n, 1);
    assert_eq!(recs[0].counts.other, 3); // X, Y, Z
}
