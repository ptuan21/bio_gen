use bio_bam::bgzf;
use bio_bam::error::BamError;
use bio_bam::{parse, read_bam};

/// Wrap raw bytes in a single valid BGZF block (CRC is left zero; the reader
/// does not verify it).
fn make_bgzf_block(data: &[u8]) -> Vec<u8> {
    let cdata = miniz_oxide::deflate::compress_to_vec(data, 6);
    let total = 12 + 6 + cdata.len() + 8;
    let bsize = (total - 1) as u16;

    let mut block = Vec::new();
    block.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x04]);
    block.extend_from_slice(&[0, 0, 0, 0]); // mtime
    block.push(0); // xfl
    block.push(0xff); // os
    block.extend_from_slice(&6u16.to_le_bytes()); // xlen
    block.extend_from_slice(b"BC");
    block.extend_from_slice(&2u16.to_le_bytes()); // slen
    block.extend_from_slice(&bsize.to_le_bytes());
    block.extend_from_slice(&cdata);
    block.extend_from_slice(&0u32.to_le_bytes()); // crc (ignored)
    block.extend_from_slice(&(data.len() as u32).to_le_bytes()); // isize
    block
}

/// A minimal BAM stream: empty header text, one reference, one aligned record.
fn build_bam() -> Vec<u8> {
    let mut b = Vec::new();
    b.extend_from_slice(b"BAM\x01");
    b.extend_from_slice(&0u32.to_le_bytes()); // l_text
    b.extend_from_slice(&1u32.to_le_bytes()); // n_ref
    b.extend_from_slice(&4u32.to_le_bytes()); // l_name
    b.extend_from_slice(b"ref\0");
    b.extend_from_slice(&100i32.to_le_bytes()); // l_ref

    let mut rec = Vec::new();
    rec.extend_from_slice(&0i32.to_le_bytes()); // ref_id
    rec.extend_from_slice(&5i32.to_le_bytes()); // pos
    rec.push(3); // l_read_name
    rec.push(30); // mapq
    rec.extend_from_slice(&0u16.to_le_bytes()); // bin
    rec.extend_from_slice(&1u16.to_le_bytes()); // n_cigar
    rec.extend_from_slice(&0u16.to_le_bytes()); // flag
    rec.extend_from_slice(&4u32.to_le_bytes()); // l_seq
    rec.extend_from_slice(&(-1i32).to_le_bytes()); // next_ref
    rec.extend_from_slice(&(-1i32).to_le_bytes()); // next_pos
    rec.extend_from_slice(&0i32.to_le_bytes()); // tlen
    rec.extend_from_slice(b"r1\0");
    rec.extend_from_slice(&64u32.to_le_bytes()); // cigar 4M
    rec.extend_from_slice(&[0x12, 0x48]); // seq ACGT
    rec.extend_from_slice(&[30, 30, 30, 30]); // qual

    b.extend_from_slice(&(rec.len() as u32).to_le_bytes());
    b.extend_from_slice(&rec);
    b
}

#[test]
fn bgzf_roundtrip_single_block() {
    let data = b"hello bioinformatics world";
    let block = make_bgzf_block(data);
    assert_eq!(bgzf::decompress(&block).unwrap(), data);
}

#[test]
fn bgzf_roundtrip_multiple_blocks() {
    let mut stream = make_bgzf_block(b"first-");
    stream.extend(make_bgzf_block(b"second"));
    assert_eq!(bgzf::decompress(&stream).unwrap(), b"first-second");
}

#[test]
fn bgzf_rejects_bad_block() {
    let bad = [0u8; 20];
    assert_eq!(bgzf::decompress(&bad), Err(BamError::BadBlock));
}

#[test]
fn parses_header_and_record() {
    let (header, records) = parse(&build_bam()).unwrap();
    assert_eq!(header.references.len(), 1);
    assert_eq!(header.references[0].name, "ref");
    assert_eq!(header.references[0].length, 100);

    assert_eq!(records.len(), 1);
    let r = &records[0];
    assert_eq!(r.name, "r1");
    assert_eq!(r.ref_name.as_deref(), Some("ref"));
    assert_eq!(r.pos, 5);
    assert_eq!(r.mapq, 30);
    assert_eq!(r.cigar, "4M");
    assert_eq!(r.seq, "ACGT");
    assert_eq!(r.qual, vec![30, 30, 30, 30]);
}

#[test]
fn read_bam_decompresses_and_parses() {
    let block = make_bgzf_block(&build_bam());
    let (header, records) = read_bam(&block).unwrap();
    assert_eq!(header.references[0].name, "ref");
    assert_eq!(records[0].seq, "ACGT");
}

#[test]
fn rejects_bad_magic() {
    let mut data = build_bam();
    data[0] = b'X';
    assert_eq!(parse(&data), Err(BamError::BadMagic));
}

#[test]
fn truncated_record_errors() {
    let mut data = build_bam();
    data.truncate(data.len() - 3);
    assert_eq!(parse(&data), Err(BamError::Truncated));
}
