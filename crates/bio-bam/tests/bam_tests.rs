use bio_bam::bgzf;
use bio_bam::error::BamError;
use bio_bam::{parse, read_bam, read_bam_region};

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

/// BAI matching `build_bam`: the single record starts at uncompressed offset 24
/// and ends at 73 (both in block 0), and covers reference bin 4681.
fn build_bai() -> Vec<u8> {
    let mut b = Vec::new();
    b.extend_from_slice(b"BAI\x01");
    b.extend_from_slice(&1u32.to_le_bytes()); // n_ref
    b.extend_from_slice(&1u32.to_le_bytes()); // n_bin
    b.extend_from_slice(&4681u32.to_le_bytes()); // bin
    b.extend_from_slice(&1u32.to_le_bytes()); // n_chunk
    b.extend_from_slice(&24u64.to_le_bytes()); // chunk beg (virtual offset)
    b.extend_from_slice(&73u64.to_le_bytes()); // chunk end
    b.extend_from_slice(&1u32.to_le_bytes()); // n_intv
    b.extend_from_slice(&24u64.to_le_bytes()); // interval 0
    b
}

#[test]
fn region_query_returns_overlapping_record() {
    let bam = make_bgzf_block(&build_bam());
    let bai = build_bai();
    let hits = read_bam_region(&bam, &bai, "ref", 0, 100).unwrap();
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].seq, "ACGT");
    assert_eq!(hits[0].pos, 5);
    assert_eq!(hits[0].ref_span, 4); // 4M
}

#[test]
fn region_query_filters_non_overlapping() {
    let bam = make_bgzf_block(&build_bam());
    let bai = build_bai();
    // Record spans [5, 9); a [50, 60) query must return nothing.
    let hits = read_bam_region(&bam, &bai, "ref", 50, 60).unwrap();
    assert!(hits.is_empty());
}

#[test]
fn region_query_unknown_reference_errors() {
    let bam = make_bgzf_block(&build_bam());
    let bai = build_bai();
    assert!(matches!(
        read_bam_region(&bam, &bai, "chrX", 0, 100),
        Err(BamError::UnknownReference(_))
    ));
}

fn header_bytes() -> Vec<u8> {
    let mut b = Vec::new();
    b.extend_from_slice(b"BAM\x01");
    b.extend_from_slice(&0u32.to_le_bytes()); // l_text
    b.extend_from_slice(&1u32.to_le_bytes()); // n_ref
    b.extend_from_slice(&4u32.to_le_bytes()); // l_name
    b.extend_from_slice(b"ref\0");
    b.extend_from_slice(&2000i32.to_le_bytes()); // l_ref
    b
}

/// One BAM record (`block_size` + body): ACGT / 4M / mapq 30 at `pos`.
fn rec_block(pos: i32, name: &str) -> Vec<u8> {
    let mut rec = Vec::new();
    rec.extend_from_slice(&0i32.to_le_bytes()); // ref_id
    rec.extend_from_slice(&pos.to_le_bytes());
    rec.push((name.len() + 1) as u8); // l_read_name
    rec.push(30); // mapq
    rec.extend_from_slice(&0u16.to_le_bytes()); // bin (ignored by reader)
    rec.extend_from_slice(&1u16.to_le_bytes()); // n_cigar
    rec.extend_from_slice(&0u16.to_le_bytes()); // flag
    rec.extend_from_slice(&4u32.to_le_bytes()); // l_seq
    rec.extend_from_slice(&(-1i32).to_le_bytes());
    rec.extend_from_slice(&(-1i32).to_le_bytes());
    rec.extend_from_slice(&0i32.to_le_bytes());
    rec.extend_from_slice(name.as_bytes());
    rec.push(0);
    rec.extend_from_slice(&64u32.to_le_bytes()); // CIGAR 4M
    rec.extend_from_slice(&[0x12, 0x48]); // seq ACGT
    rec.extend_from_slice(&[30, 30, 30, 30]); // qual
    let mut out = (rec.len() as u32).to_le_bytes().to_vec();
    out.extend_from_slice(&rec);
    out
}

#[test]
fn region_query_reads_second_block_without_inflating_first() {
    // Two BGZF blocks: header+r1 in block 1, r2 (far away) alone in block 2.
    let header = header_bytes();
    let r1 = rec_block(5, "r1");
    let r2 = rec_block(1000, "r2");
    let block1_raw = [header.clone(), r1.clone()].concat();
    let block1 = make_bgzf_block(&block1_raw);
    let block2 = make_bgzf_block(&r2);
    let coffset2 = block1.len() as u64;
    let bam = [block1, block2].concat();

    // Virtual offsets (coffset << 16 | uoffset).
    let vbeg1 = header.len() as u64;
    let vend1 = block1_raw.len() as u64;
    let vbeg2 = coffset2 << 16;
    let vend2 = (coffset2 << 16) | (r2.len() as u64);

    // BAI: both records share bin 4681 (positions < 16 kbp), one chunk each.
    let mut bai = Vec::new();
    bai.extend_from_slice(b"BAI\x01");
    bai.extend_from_slice(&1u32.to_le_bytes()); // n_ref
    bai.extend_from_slice(&1u32.to_le_bytes()); // n_bin
    bai.extend_from_slice(&4681u32.to_le_bytes());
    bai.extend_from_slice(&2u32.to_le_bytes()); // n_chunk
    for (b, e) in [(vbeg1, vend1), (vbeg2, vend2)] {
        bai.extend_from_slice(&b.to_le_bytes());
        bai.extend_from_slice(&e.to_le_bytes());
    }
    bai.extend_from_slice(&1u32.to_le_bytes()); // n_intv
    bai.extend_from_slice(&0u64.to_le_bytes()); // linear index: keep-all

    // Record living in the SECOND block resolves correctly (random access).
    let hits = read_bam_region(&bam, &bai, "ref", 1000, 1004).unwrap();
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].name, "r2");
    assert_eq!(hits[0].pos, 1000);

    // And the first block still resolves independently.
    let hits = read_bam_region(&bam, &bai, "ref", 5, 9).unwrap();
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].name, "r1");
}

#[test]
fn region_query_rejects_bad_bai() {
    let bam = make_bgzf_block(&build_bam());
    assert_eq!(
        read_bam_region(&bam, b"NOPE", "ref", 0, 100),
        Err(BamError::BadBaiMagic)
    );
}
