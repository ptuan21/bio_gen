use crate::bai;
use crate::bgzf;
use crate::error::{BamError, Result};
use crate::reader::Cursor;

const SEQ_NT: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
const CIGAR_OPS: &[u8; 9] = b"MIDNSHP=X";
/// Which CIGAR ops (by index into `CIGAR_OPS`) consume the reference: M D N = X.
const REF_CONSUMING: [bool; 9] = [true, false, true, true, false, false, false, true, true];

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reference {
    pub name: String,
    pub length: i32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BamHeader {
    pub text: String,
    pub references: Vec<Reference>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BamRecord {
    pub name: String,
    pub flag: u16,
    /// Reference name, or `None` when unmapped (`ref_id == -1`).
    pub ref_name: Option<String>,
    /// 0-based leftmost position; `-1` when unset.
    pub pos: i32,
    /// Reference bases spanned (from the CIGAR); 0 for unmapped reads.
    pub ref_span: u32,
    pub mapq: u8,
    pub cigar: String,
    pub seq: String,
    /// Phred qualities; empty when unavailable (`0xFF` fill).
    pub qual: Vec<u8>,
}

impl BamRecord {
    /// Whether this record overlaps the 0-based half-open interval `[beg, end)`.
    fn overlaps(&self, ref_name: &str, beg: i32, end: i32) -> bool {
        self.ref_name.as_deref() == Some(ref_name)
            && self.pos >= 0
            && self.pos < end
            && self.pos + self.ref_span as i32 > beg
    }
}

fn trim_nul(bytes: &[u8]) -> String {
    let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
    String::from_utf8_lossy(&bytes[..end]).into_owned()
}

fn decode_seq(bytes: &[u8], l_seq: usize) -> String {
    (0..l_seq)
        .map(|i| {
            let byte = bytes[i / 2];
            let nibble = if i % 2 == 0 { byte >> 4 } else { byte & 0x0f };
            SEQ_NT[nibble as usize] as char
        })
        .collect()
}

fn decode_cigar(cursor: &mut Cursor, n_ops: usize) -> Result<(String, u32)> {
    if n_ops == 0 {
        return Ok(("*".to_string(), 0));
    }
    let mut cigar = String::new();
    let mut ref_span = 0u32;
    for _ in 0..n_ops {
        let v = cursor.u32()?;
        let len = v >> 4;
        let op = (v & 0xf) as usize;
        cigar.push_str(&len.to_string());
        cigar.push(CIGAR_OPS[op] as char);
        if REF_CONSUMING[op] {
            ref_span += len;
        }
    }
    Ok((cigar, ref_span))
}

fn read_header(cursor: &mut Cursor) -> Result<BamHeader> {
    if cursor.take(4)? != b"BAM\x01" {
        return Err(BamError::BadMagic);
    }
    let l_text = cursor.u32()? as usize;
    let text = String::from_utf8_lossy(cursor.take(l_text)?).into_owned();

    let n_ref = cursor.u32()? as usize;
    let mut references = Vec::with_capacity(n_ref);
    for _ in 0..n_ref {
        let l_name = cursor.u32()? as usize;
        let name = trim_nul(cursor.take(l_name)?);
        let length = cursor.i32()?;
        references.push(Reference { name, length });
    }
    Ok(BamHeader { text, references })
}

fn read_record(cursor: &mut Cursor, refs: &[Reference]) -> Result<BamRecord> {
    let block_size = cursor.u32()? as usize;
    let start = cursor.position();

    let ref_id = cursor.i32()?;
    let pos = cursor.i32()?;
    let l_read_name = cursor.u8()? as usize;
    let mapq = cursor.u8()?;
    let _bin = cursor.u16()?;
    let n_cigar = cursor.u16()? as usize;
    let flag = cursor.u16()?;
    let l_seq = cursor.u32()? as usize;
    let _next_ref = cursor.i32()?;
    let _next_pos = cursor.i32()?;
    let _tlen = cursor.i32()?;

    let name = trim_nul(cursor.take(l_read_name)?);
    let (cigar, ref_span) = decode_cigar(cursor, n_cigar)?;
    let seq = decode_seq(cursor.take(l_seq.div_ceil(2))?, l_seq);
    let qual_bytes = cursor.take(l_seq)?;
    let qual = if qual_bytes.first() == Some(&0xff) {
        Vec::new()
    } else {
        qual_bytes.to_vec()
    };

    // Skip any auxiliary tag data at the end of the record.
    let consumed = cursor.position() - start;
    if consumed < block_size {
        cursor.take(block_size - consumed)?;
    }

    let ref_name = usize::try_from(ref_id)
        .ok()
        .and_then(|i| refs.get(i))
        .map(|r| r.name.clone());

    Ok(BamRecord {
        name,
        flag,
        ref_name,
        pos,
        ref_span,
        mapq,
        cigar,
        seq,
        qual,
    })
}

/// Parse already-decompressed BAM bytes into a header and all records.
pub fn parse(data: &[u8]) -> Result<(BamHeader, Vec<BamRecord>)> {
    let mut cursor = Cursor::new(data);
    let header = read_header(&mut cursor)?;
    let mut records = Vec::new();
    while !cursor.eof() {
        records.push(read_record(&mut cursor, &header.references)?);
    }
    Ok((header, records))
}

/// Parse only the BAM header (magic, SAM text and references).
pub fn parse_header(data: &[u8]) -> Result<BamHeader> {
    read_header(&mut Cursor::new(data))
}

/// Inflate just enough leading BGZF blocks to parse the header, then stop.
fn inflate_header(bgzf_bytes: &[u8]) -> Result<BamHeader> {
    let mut buf = Vec::new();
    let mut pos = 0;
    loop {
        match parse_header(&buf) {
            Ok(header) => return Ok(header),
            Err(BamError::Truncated) => {}
            Err(e) => return Err(e),
        }
        if pos >= bgzf_bytes.len() {
            return Err(BamError::Truncated);
        }
        pos += bgzf::inflate_block(bgzf_bytes, pos, &mut buf)?;
    }
}

/// Decompress a BGZF-wrapped BAM file and parse all records.
pub fn read_bam(bgzf_bytes: &[u8]) -> Result<(BamHeader, Vec<BamRecord>)> {
    parse(&bgzf::decompress(bgzf_bytes)?)
}

/// Fetch records overlapping `ref_name:[beg, end)` using a BAI index, reading
/// only the byte ranges the index points at instead of scanning every record.
pub fn read_bam_region(
    bgzf_bytes: &[u8],
    bai_bytes: &[u8],
    ref_name: &str,
    beg: i32,
    end: i32,
) -> Result<Vec<BamRecord>> {
    let index = bai::parse_bai(bai_bytes)?;
    let header = inflate_header(bgzf_bytes)?;

    let ref_id = header
        .references
        .iter()
        .position(|r| r.name == ref_name)
        .ok_or_else(|| BamError::UnknownReference(ref_name.to_string()))?;

    let mut out = Vec::new();
    for chunk in index.query_chunks(ref_id, beg, end) {
        // Inflate only the blocks this chunk spans, then slice the record range.
        let (buf, end_base) = bgzf::inflate_region(bgzf_bytes, chunk.beg >> 16, chunk.end >> 16)?;
        let from = (chunk.beg & 0xffff) as usize;
        let to = end_base + (chunk.end & 0xffff) as usize;
        if from >= to || to > buf.len() {
            continue;
        }
        let mut cursor = Cursor::new(&buf[from..to]);
        while !cursor.eof() {
            let record = read_record(&mut cursor, &header.references)?;
            if record.overlaps(ref_name, beg, end) {
                out.push(record);
            }
        }
    }
    out.sort_by_key(|r| r.pos);
    Ok(out)
}
