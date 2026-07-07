use crate::bgzf;
use crate::error::{BamError, Result};

const SEQ_NT: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
const CIGAR_OPS: &[u8; 9] = b"MIDNSHP=X";

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
    pub mapq: u8,
    pub cigar: String,
    pub seq: String,
    /// Phred qualities; empty when unavailable (`0xFF` fill).
    pub qual: Vec<u8>,
}

struct Cursor<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }

    fn eof(&self) -> bool {
        self.pos >= self.data.len()
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self.pos.checked_add(n).ok_or(BamError::Truncated)?;
        if end > self.data.len() {
            return Err(BamError::Truncated);
        }
        let slice = &self.data[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16> {
        let b = self.take(2)?;
        Ok(u16::from_le_bytes([b[0], b[1]]))
    }

    fn u32(&mut self) -> Result<u32> {
        let b = self.take(4)?;
        Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn i32(&mut self) -> Result<i32> {
        Ok(self.u32()? as i32)
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

fn decode_cigar(cursor: &mut Cursor, n_ops: usize) -> Result<String> {
    if n_ops == 0 {
        return Ok("*".to_string());
    }
    let mut cigar = String::new();
    for _ in 0..n_ops {
        let v = cursor.u32()?;
        cigar.push_str(&(v >> 4).to_string());
        cigar.push(CIGAR_OPS[(v & 0xf) as usize] as char);
    }
    Ok(cigar)
}

fn parse_records(cursor: &mut Cursor, refs: &[Reference]) -> Result<Vec<BamRecord>> {
    let mut records = Vec::new();
    while !cursor.eof() {
        let block_size = cursor.u32()? as usize;
        let start = cursor.pos;

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
        let cigar = decode_cigar(cursor, n_cigar)?;
        let seq = decode_seq(cursor.take(l_seq.div_ceil(2))?, l_seq);
        let qual_bytes = cursor.take(l_seq)?;
        let qual = if qual_bytes.first() == Some(&0xff) {
            Vec::new()
        } else {
            qual_bytes.to_vec()
        };

        // Skip any auxiliary tag data at the end of the record.
        let consumed = cursor.pos - start;
        if consumed < block_size {
            cursor.take(block_size - consumed)?;
        }

        let ref_name = usize::try_from(ref_id)
            .ok()
            .and_then(|i| refs.get(i))
            .map(|r| r.name.clone());

        records.push(BamRecord {
            name,
            flag,
            ref_name,
            pos,
            mapq,
            cigar,
            seq,
            qual,
        });
    }
    Ok(records)
}

/// Parse already-decompressed BAM bytes into a header and all records.
pub fn parse(data: &[u8]) -> Result<(BamHeader, Vec<BamRecord>)> {
    let mut cursor = Cursor::new(data);
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

    let records = parse_records(&mut cursor, &references)?;
    Ok((BamHeader { text, references }, records))
}

/// Decompress a BGZF-wrapped BAM file and parse it.
pub fn read_bam(bgzf_bytes: &[u8]) -> Result<(BamHeader, Vec<BamRecord>)> {
    parse(&bgzf::decompress(bgzf_bytes)?)
}
