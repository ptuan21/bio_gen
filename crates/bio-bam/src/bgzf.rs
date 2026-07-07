use miniz_oxide::inflate::decompress_to_vec;

use crate::error::{BamError, Result};

fn read_u16(data: &[u8], at: usize) -> u16 {
    u16::from_le_bytes([data[at], data[at + 1]])
}

/// Inflate only the blocks spanning compressed offsets `[start_coffset,
/// end_coffset]` (both block-aligned, as BAI virtual offsets always are).
///
/// This is the core of random access: for a region query we inflate just the
/// handful of blocks the index points at, never the whole file. Returns the
/// concatenated bytes and the position, within them, where the `end_coffset`
/// block begins — so the caller can slice out the exact record range.
pub fn inflate_region(input: &[u8], start_coffset: u64, end_coffset: u64) -> Result<(Vec<u8>, usize)> {
    let mut out = Vec::new();
    let mut pos = start_coffset as usize;
    let mut end_base: Option<usize> = None;

    while pos < input.len() {
        let at_end = pos as u64 == end_coffset;
        if at_end {
            end_base = Some(out.len());
        }
        pos += inflate_block(input, pos, &mut out)?;
        if at_end || pos as u64 > end_coffset {
            break;
        }
    }

    let end_base = end_base.unwrap_or(out.len());
    Ok((out, end_base))
}

/// Inflate one block starting at `pos`, append it to `out`, return its length.
pub(crate) fn inflate_block(input: &[u8], pos: usize, out: &mut Vec<u8>) -> Result<usize> {
    if pos + 12 > input.len() {
        return Err(BamError::Truncated);
    }
    if input[pos] != 0x1f || input[pos + 1] != 0x8b || input[pos + 3] & 0x04 == 0 {
        return Err(BamError::BadBlock);
    }

    let xlen = read_u16(input, pos + 10) as usize;
    let extra_start = pos + 12;
    let extra_end = extra_start + xlen;
    if extra_end > input.len() {
        return Err(BamError::Truncated);
    }

    let mut bsize: Option<usize> = None;
    let mut i = extra_start;
    while i + 4 <= extra_end {
        let slen = read_u16(input, i + 2) as usize;
        if input[i] == b'B' && input[i + 1] == b'C' && slen == 2 && i + 6 <= extra_end {
            bsize = Some(read_u16(input, i + 4) as usize);
        }
        i += 4 + slen;
    }
    let block_len = bsize.ok_or(BamError::MissingBlockSize)? + 1;
    if pos + block_len > input.len() || block_len < xlen + 26 {
        return Err(BamError::Truncated);
    }

    let cdata = &input[extra_end..pos + block_len - 8];
    let mut block = decompress_to_vec(cdata).map_err(|e| BamError::Inflate(e.to_string()))?;
    out.append(&mut block);
    Ok(block_len)
}

/// Decompress a BGZF stream (concatenated gzip blocks, as used by BAM) into the
/// raw uncompressed bytes. Blocks are located via the mandatory `BC` subfield.
pub fn decompress(input: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    let mut pos = 0;
    while pos < input.len() {
        pos += inflate_block(input, pos, &mut out)?;
    }
    Ok(out)
}
