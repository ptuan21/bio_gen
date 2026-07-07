use miniz_oxide::inflate::decompress_to_vec;

use crate::error::{BamError, Result};

fn read_u16(data: &[u8], at: usize) -> u16 {
    u16::from_le_bytes([data[at], data[at + 1]])
}

/// Decompress a BGZF stream (concatenated gzip blocks, as used by BAM) into the
/// raw uncompressed bytes. Blocks are located via the mandatory `BC` subfield.
pub fn decompress(input: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    let mut pos = 0;

    while pos < input.len() {
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
        pos += block_len;
    }

    Ok(out)
}
