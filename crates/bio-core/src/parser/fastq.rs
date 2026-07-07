use std::io::{BufRead, Lines};

use crate::error::BioError;
use crate::sequence::{split_header, SeqKind, SeqRecord, Sequence};

const PHRED_OFFSET: u8 = 33;

/// Streaming FASTQ reader for the standard 4-line record layout.
pub struct FastqReader<R: BufRead> {
    lines: std::iter::Enumerate<Lines<R>>,
    kind: SeqKind,
    finished: bool,
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R, kind: SeqKind) -> Self {
        Self {
            lines: reader.lines().enumerate(),
            kind,
            finished: false,
        }
    }

    fn fail(&mut self, err: BioError) -> Option<Result<SeqRecord, BioError>> {
        self.finished = true;
        Some(Err(err))
    }

    /// Next non-empty line, or `None` at clean end of input.
    fn next_line(&mut self) -> Result<Option<(usize, String)>, BioError> {
        loop {
            match self.lines.next() {
                None => return Ok(None),
                Some((_, Err(e))) => return Err(BioError::Io(e.to_string())),
                Some((idx, Ok(line))) => {
                    if line.trim().is_empty() {
                        continue;
                    }
                    return Ok(Some((idx + 1, line)));
                }
            }
        }
    }
}

impl<R: BufRead> Iterator for FastqReader<R> {
    type Item = Result<SeqRecord, BioError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        let (line_no, header) = match self.next_line() {
            Ok(Some(x)) => x,
            Ok(None) => {
                self.finished = true;
                return None;
            }
            Err(e) => return self.fail(e),
        };

        let header = match header.trim().strip_prefix('@') {
            Some(rest) => rest.to_string(),
            None => {
                return self.fail(BioError::MalformedFastq {
                    line: line_no,
                    reason: "expected '@' at record start".to_string(),
                })
            }
        };

        let seq_line = match self.next_line() {
            Ok(Some((_, l))) => l.trim().to_string(),
            Ok(None) => {
                return self.fail(BioError::MalformedFastq {
                    line: line_no,
                    reason: "missing sequence line".to_string(),
                })
            }
            Err(e) => return self.fail(e),
        };

        match self.next_line() {
            Ok(Some((n, l))) if !l.trim_start().starts_with('+') => {
                return self.fail(BioError::MalformedFastq {
                    line: n,
                    reason: "expected '+' separator".to_string(),
                })
            }
            Ok(None) => {
                return self.fail(BioError::MalformedFastq {
                    line: line_no,
                    reason: "missing '+' separator".to_string(),
                })
            }
            Err(e) => return self.fail(e),
            _ => {}
        }

        let qual_line = match self.next_line() {
            Ok(Some((_, l))) => l.trim_end().to_string(),
            Ok(None) => {
                return self.fail(BioError::MalformedFastq {
                    line: line_no,
                    reason: "missing quality line".to_string(),
                })
            }
            Err(e) => return self.fail(e),
        };

        if qual_line.len() != seq_line.len() {
            return self.fail(BioError::LengthMismatch {
                expected: seq_line.len(),
                found: qual_line.len(),
            });
        }

        let quality = qual_line.bytes().map(|b| b.wrapping_sub(PHRED_OFFSET)).collect();
        let (id, description) = split_header(&header);
        match Sequence::new(self.kind, seq_line.into_bytes()) {
            Ok(sequence) => Some(Ok(SeqRecord {
                id,
                description,
                sequence,
                quality: Some(quality),
            })),
            Err(e) => self.fail(e),
        }
    }
}

pub fn parse_fastq_str(input: &str, kind: SeqKind) -> Result<Vec<SeqRecord>, BioError> {
    if input.trim().is_empty() {
        return Err(BioError::EmptyInput);
    }
    FastqReader::new(input.as_bytes(), kind).collect()
}
