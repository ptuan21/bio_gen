use std::io::{BufRead, Lines};

use crate::error::BioError;
use crate::sequence::{split_header, SeqKind, SeqRecord, Sequence};

/// Streaming FASTA reader: yields one `SeqRecord` at a time so multi-gigabyte
/// files never need to be held in memory all at once.
pub struct FastaReader<R: BufRead> {
    lines: std::iter::Enumerate<Lines<R>>,
    kind: SeqKind,
    pending_header: Option<String>,
    finished: bool,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(reader: R, kind: SeqKind) -> Self {
        Self {
            lines: reader.lines().enumerate(),
            kind,
            pending_header: None,
            finished: false,
        }
    }

    fn fail(&mut self, err: BioError) -> Option<Result<SeqRecord, BioError>> {
        self.finished = true;
        Some(Err(err))
    }

    fn read_header(&mut self) -> Result<Option<String>, BioError> {
        if let Some(h) = self.pending_header.take() {
            return Ok(Some(h));
        }
        loop {
            match self.lines.next() {
                None => return Ok(None),
                Some((_, Err(e))) => return Err(BioError::Io(e.to_string())),
                Some((idx, Ok(line))) => {
                    let trimmed = line.trim();
                    if trimmed.is_empty() {
                        continue;
                    }
                    match trimmed.strip_prefix('>') {
                        Some(rest) => return Ok(Some(rest.to_string())),
                        None => {
                            return Err(BioError::MalformedFasta {
                                line: idx + 1,
                                reason: "expected '>' at record start".to_string(),
                            })
                        }
                    }
                }
            }
        }
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = Result<SeqRecord, BioError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        let header = match self.read_header() {
            Ok(Some(h)) => h,
            Ok(None) => {
                self.finished = true;
                return None;
            }
            Err(e) => return self.fail(e),
        };

        let mut data: Vec<u8> = Vec::new();
        loop {
            match self.lines.next() {
                None => {
                    self.finished = true;
                    break;
                }
                Some((_, Err(e))) => return self.fail(BioError::Io(e.to_string())),
                Some((_, Ok(line))) => {
                    let trimmed = line.trim();
                    if trimmed.is_empty() {
                        continue;
                    }
                    if let Some(rest) = trimmed.strip_prefix('>') {
                        self.pending_header = Some(rest.to_string());
                        break;
                    }
                    data.extend_from_slice(trimmed.as_bytes());
                }
            }
        }

        let (id, description) = split_header(&header);
        match Sequence::new(self.kind, data) {
            Ok(sequence) => Some(Ok(SeqRecord {
                id,
                description,
                sequence,
                quality: None,
            })),
            Err(e) => self.fail(e),
        }
    }
}

/// Parse a whole FASTA string into records (convenient for small inputs).
pub fn parse_fasta_str(input: &str, kind: SeqKind) -> Result<Vec<SeqRecord>, BioError> {
    if input.trim().is_empty() {
        return Err(BioError::EmptyInput);
    }
    FastaReader::new(input.as_bytes(), kind).collect()
}
