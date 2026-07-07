//! Push-based FASTA parser for streaming very large files chunk by chunk.
//!
//! Bytes arrive via [`FastaStreamer::push`] in arbitrary-sized chunks (e.g. from
//! `File.slice` in the browser). Only the current line and the running counts of
//! the current record are held, so memory stays flat regardless of file size.
//! Per-record summaries are emitted as each record completes.

use crate::analysis::stats::BaseCounts;
use crate::sequence::split_header;

#[derive(Debug, Clone, PartialEq)]
pub struct RecordSummary {
    pub id: String,
    pub description: String,
    pub length: usize,
    pub counts: BaseCounts,
    pub gc_content: f64,
}

struct Active {
    id: String,
    description: String,
    length: usize,
    counts: BaseCounts,
}

impl Active {
    fn into_summary(self) -> RecordSummary {
        let known = self.counts.a + self.counts.c + self.counts.g + self.counts.t;
        let gc_content = if known == 0 {
            0.0
        } else {
            (self.counts.g + self.counts.c) as f64 / known as f64
        };
        RecordSummary {
            id: self.id,
            description: self.description,
            length: self.length,
            counts: self.counts,
            gc_content,
        }
    }
}

#[derive(Default)]
pub struct FastaStreamer {
    line: Vec<u8>,
    current: Option<Active>,
}

impl FastaStreamer {
    pub fn new() -> Self {
        Self::default()
    }

    /// Feed a chunk of bytes; returns summaries for any records that completed.
    pub fn push(&mut self, chunk: &[u8]) -> Vec<RecordSummary> {
        let mut out = Vec::new();
        for &b in chunk {
            match b {
                b'\n' => {
                    let line = std::mem::take(&mut self.line);
                    self.consume_line(&line, &mut out);
                    self.line = line;
                    self.line.clear();
                }
                b'\r' => {}
                _ => self.line.push(b),
            }
        }
        out
    }

    /// Flush the final line and record after the last chunk.
    pub fn finish(&mut self) -> Vec<RecordSummary> {
        let mut out = Vec::new();
        if !self.line.is_empty() {
            let line = std::mem::take(&mut self.line);
            self.consume_line(&line, &mut out);
        }
        if let Some(active) = self.current.take() {
            out.push(active.into_summary());
        }
        out
    }

    fn consume_line(&mut self, line: &[u8], out: &mut Vec<RecordSummary>) {
        if line.is_empty() {
            return;
        }
        if line[0] == b'>' {
            if let Some(active) = self.current.take() {
                out.push(active.into_summary());
            }
            let header = String::from_utf8_lossy(&line[1..]).into_owned();
            let (id, description) = split_header(&header);
            self.current = Some(Active {
                id,
                description,
                length: 0,
                counts: BaseCounts::default(),
            });
        } else if let Some(active) = self.current.as_mut() {
            for &b in line {
                active.length += 1;
                match b.to_ascii_uppercase() {
                    b'A' => active.counts.a += 1,
                    b'C' => active.counts.c += 1,
                    b'G' => active.counts.g += 1,
                    b'T' | b'U' => active.counts.t += 1,
                    b'N' => active.counts.n += 1,
                    _ => active.counts.other += 1,
                }
            }
        }
    }
}
