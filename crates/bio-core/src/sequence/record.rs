use super::Sequence;

/// A single parsed entry from a FASTA/FASTQ file.
#[derive(Debug, Clone, PartialEq)]
pub struct SeqRecord {
    pub id: String,
    pub description: String,
    pub sequence: Sequence,
    /// Phred quality scores (FASTQ only); `None` for FASTA.
    pub quality: Option<Vec<u8>>,
}

impl SeqRecord {
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// Split a FASTA/FASTQ header into `(id, description)` at the first whitespace.
pub fn split_header(header: &str) -> (String, String) {
    let header = header.trim();
    match header.split_once(char::is_whitespace) {
        Some((id, desc)) => (id.to_string(), desc.trim().to_string()),
        None => (header.to_string(), String::new()),
    }
}
