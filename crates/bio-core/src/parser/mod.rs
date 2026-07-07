mod fasta;
mod fastq;
mod stream;

pub use fasta::{parse_fasta_str, FastaReader};
pub use fastq::{parse_fastq_str, FastqReader};
pub use stream::{FastaStreamer, RecordSummary};
