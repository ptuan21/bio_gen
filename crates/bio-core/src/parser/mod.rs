mod fasta;
mod fastq;

pub use fasta::{parse_fasta_str, FastaReader};
pub use fastq::{parse_fastq_str, FastqReader};
