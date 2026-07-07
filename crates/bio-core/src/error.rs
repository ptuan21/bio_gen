use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BioError {
    EmptyInput,
    Io(String),
    MalformedFasta { line: usize, reason: String },
    MalformedFastq { line: usize, reason: String },
    InvalidNucleotide { symbol: char, position: usize },
    LengthMismatch { expected: usize, found: usize },
    OutOfBounds { position: usize, length: usize },
}

impl fmt::Display for BioError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BioError::EmptyInput => write!(f, "input is empty"),
            BioError::Io(msg) => write!(f, "io error: {msg}"),
            BioError::MalformedFasta { line, reason } => {
                write!(f, "malformed FASTA at line {line}: {reason}")
            }
            BioError::MalformedFastq { line, reason } => {
                write!(f, "malformed FASTQ at line {line}: {reason}")
            }
            BioError::InvalidNucleotide { symbol, position } => {
                write!(f, "invalid nucleotide '{symbol}' at position {position}")
            }
            BioError::LengthMismatch { expected, found } => {
                write!(f, "length mismatch: expected {expected}, found {found}")
            }
            BioError::OutOfBounds { position, length } => {
                write!(f, "position {position} out of bounds for length {length}")
            }
        }
    }
}

impl std::error::Error for BioError {}

pub type Result<T> = std::result::Result<T, BioError>;
