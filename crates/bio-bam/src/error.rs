use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BamError {
    /// Input ended in the middle of a block or record.
    Truncated,
    /// A BGZF block header was not a valid gzip member.
    BadBlock,
    /// The BGZF `BC` extra subfield giving the block size was missing.
    MissingBlockSize,
    /// DEFLATE inflation failed.
    Inflate(String),
    /// The BAM magic bytes `BAM\1` were absent.
    BadMagic,
}

impl fmt::Display for BamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BamError::Truncated => write!(f, "unexpected end of input"),
            BamError::BadBlock => write!(f, "invalid BGZF block header"),
            BamError::MissingBlockSize => write!(f, "BGZF block missing BC subfield"),
            BamError::Inflate(msg) => write!(f, "inflate error: {msg}"),
            BamError::BadMagic => write!(f, "not a BAM stream (bad magic)"),
        }
    }
}

impl std::error::Error for BamError {}

pub type Result<T> = std::result::Result<T, BamError>;
