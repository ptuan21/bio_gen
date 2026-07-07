//! BGZF decompression and sequential BAM parsing.
//!
//! Kept in its own crate so `bio-core` stays dependency-free; only BAM support
//! pulls in a DEFLATE implementation (`miniz_oxide`). Random access via BAI is
//! not implemented yet — records are read as a single forward pass.

pub mod bam;
pub mod bgzf;
pub mod error;

pub use bam::{parse, read_bam, BamHeader, BamRecord, Reference};
pub use error::BamError;
