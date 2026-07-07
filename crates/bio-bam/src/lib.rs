//! BGZF decompression and BAM parsing, with optional BAI-indexed region access.
//!
//! Kept in its own crate so `bio-core` stays dependency-free; only BAM support
//! pulls in a DEFLATE implementation (`miniz_oxide`). Sequential reads use
//! [`read_bam`]; [`read_bam_region`] uses a BAI index to fetch just the records
//! overlapping a locus.

pub mod bai;
pub mod bam;
pub mod bgzf;
pub mod error;
pub mod pileup;
pub mod varcall;
mod reader;

pub use bai::{parse_bai, Bai, Chunk};
pub use bam::{parse, parse_header, read_bam, read_bam_region, BamHeader, BamRecord, Reference};
pub use error::BamError;
pub use pileup::{pileup, PileupColumn};
pub use varcall::{call_variants, PileupVariant};
