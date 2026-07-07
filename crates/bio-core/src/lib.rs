//! Core engine for the Interactive Genomics Viewer Lite.
//!
//! Pure Rust, no external dependencies, so it compiles natively for tests and
//! to `wasm32` for the browser. Higher layers (`bio-wasm`, UI) build on top.

pub mod analysis;
pub mod error;
pub mod parser;
pub mod sequence;

pub use error::BioError;
pub use sequence::{SeqKind, Sequence, SeqRecord};

pub mod prelude {
    pub use crate::analysis::kmer::kmer_counts;
    pub use crate::analysis::orf::{find_orfs, Orf};
    pub use crate::analysis::restriction::{digest, find_sites, Enzyme, SiteHit, ENZYMES};
    pub use crate::analysis::search::{search, Match, Strand};
    pub use crate::analysis::stats::{gc_skew, stats, BaseCounts, SeqStats};
    pub use crate::analysis::translate::{
        codon_to_aa_with, point_mutation_effect, six_frames, translate, translate_frame,
        translate_with, FrameTranslation, GeneticCode, MutationEffect,
    };
    pub use crate::analysis::variant::{call_substitutions, Variant, VariantKind};
    pub use crate::error::BioError;
    pub use crate::parser::{FastaReader, FastqReader};
    pub use crate::sequence::{SeqKind, SeqRecord, Sequence};
}
