//! Thin `wasm-bindgen` layer over `bio-core`. Each export takes plain strings
//! from JavaScript and returns JSON-friendly values via `serde-wasm-bindgen`.

use serde::Serialize;
use wasm_bindgen::prelude::*;

use bio_core::analysis::kmer::kmer_counts;
use bio_core::analysis::orf::find_orfs;
use bio_core::analysis::restriction::{digest, ENZYMES};
use bio_core::analysis::search::{search, Strand};
use bio_core::analysis::stats::{gc_skew, stats};
use bio_core::analysis::translate::{
    point_mutation_effect, six_frames, translate_with, GeneticCode, MutationEffect,
};
use bio_core::analysis::variant::call_substitutions;
use bio_core::crispr::{
    design_hdr, design_knockin, enzyme_by_name, find_guides, ENZYMES as CRISPR_ENZYMES,
};
use bio_core::vcf::{from_substitutions, write_vcf};
use bio_bam::varcall::pileup_variants_to_vcf;
use bio_core::parser::{parse_fasta_str, FastaStreamer as CoreFastaStreamer, RecordSummary};
use bio_core::sequence::{SeqKind, Sequence};

use bio_bam::bam::BamRecord;
use bio_bam::pileup::pileup;
use bio_bam::varcall::call_variants as call_pileup_variants;
use bio_bam::{read_bam, read_bam_region};

fn kind(is_rna: bool) -> SeqKind {
    if is_rna {
        SeqKind::Rna
    } else {
        SeqKind::Dna
    }
}

fn to_js_err<E: std::fmt::Display>(e: E) -> JsValue {
    JsValue::from_str(&e.to_string())
}

fn to_js<T: Serialize>(value: &T) -> Result<JsValue, JsValue> {
    serde_wasm_bindgen::to_value(value).map_err(to_js_err)
}

fn make_seq(seq: &str, is_rna: bool) -> Result<Sequence, JsValue> {
    Sequence::new(kind(is_rna), seq.as_bytes().to_vec()).map_err(to_js_err)
}

#[derive(Serialize)]
struct SummaryDto {
    id: String,
    description: String,
    length: usize,
    gc_content: f64,
}

fn summaries(records: Vec<RecordSummary>) -> Vec<SummaryDto> {
    records
        .into_iter()
        .map(|r| SummaryDto {
            id: r.id,
            description: r.description,
            length: r.length,
            gc_content: r.gc_content,
        })
        .collect()
}

/// Streaming FASTA parser for multi-gigabyte files: feed byte chunks with
/// `push`, then call `finish`. Memory stays flat; each call returns summaries
/// for the records that completed.
#[wasm_bindgen]
pub struct FastaStreamer {
    inner: CoreFastaStreamer,
}

#[wasm_bindgen]
impl FastaStreamer {
    #[wasm_bindgen(constructor)]
    #[allow(clippy::new_without_default)]
    pub fn new() -> Self {
        Self {
            inner: CoreFastaStreamer::new(),
        }
    }

    pub fn push(&mut self, chunk: &[u8]) -> Result<JsValue, JsValue> {
        to_js(&summaries(self.inner.push(chunk)))
    }

    pub fn finish(&mut self) -> Result<JsValue, JsValue> {
        to_js(&summaries(self.inner.finish()))
    }
}

#[derive(Serialize)]
struct RecordDto {
    id: String,
    description: String,
    length: usize,
    gc_content: f64,
    sequence: String,
}

#[derive(Serialize)]
struct StatsDto {
    length: usize,
    a: usize,
    c: usize,
    g: usize,
    t: usize,
    n: usize,
    other: usize,
    gc_content: f64,
}

#[derive(Serialize)]
struct MatchDto {
    start: usize,
    end: usize,
    strand: &'static str,
}

#[derive(Serialize)]
struct VariantDto {
    position: usize,
    reference: char,
    alternate: char,
    kind: &'static str,
}

#[derive(Serialize)]
struct FrameDto {
    frame: i8,
    protein: String,
}

#[derive(Serialize)]
struct OrfDto {
    start: usize,
    end: usize,
    strand: &'static str,
    frame: i8,
    protein: String,
    dna: String,
}

#[derive(Serialize)]
struct KmerDto {
    kmer: String,
    count: usize,
}

#[derive(Serialize)]
struct SiteDto {
    enzyme: &'static str,
    start: usize,
    cut: usize,
}

#[derive(Serialize)]
struct BamRecordDto {
    name: String,
    flag: u16,
    ref_name: Option<String>,
    pos: i32,
    ref_span: u32,
    mapq: u8,
    cigar: String,
    seq: String,
}

impl From<BamRecord> for BamRecordDto {
    fn from(r: BamRecord) -> Self {
        BamRecordDto {
            name: r.name,
            flag: r.flag,
            ref_name: r.ref_name,
            pos: r.pos,
            ref_span: r.ref_span,
            mapq: r.mapq,
            cigar: r.cigar,
            seq: r.seq,
        }
    }
}

#[derive(Serialize)]
struct BamResultDto {
    references: Vec<String>,
    record_count: usize,
    records: Vec<BamRecordDto>,
}

#[derive(Serialize)]
struct PileupVariantDto {
    ref_pos: i32,
    reference: char,
    alternate: char,
    depth: u32,
    alt_count: u32,
    allele_freq: f64,
    alt_fwd: u32,
    alt_rev: u32,
}

#[derive(Serialize)]
struct PileupDto {
    ref_pos: i32,
    depth: u32,
    a: u32,
    c: u32,
    g: u32,
    t: u32,
    n: u32,
    del: u32,
    consensus: Option<char>,
}

#[derive(Serialize)]
struct EnzymeDto {
    name: &'static str,
    pam: &'static str,
    spacer_len: usize,
}

#[derive(Serialize)]
struct GuideDto {
    strand: &'static str,
    start: usize,
    end: usize,
    spacer: String,
    pam: String,
    cut_site: usize,
    gc: f64,
    score: u32,
    off_targets: u32,
}

#[derive(Serialize)]
struct HdrDto {
    cut_site: usize,
    arm_len: usize,
    left_arm: String,
    insert: String,
    right_arm: String,
    donor: String,
}

#[derive(Serialize)]
struct EditDto {
    pos: usize,
    from: char,
    to: char,
}

#[derive(Serialize)]
struct KnockinDto {
    cut_site: usize,
    left_arm: String,
    insert: String,
    right_arm: String,
    donor: String,
    pam_edit: Option<EditDto>,
}

#[derive(Serialize)]
struct EffectDto {
    kind: &'static str,
    from: Option<char>,
    to: Option<char>,
    residue: Option<char>,
}

#[wasm_bindgen]
pub fn parse_fasta(input: &str, is_rna: bool) -> Result<JsValue, JsValue> {
    let records = parse_fasta_str(input, kind(is_rna)).map_err(to_js_err)?;
    let dto: Vec<RecordDto> = records
        .into_iter()
        .map(|r| {
            let gc = stats(&r.sequence).gc_content;
            RecordDto {
                id: r.id,
                description: r.description,
                length: r.sequence.len(),
                gc_content: gc,
                sequence: r.sequence.to_string(),
            }
        })
        .collect();
    to_js(&dto)
}

#[wasm_bindgen]
pub fn sequence_stats(seq: &str, is_rna: bool) -> Result<JsValue, JsValue> {
    let s = stats(&make_seq(seq, is_rna)?);
    to_js(&StatsDto {
        length: s.length,
        a: s.counts.a,
        c: s.counts.c,
        g: s.counts.g,
        t: s.counts.t,
        n: s.counts.n,
        other: s.counts.other,
        gc_content: s.gc_content,
    })
}

#[wasm_bindgen]
pub fn reverse_complement(seq: &str, is_rna: bool) -> Result<String, JsValue> {
    Ok(make_seq(seq, is_rna)?.reverse_complement().to_string())
}

#[wasm_bindgen]
pub fn transcribe(seq: &str, is_rna: bool) -> Result<String, JsValue> {
    Ok(make_seq(seq, is_rna)?.transcribe().to_string())
}

#[wasm_bindgen]
pub fn translate_seq(seq: &str, is_rna: bool, mito: bool) -> Result<String, JsValue> {
    let code = if mito {
        GeneticCode::VertebrateMito
    } else {
        GeneticCode::Standard
    };
    Ok(translate_with(&make_seq(seq, is_rna)?, code))
}

#[wasm_bindgen]
pub fn gc_skew_windows(
    seq: &str,
    window: usize,
    step: usize,
    is_rna: bool,
) -> Result<JsValue, JsValue> {
    let values = gc_skew(&make_seq(seq, is_rna)?, window, step);
    to_js(&values)
}

#[wasm_bindgen]
pub fn restriction_digest(seq: &str, is_rna: bool) -> Result<JsValue, JsValue> {
    let seq = make_seq(seq, is_rna)?;
    let dto: Vec<SiteDto> = digest(&seq, ENZYMES)
        .into_iter()
        .map(|h| SiteDto {
            enzyme: h.enzyme,
            start: h.start,
            cut: h.cut,
        })
        .collect();
    to_js(&dto)
}

#[wasm_bindgen]
pub fn search_motif(
    seq: &str,
    pattern: &str,
    is_rna: bool,
    both_strands: bool,
) -> Result<JsValue, JsValue> {
    let seq = make_seq(seq, is_rna)?;
    let dto: Vec<MatchDto> = search(&seq, pattern, both_strands)
        .into_iter()
        .map(|m| MatchDto {
            start: m.start,
            end: m.end,
            strand: match m.strand {
                Strand::Forward => "forward",
                Strand::Reverse => "reverse",
            },
        })
        .collect();
    to_js(&dto)
}

#[wasm_bindgen]
pub fn six_frame_translation(seq: &str, is_rna: bool) -> Result<JsValue, JsValue> {
    let seq = make_seq(seq, is_rna)?;
    let dto: Vec<FrameDto> = six_frames(&seq)
        .into_iter()
        .map(|f| FrameDto {
            frame: f.frame,
            protein: f.protein,
        })
        .collect();
    to_js(&dto)
}

#[wasm_bindgen]
pub fn find_open_reading_frames(
    seq: &str,
    min_aa: usize,
    is_rna: bool,
) -> Result<JsValue, JsValue> {
    let seq = make_seq(seq, is_rna)?;
    let dto: Vec<OrfDto> = find_orfs(&seq, min_aa)
        .into_iter()
        .map(|o| OrfDto {
            start: o.start,
            end: o.end,
            strand: match o.strand {
                Strand::Forward => "forward",
                Strand::Reverse => "reverse",
            },
            frame: o.frame,
            protein: o.protein,
            dna: o.dna,
        })
        .collect();
    to_js(&dto)
}

#[wasm_bindgen]
pub fn count_kmers(seq: &str, k: usize, is_rna: bool) -> Result<JsValue, JsValue> {
    let seq = make_seq(seq, is_rna)?;
    let dto: Vec<KmerDto> = kmer_counts(&seq, k)
        .into_iter()
        .map(|(kmer, count)| KmerDto { kmer, count })
        .collect();
    to_js(&dto)
}

#[wasm_bindgen]
pub fn call_variants(reference: &str, sample: &str, is_rna: bool) -> Result<JsValue, JsValue> {
    let reference = make_seq(reference, is_rna)?;
    let sample = make_seq(sample, is_rna)?;
    let variants = call_substitutions(&reference, &sample).map_err(to_js_err)?;
    let dto: Vec<VariantDto> = variants
        .into_iter()
        .map(|v| VariantDto {
            position: v.position,
            reference: v.reference as char,
            alternate: v.alternate as char,
            kind: "substitution",
        })
        .collect();
    to_js(&dto)
}

#[wasm_bindgen]
pub fn mutation_effect(seq: &str, pos: usize, alt: char, is_rna: bool) -> Result<JsValue, JsValue> {
    let seq = make_seq(seq, is_rna)?;
    let effect = point_mutation_effect(&seq, pos, alt as u8).map_err(to_js_err)?;
    let dto = match effect {
        MutationEffect::Silent { residue } => EffectDto {
            kind: "silent",
            from: None,
            to: None,
            residue: Some(residue),
        },
        MutationEffect::Missense { from, to } => EffectDto {
            kind: "missense",
            from: Some(from),
            to: Some(to),
            residue: None,
        },
        MutationEffect::Nonsense { from } => EffectDto {
            kind: "nonsense",
            from: Some(from),
            to: None,
            residue: None,
        },
        MutationEffect::StopLost { to } => EffectDto {
            kind: "stop_lost",
            from: None,
            to: Some(to),
            residue: None,
        },
    };
    to_js(&dto)
}

/// VCF from alignment-free substitution calls between `reference` and `sample`.
#[wasm_bindgen]
pub fn substitutions_vcf(
    reference: &str,
    sample: &str,
    chrom: &str,
    is_rna: bool,
) -> Result<String, JsValue> {
    let reference = make_seq(reference, is_rna)?;
    let sample = make_seq(sample, is_rna)?;
    let variants = call_substitutions(&reference, &sample).map_err(to_js_err)?;
    Ok(write_vcf(&from_substitutions(chrom, &variants)))
}

/// VCF from a BAM region pileup: region query, quality-filtered pileup, SNV
/// calling (with strand-bias filter), rendered with DP/AF/SB INFO.
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn pileup_variants_vcf(
    bam: &[u8],
    bai: &[u8],
    ref_name: &str,
    beg: i32,
    end: i32,
    reference: &str,
    ref_offset: i32,
    min_depth: u32,
    min_freq: f64,
    min_qual: u8,
    min_strand_frac: f64,
    is_rna: bool,
) -> Result<String, JsValue> {
    let records = read_bam_region(bam, bai, ref_name, beg, end).map_err(to_js_err)?;
    let columns = pileup(&records, beg, end, min_qual);
    let reference = make_seq(reference, is_rna)?;
    let variants =
        call_pileup_variants(&columns, &reference, ref_offset, min_depth, min_freq, min_strand_frac);
    Ok(pileup_variants_to_vcf(ref_name, &variants))
}

#[wasm_bindgen]
pub fn crispr_enzymes() -> Result<JsValue, JsValue> {
    let dto: Vec<EnzymeDto> = CRISPR_ENZYMES
        .iter()
        .map(|e| EnzymeDto { name: e.name, pam: e.pam, spacer_len: e.spacer_len })
        .collect();
    to_js(&dto)
}

/// Find CRISPR guides for `enzyme` (name, e.g. "SpCas9"), scored and annotated
/// with in-sequence off-target counts (Hamming distance <= `max_mismatch`).
#[wasm_bindgen]
pub fn crispr_guides(
    seq: &str,
    enzyme: &str,
    max_mismatch: usize,
    is_rna: bool,
) -> Result<JsValue, JsValue> {
    let enzyme = enzyme_by_name(enzyme).ok_or_else(|| JsValue::from_str("unknown enzyme"))?;
    let seq = make_seq(seq, is_rna)?;
    let dto: Vec<GuideDto> = find_guides(&seq, enzyme, max_mismatch)
        .into_iter()
        .map(|g| GuideDto {
            strand: match g.strand {
                Strand::Forward => "forward",
                Strand::Reverse => "reverse",
            },
            start: g.start,
            end: g.end,
            spacer: g.spacer,
            pam: g.pam,
            cut_site: g.cut_site,
            gc: g.gc,
            score: g.score,
            off_targets: g.off_targets,
        })
        .collect();
    to_js(&dto)
}

/// Design a knock-in donor around the guide at `guide_index` (0 = best). With
/// `disrupt_pam`, a single PAM base is changed in the donor to block re-cutting.
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn crispr_knockin(
    seq: &str,
    enzyme: &str,
    guide_index: usize,
    insert: &str,
    arm_len: usize,
    disrupt_pam: bool,
    is_rna: bool,
) -> Result<JsValue, JsValue> {
    let enzyme = enzyme_by_name(enzyme).ok_or_else(|| JsValue::from_str("unknown enzyme"))?;
    let seq = make_seq(seq, is_rna)?;
    let guides = find_guides(&seq, enzyme, 0);
    let guide = guides
        .get(guide_index)
        .ok_or_else(|| JsValue::from_str("no guide at that index"))?;
    let design = design_knockin(&seq, enzyme, guide, insert, arm_len, disrupt_pam).map_err(to_js_err)?;
    to_js(&KnockinDto {
        cut_site: design.template.cut_site,
        left_arm: design.template.left_arm,
        insert: design.template.insert,
        right_arm: design.template.right_arm,
        donor: design.template.donor,
        pam_edit: design.pam_edit.map(|e| EditDto { pos: e.pos, from: e.from, to: e.to }),
    })
}

/// Design an HDR donor around `cut_site`: homology arms of `arm_len` flanking
/// the cut with `insert` placed at it.
#[wasm_bindgen]
pub fn crispr_hdr(
    reference: &str,
    cut_site: usize,
    insert: &str,
    arm_len: usize,
    is_rna: bool,
) -> Result<JsValue, JsValue> {
    let reference = make_seq(reference, is_rna)?;
    let hdr = design_hdr(&reference, cut_site, insert, arm_len).map_err(to_js_err)?;
    to_js(&HdrDto {
        cut_site: hdr.cut_site,
        arm_len: hdr.arm_len,
        left_arm: hdr.left_arm,
        insert: hdr.insert,
        right_arm: hdr.right_arm,
        donor: hdr.donor,
    })
}

/// Parse a BGZF-compressed BAM file. `max_records` caps how many alignments are
/// returned to JavaScript (0 = all); `record_count` always reflects the total.
#[wasm_bindgen]
pub fn parse_bam(bytes: &[u8], max_records: usize) -> Result<JsValue, JsValue> {
    let (header, records) = read_bam(bytes).map_err(to_js_err)?;
    let limit = if max_records == 0 {
        records.len()
    } else {
        max_records.min(records.len())
    };
    let dto = BamResultDto {
        references: header.references.into_iter().map(|r| r.name).collect(),
        record_count: records.len(),
        records: records.into_iter().take(limit).map(Into::into).collect(),
    };
    to_js(&dto)
}

/// Fetch BAM records overlapping `ref_name:[beg, end)` using a BAI index.
#[wasm_bindgen]
pub fn parse_bam_region(
    bam: &[u8],
    bai: &[u8],
    ref_name: &str,
    beg: i32,
    end: i32,
) -> Result<JsValue, JsValue> {
    let records = read_bam_region(bam, bai, ref_name, beg, end).map_err(to_js_err)?;
    let dto: Vec<BamRecordDto> = records.into_iter().map(Into::into).collect();
    to_js(&dto)
}

/// Coverage pileup over `ref_name:[beg, end)`: region query plus per-position
/// depth, base counts and consensus base.
#[wasm_bindgen]
pub fn bam_pileup(
    bam: &[u8],
    bai: &[u8],
    ref_name: &str,
    beg: i32,
    end: i32,
    min_qual: u8,
) -> Result<JsValue, JsValue> {
    let records = read_bam_region(bam, bai, ref_name, beg, end).map_err(to_js_err)?;
    let dto: Vec<PileupDto> = pileup(&records, beg, end, min_qual)
        .into_iter()
        .map(|c| PileupDto {
            ref_pos: c.ref_pos,
            depth: c.depth,
            a: c.a,
            c: c.c,
            g: c.g,
            t: c.t,
            n: c.n,
            del: c.del,
            consensus: c.consensus().map(|(base, _)| base),
        })
        .collect();
    to_js(&dto)
}

/// Call SNVs in `ref_name:[beg, end)` by piling up aligned reads and comparing
/// the consensus against `reference` (whose first base is at `ref_offset`).
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn call_variants_pileup(
    bam: &[u8],
    bai: &[u8],
    ref_name: &str,
    beg: i32,
    end: i32,
    reference: &str,
    ref_offset: i32,
    min_depth: u32,
    min_freq: f64,
    min_qual: u8,
    min_strand_frac: f64,
    is_rna: bool,
) -> Result<JsValue, JsValue> {
    let records = read_bam_region(bam, bai, ref_name, beg, end).map_err(to_js_err)?;
    let columns = pileup(&records, beg, end, min_qual);
    let reference = make_seq(reference, is_rna)?;
    let dto: Vec<PileupVariantDto> =
        call_pileup_variants(&columns, &reference, ref_offset, min_depth, min_freq, min_strand_frac)
            .into_iter()
            .map(|v| PileupVariantDto {
                ref_pos: v.ref_pos,
                reference: v.reference,
                alternate: v.alternate,
                depth: v.depth,
                alt_count: v.alt_count,
                allele_freq: v.allele_freq,
                alt_fwd: v.alt_fwd,
                alt_rev: v.alt_rev,
            })
            .collect();
    to_js(&dto)
}
