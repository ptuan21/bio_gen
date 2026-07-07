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
use bio_core::parser::parse_fasta_str;
use bio_core::sequence::{SeqKind, Sequence};

use bio_bam::read_bam;

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
    mapq: u8,
    cigar: String,
    seq: String,
}

#[derive(Serialize)]
struct BamResultDto {
    references: Vec<String>,
    record_count: usize,
    records: Vec<BamRecordDto>,
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
        records: records
            .into_iter()
            .take(limit)
            .map(|r| BamRecordDto {
                name: r.name,
                flag: r.flag,
                ref_name: r.ref_name,
                pos: r.pos,
                mapq: r.mapq,
                cigar: r.cigar,
                seq: r.seq,
            })
            .collect(),
    };
    to_js(&dto)
}
