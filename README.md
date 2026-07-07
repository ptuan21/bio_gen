# bio_gen — Interactive Genomics Viewer Lite (core)

A client-side engine for parsing and analysing DNA/RNA sequences, written in
Rust and compiled to WebAssembly so large gene files can be processed in the
browser with no server and no CLI install.

This repository currently ships the **data-processing core**. The UI is
intentionally deferred — the engine exposes a clean, typed API that any
frontend (React, Svelte, vanilla) can build on later.

## Structure

```
bio_gen/
├── Cargo.toml                  workspace
├── crates/
│   ├── bio-core/               zero-dependency Rust engine (native-testable)
│   │   ├── src/
│   │   │   ├── error.rs         BioError type
│   │   │   ├── sequence/        Sequence model, IUPAC alphabet, records
│   │   │   ├── parser/          streaming FASTA / FASTQ readers
│   │   │   └── analysis/        stats, search, variants, translation,
│   │   │                        ORF finder, k-mer, restriction sites
│   │   ├── examples/pipeline.rs end-to-end demo
│   │   └── tests/               integration tests (incl. edge cases)
│   ├── bio-bam/                BGZF decompression + BAM parsing (miniz_oxide)
│   └── bio-wasm/               thin wasm-bindgen bindings over the two crates
├── web/index.html             browser demo (loads web/pkg)
└── scripts/
    ├── build-wasm.sh           build the wasm package into web/pkg
    └── serve.sh                static server for the demo
```

`bio-core` has **no external dependencies**, so it compiles fast, tests
natively, and keeps the WASM binary small. BAM support lives in `bio-bam` (the
only crate needing a DEFLATE implementation), so the core stays dependency-free.

## Features

- Streaming **FASTA** and **FASTQ** parsers (record-by-record, memory-friendly)
- DNA/RNA `Sequence` with validation, IUPAC ambiguity codes, complement,
  reverse-complement, transcription
- Base composition & **GC content**
- **Motif search** with IUPAC wildcards, both strands, forward coordinates
- **Substitution calling** between reference and sample
- Codon **translation**, **6-frame translation**, and single-base **mutation
  effect** (silent / missense / nonsense / stop-lost)
- **ORF finder** (both strands, forward-strand coordinates, protein + DNA)
- **k-mer** frequency counting
- **GC-skew** over sliding windows
- **Restriction digest** with a built-in enzyme panel (EcoRI, BamHI, …)
- **Vertebrate mitochondrial** genetic code (alongside the standard code)
- **BGZF** decompression and sequential **BAM** parsing (no BAI random access yet)

## Develop

```bash
# native build + tests
cargo test

# run the demo
cargo run -p bio-core --example pipeline

# lint
cargo clippy --all-targets -- -D warnings

# build the browser package (needs: cargo install wasm-pack)
./scripts/build-wasm.sh
```

## Browser demo

```bash
./scripts/build-wasm.sh   # generates web/pkg
./scripts/serve.sh        # http://localhost:8000/web/
```

The demo parses FASTA/BAM, searches motifs (highlighting hits on the sequence),
translates, finds ORFs, runs restriction digests and more — all client-side.

## Using it from JavaScript

After `build-wasm.sh`, import the generated module:

```js
import init, { parse_fasta, search_motif, mutation_effect } from "./crates/bio-wasm/pkg/bio_wasm.js";

await init();
const records = parse_fasta(">g1\nATGGCCTAA\n", /* is_rna */ false);
const hits = search_motif("ATGGCCTAA", "GCC", false, true);
const effect = mutation_effect("ATGGCCTAA", 1, "C", false); // { kind: "missense", from: "M", to: "T" }
```

Exports: `parse_fasta`, `sequence_stats`, `reverse_complement`, `transcribe`,
`translate_seq` (with `mito` flag), `six_frame_translation`,
`find_open_reading_frames`, `count_kmers`, `search_motif`, `gc_skew_windows`,
`restriction_digest`, `call_variants`, `mutation_effect`, `parse_bam`.

## Roadmap

- [x] BGZF decompression + sequential BAM parsing
- [ ] BAI index for BAM random access (region queries)
- [ ] Frontend viewer (canvas/WebGL track rendering) on top of the WASM API
- [ ] Proper alignment for indel-aware variant calling
