//! End-to-end demo of the core engine. Run with:
//!   cargo run -p bio-core --example pipeline

use bio_core::analysis::kmer::kmer_counts;
use bio_core::analysis::stats::stats;
use bio_core::prelude::*;

fn main() {
    let fasta = ">gene1 demo gene\nATGGCCTACGGGTAA\n>gene2\nACGTACGTAC\n";
    let records = FastaReader::new(fasta.as_bytes(), SeqKind::Dna)
        .collect::<Result<Vec<_>, _>>()
        .expect("valid FASTA");

    for rec in &records {
        let s = stats(&rec.sequence);
        println!(
            "{:<6} len={} gc={:.1}% seq={}",
            rec.id,
            s.length,
            s.gc_content * 100.0,
            rec.sequence
        );
    }

    let gene = &records[0].sequence;
    println!("\nreverse-complement: {}", gene.reverse_complement());
    println!("protein:            {}", translate(gene));

    println!("\nmotif 'GGG' (both strands):");
    for m in search(gene, "GGG", true) {
        println!("  {:?} at {}..{}", m.strand, m.start, m.end);
    }

    // Substitution at position 1 (A->C) turns start codon ATG(M) into ACG(T).
    match point_mutation_effect(gene, 1, b'C') {
        Ok(effect) => println!("\nmutation @1 A>C: {:?}", effect),
        Err(e) => eprintln!("error: {e}"),
    }

    println!("\nopen reading frames (>=2 aa):");
    for orf in find_orfs(gene, 2) {
        println!(
            "  {:?} frame {:+} {}..{} -> {}",
            orf.strand, orf.frame, orf.start, orf.end, orf.protein
        );
    }

    println!("\ntop 3-mers:");
    for (kmer, n) in kmer_counts(gene, 3).into_iter().take(3) {
        println!("  {kmer} x{n}");
    }
}
