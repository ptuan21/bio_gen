// Runs the WebAssembly engine off the main thread so large FASTA/BAM files
// never freeze the UI. The worker holds the loaded BAM/BAI bytes and exposes a
// small message API; results are posted back keyed by request id.

import init, * as bio from "./pkg/bio_wasm.js";

let bam = null, bai = null;
const ready = init().then(() => self.postMessage({ ready: true }));

self.onmessage = async (e) => {
  const { id, method, p } = e.data;
  await ready;
  try {
    let result;
    switch (method) {
      case "call": result = bio[p.fn](...p.args); break;
      case "setBam": bam = p.bam; result = bam.length; break;
      case "setBai": bai = p.bai; result = bai.length; break;
      case "parseBam": result = bio.parse_bam(bam, p.max); break;
      case "region": result = bio.parse_bam_region(bam, bai, p.ref, p.beg, p.end); break;
      case "pileup": result = bio.bam_pileup(bam, bai, p.ref, p.beg, p.end, p.minQual); break;
      case "varcall":
        result = bio.call_variants_pileup(bam, bai, p.ref, p.beg, p.end, p.reference, p.offset, p.minDepth, p.minFreq, p.minQual, p.minStrandFrac, p.rna);
        break;
      case "pileupVcf":
        result = bio.pileup_variants_vcf(bam, bai, p.ref, p.beg, p.end, p.reference, p.offset, p.minDepth, p.minFreq, p.minQual, p.minStrandFrac, p.rna);
        break;
      case "streamFasta": {
        // Read the File in slices, feeding a streaming parser. Memory stays
        // flat no matter how large the file is; only summaries are kept.
        const file = p.file;
        const chunkSize = p.chunkSize || 1 << 20; // 1 MB
        const cap = p.cap || 1000;
        const streamer = new bio.FastaStreamer();
        const records = [];
        let totalRecords = 0, totalLength = 0;
        const collect = (done) => {
          for (const r of done) {
            totalRecords++;
            totalLength += r.length;
            if (records.length < cap) records.push(r);
          }
        };
        for (let off = 0; off < file.size; off += chunkSize) {
          const buf = new Uint8Array(await file.slice(off, off + chunkSize).arrayBuffer());
          collect(streamer.push(buf));
          self.postMessage({ id, progress: Math.min(1, (off + chunkSize) / file.size) });
        }
        collect(streamer.finish());
        result = { records, totalRecords, totalLength, capped: totalRecords > records.length };
        break;
      }
      default: throw new Error("unknown method " + method);
    }
    self.postMessage({ id, ok: true, result });
  } catch (err) {
    self.postMessage({ id, ok: false, error: String((err && err.message) || err) });
  }
};
