// In-browser generator for a small, valid BAM + BAI pair, so the alignment
// features can be tested without external files. BGZF blocks use DEFLATE
// "stored" blocks, so no compression library is needed.

const SEQ_CODE = { A: 1, C: 2, G: 4, T: 8, U: 8, N: 15 };

class Bytes {
  constructor() { this.a = []; }
  u8(v) { this.a.push(v & 0xff); return this; }
  u16(v) { this.a.push(v & 0xff, (v >> 8) & 0xff); return this; }
  u32(v) { this.a.push(v & 0xff, (v >> 8) & 0xff, (v >> 16) & 0xff, (v >>> 24) & 0xff); return this; }
  i32(v) { return this.u32(v >>> 0); }
  u64(v) { return this.u32(v >>> 0).u32(Math.floor(v / 2 ** 32) >>> 0); }
  str(s) { for (const c of s) this.a.push(c.charCodeAt(0) & 0xff); return this; }
  bytes(arr) { for (const b of arr) this.a.push(b & 0xff); return this; }
  get length() { return this.a.length; }
  build() { return Uint8Array.from(this.a); }
}

// Smallest bin fully containing [beg, end) (BAM/SAM binning scheme).
function reg2bin(beg, end) {
  end -= 1;
  if (beg >> 14 === end >> 14) return 4681 + (beg >> 14);
  if (beg >> 17 === end >> 17) return 585 + (beg >> 17);
  if (beg >> 20 === end >> 20) return 73 + (beg >> 20);
  if (beg >> 23 === end >> 23) return 9 + (beg >> 23);
  if (beg >> 26 === end >> 26) return 1 + (beg >> 26);
  return 0;
}

function encodeSeq(seq) {
  const out = [];
  for (let i = 0; i < seq.length; i += 2) {
    const hi = SEQ_CODE[seq[i]] ?? 15;
    const lo = i + 1 < seq.length ? (SEQ_CODE[seq[i + 1]] ?? 15) : 0;
    out.push((hi << 4) | lo);
  }
  return out;
}

// A single DEFLATE stored block (BFINAL=1, BTYPE=00). Requires len < 65536.
function deflateStored(data) {
  const len = data.length;
  return [0x01, len & 0xff, (len >> 8) & 0xff, ~len & 0xff, (~len >> 8) & 0xff, ...data];
}

function bgzfBlock(data) {
  const cdata = deflateStored(data);
  const bsize = 12 + 6 + cdata.length + 8 - 1;
  return new Bytes()
    .bytes([0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff])
    .u16(6).str("BC").u16(2).u16(bsize)
    .bytes(cdata)
    .u32(0)            // CRC32 (ignored by the reader)
    .u32(data.length)  // ISIZE
    .build();
}

function recordBytes(index, pos, seq) {
  const name = "read" + index;
  const bin = reg2bin(pos, pos + seq.length);
  return new Bytes()
    .i32(0)                       // ref_id
    .i32(pos)
    .u8(name.length + 1)          // l_read_name
    .u8(60)                       // mapq
    .u16(bin)
    .u16(1)                       // n_cigar
    .u16(0)                       // flag
    .u32(seq.length)              // l_seq
    .i32(-1).i32(-1).i32(0)       // next_ref, next_pos, tlen
    .str(name).u8(0)
    .u32((seq.length << 4) | 0)   // CIGAR: <len>M
    .bytes(encodeSeq(seq))
    .bytes(new Array(seq.length).fill(40)) // qual
    .build();
}

// Build a sample BAM+BAI from `reference`. Half of the reads covering `mutPos`
// carry `altBase`, so the pileup shows a ~0.5-frequency SNV there.
export function buildSampleBam(reference, opts = {}) {
  const ref = reference.toUpperCase().replace(/[^ACGTUN]/g, "N");
  const refLen = ref.length;
  if (refLen < 4) throw "reference too short (need >= 4 bases)";

  const readLen = Math.min(opts.readLen || 20, refLen);
  const reads = opts.reads || 30;
  const refName = opts.refName || "ref";
  const mutPos = opts.mutPos ?? Math.floor(refLen / 2);
  const refAt = ref[mutPos] || "A";
  const altBase = refAt === "A" ? "G" : "A";

  const bam = new Bytes()
    .str("BAM\x01").u32(0).u32(1)
    .u32(refName.length + 1).str(refName).u8(0).i32(refLen);

  const span = Math.max(1, reads - 1);
  const placed = [];
  for (let i = 0; i < reads; i++) {
    const pos = Math.floor((i * (refLen - readLen)) / span);
    let seq = ref.slice(pos, pos + readLen);
    if (i % 2 === 0 && mutPos >= pos && mutPos < pos + readLen) {
      const k = mutPos - pos;
      seq = seq.slice(0, k) + altBase + seq.slice(k + 1);
    }
    const beg = bam.length;
    const rb = recordBytes(i, pos, seq);
    bam.u32(rb.length).bytes(rb);
    placed.push({ pos, len: readLen, beg, end: bam.length });
  }

  const bamRaw = bam.build();
  if (bamRaw.length >= 65536) throw "sample too large; reduce reads or reference length";

  // BAI: one chunk per record (coffset 0), linear index left zeroed (safe).
  const binMap = new Map();
  for (const r of placed) {
    const bin = reg2bin(r.pos, r.pos + r.len);
    if (!binMap.has(bin)) binMap.set(bin, []);
    binMap.get(bin).push([r.beg, r.end]);
  }
  const bai = new Bytes().str("BAI\x01").u32(1).u32(binMap.size);
  for (const [bin, cks] of binMap) {
    bai.u32(bin).u32(cks.length);
    for (const [b, e] of cks) bai.u64(b).u64(e);
  }
  const nIntv = (refLen >> 14) + 1;
  bai.u32(nIntv);
  for (let i = 0; i < nIntv; i++) bai.u64(0);

  return {
    bam: bgzfBlock(bamRaw),
    bai: bai.build(),
    refName, refLen, readLen, reads, mutPos, refAt, altBase,
  };
}
