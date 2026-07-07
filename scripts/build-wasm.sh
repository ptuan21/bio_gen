#!/usr/bin/env bash
# Build the browser-ready WebAssembly package into web/pkg.
# Prefers wasm-pack; falls back to cargo + wasm-bindgen-cli.
set -euo pipefail

cd "$(dirname "$0")/.."
OUT="web/pkg"

if command -v wasm-pack >/dev/null 2>&1; then
  wasm-pack build crates/bio-wasm --target web --release --out-dir "../../$OUT"
else
  echo "wasm-pack not found; using cargo + wasm-bindgen." >&2
  if ! command -v wasm-bindgen >/dev/null 2>&1; then
    echo "Install one of:" >&2
    echo "  cargo install wasm-pack" >&2
    echo "  cargo install wasm-bindgen-cli --version 0.2.126" >&2
    exit 1
  fi
  rustup target add wasm32-unknown-unknown >/dev/null 2>&1 || true
  cargo build -p bio-wasm --target wasm32-unknown-unknown --release
  wasm-bindgen --target web --out-dir "$OUT" \
    target/wasm32-unknown-unknown/release/bio_wasm.wasm
fi

echo "Built WebAssembly package to $OUT"
