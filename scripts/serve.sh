#!/usr/bin/env bash
# Serve the demo at http://localhost:8000/web/ (ES modules need HTTP, not file://).
set -euo pipefail
cd "$(dirname "$0")/.."
echo "Open http://localhost:8000/web/"
python3 -m http.server 8000
