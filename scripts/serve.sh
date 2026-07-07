#!/usr/bin/env bash
# Serve the demo at http://localhost:8000/ (ES modules need HTTP, not file://).
set -euo pipefail
cd "$(dirname "$0")/../web"
echo "Open http://localhost:8000/"
python3 -m http.server 8000
