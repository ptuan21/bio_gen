set -euo pipefail
cd "$(dirname "$0")/../web"
echo "Open http://localhost:8000/"
python3 -m http.server 8000
