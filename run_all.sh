#!/usr/bin/env bash
# Run the whole pipeline in order. Each step is resumable -- rerunning skips
# work whose output already exists in results/.
set -euo pipefail
cd "$(dirname "$0")/scripts"

for step in 01_extract_hypothetical 02_map_orthologs 03_reannotate \
            04_localization 05_fetch_structures 06_foldseek \
            07_epitopes 08_conservation 09_expression 10_score; do
  echo ""
  echo "=============================================================="
  echo ">>> ${step}"
  echo "=============================================================="
  python3 "${step}.py"
done

echo ""
echo "Done. Final table: results/10_FINAL_SCORECARD.tsv"
