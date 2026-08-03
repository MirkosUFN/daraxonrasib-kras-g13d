#!/usr/bin/env bash
# =============================================================================
# run_all.sh — full pipeline for every system, then the comparative analysis
#
#   bash scripts/run_all.sh [NS] [SYSTEMS...]
#     NS       production length per system in ns (default 5)
#     SYSTEMS  subset to run (default: 9BG5 9BG9 4TQA 8BLR 4OBE)
#
#   # 10 ns each on a GPU machine:
#   GPU=1 bash scripts/run_all.sh 10
#   # only the two RAS-ON systems, 5 ns:
#   bash scripts/run_all.sh 5 9BG5 9BG9
# =============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NS="${1:-5}"; shift || true
SYSTEMS=("$@"); [[ ${#SYSTEMS[@]} -eq 0 ]] && SYSTEMS=(9BG5 9BG9 4TQA 8BLR 4OBE)

echo "=================================================================="
echo " daraxonrasib / KRAS tri-complex MD — ${#SYSTEMS[@]} systems, ${NS} ns each"
echo " systems: ${SYSTEMS[*]}"
echo "=================================================================="
START=$(date +%s)
for S in "${SYSTEMS[@]}"; do
    echo; echo "################ $S ################"
    bash "$ROOT/scripts/01_build.sh"  "$S"
    bash "$ROOT/scripts/02_run_md.sh" "$S" "$NS"
done
echo; echo "################ analysis ################"
python3 "$ROOT/scripts/04_analysis.py" --root "$ROOT" --systems "${SYSTEMS[@]}"
printf "\nAll done in %s\n" "$(python3 -c "
import sys;t=int(sys.argv[1]);print(f'{t//3600}h{(t%3600)//60:02d}m')" $(( $(date +%s) - START )))"
