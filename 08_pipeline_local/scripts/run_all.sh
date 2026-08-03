#!/usr/bin/env bash
# =============================================================================
# run_all.sh — full pipeline for every system, then the comparative analysis
#
#   bash scripts/run_all.sh [NS] [SYSTEMS...]
#     NS       production length per system in ns (default 200)
#     SYSTEMS  subset to run (default: 9BG5 9BG9 4TQA 8BLR 4OBE)
#
#   # all five systems, 200 ns each (CPU: ~28 days — see MAXH below):
#   bash scripts/run_all.sh
#   # 200 ns each on a GPU machine:
#   GPU=1 bash scripts/run_all.sh 200
#   # only the two RAS-ON systems:
#   bash scripts/run_all.sh 200 9BG5 9BG9
#
# For a 200 ns campaign, running the systems SEQUENTIALLY through this script is
# rarely what you want: one system finishing after four weeks gives you nothing
# to look at in the meantime. Prefer one long-lived process per system, each
# with its own thread slice, so results accrue in parallel:
#
#   for S in 9BG5 9BG9 4TQA 8BLR 4OBE; do
#       NT=4 nohup bash scripts/02_run_md.sh $S 200 > runs/$S/run.log 2>&1 &
#   done
#
# (NT=4 x 5 systems = 20 of 24 cores. Do not oversubscribe: total NT across all
# concurrent mdrun processes must stay at or below your core count.)
# Then run the analysis whenever you want a snapshot — it works on whatever
# systems have finished, and 03_postprocess.sh can be run on a partial md.xtc.
# =============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NS="${1:-200}"; shift || true
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
