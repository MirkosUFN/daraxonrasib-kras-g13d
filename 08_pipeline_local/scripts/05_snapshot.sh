#!/usr/bin/env bash
# =============================================================================
# 05_snapshot.sh — inspect a campaign that is still running
#
#   bash scripts/05_snapshot.sh [SYSTEMS...]
#
# A 200 ns run takes days. This postprocesses and analyses whatever has been
# written so far, so you can watch the trajectory develop instead of waiting for
# the end. It is safe to run while mdrun is still writing: trjconv reads md.xtc,
# it never writes to it.
#
# Outputs land in analysis/ exactly as after a finished run, and are overwritten
# by the next snapshot. Read runs/<SYS>/pbc_check.txt each time.
#
# Env: GMX (default gmx), MOVIE_FRAMES (default 200)
# =============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SYSTEMS=("$@"); [[ ${#SYSTEMS[@]} -eq 0 ]] && SYSTEMS=(9BG5 9BG9 4TQA 8BLR 4OBE)
GMX="${GMX:-gmx}"

DONE=()
for S in "${SYSTEMS[@]}"; do
    if [[ ! -f "$ROOT/runs/$S/md.xtc" ]]; then
        echo "[$S] no md.xtc yet — skipping"
        continue
    fi
    NS=$($GMX check -f "$ROOT/runs/$S/md.xtc" 2>&1 | awk '/^Last frame/ {print $NF}')
    echo "[$S] trajectory currently reaches ${NS:-?} ps"
    bash "$ROOT/scripts/03_postprocess.sh" "$S"
    DONE+=("$S")
done

if [[ ${#DONE[@]} -eq 0 ]]; then
    echo "nothing to analyse yet"; exit 0
fi
# equil-frac 0.5 discards the first half; for an early snapshot that may be all
# you have, so drop to 0.2 while the run is short
python3 "$ROOT/scripts/04_analysis.py" --root "$ROOT" --systems "${DONE[@]}" \
        --equil-frac "${EQUIL_FRAC:-0.5}"
echo
echo "snapshot written to $ROOT/analysis  (systems: ${DONE[*]})"
