#!/usr/bin/env bash
# =============================================================================
# 01_build.sh — box, solvation, neutralisation and index groups for ONE system
#
#   bash scripts/01_build.sh <SYSTEM>
#   SYSTEM in: 9BG5 9BG9 4TQA 8BLR 4OBE
#
# Env overrides: GMX (default: gmx), BOXD (1.0 nm), CONC (0.15 M NaCl)
# Produces: runs/<SYSTEM>/{solv_ions.gro, topol.top, index.ndx, build.log}
# =============================================================================
set -euo pipefail
SYS="${1:?usage: 01_build.sh <SYSTEM>}"
GMX="${GMX:-gmx}"; export GMX_MAXBACKUP=-1
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BOXD="${BOXD:-1.0}"
CONC="${CONC:-0.15}"
IN="$ROOT/inputs/$SYS"; RUN="$ROOT/runs/$SYS"
[[ -d "$IN" ]] || { echo "ERROR: no such system: inputs/$SYS"; exit 1; }

mkdir -p "$RUN"; cd "$RUN"
cp "$IN"/*.itp "$IN"/topol.top "$IN"/system_raw.gro .
mkdir -p toppar && cp "$ROOT"/toppar/*.itp toppar/
sed -i 's|\.\./\.\./toppar/|toppar/|g' topol.top      # make includes local to runs/<SYS>
: > build.log

echo "[$SYS] 1/4  box: dodecahedron, -d ${BOXD} nm"
$GMX editconf -f system_raw.gro -o box.gro -bt dodecahedron -d "$BOXD" -c >> build.log 2>&1

echo "[$SYS] 2/4  solvate with TIP3P"
$GMX solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top >> build.log 2>&1

echo "[$SYS] 3/4  ions: ${CONC} M NaCl + neutralise"
$GMX grompp -f "$ROOT/mdp/em.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 5 >> build.log 2>&1
echo SOL | $GMX genion -s ions.tpr -o solv_ions.gro -p topol.top \
        -pname NA -nname CL -neutral -conc "$CONC" >> build.log 2>&1

echo "[$SYS] 4/4  index groups"
python3 "$ROOT/scripts/make_index.py" solv_ions.gro topol.top index.ndx | tee -a build.log

# verify the final topology is charge-neutral and grompp-clean
$GMX grompp -f "$ROOT/mdp/em.mdp" -c solv_ions.gro -p topol.top -n index.ndx \
        -o .check.tpr -maxwarn 2 >> build.log 2>&1
if grep -q "non-zero total charge" build.log; then
    echo "[$SYS] WARNING: residual net charge — inspect build.log"
fi
rm -f .check.tpr ions.tpr solv.gro box.gro \#*

echo "[$SYS] build OK: $(sed -n 2p solv_ions.gro | tr -d ' ') atoms -> $RUN/solv_ions.gro"
awk '/\[ molecules \]/,0' topol.top | tail -n +2
