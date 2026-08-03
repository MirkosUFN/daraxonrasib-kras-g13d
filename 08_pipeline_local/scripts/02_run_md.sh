#!/usr/bin/env bash
# =============================================================================
# 02_run_md.sh — EM -> NVT -> NPT -> production MD for ONE system
#
#   bash scripts/02_run_md.sh <SYSTEM> [NS]
#     NS = production length in nanoseconds (default 5)
#
# Env overrides:
#   GMX     gmx binary (default: gmx)
#   NT      OpenMP threads (default: all cores)
#   GPU     1 -> add "-nb gpu -bonded gpu -pme gpu" (needs a GPU-enabled build)
#   NVT_PS  NVT equilibration length in ps (default 100)
#   NPT_PS  NPT equilibration length in ps (default 100)
#
# Interrupted run? Just re-issue the same command: any stage whose .gro already
# exists is skipped, and an unfinished production restarts from md.cpt.
# Requires: scripts/01_build.sh already run for <SYSTEM>
# =============================================================================
set -euo pipefail
SYS="${1:?usage: 02_run_md.sh <SYSTEM> [NS]}"
NS="${2:-5}"
GMX="${GMX:-gmx}"; export GMX_MAXBACKUP=-1
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN="$ROOT/runs/$SYS"
[[ -f "$RUN/solv_ions.gro" ]] || { echo "ERROR: run 01_build.sh $SYS first"; exit 1; }
cd "$RUN"

NT="${NT:-$(nproc)}"
MDRUN=($GMX mdrun -ntomp "$NT" -pin on)
[[ "${GPU:-0}" == "1" ]] && MDRUN+=(-nb gpu -bonded gpu -pme gpu)

NSTEPS=$(python3 -c "print(int(round($NS*1e6/2)))")            # dt = 2 fs
NVT_STEPS=$(python3 -c "print(int(round(${NVT_PS:-100}/0.002)))")
NPT_STEPS=$(python3 -c "print(int(round(${NPT_PS:-100}/0.002)))")
echo "[$SYS] production: ${NS} ns = ${NSTEPS} steps | threads=$NT | GPU=${GPU:-0}"

echo "[$SYS] === 1/4 energy minimisation ==="
if [[ -f em.gro ]]; then echo "  em.gro exists — skipping"; else
  $GMX grompp -f "$ROOT/mdp/em.mdp"  -c solv_ions.gro -r solv_ions.gro -p topol.top -n index.ndx -o em.tpr  -maxwarn 2
  "${MDRUN[@]}" -deffnm em
  echo Potential | $GMX energy -f em.edr -o em_potential.xvg >/dev/null 2>&1 || true
fi

echo "[$SYS] === 2/4 NVT equilibration (${NVT_PS:-100} ps, restrained) ==="
if [[ -f nvt.gro ]]; then echo "  nvt.gro exists — skipping"; else
  sed "s/^nsteps .*/nsteps                  = ${NVT_STEPS}/" "$ROOT/mdp/nvt.mdp" > nvt_run.mdp
  $GMX grompp -f nvt_run.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2
  "${MDRUN[@]}" -deffnm nvt
  echo Temperature | $GMX energy -f nvt.edr -o nvt_temperature.xvg >/dev/null 2>&1 || true
fi

echo "[$SYS] === 3/4 NPT equilibration (${NPT_PS:-100} ps, restrained) ==="
if [[ -f npt.gro ]]; then echo "  npt.gro exists — skipping"; else
  sed "s/^nsteps .*/nsteps                  = ${NPT_STEPS}/" "$ROOT/mdp/npt.mdp" > npt_run.mdp
  $GMX grompp -f npt_run.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 2
  "${MDRUN[@]}" -deffnm npt
  printf "Pressure\nDensity\n" | $GMX energy -f npt.edr -o npt_pressure_density.xvg >/dev/null 2>&1 || true
fi

echo "[$SYS] === 4/4 production MD (${NS} ns, unrestrained) ==="
if [[ -f md.gro ]]; then
  echo "  md.gro exists — production already complete, skipping"
elif [[ -f md.cpt ]]; then
  echo "  md.cpt found — resuming from checkpoint"
  "${MDRUN[@]}" -deffnm md -cpi md.cpt
else
  sed "s/^nsteps .*/nsteps                  = ${NSTEPS}/" "$ROOT/mdp/md.mdp" > md_run.mdp
  $GMX grompp -f md_run.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr -maxwarn 2
  "${MDRUN[@]}" -deffnm md
fi
echo "[$SYS] MD finished -> $RUN/md.xtc"

bash "$ROOT/scripts/03_postprocess.sh" "$SYS"
