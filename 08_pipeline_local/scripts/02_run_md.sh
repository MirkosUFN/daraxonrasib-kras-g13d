#!/usr/bin/env bash
# =============================================================================
# 02_run_md.sh — EM -> NVT -> NPT -> production MD for ONE system
#
#   bash scripts/02_run_md.sh <SYSTEM> [NS]
#     NS = production length in nanoseconds (default 200)
#
# Env overrides:
#   GMX     gmx binary (default: gmx)
#   NT      OpenMP threads (default: all cores)
#   GPU     1 -> add "-nb gpu -bonded gpu -pme gpu" (needs a GPU-enabled build)
#   MAXH    wall-clock limit in HOURS for the production run. mdrun stops
#           cleanly just before it, writing md.cpt; re-issue the same command
#           to continue. Use this for multi-day runs you want to pause.
#   NVT_PS  NVT equilibration length in ps (default 100)
#   NPT_PS  NPT equilibration length in ps (default 100)
#
# Interrupted run? Just re-issue the same command: any stage whose .gro already
# exists is skipped, and an unfinished production restarts from md.cpt. This is
# safe to do any number of times — a 200 ns run split over ten sessions is
# bit-for-bit equivalent to one uninterrupted run.
# Requires: scripts/01_build.sh already run for <SYSTEM>
# =============================================================================
set -euo pipefail
SYS="${1:?usage: 02_run_md.sh <SYSTEM> [NS]}"
NS="${2:-200}"
GMX="${GMX:-gmx}"; export GMX_MAXBACKUP=-1
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN="$ROOT/runs/$SYS"
[[ -f "$RUN/solv_ions.gro" ]] || { echo "ERROR: run 01_build.sh $SYS first"; exit 1; }
cd "$RUN"

NT="${NT:-$(nproc)}"
MDRUN=($GMX mdrun -ntomp "$NT" -pin on)
[[ "${GPU:-0}" == "1" ]] && MDRUN+=(-nb gpu -bonded gpu -pme gpu)
PROD=("${MDRUN[@]}")
[[ -n "${MAXH:-}" ]] && PROD+=(-maxh "$MAXH")

NSTEPS=$(python3 -c "print(int(round($NS*1e6/2)))")            # dt = 2 fs
NVT_STEPS=$(python3 -c "print(int(round(${NVT_PS:-100}/0.002)))")
NPT_STEPS=$(python3 -c "print(int(round(${NPT_PS:-100}/0.002)))")
echo "[$SYS] production: ${NS} ns = ${NSTEPS} steps | threads=$NT | GPU=${GPU:-0}"
# rough resource forecast so a multi-day run is a deliberate choice
python3 - "$NS" "${GPU:-0}" "$NT" "${MAXH:-0}" <<'PY'
import sys
ns, gpu, nt, maxh = float(sys.argv[1]), sys.argv[2] == "1", int(sys.argv[3]), float(sys.argv[4])
# 32 ns/day measured on 24 CPU cores; scaling with thread count is sublinear
rate = 250.0 if gpu else 32.0 * (nt / 24.0) ** 0.85
h = ns / rate * 24
frames = ns * 1000 / 20                           # nstxout-compressed = 20 ps
print(f"  forecast: ~{h:.1f} h ({h/24:.1f} days) at ~{rate:.0f} ns/day "
      f"({'GPU' if gpu else f'CPU, {nt} threads'}), {frames:.0f} frames")
print(f"  disk: ~{frames*67000*3*2/1e9:.1f} GB md.xtc + ~{frames*5400*3*2/1e9:.2f} GB md_clean.xtc")
if maxh > 0:
    print(f"  MAXH={maxh:g} h -> this invocation covers ~{maxh*rate/24:.1f} ns; "
          f"re-issue the same command ~{max(1, -(-h//maxh)):.0f}x in total to finish.")
elif h > 24 and not gpu:
    print("  NOTE: over a day per system. Set MAXH to pause/resume, or GPU=1 if available.")
PY

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
  "${PROD[@]}" -deffnm md -cpi md.cpt
else
  sed "s/^nsteps .*/nsteps                  = ${NSTEPS}/" "$ROOT/mdp/md.mdp" > md_run.mdp
  $GMX grompp -f md_run.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr -maxwarn 2
  "${PROD[@]}" -deffnm md
fi

# with MAXH the run stops early on purpose — do not postprocess a partial run
if [[ ! -f md.gro ]]; then
  echo "[$SYS] production not finished yet (MAXH reached or interrupted)."
  echo "[$SYS] re-run the same command to continue from md.cpt."
  exit 0
fi
echo "[$SYS] MD finished -> $RUN/md.xtc"

bash "$ROOT/scripts/03_postprocess.sh" "$SYS"
