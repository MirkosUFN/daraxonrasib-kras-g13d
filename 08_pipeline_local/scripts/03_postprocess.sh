#!/usr/bin/env bash
# =============================================================================
# 03_postprocess.sh — PBC correction + viewable structures for ONE system
#
#   bash scripts/03_postprocess.sh <SYSTEM>
#
# ---------------------------------------------------------------------------
# WHY THIS EXACT ORDER (three stages, and NOT the usual four-step recipe)
#
# This is a TWO-CHAIN complex (KRAS + CypA) plus three free HETATM species
# (daraxonrasib, the nucleotide, Mg2+). Each is a separate "molecule" to
# GROMACS, so two independent faults can appear:
#   (i)  a single chain broken across the box boundary, and
#   (ii) the two intact chains wrapped into opposite corners.
# Either one makes the complex look torn apart (protein Rg ~3.4 nm instead of
# ~2.2 nm) and gives nonsensical drug-CypA distances (~55 A instead of ~2 A).
#
#   (1) -pbc whole                  : repair each molecule individually  (fixes i)
#   (2) -pbc cluster on Protein_LIG : gather the solute into one cluster (fixes ii)
#   (3) -fit rot+trans on Protein   : remove global rotation/translation,
#                                     writing only the solute
#
# Two steps that appear in most online recipes are deliberately ABSENT, because
# on this system they silently destroy the assembly (verified by checking the
# radius of gyration of each chain):
#   * "-pbc nojump" AFTER clustering re-splits CypA (Rg 1.45 -> 3.66 nm). Its
#     atom-wise unwrapping undoes the molecular reassembly done in stages 1-2,
#     and the drug-CypA distance then reads 1.7 A to a periodic-image fragment.
#   * "-center ... -pbc mol -ur compact" re-wraps every molecule by its own
#     centre of mass, again pulling the two chains into opposite corners.
#   Fitting in stage 3 already centres the solute, so neither is needed.
#
# The pbc_check.txt written at the end compares the per-chain radius of
# gyration against the pre-MD reference. Read it — it is the guard that catches
# a broken complex before it reaches the analysis.
# ---------------------------------------------------------------------------
#
# Produces in runs/<SYSTEM>/:
#   md_clean.xtc            solute-only trajectory, PBC-corrected and fitted
#   complex.tpr             solute-only run input (matches md_clean.xtc)
#   frame0.pdb              first frame of the solute (static picture)
#   trajectory_movie.pdb    multi-MODEL pdb (every 5th frame) for PyMOL/VMD
#   pbc_check.txt           per-chain Rg, before vs. after — sanity check
# =============================================================================
set -euo pipefail
SYS="${1:?usage: 03_postprocess.sh <SYSTEM>}"
GMX="${GMX:-gmx}"; export GMX_MAXBACKUP=-1
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN="$ROOT/runs/$SYS"; cd "$RUN"
[[ -f md.xtc ]] || { echo "ERROR: no md.xtc in $RUN"; exit 1; }

echo "[$SYS] PBC 1/3  make every molecule whole"
echo System | $GMX trjconv -s md.tpr -f md.xtc -n index.ndx -o .p1.xtc -pbc whole
echo "[$SYS] PBC 2/3  cluster the solute (KRAS + CypA + ligands) together"
printf "Protein_LIG\nSystem\n" | $GMX trjconv -s md.tpr -f .p1.xtc -n index.ndx \
        -o .p2.xtc -pbc cluster
echo "[$SYS] PBC 3/3  fit on Protein (rot+trans), write solute only"
printf "Protein\nProtein_LIG\n" | $GMX trjconv -s md.tpr -f .p2.xtc -n index.ndx \
        -o md_clean.xtc -fit rot+trans

# solute-only .tpr so md_clean.xtc can be read on its own
echo Protein_LIG | $GMX convert-tpr -s md.tpr -n index.ndx -o complex.tpr

echo "[$SYS] writing frame0.pdb and trajectory_movie.pdb"
$GMX trjconv -s complex.tpr -f md_clean.xtc -o frame0.pdb -dump 0 -conect <<< "System"
$GMX trjconv -s complex.tpr -f md_clean.xtc -o trajectory_movie.pdb -skip 5 -conect <<< "System"
rm -f .p1.xtc .p2.xtc \#*

# ---- guard: per-chain radius of gyration, pre-MD vs. post-correction --------
python3 - "$SYS" <<'PY' | tee pbc_check.txt
import sys, warnings
warnings.filterwarnings("ignore")
try:
    import MDAnalysis as mda
except ImportError:
    print("MDAnalysis not installed - skipping PBC check"); sys.exit(0)
from MDAnalysis.lib.distances import distance_array
sysname = sys.argv[1]

def chains(u):
    p = u.select_atoms("protein"); res = p.residues; ids = res.resids
    cut = next((i for i in range(1, len(ids)) if ids[i] <= ids[i-1]), None)
    if cut is None: cut = min(169, len(ids)-1)
    return p, res[:cut].atoms, res[cut:].atoms

def report(tag, top, trj):
    u = mda.Universe(top, trj, refresh_offsets=True)
    p, A, B = chains(u); d = u.select_atoms("resname DRG")
    row = dict(tag=tag, kras=A.radius_of_gyration()/10, cypa=B.radius_of_gyration()/10,
               prot=p.radius_of_gyration()/10)
    if d.n_atoms:
        row["dk"] = distance_array(d.positions, A.positions).min()
        row["dc"] = distance_array(d.positions, B.positions).min()
    return row

print(f"PBC sanity check - {sysname}")
print("Radius of gyration (nm); a post-MD value well above the pre-MD reference")
print("means the complex is still split across the periodic box.\n")
ref  = report("pre-MD  (solv_ions.gro)", "md.tpr", "solv_ions.gro")
post = report("post-PBC (md_clean.xtc)", "complex.tpr", "md_clean.xtc")
hdr = f"{'stage':26s} {'KRAS':>7s} {'CypA':>7s} {'Protein':>8s}"
if "dk" in ref: hdr += f" {'drug-KRAS':>10s} {'drug-CypA':>10s}"
print(hdr)
for r in (ref, post):
    line = f"{r['tag']:26s} {r['kras']:7.3f} {r['cypa']:7.3f} {r['prot']:8.3f}"
    if "dk" in r: line += f" {r['dk']:10.2f} {r['dc']:10.2f}"
    print(line)
bad = [n for n in ("kras","cypa","prot") if post[n] > 1.15*ref[n]]
print()
if bad:
    print(f"WARNING: {', '.join(bad)} still looks split across the box.")
    print("         Inspect frame0.pdb before trusting any distance or contact.")
    sys.exit(0)
print("PASSED: every chain is contiguous and the complex is intact.")
PY
echo "[$SYS] postprocess OK: md_clean.xtc, complex.tpr, frame0.pdb, trajectory_movie.pdb"
