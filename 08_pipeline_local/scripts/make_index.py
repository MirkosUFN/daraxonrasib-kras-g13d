#!/usr/bin/env python3
"""Write a GROMACS index file deterministically from solv_ions.gro + topol.top.

Chain boundaries are taken from the [ molecules ] block of topol.top and the
atom counts of each topol_Protein_chain_*.itp, so the split is exact and does
not depend on residue numbering (which restarts per chain in .gro files).

Groups written:
  System, Protein, KRAS, CypA, DRG, LIG, Protein_LIG, Water, Ions, Water_and_ions

  LIG          = daraxonrasib (DRG) + nucleotide (GNP/GDP) + Mg2+
  Protein_LIG  = Protein + LIG  -> thermostat group and analysis/centring group

Usage: python3 make_index.py <solv_ions.gro> <topol.top> <index.ndx>
"""
import sys, os, re, collections

SOLV = {"SOL", "WAT", "HOH", "TIP3"}
IONS = {"NA", "CL", "K", "BR", "ZN"}
COFAC = {"DRG", "GNP", "GDP", "MG"}

def itp_natoms(path):
    txt = open(path).read()
    m = re.search(r"\[\s*atoms\s*\](.*?)(?=\n\s*\[)", txt, re.S)
    if not m:
        return 0
    return sum(1 for l in m.group(1).splitlines()
               if l.strip() and not l.strip().startswith(";"))

def molecules(top):
    txt = open(top).read()
    m = re.search(r"\[\s*molecules\s*\](.*)", txt, re.S)
    out = []
    for l in m.group(1).splitlines():
        l = l.split(";")[0].strip()
        if not l:
            continue
        f = l.split()
        if len(f) >= 2 and f[1].isdigit():
            out.append((f[0], int(f[1])))
    return out

def parse_gro(path):
    lines = open(path).read().splitlines()
    n = int(lines[1])
    return n, [l[5:10].strip() for l in lines[2:2 + n]]

def main():
    gro, top, out = sys.argv[1], sys.argv[2], sys.argv[3]
    topdir = os.path.dirname(os.path.abspath(top))
    natoms, resn = parse_gro(gro)

    # ---- protein chain sizes from the chain .itp files, in [ molecules ] order
    chain_sizes = []
    for name, nmol in molecules(top):
        if name.startswith("Protein_chain") or name.startswith("Protein"):
            p = os.path.join(topdir, f"topol_{name}.itp")
            if os.path.exists(p):
                chain_sizes += [itp_natoms(p)] * nmol
    g = collections.OrderedDict()
    g["System"] = list(range(1, natoms + 1))

    prot_n = sum(chain_sizes)
    g["Protein"] = list(range(1, prot_n + 1))
    if len(chain_sizes) >= 2:
        g["KRAS"] = list(range(1, chain_sizes[0] + 1))
        g["CypA"] = list(range(chain_sizes[0] + 1, chain_sizes[0] + chain_sizes[1] + 1))
    elif chain_sizes:
        g["KRAS"] = list(g["Protein"])

    drg, lig, wat, ion = [], [], [], []
    for i in range(prot_n + 1, natoms + 1):
        r = resn[i - 1]
        if r in COFAC:
            lig.append(i)
            if r == "DRG":
                drg.append(i)
        elif r in SOLV:
            wat.append(i)
        elif r in IONS:
            ion.append(i)
        else:
            lig.append(i)
    g["DRG"] = drg
    g["LIG"] = lig
    g["Protein_LIG"] = g["Protein"] + lig
    g["Water"] = wat
    g["Ions"] = ion
    g["Water_and_ions"] = sorted(wat + ion)

    with open(out, "w") as fh:
        for name, ids in g.items():
            if not ids:
                continue
            fh.write(f"[ {name} ]\n")
            for k in range(0, len(ids), 15):
                fh.write(" ".join(f"{x:6d}" for x in ids[k:k + 15]) + "\n")
    print("index:", {k: len(v) for k, v in g.items() if v})
    assert len(g["Protein_LIG"]) + len(g["Water_and_ions"]) == natoms, \
        "group partition does not cover the system"

if __name__ == "__main__":
    main()
