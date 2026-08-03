#!/usr/bin/env python3
# =============================================================================
# 04_analysis.py — trajectory analysis for the daraxonrasib / KRAS tri-complex
#
#   python3 scripts/04_analysis.py --root . --systems 9BG5 9BG9 4TQA 8BLR 4OBE
#   python3 scripts/04_analysis.py --root . --equil-frac 0.5     # drop first 50%
#
# Reads  runs/<SYS>/{complex.tpr, md_clean.xtc}   (written by 03_postprocess.sh)
# Writes analysis/<SYS>/*.png|csv  and  analysis/00_comparison/*
#
# All figure text is in ENGLISH.
#
# Metrics per system
#   - backbone RMSD of KRAS and of CypA vs. frame 0
#   - per-residue C-alpha RMSF (KRAS and CypA), switch-I / switch-II shaded
#   - radius of gyration of the protein
#   - minimum distance drug--KRAS and drug--CypA
#   - number of drug heavy-atom contacts < 4.0 A with each partner
#   - fraction of frames with the drug engaged (min. distance < 4.5 A)
#   - drug RMSD in the KRAS frame (does the glue stay in its pocket?)
# Comparative figures
#   - overlaid RMSD / RMSF / min. distance / contacts for all systems
#   - bar chart of drug engagement and a summary CSV table
# =============================================================================
import argparse, warnings
from pathlib import Path

import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.lib.distances import distance_array

warnings.filterwarnings("ignore")

LIG = "DRG"
SWITCH_I  = (30, 40)
SWITCH_II = (60, 76)
CONTACT_CUT = 4.0        # A, heavy-atom contact
ENGAGE_CUT  = 4.5        # A, "drug engaged with partner"

# Colour per system: RAS-ON = warm, RAS-OFF = cool
SYS_META = {
    "9BG5": dict(label="KRAS G13D + CypA (RAS-ON, GppNHp)", colour="#d62728", state="ON"),
    "9BG9": dict(label="KRAS WT + CypA (RAS-ON, GppNHp)",   colour="#ff7f0e", state="ON"),
    "4TQA": dict(label="KRAS G13D + CypA (RAS-OFF, GDP)",   colour="#1f77b4", state="OFF"),
    "8BLR": dict(label="KRAS G13D open switch-I (GDP)",     colour="#17becf", state="OFF"),
    "4OBE": dict(label="KRAS WT + CypA (RAS-OFF, GDP)",     colour="#9467bd", state="OFF"),
}
GREY = "#4d4d4d"

def style():
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
        "legend.frameon": False, "legend.fontsize": 8,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
    })

def shade_switches(ax, label=False):
    """Shade switch-I / switch-II; labels are drawn vertically inside the band
    so they never collide with each other or with the traces."""
    for (lo, hi), nm in ((SWITCH_I, "switch-I"), (SWITCH_II, "switch-II")):
        ax.axvspan(lo, hi, color="#ffcc66", alpha=0.35, lw=0, zorder=0)
        if label:
            ax.text((lo + hi) / 2, 0.985, nm, transform=ax.get_xaxis_transform(),
                    ha="center", va="top", rotation=90, fontsize=6.5, color=GREY,
                    zorder=1)

def panel(ax, letter):
    ax.text(-0.13, 1.06, letter, transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top", ha="left")

def native_numbering(root, sysname):
    """Recover the crystallographic residue numbering for each chain.

    pdb2gmx renumbers every chain from 1 and the two chains are concatenated in
    the .tpr, so a raw MDAnalysis resid is meaningless for a reader: CypA Arg55
    shows up as residue 223. inputs/<SYS>/protein_only.pdb still carries the
    original numbering, so build concat_resid -> (native_resid, resname) from it.
    Returns (map_dict, ok) — ok is False when the file is missing, in which case
    the analysis falls back to the concatenated numbering.
    """
    pdb = root / "inputs" / sysname / "protein_only.pdb"
    if not pdb.exists():
        return {}, False
    chains = {}
    order = []
    for line in pdb.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        ch, ri, rn = line[21], int(line[22:26]), line[17:20].strip()
        lst = chains.setdefault(ch, [])
        if ch not in order:
            order.append(ch)
        if not lst or lst[-1][0] != ri:
            lst.append((ri, rn))
    out, n = {}, 0
    for ch in order:
        for ri, rn in chains[ch]:
            n += 1
            out[n] = (ri, rn, ch)
    return out, bool(out)

def protein_chains(u):
    """Split the protein into KRAS (first chain) and CypA (second).

    A .tpr written by convert-tpr carries a single segid whose name contains
    spaces, so segid-based selections are unusable; the two chains are instead
    separated at the point where residue numbering restarts (pdb2gmx numbers
    each chain from 1). Falls back to the KRAS length (~169 residues) if the
    numbering happens to be continuous.
    """
    prot = u.select_atoms("protein")
    res = prot.residues
    ids = res.resids
    cut = next((i for i in range(1, len(ids)) if ids[i] <= ids[i - 1]), None)
    if cut is None:
        cut = min(169, len(ids) - 1)
    return res[:cut].atoms, res[cut:].atoms

def analyse(sysname, root, equil_frac, stride, rmsf_max_frames=2000,
            n_blocks=5):  # noqa: C901
    run = root / "runs" / sysname
    tpr, xtc = run / "complex.tpr", run / "md_clean.xtc"
    if not (tpr.exists() and xtc.exists()):
        print(f"  [{sysname}] SKIP — complex.tpr / md_clean.xtc not found")
        return None
    u = mda.Universe(str(tpr), str(xtc))
    kras, cypa = protein_chains(u)
    drug = u.select_atoms(f"resname {LIG}")
    drug_heavy = drug.select_atoms("not name H*")
    kras_heavy = kras.select_atoms("not name H*")
    cypa_heavy = cypa.select_atoms("not name H*")
    out = dict(name=sysname, n_frames=len(u.trajectory),
               kras_res=len(kras.residues), cypa_res=len(cypa.residues),
               n_drug=drug.n_atoms)

    natmap, nat_ok = native_numbering(root, sysname)
    out["native_ok"] = nat_ok

    def to_native(resids):
        if not nat_ok:
            return np.asarray(resids)
        return np.array([natmap.get(int(r), (int(r), "", ""))[0] for r in resids])

    ref = mda.Universe(str(tpr), str(xtc)); ref.trajectory[0]
    kras_bb = f"index {kras.atoms.ix[0]}:{kras.atoms.ix[-1]} and backbone"
    cypa_bb = f"index {cypa.atoms.ix[0]}:{cypa.atoms.ix[-1]} and backbone"

    # ---- RMSD (KRAS backbone as the fitting frame) + drug RMSD in that frame
    groups = [f"resname {LIG}"] if drug.n_atoms else []
    R = rms.RMSD(u, ref, select=kras_bb, groupselections=groups + [cypa_bb],
                 ref_frame=0).run(step=stride)
    res = R.results.rmsd
    t_ns = res[:, 1] / 1000.0
    out["t_ns"] = t_ns
    out["rmsd_kras"] = res[:, 2]
    if drug.n_atoms:
        out["rmsd_drug"] = res[:, 3]
        out["rmsd_cypa"] = res[:, 4]
    else:
        out["rmsd_drug"] = None
        out["rmsd_cypa"] = res[:, 3]

    # ---- per-frame geometry: min distances, contacts, Rg, per-residue contacts
    dmin_k, dmin_c, nc_k, nc_c, rg = [], [], [], [], []
    prot = u.select_atoms("protein")
    # residue-resolved contact frequency (fraction of frames with any heavy atom
    # of the drug within CONTACT_CUT of that residue)
    occ = {}
    for tag, grp in (("kras", kras_heavy), ("cypa", cypa_heavy)):
        occ[tag] = dict(resids=grp.residues.resids,
                        names=grp.residues.resnames,
                        hits=np.zeros(grp.n_residues),
                        rindex=grp.resindices - grp.residues.resindices[0])
    for ts in u.trajectory[::stride]:
        rg.append(prot.radius_of_gyration())
        if drug_heavy.n_atoms:
            dk = distance_array(drug_heavy.positions, kras_heavy.positions)
            dc = distance_array(drug_heavy.positions, cypa_heavy.positions)
            dmin_k.append(dk.min()); dmin_c.append(dc.min())
            nc_k.append(int((dk < CONTACT_CUT).sum()))
            nc_c.append(int((dc < CONTACT_CUT).sum()))
            for tag, dmat, grp in (("kras", dk, kras_heavy), ("cypa", dc, cypa_heavy)):
                close = dmat.min(axis=0) < CONTACT_CUT          # per protein atom
                ridx = occ[tag]["rindex"][close]
                if ridx.size:
                    occ[tag]["hits"][np.unique(ridx)] += 1
    out["rg"] = np.array(rg)
    nfr = max(len(rg), 1)
    for tag in ("kras", "cypa"):
        occ[tag]["freq"] = occ[tag]["hits"] / nfr
        occ[tag]["native"] = to_native(occ[tag]["resids"])
    out["occ"] = occ
    for k, v in (("dmin_kras", dmin_k), ("dmin_cypa", dmin_c),
                 ("nc_kras", nc_k), ("nc_cypa", nc_c)):
        out[k] = np.array(v) if v else None

    # ---- RMSF on the equilibrated part, computed on an aligned copy
    #
    # RMSF needs the trajectory in memory to superpose it, which costs
    # n_frames * n_atoms * 12 bytes. A 200 ns run at 20 ps/frame is 10,000
    # frames, so the equilibrated half is subsampled to at most
    # rmsf_max_frames before the transfer — RMSF is an average over frames and
    # is insensitive to the sampling interval, whereas an unbounded transfer
    # would need several GB.
    n_eq = int(len(t_ns) * equil_frac)
    v = mda.Universe(str(tpr), str(xtc))
    n_avail = len(v.trajectory) - n_eq
    rstep = max(1, -(-n_avail // rmsf_max_frames)) if n_avail > 0 else 1
    v.transfer_to_memory(start=n_eq, step=rstep)
    if rstep > 1:
        print(f"  [{sysname}] RMSF: {n_avail} equilibrated frames subsampled "
              f"every {rstep} -> {len(v.trajectory)} frames in memory")
    align.AlignTraj(v, v, select=kras_bb, in_memory=True).run()
    vk, vc = protein_chains(v)
    for tag, grp in (("kras", vk), ("cypa", vc)):
        ca = grp.select_atoms("name CA")
        if not ca.n_atoms:
            out[f"rmsf_{tag}"] = out[f"resid_{tag}"] = None
            continue
        r = rms.RMSF(ca).run()          # v already starts at the equilibration point
        out[f"rmsf_{tag}"] = r.results.rmsf
        out[f"resid_{tag}"] = to_native(ca.resids)
    del v
    out["n_eq"] = n_eq

    def mstd(a):
        return (float(np.mean(a[n_eq:])), float(np.std(a[n_eq:]))) if a is not None else (np.nan, np.nan)
    out["stats"] = {
        "rmsd_kras": mstd(out["rmsd_kras"]), "rmsd_cypa": mstd(out["rmsd_cypa"]),
        "rmsd_drug": mstd(out["rmsd_drug"]), "rg": mstd(out["rg"]),
        "dmin_kras": mstd(out["dmin_kras"]), "dmin_cypa": mstd(out["dmin_cypa"]),
        "nc_kras": mstd(out["nc_kras"]), "nc_cypa": mstd(out["nc_cypa"]),
    }
    for tag in ("kras", "cypa"):
        f = out[f"rmsf_{tag}"]
        out["stats"][f"rmsf_{tag}"] = (float(np.mean(f)), float(np.max(f))) if f is not None else (np.nan, np.nan)
    if out["dmin_kras"] is not None:
        out["engage_kras"] = float(np.mean(out["dmin_kras"][n_eq:] < ENGAGE_CUT))
        out["engage_cypa"] = float(np.mean(out["dmin_cypa"][n_eq:] < ENGAGE_CUT))
    else:
        out["engage_kras"] = out["engage_cypa"] = np.nan

    # ---- block averages: the convergence test that only a long run allows
    # Splitting the equilibrated part into independent blocks and comparing them
    # is what distinguishes "the metric has converged" from "the metric happened
    # to sit there during the window I looked at". On a 5 ns run the blocks are
    # too short to be independent; at 200 ns they are informative.
    nb = min(n_blocks, max(1, len(t_ns) - n_eq))
    edges = np.linspace(n_eq, len(t_ns), nb + 1).astype(int)
    blocks = []
    for i in range(nb):
        s, e = edges[i], edges[i + 1]
        if e <= s:
            continue
        b = {"block": i + 1, "t_start_ns": float(t_ns[s]), "t_end_ns": float(t_ns[e - 1])}
        b["rmsd_kras"] = float(np.mean(out["rmsd_kras"][s:e]))
        b["rmsd_cypa"] = float(np.mean(out["rmsd_cypa"][s:e]))
        b["rg"] = float(np.mean(out["rg"][s:e]))
        if out["dmin_kras"] is not None:
            b["dmin_kras"] = float(np.mean(out["dmin_kras"][s:e]))
            b["dmin_cypa"] = float(np.mean(out["dmin_cypa"][s:e]))
            b["engage_kras"] = float(np.mean(out["dmin_kras"][s:e] < ENGAGE_CUT))
            b["engage_cypa"] = float(np.mean(out["dmin_cypa"][s:e] < ENGAGE_CUT))
        blocks.append(b)
    out["blocks"] = blocks
    # block-to-block spread of the headline metric = an honest uncertainty
    if len(blocks) > 1:
        out["rmsd_kras_block_sd"] = float(np.std([b["rmsd_kras"] for b in blocks]))
        out["engage_kras_block_sd"] = float(np.std([b.get("engage_kras", np.nan) for b in blocks]))
    else:
        out["rmsd_kras_block_sd"] = out["engage_kras_block_sd"] = np.nan
    print(f"  [{sysname}] {out['n_frames']} frames | KRAS RMSD "
          f"{out['stats']['rmsd_kras'][0]:.2f} A | drug-KRAS {out['stats']['dmin_kras'][0]:.2f} A "
          f"| engaged {out['engage_kras']*100:.0f}% KRAS / {out['engage_cypa']*100:.0f}% CypA")
    return out

# ----------------------------------------------------------------- per system
def figure_system(m, outdir):
    meta = SYS_META.get(m["name"], dict(label=m["name"], colour="#333333"))
    fig, axs = plt.subplots(2, 2, figsize=(9.5, 6.6))
    t = m["t_ns"]; eq = t[m["n_eq"]] if m["n_eq"] < len(t) else t[0]

    ax = axs[0, 0]
    ax.plot(t, m["rmsd_kras"], color="#1f77b4", lw=1.0, label="KRAS backbone")
    if m["rmsd_cypa"] is not None:
        ax.plot(t, m["rmsd_cypa"], color="#d62728", lw=1.0, label="CypA backbone")
    if m["rmsd_drug"] is not None:
        ax.plot(t, m["rmsd_drug"], color="#2ca02c", lw=1.0, label="daraxonrasib")
    ax.axvline(eq, color=GREY, ls=":", lw=0.8)
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("RMSD (Å)")
    ax.set_title("Structural drift vs. first frame"); ax.legend(loc="upper left")
    panel(ax, "a")

    ax = axs[0, 1]
    for tag, col, nm in (("kras", "#1f77b4", "KRAS"), ("cypa", "#d62728", "CypA")):
        if m[f"rmsf_{tag}"] is not None:
            ax.plot(m[f"resid_{tag}"], m[f"rmsf_{tag}"], color=col, lw=1.0, label=nm)
    shade_switches(ax, label=True)
    ax.set_xlabel("Residue number (crystallographic numbering, both chains)")
    ax.set_ylabel("Cα RMSF (Å)")
    ax.set_title("Per-residue flexibility (equilibrated part)"); ax.legend(loc="upper right")
    panel(ax, "b")

    ax = axs[1, 0]
    if m["dmin_kras"] is not None:
        ax.plot(t, m["dmin_kras"], color="#1f77b4", lw=1.0, label="drug–KRAS")
        ax.plot(t, m["dmin_cypa"], color="#d62728", lw=1.0, label="drug–CypA")
        ax.axhline(ENGAGE_CUT, color=GREY, ls="--", lw=0.8)
        ax.text(t[-1], ENGAGE_CUT, f" {ENGAGE_CUT} Å ", ha="right", va="bottom",
                fontsize=7, color=GREY)
        ax.set_ylabel("Minimum heavy-atom distance (Å)")
        ax.legend(loc="upper left")
    ax.set_xlabel("Time (ns)"); ax.set_title("Drug contact with each partner")
    panel(ax, "c")

    ax = axs[1, 1]
    if m["nc_kras"] is not None:
        ax.plot(t, m["nc_kras"], color="#1f77b4", lw=1.0, label="drug–KRAS")
        ax.plot(t, m["nc_cypa"], color="#d62728", lw=1.0, label="drug–CypA")
        ax.legend(loc="upper left")
    ax.set_xlabel("Time (ns)"); ax.set_ylabel(f"Heavy-atom contacts < {CONTACT_CUT} Å")
    ax.set_title("Interface contact count")
    panel(ax, "d")

    fig.suptitle(f"{m['name']} — {meta['label']}", fontsize=11, y=1.00)
    fig.tight_layout()
    p = outdir / f"md_analysis_{m['name']}.png"
    fig.savefig(p); plt.close(fig)
    return p

def figure_contacts(m, outdir, top_n=15):
    """Which residues does the glue actually touch, and how persistently."""
    occ = m.get("occ")
    if not occ:
        return None
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.9))
    for ax, tag, nm, col in ((axs[0], "kras", "KRAS", "#1f77b4"),
                             (axs[1], "cypa", "CypA", "#d62728")):
        o = occ[tag]
        order = np.argsort(-o["freq"])[:top_n]
        order = order[o["freq"][order] > 0]
        if not len(order):
            ax.text(0.5, 0.5, f"no drug contacts with {nm}", transform=ax.transAxes,
                    ha="center", va="center", color=GREY, fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
        else:
            key = "native" if "native" in o else "resids"
            lbl = [f"{o['names'][i].capitalize()}{o[key][i]}" for i in order]
            y = np.arange(len(order))[::-1]
            ax.barh(y, o["freq"][order] * 100, color=col, height=0.72)
            ax.set_yticks(y); ax.set_yticklabels(lbl, fontsize=7)
            ax.set_xlim(0, 105); ax.set_xlabel("Frames in contact (%)")
            for yi, v in zip(y, o["freq"][order] * 100):
                ax.text(v + 1.5, yi, f"{v:.0f}", va="center", fontsize=6.5, color=GREY)
        ax.set_title(f"daraxonrasib – {nm} contact residues (< {CONTACT_CUT} Å)")
        ax.grid(axis="y", alpha=0)
    panel(axs[0], "a"); panel(axs[1], "b")
    fig.suptitle(f"{m['name']} — residue-resolved interaction map", fontsize=11, y=1.02)
    fig.tight_layout()
    p = outdir / f"contacts_{m['name']}.png"
    fig.savefig(p); plt.close(fig)
    return p

def csv_contacts(m, outdir):
    import csv
    occ = m.get("occ")
    if not occ:
        return None
    p = outdir / f"contact_residues_{m['name']}.csv"
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chain", "resid", "resname", "contact_frequency_pct"])
        for tag, nm in (("kras", "KRAS"), ("cypa", "CypA")):
            o = occ[tag]
            key = "native" if "native" in o else "resids"
            for i in np.argsort(-o["freq"]):
                if o["freq"][i] <= 0:
                    break
                w.writerow([nm, int(o[key][i]), o["names"][i],
                            f"{o['freq'][i]*100:.1f}"])
    return p

def figure_convergence(m, outdir):
    """Has the run converged? Cumulative averages + independent block averages.

    Only meaningful for long trajectories: with a few hundred ps the blocks are
    correlated and a flat cumulative mean says nothing.
    """
    nm = m["name"]
    t, n_eq = m["t_ns"], m["n_eq"]
    blocks = m.get("blocks") or []
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.5))
    fig.suptitle(f"{nm} — convergence diagnostics", y=0.98)

    # (a) cumulative (running) mean of KRAS RMSD: flattens when converged
    ax = axs[0, 0]
    for key, col, lab in (("rmsd_kras", "#1f77b4", "KRAS backbone"),
                          ("rmsd_cypa", "#d62728", "CypA backbone")):
        y = m.get(key)
        if y is None:
            continue
        cum = np.cumsum(y) / np.arange(1, len(y) + 1)
        ax.plot(t, cum, color=col, lw=1.4, label=lab)
    ax.axvline(t[n_eq] if n_eq < len(t) else t[-1], color=GREY, ls=":", lw=1.0)
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("Cumulative mean RMSD (Å)")
    ax.set_title("Running mean RMSD (flat = converged)")
    ax.legend(frameon=False, fontsize=8)
    panel(ax, "a")

    # (b) per-block mean RMSD with the block spread as the error bar
    ax = axs[0, 1]
    if blocks:
        xb = [b["block"] for b in blocks]
        ax.plot(xb, [b["rmsd_kras"] for b in blocks], "o-", color="#1f77b4",
                lw=1.4, ms=5, label="KRAS")
        ax.plot(xb, [b["rmsd_cypa"] for b in blocks], "s-", color="#d62728",
                lw=1.4, ms=5, label="CypA")
        ax.set_xticks(xb)
        ax.set_xticklabels([f"{b['t_start_ns']:.0f}–{b['t_end_ns']:.0f}" for b in blocks],
                           rotation=30, ha="right", fontsize=7)
        ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel("Block (ns)"); ax.set_ylabel("Block mean RMSD (Å)")
    ax.set_title("Independent block averages")
    panel(ax, "b")

    # (c) engagement per block — the mechanistic readout, block by block
    ax = axs[1, 0]
    if blocks and "engage_kras" in blocks[0]:
        xb = np.arange(len(blocks)); w = 0.38
        ax.bar(xb - w / 2, [b["engage_kras"] * 100 for b in blocks], w,
               color="#1f77b4", label="drug–KRAS")
        ax.bar(xb + w / 2, [b["engage_cypa"] * 100 for b in blocks], w,
               color="#d62728", label="drug–CypA")
        ax.set_xticks(xb)
        ax.set_xticklabels([f"{b['t_start_ns']:.0f}–{b['t_end_ns']:.0f}" for b in blocks],
                           rotation=30, ha="right", fontsize=7)
        ax.set_ylim(0, 128)
        ax.legend(frameon=False, fontsize=8, loc="upper center", ncol=2)
    ax.set_xlabel("Block (ns)")
    ax.set_ylabel(f"Frames with contact < {ENGAGE_CUT} Å (%)")
    ax.set_title("Tri-complex engagement per block")
    panel(ax, "c")

    # (d) radius of gyration — a slow global coordinate; drift here means the
    #     run is still relaxing even if RMSD looks flat
    ax = axs[1, 1]
    ax.plot(t, m["rg"], color=GREY, lw=0.8)
    if len(t) > 50:
        w = max(5, len(t) // 50)
        k = np.ones(w) / w
        ax.plot(t[w - 1:], np.convolve(m["rg"], k, mode="valid"),
                color="#1f77b4", lw=1.6, label=f"{w}-frame mean")
        ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("Radius of gyration (nm)")
    ax.set_title("Complex compactness over time")
    panel(ax, "d")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    f = outdir / f"convergence_{nm}.png"
    fig.savefig(f, dpi=200); plt.close(fig)
    return f


def csv_blocks(m, outdir):
    """Per-block table — the numbers behind the convergence figure."""
    blocks = m.get("blocks") or []
    if not blocks:
        return None
    cols = ["block", "t_start_ns", "t_end_ns", "rmsd_kras", "rmsd_cypa", "rg",
            "dmin_kras", "dmin_cypa", "engage_kras", "engage_cypa"]
    f = outdir / f"blocks_{m['name']}.csv"
    with open(f, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for b in blocks:
            w.writerow([f"{b[c]:.3f}" if isinstance(b.get(c), float) else b.get(c, "")
                        for c in cols])
    return f


def csv_system(m, outdir):
    import csv
    p = outdir / f"timeseries_{m['name']}.csv"
    cols = ["time_ns", "rmsd_kras_A", "rmsd_cypa_A", "rmsd_drug_A", "rg_A",
            "dmin_drug_kras_A", "dmin_drug_cypa_A", "contacts_kras", "contacts_cypa"]
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols)
        n = len(m["t_ns"])
        def col(k): 
            v = m.get(k)
            return v if v is not None else [""] * n
        for row in zip(m["t_ns"], m["rmsd_kras"], col("rmsd_cypa"), col("rmsd_drug"),
                       m["rg"], col("dmin_kras"), col("dmin_cypa"),
                       col("nc_kras"), col("nc_cypa")):
            w.writerow([f"{x:.4f}" if isinstance(x, (int, float, np.floating)) else x for x in row])
    return p

# ------------------------------------------------------------- comparative
def figure_comparison(mets, outdir):
    fig, axs = plt.subplots(2, 2, figsize=(9.8, 6.8))

    ax = axs[0, 0]
    for m in mets:
        c = SYS_META.get(m["name"], {}).get("colour", "#333")
        ax.plot(m["t_ns"], m["rmsd_kras"], color=c, lw=1.0, label=m["name"])
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("KRAS backbone RMSD (Å)")
    ax.set_title("KRAS stability across systems"); ax.legend(loc="upper left", ncol=2)
    panel(ax, "a")

    ax = axs[0, 1]
    for m in mets:
        if m["rmsf_kras"] is None: continue
        c = SYS_META.get(m["name"], {}).get("colour", "#333")
        ax.plot(m["resid_kras"], m["rmsf_kras"], color=c, lw=1.0, label=m["name"])
    shade_switches(ax, label=True)
    ax.set_xlabel("KRAS residue number"); ax.set_ylabel("Cα RMSF (Å)")
    ax.set_title("KRAS flexibility (switch regions shaded)")
    panel(ax, "b")

    ax = axs[1, 0]
    for m in mets:
        if m["dmin_kras"] is None: continue
        c = SYS_META.get(m["name"], {}).get("colour", "#333")
        ax.plot(m["t_ns"], m["dmin_kras"], color=c, lw=1.0, label=m["name"])
    ax.axhline(ENGAGE_CUT, color=GREY, ls="--", lw=0.8)
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("Minimum drug–KRAS distance (Å)")
    ax.set_title("Minimum drug–KRAS distance over time")
    panel(ax, "c")

    ax = axs[1, 1]
    names = [m["name"] for m in mets]
    x = np.arange(len(names)); w = 0.38
    ek = [m["engage_kras"] * 100 for m in mets]
    ec = [m["engage_cypa"] * 100 for m in mets]
    ax.bar(x - w/2, ek, w, color="#1f77b4", label="drug–KRAS")
    ax.bar(x + w/2, ec, w, color="#d62728", label="drug–CypA")
    ax.set_xticks(x); ax.set_xticklabels(names)
    ax.set_ylabel(f"Frames with contact < {ENGAGE_CUT} Å (%)")
    ax.set_xlim(-0.6, len(names) - 0.4)
    ax.set_ylim(0, 125)
    ax.set_title("Tri-complex engagement")
    ax.legend(loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.0))
    for xi, (a_, b_) in enumerate(zip(ek, ec)):
        ax.text(xi - w/2, a_ + 2, f"{a_:.0f}", ha="center", va="bottom", fontsize=7)
        ax.text(xi + w/2, b_ + 2, f"{b_:.0f}", ha="center", va="bottom", fontsize=7)
    panel(ax, "d")

    fig.suptitle("daraxonrasib / KRAS / CypA tri-complex — comparison across systems",
                 fontsize=11, y=1.00)
    fig.tight_layout()
    p = outdir / "comparison_all_systems.png"
    fig.savefig(p); plt.close(fig)
    return p

def csv_summary(mets, outdir):
    import csv
    p = outdir / "summary_all_systems.csv"
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["system", "label", "ras_state", "frames",
                    "rmsd_kras_mean_A", "rmsd_kras_sd_A",
                    "rmsd_cypa_mean_A", "rmsd_drug_mean_A",
                    "rmsf_kras_mean_A", "rmsf_kras_max_A",
                    "rg_mean_A", "dmin_drug_kras_mean_A", "dmin_drug_cypa_mean_A",
                    "contacts_kras_mean", "contacts_cypa_mean",
                    "engaged_kras_pct", "engaged_cypa_pct",
                    "rmsd_kras_block_sd_A", "engaged_kras_block_sd_pct",
                    "n_blocks", "sim_time_ns"])
        for m in mets:
            s = m["stats"]; meta = SYS_META.get(m["name"], {})
            w.writerow([m["name"], meta.get("label", ""), meta.get("state", ""), m["n_frames"],
                        f"{s['rmsd_kras'][0]:.3f}", f"{s['rmsd_kras'][1]:.3f}",
                        f"{s['rmsd_cypa'][0]:.3f}", f"{s['rmsd_drug'][0]:.3f}",
                        f"{s['rmsf_kras'][0]:.3f}", f"{s['rmsf_kras'][1]:.3f}",
                        f"{s['rg'][0]:.3f}", f"{s['dmin_kras'][0]:.3f}", f"{s['dmin_cypa'][0]:.3f}",
                        f"{s['nc_kras'][0]:.1f}", f"{s['nc_cypa'][0]:.1f}",
                        f"{m['engage_kras']*100:.1f}", f"{m['engage_cypa']*100:.1f}",
                        f"{m.get('rmsd_kras_block_sd', float('nan')):.3f}",
                        f"{m.get('engage_kras_block_sd', float('nan'))*100:.1f}",
                        len(m.get("blocks") or []), f"{m['t_ns'][-1]:.1f}"])
    return p

def main():
    ap = argparse.ArgumentParser(description="Analyse daraxonrasib/KRAS MD runs")
    ap.add_argument("--root", default=".", type=Path)
    ap.add_argument("--systems", nargs="*", default=["9BG5", "9BG9", "4TQA", "8BLR", "4OBE"])
    ap.add_argument("--equil-frac", type=float, default=0.5,
                    help="fraction of the trajectory treated as equilibration (default 0.5)")
    ap.add_argument("--stride", type=int, default=1,
                    help="read every Nth frame of the trajectory (default 1)")
    ap.add_argument("--n-blocks", type=int, default=5,
                    help="number of independent blocks for the convergence test "
                         "(default 5; use 1 to disable)")
    ap.add_argument("--rmsf-max-frames", type=int, default=2000,
                    help="cap on frames held in memory for the RMSF superposition "
                         "(default 2000; raise it only if you have the RAM)")
    a = ap.parse_args()
    style()
    root = a.root.resolve()
    outroot = root / "analysis"; outroot.mkdir(exist_ok=True)
    mets = []
    print("Analysing trajectories:")
    for s in a.systems:
        m = analyse(s, root, a.equil_frac, a.stride, a.rmsf_max_frames,
                    a.n_blocks)
        if m is None: continue
        d = outroot / s; d.mkdir(parents=True, exist_ok=True)
        figure_system(m, d); csv_system(m, d)
        figure_contacts(m, d); csv_contacts(m, d)
        figure_convergence(m, d); csv_blocks(m, d)
        mets.append(m)
    if not mets:
        print("No trajectories found — run scripts/run_all.sh first."); return
    cmp_dir = outroot / "00_comparison"; cmp_dir.mkdir(exist_ok=True)
    figure_comparison(mets, cmp_dir)
    p = csv_summary(mets, cmp_dir)
    print(f"\nWrote {len(mets)} per-system reports + comparison in {outroot}")
    print(f"Summary table: {p}")

if __name__ == "__main__":
    main()
