# daraxonrasib / KRAS / CypA — complete GROMACS MD pipeline

Self-contained package to run molecular dynamics of the **daraxonrasib (RMC-6236)
tri-complex** with KRAS and cyclophilin A (CypA), from system build through
production MD to publication-ready analysis figures.

Everything runs **locally**. No internet access is required at any stage.

---

## 1. Mechanistic background

Daraxonrasib is a **non-covalent tri-complex RAS(ON) inhibitor** ("molecular
glue"). It does not bind RAS on its own: it binds **cyclophilin A (PPIA,
UniProt P62937)** first, and the resulting binary CypA-drug complex then engages
the **switch-I and switch-II** regions of **active (GTP-bound) RAS**. Any
simulation that omits CypA cannot reproduce the pharmacology.

Five systems are provided, chosen to test that mechanism rather than to repeat it:

| ID   | KRAS   | State   | Nucleotide | Mg2+ | Drug pose | Purpose |
|------|--------|---------|------------|------|-----------|---------|
| 9BG5 | G13D   | RAS-ON  | GNP (GppNHp) | yes | crystallographic | primary system of interest |
| 9BG9 | WT     | RAS-ON  | GNP          | yes | crystallographic | mutation control (G13D vs WT) |
| 4TQA | G13D   | RAS-OFF | GDP          | yes | docked           | nucleotide-state control |
| 8BLR | G13D   | RAS-OFF | GDP          | no  | docked           | open switch-I, no Mg2+ |
| 4OBE | WT     | RAS-OFF | GDP          | yes | docked           | double control (WT + OFF) |

The RAS-ON systems (9BG5, 9BG9) carry the drug in its **experimentally observed**
pose. The RAS-OFF systems carry a **docked** pose and are expected to be less
stable — that contrast is the scientific readout, not a defect.

**GNP/GDP and Mg2+ are part of the model.** They are retained in every system.
A pipeline that strips HETATM records (e.g. `grep -v HETATM`) silently deletes
the nucleotide and the magnesium ion, which destroys the P-loop and the switch
conformations. Do not do that.

---

## 2. Requirements

| Software | Version | Notes |
|----------|---------|-------|
| GROMACS  | >= 2021 | `gmx` on `PATH`, or set `GMX=/path/to/gmx` |
| Python   | >= 3.9  | analysis only |
| numpy, matplotlib, MDAnalysis | any recent | `pip install numpy matplotlib MDAnalysis` |

Force field: **AMBER99SB-ILDN** (protein) + **GAFF2** (drug, GNP, GDP) + **TIP3P**
water. The GAFF2 ligand topologies are pre-generated in `toppar/` — AmberTools /
ACPYPE are **not** needed to run this package.

Check your installation:

```bash
gmx --version          # or: $GMX --version
python3 -c "import MDAnalysis, numpy, matplotlib; print('analysis deps OK')"
```

---

## 3. Directory layout

```
kras_md_pipeline/
├── README.md
├── toppar/                     # force-field fragments, shared by all systems
│   ├── ligand_atomtypes.itp    # unified GAFF2 [atomtypes] (included once)
│   ├── DRG.itp                 # daraxonrasib, 116 atoms, net charge 0
│   ├── GNP.itp                 # GppNHp, 45 atoms, net charge -4
│   └── GDP.itp                 # GDP, 40 atoms, net charge -3
├── inputs/<SYS>/               # one directory per system, ready for grompp
│   ├── protein_only.pdb        # crystallographic numbering (used by analysis)
│   ├── system_raw.gro          # protein + drug + nucleotide + Mg2+
│   ├── topol.top
│   ├── topol_Protein_chain_A.itp / _B.itp
│   └── posre_*.itp             # position restraints (protein, drug, nucleotide)
├── mdp/                        # em / nvt / npt / md
├── scripts/
│   ├── 01_build.sh             # box, solvate, ions, grompp check
│   ├── 02_run_md.sh            # EM -> NVT -> NPT -> production
│   ├── 03_postprocess.sh       # PBC treatment, frame0.pdb, movie
│   ├── 04_analysis.py          # figures + CSV tables (English text)
│   ├── make_index.py           # builds index.ndx (Protein_LIG, Water_and_ions)
│   └── run_all.sh              # everything, all systems, one command
├── runs/<SYS>/                 # created by the scripts
└── analysis/                   # created by 04_analysis.py
```

---

## 4. Quick start

```bash
cd kras_md_pipeline
bash scripts/run_all.sh                       # all five systems, 5 ns each
```

Or one system at a time:

```bash
bash scripts/01_build.sh      9BG5             # solvate + ionise + grompp check
bash scripts/02_run_md.sh     9BG5             # EM, NVT, NPT, production
bash scripts/03_postprocess.sh 9BG5            # PBC fix, frame0.pdb, movie
python3 scripts/04_analysis.py --root . --systems 9BG5
```

Comparison figure across everything you have run:

```bash
python3 scripts/04_analysis.py --root . --systems 9BG5 9BG9 4TQA 8BLR 4OBE
```

### Smoke test first (recommended)

Run a 20 ps production with short equilibration before committing hours of CPU:

```bash
NVT_PS=10 NPT_PS=10 bash scripts/02_run_md.sh 4TQA 0.02
bash scripts/03_postprocess.sh 4TQA
python3 scripts/04_analysis.py --root . --systems 4TQA --equil-frac 0.0
```

The second argument of `02_run_md.sh` is the production length **in ns**.
If it finishes and the analysis produces figures, the full run will too.

---

## 5. Environment variables

| Variable | Default | Meaning |
|----------|---------|---------|
| `GMX`     | `gmx`  | GROMACS binary |
| `NT`      | all cores | threads for `mdrun` |
| `GPU`     | `0`    | `GPU=1` adds `-nb gpu -bonded gpu -pme gpu` |
| `NVT_PS`  | `100`  | NVT equilibration length (ps) |
| `NPT_PS`  | `100`  | NPT equilibration length (ps) |
| `BOXD`    | `1.2`  | minimum solute-to-box distance (nm) |
| `CONC`    | `0.15` | NaCl concentration (M) |

Example — 10 ns on a GPU with 16 threads:

```bash
GPU=1 NT=16 bash scripts/02_run_md.sh 9BG5 10
```

`02_run_md.sh` **skips any stage whose output already exists** and resumes an
interrupted production run from `md.cpt`, so re-running it after a crash or a
`Ctrl-C` is safe and cheap.

---

## 6. Expected runtimes

Measured on 24 CPU cores (no GPU), ~59,000-67,000 atoms depending on the system:

| Stage | Wall time |
|-------|-----------|
| Build (`01_build.sh`) | < 1 min |
| EM | ~1-2 min |
| NVT 100 ps | ~5 min |
| NPT 100 ps | ~5 min |
| Production, per ns | ~40 min (29-36 ns/day) |
| **5 ns production** | **~3.4 h** |
| All five systems, 5 ns each | ~18 h |

A single consumer GPU typically gives a 5-15x speed-up on the production stage.
Systems are independent — run them in parallel if you have the cores, but give
each `mdrun` its own `-nt` slice rather than oversubscribing.

---

## 7. Periodic-boundary treatment (read this before trusting a distance)

Both KRAS and CypA are **separate molecules**. Under periodic boundary
conditions two independent failures are possible: a single chain can be split
across the box edge, and the two chains can be placed in opposite corners of
the box. Either one produces nonsense interface distances — a drug-CypA minimum
distance of ~56 A when the true value is ~2 A.

`03_postprocess.sh` applies **exactly three** `trjconv` stages, in this order:

1. `-pbc whole` — reassemble each molecule that was split across a boundary.
2. `-pbc cluster` (centred on `Protein_LIG`) — bring KRAS, CypA and the drug
   into the same periodic image as one cluster.
3. `-fit rot+trans` (fit on `Protein`) — remove global rotation and translation.

Two steps that seem natural are **deliberately omitted**:

* **`-pbc nojump` after clustering** re-splits CypA. Verified by radius of
  gyration: CypA Rg jumps from 1.45 nm to 3.66 nm, and the drug-CypA distance
  then reads 1.74 A against a *periodic-image fragment* — a false positive that
  looks like a good result.
* **`-center ... -pbc mol -ur compact`** re-wraps every molecule about its own
  centre of mass, which throws the two chains back into opposite corners.
  Stage 3's `-fit` already centres the system.

### The `pbc_check.txt` guard

`03_postprocess.sh` writes `runs/<SYS>/pbc_check.txt`, which compares the radius
of gyration of KRAS, CypA and the whole protein **before** MD (from
`solv_ions.gro`) against the **post-PBC** trajectory, and prints the drug-KRAS
and drug-CypA distances. A chain whose Rg grew by more than **1.15x** is flagged
`WARNING`. **Always read this file before interpreting any figure.**

Reference values from the validated 4TQA run:

```
KRAS    1.537 nm (pre-MD 1.537)     drug-KRAS ~2.1 A
CypA    1.447 nm (pre-MD 1.447)     drug-CypA ~2.2 A
Protein 2.217 nm (pre-MD 2.226)     -> PASSED
```

A `Protein` Rg near 3.4 nm instead of ~2.2 nm means the PBC treatment failed.

---

## 8. Outputs

Per system, in `runs/<SYS>/`:

| File | Contents |
|------|----------|
| `md.xtc`, `md.tpr`, `md.edr`, `md.log` | raw production trajectory |
| `md_clean.xtc` | PBC-corrected, fitted trajectory (**use this**) |
| `complex.tpr` | protein + drug + nucleotide only (matches `md_clean.xtc`) |
| `frame0.pdb` | first frame, with CONECT records |
| `trajectory_movie.pdb` | multi-model PDB (every 5th frame) for PyMOL/VMD/ChimeraX |
| `pbc_check.txt` | PBC sanity report |

Per system, in `analysis/<SYS>/`:

| File | Contents |
|------|----------|
| `md_analysis_<SYS>.png` | 4 panels: RMSD, per-residue RMSF, drug-partner distances, contact counts |
| `contacts_<SYS>.png` | residue-resolved interaction map, top 15 per chain |
| `timeseries_<SYS>.csv` | every time series, one row per frame |
| `contact_residues_<SYS>.csv` | contact frequency per residue (crystallographic numbering) |

Across systems, in `analysis/00_comparison/`:

| File | Contents |
|------|----------|
| `comparison_all_systems.png` | RMSD, RMSF, drug-KRAS distance, engagement bars |
| `summary_all_systems.csv` | one row per system: means, SDs, engagement percentages |

All figure and axis text is in **English**.

### Residue numbering

`pdb2gmx` renumbers each chain from 1 and the `.tpr` concatenates them, so a raw
MDAnalysis resid is not interpretable (CypA Arg55 appears as residue 223).
`04_analysis.py` reads `inputs/<SYS>/protein_only.pdb` and restores the
**crystallographic numbering** in every figure and CSV. Keep that file in place.

Validation: for 4TQA the contact map returns CypA **Arg55, Phe60, Met61, Asn102,
Ala103, Trp121, Lys125, His126, Asn149** — the canonical cyclophilin active site
— and KRAS **Thr35, Glu37** (switch-I) with **Ala59, Gln61, Glu62, Glu63**
(switch-II). That is the expected molecular-glue interface, and it is a good
first thing to confirm in your own runs.

---

## 9. Analysis options

```bash
python3 scripts/04_analysis.py --root . \
    --systems 9BG5 9BG9 4TQA 8BLR 4OBE \
    --equil-frac 0.5 \
    --stride 1
```

| Flag | Default | Meaning |
|------|---------|---------|
| `--root` | `.` | package root |
| `--systems` | all five | subset to analyse |
| `--equil-frac` | `0.5` | leading fraction of the trajectory discarded as equilibration; RMSF and all averages use the remainder |
| `--stride` | `1` | read every Nth frame (use 2-5 for long trajectories) |

Missing systems are skipped with a message, so the comparison figure works with
however many runs you have finished.

---

## 10. Interpreting the results

* **Engagement percentages** (panel d of the comparison figure) are the primary
  readout: the fraction of frames with any heavy-atom contact below 4.5 A to
  KRAS and to CypA. A genuine tri-complex holds **both** near 100%. Losing the
  CypA arm means the molecular-glue mode has broken down.
* **RMSD** measures drift from the starting structure, not correctness. A docked
  pose is expected to drift more than a crystallographic one; report the value
  rather than treating a rise as failure.
* **RMSF with the switch regions shaded** shows where the mutation acts.
  Switch-I is residues 30-40 and switch-II 60-76, both shaded in the figures.
* **Contact-frequency maps** are the most transferable result — they can be
  compared directly against the crystallographic interface.

### Caveats

Single replica per system, 5 ns default. This is enough to check that a complex
holds together and to see the interface, and **not** enough for free energies,
converged conformational sampling, or a claim about relative affinity. For
publication: run 3 independent replicas with different velocity seeds
(`gen-seed` in `nvt.mdp`), extend production to at least 100 ns, and add
MM/PBSA or an alchemical free-energy calculation if you need energetics.
