# daraxonrasib / KRAS / CypA — complete GROMACS MD pipeline

Self-contained package to run molecular dynamics of the **daraxonrasib (RMC-6236)
tri-complex** with KRAS and cyclophilin A (CypA), from system build through
production MD to publication-ready analysis figures.

Configured for **200 ns of production per system** across five systems, with
resumable runs, mid-campaign snapshots, and block-averaged convergence
diagnostics.

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
│   ├── 02_run_md.sh            # EM -> NVT -> NPT -> production (MAXH-resumable)
│   ├── 03_postprocess.sh       # PBC treatment, frame0.pdb, movie
│   ├── 04_analysis.py          # figures + CSV tables (English text)
│   ├── 05_snapshot.sh          # analyse a run that is still in progress
│   ├── make_index.py           # builds index.ndx (Protein_LIG, Water_and_ions)
│   └── run_all.sh              # everything, all systems, one command
├── runs/<SYS>/                 # created by the scripts
└── analysis/                   # created by 04_analysis.py
```

---

## 4. Quick start

The default production length is **200 ns per system**. On CPU that is roughly
**6 days per system** — read section 6 before launching, and do the smoke test
first.

### Step 0 — smoke test (do this first, it costs 2 minutes)

Run a 20 ps production with short equilibration. It exercises every stage of the
pipeline, so if it finishes and produces figures, the 200 ns run will too.

```bash
cd kras_md_pipeline
NVT_PS=10 NPT_PS=10 bash scripts/01_build.sh 4TQA
NVT_PS=10 NPT_PS=10 bash scripts/02_run_md.sh 4TQA 0.02
python3 scripts/04_analysis.py --root . --systems 4TQA --equil-frac 0.0
```

Then clear the test before the real run:

```bash
rm -rf runs/4TQA analysis
```

### Step 1 — build all five systems (~5 minutes total)

```bash
for S in 9BG5 9BG9 4TQA 8BLR 4OBE; do bash scripts/01_build.sh $S; done
```

### Step 2 — launch production

**Recommended: one process per system, in parallel.** Five systems run
concurrently on 20 of your 24 cores; results accrue together instead of the
fifth system finishing a month after the first.

```bash
for S in 9BG5 9BG9 4TQA 8BLR 4OBE; do
    NT=4 nohup bash scripts/02_run_md.sh $S 200 > runs/$S/run.log 2>&1 &
done
```

The total of `NT` across concurrent `mdrun` processes must stay **at or below**
your core count — oversubscribing makes everything slower. With 24 cores, `NT=4`
x 5 systems leaves 4 cores for the machine.

Sequentially instead (simpler, but ~28 days before the last result):

```bash
bash scripts/run_all.sh 200
```

On a GPU, run one system at a time with all cores — GPU offload does not
parallelise across processes well:

```bash
for S in 9BG5 9BG9 4TQA 8BLR 4OBE; do GPU=1 bash scripts/02_run_md.sh $S 200; done
```

### Step 3 — watch it develop

`02_run_md.sh` calls postprocessing and stops there; the analysis is a separate
step you can run at any time on a **partial** trajectory:

```bash
bash scripts/05_snapshot.sh                    # all systems, whatever exists
bash scripts/05_snapshot.sh 9BG5 9BG9          # just these two
```

This is safe while `mdrun` is still writing. Check progress from the log:

```bash
tail -2 runs/9BG5/md.log
grep -c "^" runs/9BG5/run.log
```

### Step 4 — final analysis

```bash
python3 scripts/04_analysis.py --root . --systems 9BG5 9BG9 4TQA 8BLR 4OBE
```

### Pausing a multi-day run

`MAXH` gives `mdrun` a wall-clock budget in hours; it stops cleanly just before
it, writing `md.cpt`. Re-issue the same command to continue. A 200 ns run split
over ten sessions is equivalent to one uninterrupted run.

```bash
MAXH=12 bash scripts/02_run_md.sh 9BG5 200      # today
MAXH=12 bash scripts/02_run_md.sh 9BG5 200      # tomorrow, continues
```

Re-running any command is always safe: finished stages are skipped, an
unfinished production resumes from its checkpoint.

---

## 5. Environment variables

| Variable | Default | Meaning |
|----------|---------|---------|
| `GMX`     | `gmx`  | GROMACS binary |
| `NT`      | all cores | threads for `mdrun` |
| `GPU`     | `0`    | `GPU=1` adds `-nb gpu -bonded gpu -pme gpu` |
| `NVT_PS`  | `100`  | NVT equilibration length (ps) |
| `NPT_PS`  | `100`  | NPT equilibration length (ps) |
| `MAXH`    | unset  | wall-clock budget in **hours** for production; stops cleanly and writes `md.cpt` |
| `MOVIE_FRAMES` | `200` | target number of models in `trajectory_movie.pdb` |
| `EQUIL_FRAC` | `0.5` | passed to the analysis by `05_snapshot.sh` |
| `BOXD`    | `1.2`  | minimum solute-to-box distance (nm) |
| `CONC`    | `0.15` | NaCl concentration (M) |

Example — 200 ns on a GPU with 16 threads:

```bash
GPU=1 NT=16 bash scripts/02_run_md.sh 9BG5 200
```

`02_run_md.sh` **skips any stage whose output already exists** and resumes an
interrupted production run from `md.cpt`, so re-running it after a crash or a
`Ctrl-C` is safe and cheap.

---

## 6. Expected runtimes

Measured on 24 CPU cores (no GPU), ~59,000-67,000 atoms depending on the system.
`02_run_md.sh` prints a forecast for your requested length before it starts.

| Stage | Wall time |
|-------|-----------|
| Build (`01_build.sh`) | < 1 min |
| EM | ~1-2 min |
| NVT 100 ps | ~5 min |
| NPT 100 ps | ~5 min |
| Production, per ns | ~40 min (29-36 ns/day, all 24 cores) |

**200 ns campaign:**

| Layout | Wall time |
|--------|-----------|
| One system, all 24 cores | **~6 days** |
| Five systems sequentially | ~28-30 days |
| Five systems in parallel, `NT=4` each | **~10-12 days** (all five finish together) |
| Five systems, single consumer GPU, one at a time | ~4-8 days total |

Per-process throughput drops with fewer threads, so five parallel `NT=4` runs
are not five times faster than sequential — but every system advances at once,
which is usually worth more than the total-throughput difference.

A GPU is the single biggest win here: production typically runs 5-15x faster.
Your machine reports an **NVIDIA RTX A4000**, which should give roughly
150-300 ns/day on a system this size — about **1 day per system** instead of six.
Confirm the build supports it:

```bash
gmx --version | grep -i -E "GPU|CUDA"
nvidia-smi
```

If that shows CUDA support, use `GPU=1` and run the systems one at a time.

### Disk space

At 20 ps per frame, 200 ns is 10,000 frames:

| File | Per system | Five systems |
|------|-----------|--------------|
| `md.xtc` (all atoms, compressed) | ~4 GB | ~20 GB |
| `md_clean.xtc` (solute only) | ~0.35 GB | ~1.7 GB |
| `md.edr`, `md.log`, checkpoints | ~0.5 GB | ~2.5 GB |
| **Total** | **~5 GB** | **~25 GB** |

Budget ~40 GB to be comfortable. `02_run_md.sh` prints the same estimate for the
length you actually request, before it starts. If space is tight, raise
`nstxout-compressed` in `mdp/md.mdp` (40000 = 80 ps per frame quarters the
trajectory size; 5,000 frames is still ample for every analysis here).
`md.xtc` is the only file you cannot regenerate — `md_clean.xtc`, the figures and
the CSVs all rebuild from it in minutes.

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
| `trajectory_movie.pdb` | multi-model PDB, ~200 models regardless of run length (`MOVIE_FRAMES`) |
| `pbc_check.txt` | PBC sanity report |

Per system, in `analysis/<SYS>/`:

| File | Contents |
|------|----------|
| `md_analysis_<SYS>.png` | 4 panels: RMSD, per-residue RMSF, drug-partner distances, contact counts |
| `contacts_<SYS>.png` | residue-resolved interaction map, top 15 per chain |
| `convergence_<SYS>.png` | running means, independent block averages, engagement per block, Rg |
| `timeseries_<SYS>.csv` | every time series, one row per frame |
| `contact_residues_<SYS>.csv` | contact frequency per residue (crystallographic numbering) |
| `blocks_<SYS>.csv` | per-block means — the numbers behind the convergence figure |

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
    --stride 1 \
    --n-blocks 5
```

| Flag | Default | Meaning |
|------|---------|---------|
| `--root` | `.` | package root |
| `--systems` | all five | subset to analyse |
| `--equil-frac` | `0.5` | leading fraction of the trajectory discarded as equilibration; RMSF and all averages use the remainder |
| `--stride` | `1` | read every Nth frame (use 2-5 for long trajectories) |
| `--n-blocks` | `5` | independent blocks for the convergence test (`1` disables it) |
| `--rmsf-max-frames` | `2000` | cap on frames held in memory for the RMSF superposition |

Missing systems are skipped with a message, so the comparison figure works with
however many runs you have finished.

RMSF requires the trajectory in memory to superpose it, at
`n_frames x n_atoms x 12` bytes. For 200 ns the equilibrated half is subsampled
to at most `--rmsf-max-frames` before the transfer — RMSF is an average over
frames and is insensitive to the sampling interval, whereas an unbounded
transfer would need several GB. The script reports the subsampling when it
applies. Every other metric still uses **every** frame.

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
* **Convergence, at 200 ns, is testable** — and `convergence_<SYS>.png` is where
  to look before quoting any average:
  * panel (a) the running mean flattens when a metric has converged; if it is
    still climbing at the end, your averages are window-dependent;
  * panel (b) the five independent block averages. Their **spread is your
    uncertainty** — it is reported as `rmsd_kras_block_sd_A` in the summary
    table. Quote means as `value +/- block SD`, not as bare numbers;
  * panel (c) engagement block by block. A tri-complex that holds should be near
    100% in *every* block; a block where it drops is a dissociation event worth
    inspecting in the movie;
  * panel (d) radius of gyration is a slow global coordinate. Drift here means
    the system is still relaxing even when RMSD looks flat.

### What 200 ns does and does not buy you

200 ns per system is enough to state that a complex is stable on that timescale,
to resolve which interface contacts are persistent versus transient, to see
switch-I/switch-II reorganisation, and — via the block analysis — to attach an
uncertainty to each average rather than quoting a single number.

It is still **one replica per system**. A single trajectory samples one basin
well and tells you nothing about how reproducible the behaviour is: an event
seen once at 150 ns could be routine or could be an artefact of this particular
set of initial velocities. The block SD measures sampling *within* a trajectory,
which is a weaker statement than agreement *between* trajectories.

Recommended next steps, in order of value per CPU-hour:

1. **Three replicas x 200 ns per system.** Copy the built system and change
   `gen-seed` in `mdp/nvt.mdp` (or set `gen_vel = yes` with a fresh seed) so each
   replica gets independent starting velocities. Agreement across replicas is
   what makes the RAS-ON vs RAS-OFF contrast publishable.
2. **MM/PBSA or alchemical free energy** if you need binding energetics — no
   amount of plain MD gives you an affinity.
3. **Hydrogen-bond analysis at the interface**, which the contact CSVs already
   set up: `gmx hbond -s complex.tpr -f md_clean.xtc -n index.ndx`.
4. **Cluster the interface conformations** to identify the dominant binding modes
   over the 200 ns rather than describing only the average.
