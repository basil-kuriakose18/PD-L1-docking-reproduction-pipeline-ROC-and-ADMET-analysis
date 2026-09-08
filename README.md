# PD-L1 Flavonoid Docking Reproduction

This project reproduces and evaluates the docking, ROC validation, and ADMET analysis from:

> Jiang, D.; Kwon, H.-K.; Kwon, O.W.; Choi, Y. *A Comparative Molecular Dynamics Study of
> Food-Derived Compounds as PD-L1 Inhibitors: Insights Across Six Flavonoid Subgroups.*
> Molecules 2025, 30(4), 907. https://doi.org/10.3390/molecules30040907
> (open access, CC BY 4.0)

All 62 compounds from the paper's Table 1 — 60 flavonoids plus the two reference PD-L1
inhibitors BMS-202 and BMS-1166 — were docked against PDB 7DY7 with AutoDock Vina. The
results are compared statistically against the published affinities, validated with an ROC
analysis against real DUD-E decoys, and checked against the paper's ADMET predictions for
its six highlighted compounds.

Molecular dynamics is not part of this project. I didn't have access to a machine that could
run GROMACS, so the scope here is docking, statistical comparison, ROC validation, and ADMET
only.

## Results

Redocking the native co-crystallized ligand ("Compound 17") into 7DY7 using the paper's grid
box gave a heavy-atom RMSD of 1.63 Å to the crystal pose — under the ~2 Å threshold generally
used to call a docking setup validated.

Across the full 62-ligand set, my docking scores correlate with the published affinities at
r = 0.67 (Pearson), ρ = 0.54 (Spearman), with an RMSE of 1.09 kcal/mol. Isoflavones stood out
as the subgroup with the largest systematic deviation — my scores ran about 1.11 kcal/mol
weaker on average across all 10 isoflavones — though I didn't track down why.

The ROC analysis is a mixed result, and the two numbers mean different things. Against the
full 62-compound set and 99 successfully docked DUD-E decoys, AUC came out to 0.53 —
essentially no discrimination. But that test isn't really comparable to the paper's own
design (more on this below), and a cleaner version limited to the two experimentally
confirmed inhibitors against their correctly matched decoys gave AUC = 0.97. That's a strong
number, but it's only two positive cases, so treat it as suggestive rather than conclusive.

ADMET came out well. Reproducing the paper's Table 3 via SwissADME and pkCSM matched exactly
for four of the six highlighted compounds (Diosmin, Formononetin, Idaein, Neohesperidin),
several values agreeing to two or three decimal places. Ginkgetin and Theaflavin each had
one or two values that didn't match — details below.

## Repository layout

```
data/
  published/    Table 1 transcribed from the source paper (CSV)
  raw/          PDB structures, ligand lookups, DUD-E decoy files, other raw inputs
docking/
  config/       Vina grid box configuration
  receptor_pdbqt/
  ligands_pdbqt/  62 ligand + 100 decoy PDBQT files
R/              Analysis scripts
results/
  figures/      Saved plots
  *.csv         Docking scores, affinity comparison, ADMET comparison
renv.lock       Locked R package versions
```

## Requirements

- R 4.6+ with `readr`, `dplyr`, `bio3d`, `ggplot2`, `ggrepel`, `pROC`, `webchem`, and `renv`
  (run `renv::restore()` for exact package versions)
- AutoDock Vina (this project used the standalone Scripps Research build)
- OpenBabel 3.1.1 — Windows users, see the note at the bottom before running `--gen3d`

## Where this diverges from the paper

I wanted this section to be honest rather than tidy, so here's everywhere this reproduction
doesn't match the original study, and why.

No MD. I didn't have a GROMACS-capable machine for this project, so the paper's molecular
dynamics work (AMBER, 100 ns × 3 replicates per compound) isn't reproduced here at all.

The paper doesn't state what exhaustiveness value it used for Vina, so I used the default
(8) and wrote it into the config explicitly rather than leaving it implicit.

Receptor preparation used different tools than the paper. The paper cleaned its structure in
PyMOL and assigned receptor charges with MGLTools 1.5.7. I used `bio3d` for cleaning and
OpenBabel for the whole receptor PDBQT conversion, which threw a "failed to kekulize aromatic
bonds" warning I never fully chased down — it's possible this adds some scoring noise, but I
didn't isolate it as the cause of anything specific.

The ROC decoys aren't set up the way the paper's were. Their methods describe using all 60
flavonoids as actives against 60 DUD-E decoys, roughly a 1:1 ratio. I only submitted the two
confirmed reference inhibitors to DUD-E, mainly to avoid the ~3,000 decoys that would come
from running all 60 actives at DUD-E's default 50-per-active rate. That means my two ROC
numbers aren't measuring quite the same thing as the paper's: the full-set AUC (0.53) is
being tested against decoys that were never property-matched to the flavonoids' chemistry,
while the BMS-only AUC (0.97) is the cleaner test but rests on just two actives.

Two separate data problems turned up in the anthocyanidins while building this. Four
compounds — Delphinidin, Malvidin, Pelargonidin, Rosinidin — failed PubChem's CID lookup by
CAS number and had to be resolved by name instead. A different four — Peonidin, Petunidin,
Callistephin, Idaein — came back from PubChem as salt-form SMILES (with a trailing
`.[Cl-]`), which silently produced broken single-atom PDBQT files until I stripped the
counterion out.

On the ADMET side, pkCSM's own site currently shows a banner for a newer prediction engine
("Deep-PK"), which suggests the underlying models may have changed since the paper was
submitted in December 2024. That's my best guess for why Ginkgetin (CYP3A4 substrate call,
Lipinski count) and Theaflavin (BBB permeability) didn't match — I haven't confirmed it, just
noting it as the likely explanation.

Finally, exact numerical agreement with Vina was never really on the table to begin with —
it's a stochastic algorithm, so some spread between runs is expected regardless of how
carefully the setup matches.

## Running this yourself

1. `renv::restore()` to get the locked package versions.
2. Run `R/00_explore_data.R` from the top — it covers structure retrieval, PDBQT conversion,
   redocking validation, the full docking batch, the statistical comparison, ROC, and ADMET,
   in that order.
3. On Windows, if `--gen3d` fails with an "all zero coordinates" error, OpenBabel's
   `rigid-fragments.txt` (usually under `%APPDATA%\OpenBabel-3.1.1\data\`) has a line-ending
   bug. Rewrite it with Unix line endings and the error goes away.
