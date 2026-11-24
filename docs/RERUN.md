# How to run this again

Merged and corrected from the two setup notes that came with the original project.
Both had drifted from what the code does: the dependency list named packages
nothing imports, and the structure step pointed at the wrong web server.

Stages 1 to 3 are pure Python and take a couple of minutes. Stage 4 needs two
external programs and a lot of manual clicking, which is the real bottleneck.

## What you need

```
pip install pandas numpy matplotlib jupyter
```

That is everything the notebooks import. The original note also listed meeko,
MDAnalysis, prolif, RDKit and gemmi, none of which appear anywhere in the code.
Don't install them unless you are adding that analysis yourself.

Also needed, separately:

- **Open Babel**, the `obabel` command-line tool, on your PATH. The notebook calls
  it through `os.system`.
- **AutoDock Vina 1.2.0**, from
  https://github.com/ccsb-scripps/AutoDock-Vina/releases. Not bundled here. Use
  1.2.0 specifically, since 1.1.2 is the version most papers cite and it scores
  differently, so the numbers here will not reproduce against it. You do not need
  `vina_split`.
- **Avogadro**, for adding hydrogens.

## Working folder

The notebooks expect one flat folder, which is not how this repo is laid out. The
folder names under `data/` match what they want, so you can copy them across:

```
work/
    vina.exe                  (or vina)
    04_docking.ipynb
    template/template.pdb     the IQS ligand
    dna/                      the 100 .pdb files, hydrogens added
    pdbqt/                    empty
    output/                   empty
    generated_sequences/      stage 1-3 notebooks read and write here
```

## Stages 1 to 3

Run the notebooks in order:

1. `01_sequence_generation.ipynb` writes `stage1_candidates.csv`, 10,000 sequences
2. `02_prefiltering.ipynb` writes `stage2_filtered.csv`, 2,000
3. `03_sequence_filtering.ipynb` writes `stage3_candidates.csv` (500),
   `stage3_for_docking.csv` (100) and `sequences_for_structure.txt`

`RANDOM_SEED = 41` in notebook 01 fixes the pool. Change it for a different pool,
leave it to reproduce this one exactly.

## Building the structures

This is 100 separate downloads and there is no way around it short of scripting
the server.

1. Open `sequences_for_structure.txt`, which lists the 100 sequences in rank order.
2. Go to https://web.x3dna-dssr.org and use the fiber model builder. The original
   note said web3dna.rutgers.edu, which is wrong: every PDB here carries a
   `REMARK DSSR` header, so they came from the DSSR platform, the current successor
   to the old w3DNA server.
3. Choose the B-DNA fiber model, paste one sequence, generate, download.
4. Repeat in rank order. Downloading one after another is what produces the
   `model (1).pdb` to `model (100).pdb` naming, and it is also what makes rank N
   equal file N. **That ordering is the only thing linking a score back to a
   sequence**, so do not reorder, re-download or rename anything midway. If you
   think the order slipped, `scripts/verify_pipeline.py` will tell you.
5. Open each PDB in Avogadro, add hydrogens, save.
6. Put the finished files in `dna/`.

There is an HTTP API at `http://api.x3dna.org/fiber` that would let you script
steps 2 to 4 and name files by sequence instead of by download order. If this
pipeline gets used again, do that. It removes the single most fragile step.

## Stage 4, docking

Run `04_docking.ipynb`. It converts the ligand to PDBQT with Open Babel (skipped if
`template.pdbqt` already exists), converts every `dna/*.pdb` with Gasteiger
charges, strips the torsion tree records (`ROOT`, `BRANCH`, `TORSDOF` and friends)
from the DNA PDBQT files, computes each receptor's centroid and docks with a 40 Å
cube centred there at exhaustiveness 32, 30 modes, 12 threads, then parses the best
`REMARK VINA RESULT` from each output into a ranking table.

The strip step is not optional. Open Babel writes PDBQT in ligand form and Vina
has to read these files as rigid receptors.

Expect a few hours for 100 receptors on 12 threads.

### Before a second run

Vina output accumulates rather than overwriting cleanly, so:

- rename or move the existing ranking file, otherwise it is overwritten and you
  lose the first run. This is where the A/B labelling in `data/rankings/` came from.
- **empty `output/`**. Stale pose files from a previous run get silently re-parsed
  and you end up with two runs mixed in one table.
- `pdbqt/` can stay, since conversion is skipped for files that already exist.

### Two things worth changing

Set a seed. Add `"--seed", "42"` to the command list. Without it the search is
stochastic and two runs disagree on the ranking, which is the problem in
LIMITATIONS item 5. A fixed seed makes one run reproducible, but it does not make
the scoring more accurate, so you still want replicates.

Run five replicates, not two, and select on mean and SD across them with
consistency as a criterion, the way the earlier PQS docking work was done. Two runs
is enough to show the ranking is unstable but not enough to average it out.

## Checking a run

```
python scripts/verify_pipeline.py
```

No dependencies, run from the repo root, about ten seconds. It checks the
sequence-to-structure mapping, recomputes every CSV value from the published
formulas, re-reads every Vina score from the pose files, and writes
`results/candidates_master.csv`, `results/top20_candidates.csv` and
`results/verification_report.txt`. Run it after any re-run. If the download order
slipped, check 1 fails and names the files.

## Line endings

Every data file here was written on Windows and has CRLF endings. The
`.gitattributes` tells Git to leave them alone, because the verification script
checks values against these files as they are. If you regenerate the data on Linux
or macOS you will get LF files, which is fine, but don't let Git normalise the
existing ones.
