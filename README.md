# dsDNA candidate screening for IQS binding

Screening pipeline that takes 10,000 random DNA sequences down to a shortlist of
20 double-stranded candidates for binding **IQS**
(2-(2-hydroxyphenyl)-thiazole-4-carbaldehyde, also called aeruginaldehyde), one of
the quorum sensing signals of *Pseudomonas aeruginosa*. The intended use is a
dsDNA-based electrochemical sensor, with the duplex sitting on the electrode and
the small molecule binding in the groove.

Three Python filtering stages, then AutoDock Vina docking of the 100 survivors.

**The output is a shortlist, not a ranking.** All 100 docked candidates fall
inside 0.478 kcal/mol of each other, and two independent Vina runs disagree about
the order: they share 13 of their top 20 and pick different winners. Claiming
candidate 1 binds better than candidate 15 is not supported by this data.
[docs/LIMITATIONS.md](docs/LIMITATIONS.md) covers that and the rest, and is the
file to read first if you want to build on any of this.

## Stages

| Stage | What it does | Tool | In | Out |
|---|---|---|---|---|
| 1 | Random generation, five hard filters (length, GC, homopolymer, base dominance, palindrome) | Python | 35,570 draws | 10,000 |
| 2 | Score: Wallace Tm (40%) + IQS interaction heuristic (60%) | Python | 10,000 | 2,000 |
| 3.1 | SantaLucia nearest-neighbour ΔG filter | Python | 2,000 | 1,654 |
| 3.2 | Levenshtein diversity filter, edit distance ≥ 4 | Python | 1,654 | 500 |
| 3.3 | Top 100 by stage-2 score, for structure building | Python | 500 | 100 |
| 4 | Docking, two independent runs | Vina 1.2.0 | 100 | 20 |

Stage 1 rejections were almost all GC content: 23,349 out of 25,570 rejected
draws. The other four filters removed about 2,200 between them.

## Shortlist

Top five by the mean of the two runs. Full 20 in
[`results/top20_candidates.csv`](results/top20_candidates.csv), all 100 with every
property in [`results/candidates_master.csv`](results/candidates_master.csv).

| Sequence (5'→3') | bp | GC | ΔG | Vina mean | Top 20 in both runs |
|---|---|---|---|---|---|
| GAGTGCCAGCGTGAGAGGGAA | 21 | 62% | −34.13 | −6.939 | yes |
| GACGGACACGGTGGGAACAAC | 21 | 62% | −33.30 | −6.929 | yes |
| CTAGAGGTCGACGGGATGAC | 20 | 60% | −30.33 | −6.884 | yes |
| CTGGCCCAAATAACGGGAAGGG | 22 | 59% | −34.37 | −6.878 | yes |
| CGGCGACGAAGGGAGTTCCT | 20 | 65% | −33.01 | −6.872 | yes |

13 of the 20 keep a top-20 place in both runs. Those are the ones worth taking to
the bench. The other 7 sit inside the run-to-run noise.

## Checking it

The 100 structures were downloaded one at a time and ended up named
`model (1).pdb` to `model (100).pdb` in download order. Nothing in a filename says
which sequence is inside, so the whole chain from sequence to score rests on that
order being right. Rather than ask anyone to take that on trust:

```
python scripts/verify_pipeline.py
```

No dependencies, about ten seconds. It reads the bases back out of every PDB atom
record, recomputes every number in the stage CSVs from the sequence alone,
re-reads every Vina score from the pose files, and writes the two result tables
plus `results/verification_report.txt`.

All 10 checks pass. What that confirms: `model (N).pdb` really is stage-3 rank N
for all 100, chain B really is the exact reverse complement of chain A in all 100,
and every stored value recomputes from the sequence (ΔG and raw IQS score to zero
deviation).

## Layout

```
data/sequences/   stage 1-3 CSVs and the sequence list used for structure building
data/dna/         100 B-form duplexes, hydrogens added
data/pdbqt/       the same, converted to Vina receptors
data/template/    the IQS ligand
data/output/      100 Vina pose files, 28-30 poses each
data/rankings/    the two run tables, and PROVENANCE.txt on which is which
notebooks/        stages 1-3 and the docking driver, as they were run
results/          master table, shortlist, verification report
scripts/          verify_pipeline.py
figures/          stage plots, and the Origin export notebooks
docs/             METHODS, LIMITATIONS, RERUN
```

Folder names under `data/` match what the notebooks expect, so you can copy them
into a flat working folder and re-run without editing paths. Notebooks keep their
stored outputs, so you can read what each stage produced without running anything.
Code cells are untouched from the original run; only markdown text was fixed where
it contradicted the code.

## Software

Not bundled, install yourself:

- **AutoDock Vina 1.2.0**, from https://github.com/ccsb-scripps/AutoDock-Vina/releases.
  Use 1.2.0 specifically. Most papers cite 1.1.2 and it scores differently.
- **Open Babel**, for PDB to PDBQT conversion
- **Avogadro**, for adding hydrogens
- **Python 3** with pandas, numpy, matplotlib, jupyter

Structures came from the fiber model builder at https://web.x3dna-dssr.org.
Procedure in [docs/RERUN.md](docs/RERUN.md).

## References

- Lee et al. 2013, *Nat Chem Biol* 9:339-343. IQS identified as a QS signal.
  [doi:10.1038/nchembio.1225](https://doi.org/10.1038/nchembio.1225)
- Raya, Montagut & Marco 2025, *Front Cell Infect Microbiol* 15:1575421. IQS and
  aeruginaldehyde, and the dispute over where IQS comes from.
  [doi:10.3389/fcimb.2025.1575421](https://doi.org/10.3389/fcimb.2025.1575421)
- Eberhardt et al. 2021, *J Chem Inf Model* 61:3891-3898. Vina 1.2.0.
  [doi:10.1021/acs.jcim.1c00203](https://doi.org/10.1021/acs.jcim.1c00203)
- Trott & Olson 2010, *J Comput Chem* 31:455-461. AutoDock Vina.
  [doi:10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334)
- SantaLucia 1998, *PNAS* 95:1460-1465. Nearest-neighbour duplex thermodynamics.
  [doi:10.1073/pnas.95.4.1460](https://doi.org/10.1073/pnas.95.4.1460)
- Wallace et al. 1979, *Nucleic Acids Res* 6:3543-3557. Tm rule.
  [doi:10.1093/nar/6.11.3543](https://doi.org/10.1093/nar/6.11.3543)
- Li, Olson & Lu 2019, *Nucleic Acids Res* 47:W26-W34. Web 3DNA 2.0, fiber models.
  [doi:10.1093/nar/gkz394](https://doi.org/10.1093/nar/gkz394)
- O'Boyle et al. 2011, *J Cheminform* 3:33. Open Babel.
  [doi:10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33)
- Hanwell et al. 2012, *J Cheminform* 4:17. Avogadro.
  [doi:10.1186/1758-2946-4-17](https://doi.org/10.1186/1758-2946-4-17)

## Status

Computational work finished, nothing tested experimentally. No sequence from this
shortlist has been synthesised. The next step is wet lab, not more docking, and
before that read LIMITATIONS, because one entry there changes how much the stage-2
ranking is worth.

## 🧑‍💻 Author

Developed by: Vu Bao Chau Nguyen, Ph.D.

Keywords: dsDNA screening, molecular docking, AutoDock Vina, IQS (aeruginaldehyde),
quorum sensing, Pseudomonas aeruginosa, electrochemical biosensor,
nearest-neighbour thermodynamics.
