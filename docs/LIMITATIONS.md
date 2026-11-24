# Limitations and known problems

Notes for whoever picks this up, which will probably be me. Roughly worst first.
Nothing was fixed retroactively: the code and numbers in this repo are what was
actually run, so where something is wrong it is flagged rather than quietly
corrected.

## 1. Stage 2 was built around the wrong molecule

This is the one that matters.

The molecule that was docked is correct. `data/template/template.pdb` is
C10H7NO2S, a thiazole ring with a phenol OH and an aldehyde, which is
2-(2-hydroxyphenyl)-thiazole-4-carbaldehyde, real IQS, MW about 205.2. Every
docking result here used it.

The problem is the justification for the stage-2 interaction score. The original
write-up described the target as "IQS (2-amino acetophenone, MW 135.16)" and built
the scoring around a benzene ring, a carbonyl and an amine. That is
2-aminoacetophenone, a different *P. aeruginosa* metabolite (C8H9NO), the volatile
responsible for the grape-like smell of cultures, regulated by MvfR. The MW 135.16
is 2-AA's, not IQS's. IQS has no free amine, its donor is a phenol OH, and it
carries an aldehyde plus a thiazole N and S, so the groove complementarity
argument does not carry over.

What this damages is only the pre-filter. The 100 sequences that reached docking
were picked using a heuristic designed for the wrong molecule, so there is no
reason to think they are the best 100 of the 10,000 for IQS. The docking and the
shortlist are unaffected, but the shortlist is the best of a pool chosen on a bad
prior. Consistent with that, the stage-2 score correlates +0.307 with the mean
Vina affinity, i.e. a higher pre-filter score goes with a slightly *weaker*
predicted binding.

If this is revisited, rebuild the interaction score on IQS's real groups and
re-screen the 10,000. Patching the weights is not enough.

Refs: Lee et al. 2013 (doi:10.1038/nchembio.1225) for IQS as a signal; Raya,
Montagut & Marco 2025 (doi:10.3389/fcimb.2025.1575421) for IQS as aeruginaldehyde.
Note from the second one that whether IQS is a real signal molecule or a byproduct
of pyochelin metabolism is still disputed. That does not affect the sensing
application but it does affect how the motivation should be worded.

## 2. The Tm part of the score does nothing

`total = 0.40 x score_tm + 0.60 x score_iqs` looks like two components. It isn't.
`score_tm` hits its ceiling of 1.0 for 9,841 of the 10,000 stage-1 sequences, and
after the ΔG filter it is exactly 1.0 for all 500 survivors and all 100 docked
sequences. Zero variance. For every docked candidate
`total = 0.4 + 0.6 x score_iqs` holds exactly. The ranking is a single-score
ranking with a constant added.

Knock-on: the weight sensitivity test in notebook 02 reports an identical top 10
under 40/60, 30/70 and 50/50 weights and concludes the result is robust. It isn't
a robustness result, it just cannot fail, because reweighting a constant changes
nothing. Cause is the wide 30-72 °C window combined with the GC ≥ 55% filter from
stage 1. Either widen the GC range or drop the Tm term.

## 3. The diversity filter is nearly a no-op

Stage 3.2 is meant to stop near-duplicates wasting docking runs. In practice it
rejected 2 sequences out of 1,654, then hit its 500 cap and stopped after looking
at about 502. The remaining 1,150 were never examined. So stage 3.2 is really
"take the top 500 by score". That happens to be fine, because random 12-22mers are
rarely within edit distance 4 of each other, but the filter is not doing the work
the method describes and the threshold of 4 was never tested against anything.

## 4. The ΔG window only has one side

`DG_MIN = -35.0` rejected 346 sequences. `DG_MAX = -5.0` rejected none. Describe it
as "reject duplexes that are too stable", because that is all it ever did. A
GC-rich 12-22mer cannot reach -5 kcal/mol, so that threshold is decorative.

## 5. Docking cannot separate the candidates it is ranking

All 100 candidates fall between -6.934 and -6.456 kcal/mol, a spread of 0.478. The
difference between the two runs for the same sequence has sd 0.053 and reaches
0.228. The two runs share 13 of their top 20, 8 of their top 10, and disagree on
first place (model 87 vs model 59). Model 23 moves 0.202 between runs, which is
42% of the entire spread across all 100 candidates.

So the gap between neighbouring candidates is smaller than the scatter of the
method on the same candidate. The order is noise, the top band is not. This is
before considering that Vina's absolute error against experiment is far larger
than 0.478 kcal/mol.

The earlier PQS docking work on this bench handled this properly: five independent
runs per sequence, mean and SD, selection on consistency across replicates rather
than a single best score. This project did two runs and the saved ranking used one
of them. Redo with five or more runs and a fixed `--seed` if the shortlist ever
needs defending.

## 6. The grid box does not cover the duplex

The box was 40 x 40 x 40 Å on the DNA centroid. The duplexes are 53.2 to 76.8 Å
along the helix axis, mean 71.0, so all 100 receptors are longer than the box.
Between 6.6 and 18.4 Å is excluded at each end, mean 15.5 Å, which is roughly 2 to
5 base pairs per end that the ligand could never visit. Sampling was confined to
about the central 12 bp.

That is defensible as a choice and it lines up with the stage-2 position weighting
that favours central bases, but it was never stated as intentional, and the
original claim that the box "allows sufficient coverage of the DNA groove region"
is too generous. If end effects matter, the box has to scale with sequence length.

## 7. Physics left out

Standard for this kind of screen, but it all limits how far the numbers can be
pushed. Receptor rigid and ligand flexible, so no DNA breathing, no groove
widening, no intercalation. No explicit solvent, no counter-ions, no salt or pH
dependence, while the real sensor works in buffer on a charged surface. The
structures are ideal B-form fiber models, so every duplex has identical backbone
geometry and all sequence-dependent structure is absent, which is a real problem
for a groove-binding question since groove geometry is most of what distinguishes
one sequence from another. And the duplex is docked free in solution, not tethered
to a surface.

## 8. No controls were docked

Everything in the screen is a random sequence. Nothing with known behaviour went
through the same protocol: no established groove binder, no scrambled or poly-AT
negative, no published aptamer. So the 0.478 kcal/mol spread has no scale attached
and there is no way to say whether -6.9 is good.

Cheapest useful thing available: push a few ligand-DNA pairs of known affinity
through this exact pipeline and see where they land. Until then "strongest binder"
only means "least weak of 100 arbitrary sequences".

## 9. G-quadruplex risk the filters miss

`MAX_HOMOPOLYMER = 3` allows GGG, and the scoring rewards G and GG steps, so the
surviving pool is G-rich: 73 of the 100 docked sequences contain at least one GGG
run. On a crude canonical G4 test (four or more G-tracts of length ≥ 2), 17 of the
100 match, though none has four GGG tracts, which would be the strong case. Only 1
of the 20 shortlisted sequences matches.

The modelling is blind to this, since every structure was built as an ideal B-form
duplex regardless of what the sequence would really do. Run anything headed for
synthesis through a G4 predictor first. The risk here is modest but a quadruplex on
the electrode would destroy the sensor response rather than just weaken it.

## 10. Smaller things

The two ranking files are labelled A and B, not 1 and 2, because the re-run
procedure renamed the old file and nothing records which came first. What is
established: all 100 archived pose files match run A exactly and run B never, and
since re-running emptied the pose folder, run A was the later run. See
`data/rankings/PROVENANCE.txt`.

No random seed was set, which is the direct cause of item 5 being unquantifiable
beyond two samples.

I tried adding Vina's intermolecular (INTER) energy term to the results table,
expecting it to discriminate better than the total affinity. It correlates with
affinity at +1.000, because the ligand is the same molecule in all 100 dockings so
the intra-ligand and torsion terms are constant. Dropped the column. Also dropped
a standard deviation column, since an SD over two runs is not worth reporting.

Markdown in the notebooks referred to 20,000 and 5,000 sequences where the code
generated 10,000 and 2,000, to a Tm window of 40-72 °C where the code used 30-72,
and notebook 02 forward-referenced a stage 3 that builds 3D structures and applies
"groove geometry filters", which stage 3 does not do. Fourteen such statements were
corrected. No executable code or stored output was touched, confirmed by diffing
every code cell against the originals.

Because code was left alone, some wrong text survives where it sits inside code
rather than markdown: two plot titles still say "all 20,000 sequences" and "all
5,000 sequences", a progress message still says 20,000, a column header still reads
"All 5k", and notebook 03 still writes `web3dna.rutgers.edu` into
`sequences_for_structure.txt`. The structures actually came from
`web.x3dna-dssr.org`, which every PDB confirms with a `REMARK DSSR` header. Fixed
in RERUN.md, left alone in the code so each cell still matches its output.

The stage-1 sanity check block labels `GCATGCGCATGC` a valid candidate, but it
fails the palindrome filter (3 of 7 windows self-complementary against a 35%
threshold). All six test cases in that block report FAIL, so it demonstrates
nothing. The label is a string in a code cell, so changing it would desynchronise
the cell from its stored output; a markdown note was added underneath instead. The
filters are correct and are checked by `verify_pipeline.py`.

`vina_split` was never used. Scores are read by parsing `REMARK VINA RESULT`
straight out of the multi-model pose files.

The original setup note listed meeko, MDAnalysis, prolif, RDKit and gemmi. Nothing
here imports them, probably aimed at an interaction fingerprint analysis that was
never written. Real requirements are in `requirements.txt`.

Everything is CRLF, since all data files were written on Windows. There is a
`.gitattributes` to stop Git rewriting them, because the byte-level checks in the
verification script run against them as they are.

## Short version

What holds up: the pipeline is fully traceable, every number re-derives from the
raw files, the duplexes are built correctly, the docking settings are sensible, and
the 13 cross-run-stable sequences are a reasonable shortlist for the bench.

What needs redoing before publishing any of it: stage 2 was built around the wrong
molecule, so the 100 that reached docking were not selected for IQS. Rebuild the
interaction score, re-screen, then dock the new top 100 with five runs and a fixed
seed. Maybe two days of compute, and it would make the result defensible.
