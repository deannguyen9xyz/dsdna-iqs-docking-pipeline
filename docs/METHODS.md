# Methods

Written as methods prose so it can be lifted into a manuscript with light
editing. It describes what was actually run, with the numbers from this repo
rather than the planned ones. Where a choice has a known weakness there is a
pointer to [LIMITATIONS.md](LIMITATIONS.md) instead of a defence of it.

Target molecule throughout: IQS, 2-(2-hydroxyphenyl)-thiazole-4-carbaldehyde
(aeruginaldehyde), C10H7NO2S, MW about 205.2.

## Sequence generation

Single-strand DNA sequences were generated at random from a uniform base
distribution, with length drawn uniformly from 12 to 22 nucleotides
(`random.seed(41)`, so the pool is reproducible). Only the sense strand was
generated; the complement follows from Watson-Crick pairing.

Each sequence had to pass five filters. Length 12 to 22 nt, since below 12 there
is too little groove to hold a small molecule and above 22 the synthesis cost and
unwanted flexibility rise. GC content 55 to 80%, since below 55% the duplex is too
weak at room temperature and above 80% it becomes rigid and the groove closes up.
No more than 3 identical bases in a row, because G runs of four or more favour
quadruplex over duplex and long C runs self-associate. No single base above 60% of
the sequence, since a sequence can sit inside the GC window and still behave like
a poly-G strand. And no more than 35% of sliding 6-nt windows self-complementary,
because self-complementary stretches fold back into hairpins instead of pairing
with the complement, which gives unreliable sensor signals.

35,570 random draws gave 10,000 passing sequences, a 28.1% pass rate. Rejections
were dominated by GC content (23,349 of 25,570), with homopolymer runs accounting
for 2,194 and base dominance and palindrome fraction for 27 between them.

## Scoring and pre-filtering

Each sequence received a composite score, `0.40 x score_Tm + 0.60 x score_IQS`,
and the top 2,000 were kept.

Melting temperature came from the Wallace rule, Tm = 2(A+T) + 4(G+C) °C (Wallace
et al. 1979), which is valid for short oligonucleotides and so suits the 12 to
22 nt range. Sequences with Tm between 30 and 72 °C scored 1.0, and outside that
window the score fell linearly to zero over 20 °C. The lower bound was set at 30
rather than 40 °C deliberately, to avoid systematically discarding the shorter
sequences.

The interaction term was a heuristic meant to capture groove complementarity. It
combined pi-stacking from dinucleotide steps (GG = 3.0, GA or AG = 2.0, AA = 1.0,
everything else 0), hydrogen bond potential per base (G = 2.0, C = 1.5, A = 1.0,
T = 0.5), and a position weight of 1.5 for bases in the central 50% of the
sequence against 1.0 for the outer 25% at each end, on the basis that the groove
is better defined away from the duplex ends. These were combined as
`0.60 x pi-stacking + 0.40 x H-bond` and then min-max normalised across the full
10,000-sequence pool, whose raw range was 6.600 to 40.100.

Two caveats have to be carried into any write-up. The functional group reasoning
behind these weights was derived from 2-aminoacetophenone rather than IQS, a
different *P. aeruginosa* metabolite that was conflated with IQS during planning;
the docking used correct IQS throughout but the pre-filter that chose which 100
sequences reached docking did not (LIMITATIONS item 1). And `score_Tm` carries no
variance past stage 3, so this composite is effectively a single-score ranking
(item 2).

## Thermodynamic and diversity filtering

ΔG of duplex formation was computed with the nearest-neighbour model (SantaLucia
1998) as ΔG = ΔH - T·ΔS at T = 298.15 K, summing published dinucleotide increments
and applying terminal initiation corrections for A/T against G/C ends. This
replaces the Wallace estimate with a sequence-context-aware value. Sequences
outside -35.0 to -5.0 kcal/mol were rejected: 346 of 2,000 fell below the lower
bound and none above the upper one, so in practice the filter removed over-stable
duplexes only.

A greedy Levenshtein filter then kept a sequence only if it differed from every
already-accepted sequence by at least 4 edits, working down from the highest
scoring. Edit distance rather than Hamming distance, because the sequences differ
in length. Two near-duplicates were removed before the 500-sequence cap was
reached (item 3). The top 100 of the 500 by stage-2 score were taken forward.

## Structure preparation

Ideal B-form duplexes for each of the 100 sequences were built with the fiber
model builder of the 3DNA/DSSR web platform (https://web.x3dna-dssr.org, B-DNA
model; Li, Olson & Lu 2019). Hydrogens were added in Avogadro. The files retain
the `REMARK DSSR` provenance header.

Because structures were generated and downloaded one at a time, they are named
`model (1).pdb` to `model (100).pdb` in download order, which corresponds to
descending stage-2 score. This correspondence was verified rather than assumed:
chain A of `model (N).pdb` was read back from the atom records and matched to
stage-3 rank N for all 100 files, and chain B was confirmed as the exact reverse
complement of chain A in every case (`scripts/verify_pipeline.py`).

Receptors and ligand were converted to PDBQT with Open Babel using Gasteiger
partial charges, which merges non-polar hydrogens and keeps polar ones: 1,373 to
992 atoms for a 21 bp duplex, and 21 to 15 for the ligand. Because Open Babel
writes PDBQT in ligand form, the torsion tree records (`ROOT`, `ENDROOT`,
`BRANCH`, `ENDBRANCH`, `TORSDOF`) were stripped from the DNA files so Vina would
read them as rigid receptors. The ligand kept 3 active torsions. Receptor atom
types after conversion were A, C, HD, N, NA, OA and P, all valid AutoDock types.

## Docking

Docking used AutoDock Vina 1.2.0 (Eberhardt et al. 2021; Trott & Olson 2010) with
the receptor rigid and the ligand flexible. For each duplex the search box was a
40 x 40 x 40 Å cube centred on the geometric centre of the DNA, with exhaustiveness
32, up to 30 binding modes per receptor and 12 CPU threads. 28 to 30 poses were
returned per receptor.

Vina approximates the binding free energy as a weighted sum of steric, hydrogen
bonding and hydrophobic terms plus a penalty for ligand conformational entropy,
and searches with stochastic global optimisation refined by local BFGS
minimisation. No random seed was fixed, so repeat runs differ.

Note that the 40 Å box is shorter than every receptor, which ran 53.2 to 76.8 Å
along the helix axis, so sampling was confined to roughly the central 12 base pairs
and the terminal 2 to 5 base pairs at each end were not accessible to the ligand
(item 6).

## Analysis and selection

Two independent docking runs were carried out over all 100 receptors. For each run
the best-scoring pose per receptor was taken as that sequence's affinity.
Candidates were ranked by the mean of the two runs and the best 20 taken as the
shortlist.

Because Vina's search is stochastic, agreement between runs was used as the
stability criterion rather than the score alone. Of the 20 shortlisted sequences,
13 occupy a top-20 position in both runs independently, and those 13 are reported
as the robust subset. Across all 100 candidates the affinities span only
0.478 kcal/mol while the between-run standard deviation is 0.053 kcal/mol, the
largest single-sequence shift is 0.228 kcal/mol, and the two runs do not agree on
the highest-ranked sequence. The results are therefore reported as a shortlist
rather than a ranking, and absolute affinities are not interpreted as predicted
binding constants.

Per-candidate results for both runs with the cross-run stability flag are in
`results/candidates_master.csv`.

## Limitations

Receptors were treated as rigid ideal B-form duplexes, so sequence-dependent groove
geometry, arguably the main determinant of groove-binding selectivity, is absent
from the model. Solvent, counter-ions and the electrode surface were not
represented. No reference ligand-DNA complex of known affinity was docked through
the same protocol, so the affinity range has no external calibration. Docking of
this kind supports relative shortlisting only and needs experimental validation.
Full list in [LIMITATIONS.md](LIMITATIONS.md).
