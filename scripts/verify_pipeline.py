#!/usr/bin/env python3
"""
Re-checks the dsDNA/IQS pipeline from the raw files in data/ and writes the
result tables.

The 100 structures were downloaded one at a time and ended up named
"model (1).pdb" ... "model (100).pdb" in download order. Nothing in a filename
says which sequence is inside, so the link from a Vina score back to a sequence
depends on that order being right. This script doesn't trust it: it reads the
bases back out of the PDB atom records and checks. It then recomputes every
number in the stage CSVs from the sequence alone, and re-reads every Vina score
out of the pose files.

Standard library only. Run from the repo root:  python scripts/verify_pipeline.py
"""

import csv
import os
import re
import statistics
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEQ = os.path.join(ROOT, "data", "sequences")
DNA = os.path.join(ROOT, "data", "dna")
POSES = os.path.join(ROOT, "data", "output")
RANKS = os.path.join(ROOT, "data", "rankings")
RESULTS = os.path.join(ROOT, "results")

STAGE1 = os.path.join(SEQ, "stage1_candidates.csv")
STAGE2 = os.path.join(SEQ, "stage2_filtered.csv")
STAGE3 = os.path.join(SEQ, "stage3_for_docking.csv")
RUN_A = os.path.join(RANKS, "docking_ranking_runA.txt")
RUN_B = os.path.join(RANKS, "docking_ranking_runB.txt")

N = 100
TOL = 1e-4
log = []


def say(line=""):
    print(line)
    log.append(line)


# --- scoring functions, as used by the notebooks ---------------------------

# Wallace et al. 1979, doi:10.1093/nar/6.11.3543
def wallace_tm(seq):
    return 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))


def tm_score(tm):
    if 30 <= tm <= 72:
        return 1.0
    off = 30 - tm if tm < 30 else tm - 72
    return max(0.0, 1.0 - off / 20)


# IQS interaction heuristic (stage 2)
STEP = {"GG": 3.0, "GA": 2.0, "AG": 2.0, "AA": 1.0}
HBOND = {"G": 2.0, "C": 1.5, "A": 1.0, "T": 0.5}


def pos_weight(i, length):
    return 1.5 if length / 4 <= i < 3 * length / 4 else 1.0


def iqs_raw(seq):
    n = len(seq)
    stack = sum(STEP[seq[i:i + 2]] * (pos_weight(i, n) + pos_weight(i + 1, n)) / 2
                for i in range(n - 1) if seq[i:i + 2] in STEP)
    hbond = sum(HBOND[b] * pos_weight(i, n) for i, b in enumerate(seq))
    return 0.60 * stack + 0.40 * hbond


# SantaLucia 1998, doi:10.1073/pnas.95.4.1460
NN = {"AA": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
      "CA": (-8.5, -22.7), "GT": (-8.4, -22.4), "CT": (-7.8, -21.0),
      "GA": (-8.2, -22.2), "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
      "GG": (-8.0, -19.9), "TT": (-7.9, -22.2), "AC": (-7.8, -21.0),
      "TC": (-8.2, -22.2), "AG": (-7.8, -21.0), "TG": (-8.5, -22.7),
      "CC": (-8.0, -19.9)}


def delta_g(seq):
    dH = dS = 0.0
    for i in range(len(seq) - 1):
        if seq[i:i + 2] in NN:
            h, s = NN[seq[i:i + 2]]
            dH, dS = dH + h, dS + s
    for end in (seq[0], seq[-1]):
        h, s = (2.3, 4.1) if end in "AT" else (0.1, -2.8)
        dH, dS = dH + h, dS + s
    return round(dH - 298.15 * dS / 1000.0, 3)


def rev_comp(seq):
    return "".join({"A": "T", "T": "A", "G": "C", "C": "G"}[b] for b in reversed(seq))


# --- readers ---------------------------------------------------------------

def chain_seq(path, chain):
    """Bases of one strand, from the atom records.

    Avogadro appends added hydrogens after the heavy atoms, so residues
    reappear later in the file. Keying on residue number recovers 5'->3'.
    """
    res = {}
    for line in open(path, errors="replace"):
        if line.startswith(("ATOM", "HETATM")) and line[21] == chain:
            res[int(line[22:26])] = line[17:20].strip()
    return "".join(res[k] for k in sorted(res))


def poses_in(path):
    """[(affinity, inter), ...] for every pose in a Vina output file."""
    out, aff = [], None
    for line in open(path, errors="replace"):
        if line.startswith("REMARK VINA RESULT"):
            aff = float(line.split()[3])
        elif line.startswith("REMARK INTER:") and aff is not None:
            out.append((aff, float(line.split()[-1])))
            aff = None
    return out


def ranking(path):
    scores = {}
    for line in open(path, errors="replace"):
        m = re.match(r"\s*\d+\s+model \((\d+)\)\s+(-?\d+\.?\d*)", line)
        if m:
            scores[int(m.group(1))] = float(m.group(2))
    return scores


def load(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


# --- checks ----------------------------------------------------------------

def main():
    failed = []

    def check(name, ok, detail=""):
        say("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, "  " + detail if detail else ""))
        if not ok:
            failed.append(name)

    for p in (STAGE1, STAGE2, STAGE3, RUN_A, RUN_B):
        if not os.path.exists(p):
            sys.exit("Missing %s. Run this from the repo root." % p)

    s1, s2, s3 = load(STAGE1), load(STAGE2), load(STAGE3)
    say("Pipeline check. Pools: %d -> %d -> %d docked" % (len(s1), len(s2), len(s3)))

    say()
    say("1. Structure files map to the ranked sequences")
    csv_seqs = [r["sequence"] for r in s3]
    seqs, wrong, bad_bp, missing = {}, [], [], []
    for n in range(1, N + 1):
        path = os.path.join(DNA, "model (%d).pdb" % n)
        if not os.path.exists(path):
            missing.append(n)
            continue
        a, b = chain_seq(path, "A"), chain_seq(path, "B")
        seqs[n] = a
        if a != csv_seqs[n - 1]:
            wrong.append(n)
        if b != rev_comp(a):
            bad_bp.append(n)
    check("all %d PDB files present" % N, not missing, str(missing) if missing else "")
    check("chain A of model (N) is stage-3 rank N", not wrong and len(seqs) == N,
          str(wrong) if wrong else "%d/%d" % (len(seqs), N))
    check("chain B is the reverse complement of chain A", not bad_bp,
          str(bad_bp) if bad_bp else "all %d duplexes correct" % len(seqs))
    check("the 100 docked sequences are all different", len(set(seqs.values())) == len(seqs))

    say()
    say("2. Stored values recompute from the sequence")
    pool = [iqs_raw(r["sequence"]) for r in s1]
    lo, hi = min(pool), max(pool)
    say("     iqs_raw range over the stage-1 pool: %.3f to %.3f" % (lo, hi))

    def recompute(rows, label, cols):
        worst, bad = 0.0, 0
        for r in rows:
            q = r["sequence"]
            got = {"gc_content": (q.count("G") + q.count("C")) / len(q),
                   "tm_celsius": wallace_tm(q), "score_tm": tm_score(wallace_tm(q)),
                   "iqs_raw": iqs_raw(q), "score_iqs": (iqs_raw(q) - lo) / (hi - lo),
                   "delta_g": delta_g(q)}
            got["total_score"] = 0.40 * got["score_tm"] + 0.60 * got["score_iqs"]
            for c in cols:
                d = abs(got[c] - float(r[c]))
                worst = max(worst, d)
                bad += d > TOL
        check("%s: %d columns x %d rows recompute" % (label, len(cols), len(rows)),
              bad == 0, "max deviation %.1e%s" % (worst, "" if not bad else ", %d bad" % bad))

    recompute(s2, "stage 2", ["gc_content", "tm_celsius", "score_tm", "iqs_raw",
                              "score_iqs", "total_score"])
    recompute(s3, "stage 3", ["gc_content", "tm_celsius", "score_tm", "iqs_raw",
                              "score_iqs", "total_score", "delta_g"])

    say()
    say("3. Vina scores in the pose files match the ranking tables")
    a_rank, b_rank = ranking(RUN_A), ranking(RUN_B)
    pose, odd, counts = {}, [], []
    for n in range(1, N + 1):
        path = os.path.join(POSES, "model (%d)_out.pdbqt" % n)
        if not os.path.exists(path):
            odd.append("model (%d) missing" % n)
            continue
        p = poses_in(path)
        if not p:
            odd.append("model (%d) unparseable" % n)
            continue
        pose[n] = p
        counts.append(len(p))
        best = min(s for s, _ in p)
        if min(abs(best - a_rank.get(n, 9e9)), abs(best - b_rank.get(n, 9e9))) > 1e-3:
            odd.append("model (%d) best pose %.3f in neither table" % (n, best))
    check("pose file readable for all %d models" % N, len(pose) == N, "%d read" % len(pose))
    check("best pose matches one of the two runs", not odd,
          "; ".join(odd[:3]) if odd else "%d-%d poses per file" % (min(counts), max(counts)))
    check("both ranking tables cover all %d models" % N,
          len(a_rank) == N and len(b_rank) == N)
    check("archived poses all belong to run A",
          all(abs(min(s for s, _ in pose[n]) - a_rank[n]) < 1e-9 for n in pose))

    say()
    say("4. Agreement between the two docking runs")
    diff = [b_rank[n] - a_rank[n] for n in a_rank]
    vals = list(a_rank.values())
    ord_a = sorted(a_rank, key=lambda n: a_rank[n])
    ord_b = sorted(b_rank, key=lambda n: b_rank[n])
    say("     affinity spread over 100 candidates : %.3f to %.3f kcal/mol (range %.3f)"
        % (min(vals), max(vals), max(vals) - min(vals)))
    say("     between-run difference              : sd %.3f, largest %.3f kcal/mol"
        % (statistics.stdev(diff), max(abs(d) for d in diff)))
    say("     top-20 shared / top-10 shared       : %d of 20 / %d of 10"
        % (len(set(ord_a[:20]) & set(ord_b[:20])), len(set(ord_a[:10]) & set(ord_b[:10]))))
    say("     best candidate                      : run A model (%d), run B model (%d)"
        % (ord_a[0], ord_b[0]))
    say()
    say("     The scatter between runs is about the size of the gap between")
    say("     neighbouring candidates, so the order is not reproducible but the")
    say("     top band is. Use the output as a shortlist. See docs/LIMITATIONS.md.")

    # --- tables ---
    say()
    say("Writing tables")
    os.makedirs(RESULTS, exist_ok=True)
    top20_a, top20_b = set(ord_a[:20]), set(ord_b[:20])

    rows = []
    for n in range(1, N + 1):
        r = s3[n - 1]
        a, b = a_rank[n], b_rank[n]
        rows.append({
            "model_file": "model (%d)" % n,
            "prefilter_rank": n,
            "sequence": seqs.get(n, r["sequence"]),
            "reverse_complement": r["reverse_comp"],
            "length_bp": r["length"],
            "gc_fraction": "%.4f" % float(r["gc_content"]),
            "tm_wallace_c": r["tm_celsius"],
            "delta_g_kcal_mol": "%.3f" % float(r["delta_g"]),
            "gg_steps": r["gg_steps"],
            "purine_fraction": r["purine_fraction"],
            "prefilter_score_iqs": "%.4f" % float(r["score_iqs"]),
            "prefilter_total_score": "%.4f" % float(r["total_score"]),
            "vina_runA": "%.3f" % a,
            "vina_runB": "%.3f" % b,
            "vina_mean": "%.3f" % ((a + b) / 2),
            "n_poses": len(pose.get(n, [])),
            "in_top20_both_runs": "yes" if (n in top20_a and n in top20_b) else "no",
        })

    order = sorted(rows, key=lambda r: float(r["vina_mean"]))
    for i, r in enumerate(order, 1):
        r["docking_rank"] = i
    cols = ["docking_rank"] + [c for c in rows[0] if c != "docking_rank"]

    def write(path, data):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(data)

    write(os.path.join(RESULTS, "candidates_master.csv"), order)
    write(os.path.join(RESULTS, "top20_candidates.csv"), order[:20])
    stable = sum(1 for r in order[:20] if r["in_top20_both_runs"] == "yes")
    say("  results/candidates_master.csv  %d rows" % len(order))
    say("  results/top20_candidates.csv   20 rows, %d of them in the top 20 of both runs" % stable)

    say()
    if failed:
        say("%d check(s) FAILED: %s" % (len(failed), ", ".join(failed)))
    else:
        say("All checks passed.")

    with open(os.path.join(RESULTS, "verification_report.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
