#!/usr/bin/env python3

import argparse
import sys
from collections import Counter

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
AMBIG = ({"A", "T"}, {"C", "G"})


def load(bim_path):
    rows, pos_count = [], Counter()
    with open(bim_path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 6:
                continue
            vid = f"{f[0]}:{f[3]}"
            rows.append((vid, f))
            pos_count[vid] += 1

    dups = {p for p, n in pos_count.items() if n > 1}
    usable, n_amb, n_nonacgt = {}, 0, 0
    for vid, f in rows:
        a1, a2 = f[4], f[5]
        if vid in dups:
            continue
        if a1 not in COMP or a2 not in COMP or a1 == a2:
            n_nonacgt += 1
            continue
        if {a1, a2} in AMBIG:
            n_amb += 1
            continue
        usable[vid] = (a1, a2)
    return rows, usable, dups, n_amb, n_nonacgt


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref-bim", required=True)
    ap.add_argument("--target-bim", required=True)
    ap.add_argument("--out-ref-bim", required=True)
    ap.add_argument("--out-target-bim", required=True)
    ap.add_argument("--out-shared", required=True)
    ap.add_argument("--min-shared", type=int, default=20000)
    a = ap.parse_args()

    ref_rows, ref_ok, ref_dups, r_amb, r_bad = load(a.ref_bim)
    tgt_rows, tgt_ok, tgt_dups, t_amb, t_bad = load(a.target_bim)

    print(f"  reference: {len(ref_rows):>8} variants, {len(ref_ok):>8} usable "
          f"({len(ref_dups)} duplicate positions, {r_amb} ambiguous, {r_bad} non-ACGT)")
    print(f"  target:    {len(tgt_rows):>8} variants, {len(tgt_ok):>8} usable "
          f"({len(tgt_dups)} duplicate positions, {t_amb} ambiguous, {t_bad} non-ACGT)")

    # reference .bim: only the ID column changes
    with open(a.out_ref_bim, "w") as out:
        for vid, f in ref_rows:
            out.write(f"{f[0]}\t{vid}\t{f[2]}\t{f[3]}\t{f[4]}\t{f[5]}\n")

    # target .bim: ID column, plus allele columns where a strand flip is needed
    shared, n_flip, n_incompat = [], 0, 0
    with open(a.out_target_bim, "w") as out:
        for vid, f in tgt_rows:
            a1, a2 = f[4], f[5]
            r = ref_ok.get(vid)
            if r is not None and vid in tgt_ok:
                ra1, ra2 = r
                if {a1, a2} == {ra1, ra2}:
                    shared.append(vid)
                elif {COMP[a1], COMP[a2]} == {ra1, ra2}:
                    a1, a2 = COMP[a1], COMP[a2]
                    n_flip += 1
                    shared.append(vid)
                else:
                    # same position, genuinely different alleles -- a triallelic
                    # site split across datasets, or a position mismatch
                    n_incompat += 1
            out.write(f"{f[0]}\t{vid}\t{f[2]}\t{f[3]}\t{a1}\t{a2}\n")

    with open(a.out_shared, "w") as fh:
        fh.write("\n".join(shared) + "\n")

    print(f"  shared:    {len(shared):>8} SNPs "
          f"({n_flip} strand-flipped, {n_incompat} allele-incompatible)")

    if len(shared) < a.min_shared:
        sys.exit(
            f"\nERROR: only {len(shared)} shared SNPs, expected at least "
            f"{a.min_shared}.\n"
            "Nearly always a genome-build mismatch -- the two datasets must both\n"
            "be GRCh37/b37 or both GRCh38/b38, because the IDs here are chr:pos.\n"
            "Check a known SNP: rs3131972 is at 1:752721 in b37, 1:817341 in b38.\n"
            "  grep -P '^1\\t' target.bim | head -3\n"
            "If your data is b38, use the GRCh38 reference panel instead, or\n"
            "liftOver your cohort to b37 first.")


if __name__ == "__main__":
    main()
