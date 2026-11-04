#!/usr/bin/env python3

import argparse
import os
import sys
from collections import Counter

SUPERPOPS = {"AFR", "AMR", "EAS", "EUR", "SAS"}


def from_psam(path):
    with open(path) as fh:
        header = fh.readline()
        if not header.startswith("#"):
            return []
        cols = [c.lower() for c in header.lstrip("#").split()]
        if "superpop" not in cols or "population" not in cols:
            return []
        i_iid = cols.index("iid") if "iid" in cols else 0
        i_sp, i_pop = cols.index("superpop"), cols.index("population")
        rows = []
        for line in fh:
            f = line.split()
            if len(f) > max(i_iid, i_sp, i_pop):
                rows.append((f[i_iid], f[i_pop], f[i_sp]))
        return rows


def from_panel(path):
    rows = []
    with open(path) as fh:
        first = fh.readline().split()
        # the official panel file has a header line: sample pop super_pop gender
        if first and first[0].lower() not in ("sample", "iid"):
            fh.seek(0)
        for line in fh:
            f = line.split()
            if len(f) >= 3:
                rows.append((f[0], f[1], f[2]))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refdir", default="1000G")
    ap.add_argument("--out", default="ancestry_out/ref_labels.tsv")
    a = ap.parse_args()

    candidates = [
        (os.path.join(a.refdir, "all_phase3.psam"), from_psam),
        (os.path.join(a.refdir, "phase3_corrected.psam"), from_psam),
        (os.path.join(a.refdir, "integrated_call_samples_v3.20130502.ALL.panel"),
         from_panel),
    ]

    rows, used = [], None
    for path, reader in candidates:
        if os.path.exists(path):
            rows = reader(path)
            if rows:
                used = path
                break

    if not rows:
        sys.exit(
            "ERROR: no usable population labels found in " + a.refdir + "\n\n"
            "Put ONE of these there:\n"
            "  all_phase3.psam   -- rename phase3_corrected.psam, from\n"
            "                       https://www.cog-genomics.org/plink/2.0/resources\n"
            "  integrated_call_samples_v3.20130502.ALL.panel -- from\n"
            "    https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/\n\n"
            "If you do have a .psam, check it has SuperPop and Population\n"
            "columns -- phase3_orig.psam does not always carry them:\n"
            "  head -2 " + a.refdir + "/all_phase3.psam")

    bad = {sp for _, _, sp in rows} - SUPERPOPS
    if bad:
        print(f"note: unexpected superpopulation codes present: {sorted(bad)}")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("IID\tpopulation\tsuperpopulation\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    counts = Counter(sp for _, _, sp in rows)
    print(f"read {len(rows)} samples from {used}")
    for sp, n in sorted(counts.items()):
        print(f"  {sp}  {n:>5}")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
