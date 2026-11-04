#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SUPERPOP_NAMES = {
    "AFR": "African", "AMR": "Admixed American", "EAS": "East Asian",
    "EUR": "European", "SAS": "South Asian",
}
COLOURS = {"AFR": "#E69F00", "AMR": "#009E73", "EAS": "#0072B2",
           "EUR": "#CC79A7", "SAS": "#56B4E9"}

MISS_MAX = 0.02      # call rate below 98% is a bad sample
HET_Z_MAX = 3.0      # heterozygosity this far from its ancestry group's mean
INNER_Q = 0.95       # 'inside' = within this quantile of the group's own samples
OUTER_Q = 0.999      # 'outside' = beyond this quantile times MARGIN
MARGIN = 1.5
UNCERTAIN = "Admixed/Uncertain"


# readers 
def read_eigenvec(path: Path) -> pd.DataFrame:
    """PLINK 1.9 .eigenvec -> DataFrame indexed by IID, columns PC1..PCn."""
    with open(path) as fh:
        first = fh.readline().split()
    if not first:
        sys.exit(f"{path} is empty")
    if first[0].lstrip("#").upper() in ("FID", "IID"):
        df = pd.read_csv(path, sep=r"\s+")
        df.columns = [c.lstrip("#") for c in df.columns]
    else:
        names = ["FID", "IID"] + [f"PC{i + 1}" for i in range(len(first) - 2)]
        df = pd.read_csv(path, sep=r"\s+", header=None, names=names)
    if "IID" not in df.columns:
        sys.exit(f"{path}: no IID column (found {list(df.columns)})")
    pc = sorted([c for c in df.columns if c.upper().startswith("PC")],
                key=lambda c: int(c[2:]))
    if not pc:
        sys.exit(f"{path}: no PC columns")
    out = df[pc].copy()
    out.columns = [f"PC{i + 1}" for i in range(len(pc))]
    out.index = df["IID"].astype(str)
    out.index.name = "IID"
    return out


def read_qc(prefix: str) -> pd.DataFrame:
    out = pd.DataFrame()
    imiss, het = Path(f"{prefix}.imiss"), Path(f"{prefix}.het")
    if imiss.exists():
        d = pd.read_csv(imiss, sep=r"\s+")
        out["missing_rate"] = pd.Series(d["F_MISS"].values, index=d["IID"].astype(str))
    if het.exists():
        d = pd.read_csv(het, sep=r"\s+")
        out["het_F"] = pd.Series(d["F"].values, index=d["IID"].astype(str))
    out.index.name = "IID"
    return out


# distance geometry 
def _mahal(A, mu, inv):
    d = A - mu
    return np.sqrt(np.einsum("ij,jk,ik->i", d, inv, d))


def group_stats(X, y, classes, npc):
    st = {}
    for cls in classes:
        G = X[y == cls]
        mu = G.mean(axis=0)
        inv = np.linalg.pinv(np.cov(G, rowvar=False) + np.eye(npc) * 1e-10)
        d = _mahal(G, mu, inv)
        st[cls] = (mu, inv, np.quantile(d, INNER_Q), np.quantile(d, OUTER_Q) * MARGIN)
    return st


def distances(Xt, classes, st):
    D = np.empty((len(Xt), len(classes)))
    for k, cls in enumerate(classes):
        mu, inv, _, _ = st[cls]
        D[:, k] = _mahal(Xt, mu, inv)
    inner = np.array([st[c][2] for c in classes])
    outer = np.array([st[c][3] for c in classes])
    return D, inner, outer


def assign(proba, classes, D, inner, outer, min_prob):
    classes = np.asarray(classes)
    best, top_p = proba.argmax(axis=1), proba.max(axis=1)
    is_inside = D <= inner[None, :]
    is_within = D <= outer[None, :]
    n_inside = is_inside.sum(axis=1)
    outside_all = ~is_within.any(axis=1)
    rel = D / outer[None, :]                      # <1 means within that group
    nearest = rel.argmin(axis=1)
    rows = np.arange(len(D))
    _ = is_inside[rows, best]      # kept for readability of the branches below

    fit = np.where(outside_all, "outside", np.where(n_inside > 0, "inside", "edge"))
    call = np.empty(len(D), dtype=object)
    conflict = np.array([""] * len(D), dtype=object)

    for i in range(len(D)):
        if outside_all[i]:
            call[i] = "Outlier"
        elif n_inside[i] == 1:
            j = int(np.flatnonzero(is_inside[i])[0])
            call[i] = classes[j]
            if j != best[i]:
                conflict[i] = f"forest said {classes[best[i]]} (p={top_p[i]:.2f})"
        elif n_inside[i] > 1:
            call[i] = classes[best[i]] if top_p[i] >= min_prob else UNCERTAIN
            conflict[i] = "inside " + "+".join(classes[np.flatnonzero(is_inside[i])])
        elif top_p[i] >= min_prob:
            call[i] = classes[best[i]]
        else:
            call[i] = UNCERTAIN
    return dict(call=call, fit=fit, best=best, top_p=top_p, nearest=nearest,
                rel=rel, D=D, conflict=conflict, outside_all=outside_all)


def make_forest(seed):
    from sklearn.ensemble import RandomForestClassifier
    return RandomForestClassifier(n_estimators=1000, min_samples_leaf=2,
                                  random_state=seed, n_jobs=-1)


def cv_full_rule(X, y, npc, min_prob, seed):
    from sklearn.model_selection import StratifiedKFold
    out = np.empty(len(y), dtype=object)
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed).split(X, y):
        clf = make_forest(seed).fit(X[tr], y[tr])
        st = group_stats(X[tr], y[tr], clf.classes_, npc)
        D, inner, outer = distances(X[te], clf.classes_, st)
        out[te] = assign(clf.predict_proba(X[te]), clf.classes_, D, inner, outer,
                         min_prob)["call"]
    return out


def classify(ref, tgt_pcs, pcs, npc, level, include_amr, min_prob, seed, verbose=True):
    train = ref if include_amr else ref[ref["superpopulation"] != "AMR"]
    y, X = train[level].values, train[pcs].values
    Xt = tgt_pcs[pcs].values

    clf = make_forest(seed).fit(X, y)
    st = group_stats(X, y, clf.classes_, npc)
    D, inner, outer = distances(Xt, clf.classes_, st)
    res = assign(clf.predict_proba(Xt), clf.classes_, D, inner, outer, min_prob)
    res["classes"] = clf.classes_
    res["proba"] = clf.predict_proba(Xt)

    if verbose:
        cv_calls = cv_full_rule(X, y, npc, min_prob, seed)
        agree = (cv_calls == y).mean()
        print(f"5-fold CV of the complete rule, on reference samples: {agree:.3f}")
        bad = pd.Series(cv_calls)[cv_calls != y].value_counts()
        if len(bad):
            print("  reference samples NOT given their own label:")
            for lab, n in bad.items():
                print(f"    {lab:<20} {n:>4}  ({n / len(y):.1%})")
        print("  This is 1000 Genomes evaluating itself: sequenced data, the same\n"
              "  marker set, the same cohort the thresholds came from. It is a\n"
              "  self-consistency check, NOT an accuracy figure for your cohort,\n"
              "  which is a different platform with different SNPs and missingness.")
        res["cv_calls"] = cv_calls
    return res


# main 
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="ancestry_out")
    ap.add_argument("--eigenvec", required=True)
    ap.add_argument("--qc", default=None)
    ap.add_argument("--npc", type=int, default=8)
    ap.add_argument("--level", choices=["superpopulation", "population"],
                    default="superpopulation")
    ap.add_argument("--min-prob", type=float, default=0.80)
    ap.add_argument("--include-amr", action="store_true",
                    help="train on AMR as well as the four other groups")
    ap.add_argument("--compare-amr", action="store_true",
                    help="classify both ways and report every sample whose call "
                         "changes, instead of assuming one setting is better")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = Path(args.outdir)
    labels = pd.read_csv(out / "ref_labels.tsv", sep="\t", dtype=str).set_index("IID")

    allpcs = read_eigenvec(Path(args.eigenvec))
    is_ref = allpcs.index.isin(labels.index)
    ref_pcs, tgt_pcs = allpcs[is_ref], allpcs[~is_ref]
    print(f"{len(allpcs)} samples in {args.eigenvec}")
    if tgt_pcs.empty:
        sys.exit("Every sample matched a reference label -- your IDs collide "
                 "with the 1000 Genomes IDs.")
    if ref_pcs.empty:
        sys.exit("No sample matched ref_labels.tsv. Compare the ID formats:\n"
                 f"  head -2 {args.eigenvec}\n  head -2 {out / 'ref_labels.tsv'}")

    npc = min(args.npc, ref_pcs.shape[1])
    pcs = [f"PC{i + 1}" for i in range(npc)]
    ref = ref_pcs.join(labels, how="inner")
    print(f"  reference: {len(ref)}    cohort: {len(tgt_pcs)}    PCs used: {npc}\n")

    # AMR
    if args.compare_amr:
        print("=== classifying WITHOUT AMR ===")
        a = classify(ref, tgt_pcs, pcs, npc, args.level, False, args.min_prob, args.seed)
        print("\n=== classifying WITH AMR ===")
        b = classify(ref, tgt_pcs, pcs, npc, args.level, True, args.min_prob, args.seed)
        changed = a["call"] != b["call"]
        print(f"\n=== {changed.sum()} of {len(tgt_pcs)} samples change call ===")
        if changed.any():
            cmp = pd.DataFrame({
                "IID": tgt_pcs.index[changed],
                "without_AMR": a["call"][changed], "with_AMR": b["call"][changed],
                "fit_without": a["fit"][changed], "fit_with": b["fit"][changed],
            })
            print(cmp.to_string(index=False))
            cmp.to_csv(out / "amr_comparison.tsv", sep="\t", index=False)
            print(f"\nwritten: {out / 'amr_comparison.tsv'}")
        print("\nDecide from these samples, not from a default. AMR is admixed and\n"
              "overlaps EUR, so including it can pull genuine EUR samples into AMR\n"
              "or into the uncertain band; excluding it does NOT guarantee that\n"
              "Latin American samples become Admixed/Uncertain -- they may be\n"
              "assigned elsewhere or flagged Outlier. Only the table above says\n"
              "what actually happens to your data.\n")
        res = b if args.include_amr else a
    else:
        if not args.include_amr:
            n_amr = (ref["superpopulation"] == "AMR").sum()
            if n_amr:
                print(f"AMR excluded from training ({n_amr} samples). This is a "
                      f"default, not a finding --\n  run --compare-amr to see what "
                      f"it does to your samples.\n")
        res = classify(ref, tgt_pcs, pcs, npc, args.level, args.include_amr,
                       args.min_prob, args.seed)

    classes = res["classes"]

    # assemble 
    tab = pd.DataFrame({
        "IID": tgt_pcs.index,
        f"best_{args.level}": classes[res["best"]],
        "probability": res["top_p"].round(4),
        "nearest_by_distance": classes[res["nearest"]],
        "call": res["call"],
        "group_fit": res["fit"],
        "note": res["conflict"],
    })
    for k, cls in enumerate(classes):
        tab[f"d_{cls}"] = res["D"][:, k].round(2)
    for k, cls in enumerate(classes):
        tab[f"rel_{cls}"] = res["rel"][:, k].round(3)   # <1 = inside that group's cut-off
    for k, cls in enumerate(classes):
        tab[f"p_{cls}"] = res["proba"][:, k].round(4)

    # sample QC 
    qc_notes = []
    if args.qc:
        qc = read_qc(args.qc)
        tab = tab.merge(qc, left_on="IID", right_index=True, how="left")
        if "missing_rate" in tab:
            bad = tab["missing_rate"] > MISS_MAX
            tab["low_callrate"] = np.where(bad.fillna(False), "YES", "")
            if bad.sum():
                qc_notes.append(f"{int(bad.sum())} sample(s) with missing rate > "
                                f"{MISS_MAX:.0%} -- unambiguously poor data")
        if "het_F" in tab:
            # F is computed against cohort-wide allele frequencies, so a sample of
            # minority ancestry gets an inflated |F| purely from being different.
            # Standardise inside each ancestry call instead.
            z = pd.Series(np.nan, index=tab.index)
            for _, rows in tab.groupby("call"):
                f = rows["het_F"]
                if len(rows) >= 20 and f.std(ddof=0) > 0:
                    z.loc[rows.index] = (f - f.mean()) / f.std(ddof=0)
            miss = z.isna()
            if miss.any() and tab["het_F"].std(ddof=0) > 0:
                f = tab["het_F"]
                z.loc[miss] = ((f - f.mean()) / f.std(ddof=0)).loc[miss]
            tab["het_z"] = z.round(2)
            badh = z.abs() > HET_Z_MAX
            tab["het_outlier"] = np.where(badh.fillna(False), "YES", "")
            if badh.sum():
                qc_notes.append(f"{int(badh.sum())} sample(s) more than {HET_Z_MAX:g} "
                                f"sd from their group's heterozygosity")

    tab = pd.concat([tab, tgt_pcs[pcs].round(5).reset_index(drop=True)], axis=1)
    res_path = out / "target_ancestry.tsv"
    tab.to_csv(res_path, sep="\t", index=False)

    # report 
    print("group fit: " + ", ".join(f"{k} {v}" for k, v in
                                    pd.Series(res["fit"]).value_counts().items()))
    print("\ncalls:")
    for lab, n in tab["call"].value_counts().items():
        print(f"  {lab:<18} {n:>5}" +
              (f"   ({SUPERPOP_NAMES[lab]})" if lab in SUPERPOP_NAMES else ""))
    n_conf = int((res["conflict"] != "").sum())
    if n_conf:
        print(f"\n{n_conf} sample(s) where the forest and the distances disagreed "
              f"-- see the 'note' column.")
    n_out = int(res["outside_all"].sum())
    if n_out:
        print(f"\n{n_out} sample(s) beyond the cut-off of EVERY reference group. "
              f"That covers\nboth genuine admixture and technical failure -- read "
              f"the QC columns before dropping any.")
    for note in qc_notes:
        print(f"  {note}")
    print(f"\nwritten: {res_path}")

    # plot 
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        pairs = [(1, 2), (3, 4)] if npc >= 4 else [(1, 2)]
        fig, axes = plt.subplots(1, len(pairs), figsize=(7 * len(pairs), 6.2))
        axes = np.atleast_1d(axes)
        outl = res["outside_all"]
        for ax, (i, j) in zip(axes, pairs):
            for cls in sorted(ref["superpopulation"].unique()):
                m = ref["superpopulation"] == cls
                ax.scatter(ref.loc[m, f"PC{i}"], ref.loc[m, f"PC{j}"], s=9,
                           alpha=0.45, linewidths=0, c=COLOURS.get(cls, "#999999"),
                           label=f"{cls} – {SUPERPOP_NAMES.get(cls, cls)}")
            ax.scatter(tgt_pcs[f"PC{i}"], tgt_pcs[f"PC{j}"], s=28, marker="x",
                       c="black", linewidths=1.0, label="your cohort")
            if n_out:
                ax.scatter(tgt_pcs[f"PC{i}"][outl], tgt_pcs[f"PC{j}"][outl],
                           s=130, facecolors="none", edgecolors="red",
                           linewidths=1.3, label="outside every group")
            ax.set_xlabel(f"PC{i}"); ax.set_ylabel(f"PC{j}")
            ax.spines[["top", "right"]].set_visible(False)
        axes[0].legend(frameon=False, fontsize=9, loc="best")
        fig.suptitle("Joint PCA of the cohort and 1000 Genomes phase 3", fontsize=13)
        fig.tight_layout()
        png = out / "ancestry_pca.png"
        fig.savefig(png, dpi=160)
        print(f"written: {png}")
    except ImportError:
        print("(matplotlib not installed -- plot skipped)")


if __name__ == "__main__":
    main()
