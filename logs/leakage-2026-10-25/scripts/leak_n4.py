#!/usr/bin/env python3
"""leak_n4.py <tiles.tsv> - leakage by locus on the four-assembly set, from the window-to-reference table that
`locate tilemap` writes.

A training window and a test window that answer to the same stretch of the reference are not independent samples: a
model that memorised one has seen the other. This measures how often that happens under two ways of splitting a
pangenome into train and test.

A window is used only if both of its ends have a reference coordinate and the hull of the two is within 10 % of the
window length; a window that crosses a rearrangement or a long insertion is not one locus and is counted separately.
Overlap is the share of a test window's reference interval that the union of the training windows' intervals covers.
"""
import sys

import numpy as np

W = 8192
path = sys.argv[1] if len(sys.argv) > 1 else "tiles8k.tsv"

asm, lo, hi = [], [], []
with open(path) as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")          # '#' appears in contig names: do not let a csv reader treat it as a comment
        asm.append(p[0]); lo.append(int(p[3])); hi.append(int(p[4]))
asm = np.array(asm); lo = np.array(lo, dtype=np.int64); hi = np.array(hi, dtype=np.int64)
ok = hi > lo
good = ok & (np.abs(hi - lo - W) <= W // 10)
names = sorted(set(asm.tolist()))
T2T = 3117292070

print(f"assemblies {names}")
print(f"windows {len(asm):,}; with a reference interval {ok.sum():,} ({100*ok.mean():.3f} %); "
      f"of those, hull within 10 % of W {good.sum():,} ({100*good.sum()/ok.sum():.2f} %)")


def merged(a, b):
    o = np.argsort(a, kind="stable"); a, b = a[o], b[o]
    S, E = [a[0]], [b[0]]
    for x, y in zip(a[1:], b[1:]):
        if x <= E[-1]:
            E[-1] = max(E[-1], y)
        else:
            S.append(x); E.append(y)
    return np.array(S, dtype=np.int64), np.array(E, dtype=np.int64)


def coverer(S, E):
    CL = np.concatenate(([0], np.cumsum(E - S)))

    def F(x):
        i = np.searchsorted(S, x, "right") - 1
        j = np.clip(i, 0, len(S) - 1)
        return np.where(i >= 0, CL[j] + np.clip(np.minimum(x, E[j]) - S[j], 0, None), 0)
    return F


def leak(trm, tem, label):
    tr, te = trm & good, tem & good
    S, E = merged(lo[tr], hi[tr]); F = coverer(S, E)
    cov = (F(hi[te]) - F(lo[te])) / (hi[te] - lo[te])
    print(f"{label}")
    print(f"  train {tr.sum():>9,} windows   test {te.sum():>9,} windows   merged train intervals {len(S):>7,}"
          f"   train covers {int((E-S).sum()):,} reference bases ({100*(E-S).sum()/T2T:.1f} % of T2T)")
    print(f"  test windows overlapping a train window by >= 50 %: {100*(cov>=0.5).mean():7.3f} %")
    print(f"  test windows overlapping a train window by >= 90 %: {100*(cov>=0.9).mean():7.3f} %")
    print(f"  median overlap {100*np.median(cov):.2f} %, zero overlap {100*(cov==0).mean():.3f} %")


rng = np.random.default_rng(20261003)
tr = set(rng.choice(names, 2, replace=False).tolist())
print(f"\n--- random split by sample: train {sorted(tr)}, test {sorted(set(names)-tr)}")
leak(np.isin(asm, list(tr)), ~np.isin(asm, list(tr)), "random by sample")

tr2 = {"y1_HG00438.1", "y1_HG00438.2"}
print(f"\n--- both haplotypes of one individual in train: {sorted(tr2)}")
leak(np.isin(asm, list(tr2)), ~np.isin(asm, list(tr2)), "by individual")

for CH in (10_000_000, 1_000_000):
    chunk = lo // CH
    cid = np.unique(chunk[good])
    r2 = np.random.default_rng(20261003)
    trc = set(r2.choice(cid, len(cid) // 2, replace=False).tolist())
    m = np.isin(chunk, list(trc))
    print(f"\n--- split by locus: T2T cut into {CH//10**6} Mb chunks, {len(trc)} of {len(cid)} to train")
    leak(m, ~m, f"by locus ({CH//10**6} Mb)")
