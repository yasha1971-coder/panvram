#!/usr/bin/env python3
"""bench_c.py <bgzip dir> <names file> <coords prefix> <out prefix> - pyfaidx over bgzip-compressed assembly FASTA.
Replays exactly the windows of bench_a.py (same (assembly, contig index, start)), same batch, same warm-up and runs.
Contig index -> name from the .fai of each file (same order as the archive's records).
"""
import json
import sys
import time

import numpy as np
import torch
from pyfaidx import Fasta

BATCH, RUNS, WS = 1024, 9, (8192, 32768)
bdir, namesf, cpre, outp = sys.argv[1:5]
names = [l.strip() for l in open(namesf) if l.strip()]

t0 = time.perf_counter()
fa, ctg = [], []
for n in names:
    p = f"{bdir}/{n}.fa.gz"
    fa.append(Fasta(p, as_raw=True, sequence_always_upper=False, rebuild=False))
    ctg.append([l.split("\t")[0] for l in open(p + ".fai")])
topen = time.perf_counter() - t0
print(f"OPEN\t{len(fa)} assemblies\t{topen:.3f} s", flush=True)

out = {"tool": f"pyfaidx {__import__('pyfaidx').__version__} on bgzip FASTA", "device": "cpu", "open_s": topen,
       "assemblies": names, "batch": BATCH, "runs": RUNS, "rows": []}
for W in WS:
    co = json.load(open(f"{cpre}.coords.W{W}.json"))
    assert len(co) == BATCH

    def batch():
        x = np.empty((BATCH, W), dtype=np.uint8)
        for i, (a, c, s) in enumerate(co):
            x[i] = np.frombuffer(str(fa[a][ctg[a][c]][s:s + W]).encode(), dtype=np.uint8)
        return torch.from_numpy(x)

    for _ in range(3):
        batch()
    ts = []
    for _ in range(RUNS):
        t = time.perf_counter()
        y = batch()
        ts.append(time.perf_counter() - t)
        assert y.shape == (BATCH, W)
    ts.sort()
    med = ts[len(ts) // 2]
    print(f"C\tW {W}\tbatch {BATCH}\tmedian {med*1e3:.2f} ms\t{BATCH/med:.0f} windows/s\t{BATCH*W/med/1e9:.3f} GB/s\tmin {min(ts)*1e3:.2f}\tmax {max(ts)*1e3:.2f}", flush=True)
    out["rows"].append({"W": W, "median_s": med, "min_s": min(ts), "max_s": max(ts), "windows_per_s": BATCH / med, "gb_per_s": BATCH * W / med / 1e9})
json.dump(out, open(f"{outp}.json", "w"), indent=1)
