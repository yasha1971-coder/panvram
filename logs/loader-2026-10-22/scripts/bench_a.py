#!/usr/bin/env python3
"""bench_a.py <cohort dir> <out prefix> - panvram v1 q16k window loader, CPU path.
Batches of BATCH windows of W bases, 3 warm-up batches + RUNS timed batches, median.
Writes the request coordinates of the first timed batch so the other loaders replay the same windows.
"""
import json
import sys
import time

import torch

import panvram

BATCH, RUNS, WS = 1024, 9, (8192, 32768)
cdir, outp = sys.argv[1], sys.argv[2]

t0 = time.perf_counter()
THREADS = int(__import__("os").environ.get("PV_THREADS", "16"))
c = panvram.Cohort.open(cdir, device="cpu", dataset="q16k", threads=THREADS)
topen = time.perf_counter() - t0
print(f"OPEN\t{len(c.names)} assemblies\t{topen:.3f} s\tthreads {THREADS}", flush=True)

out = {"tool": "panvram-v1-q16k", "device": "cpu", "open_s": topen, "assemblies": list(c.names), "batch": BATCH, "runs": RUNS, "threads": THREADS, "rows": []}
for W in WS:
    g = torch.Generator().manual_seed(20261003)
    for _ in range(3):
        c.sample(BATCH, W, generator=g)
    g = torch.Generator().manual_seed(20261003)
    x, coords = c.sample(BATCH, W, generator=g, return_coords=True)
    json.dump([[int(r[0]), int(r[1]), int(r[2])] for r in coords], open(f"{outp}.coords.W{W}.json", "w"))
    ts = []
    for _ in range(RUNS):
        gg = torch.Generator().manual_seed(20261003)
        a = time.perf_counter()
        y = c.sample(BATCH, W, generator=gg)
        ts.append(time.perf_counter() - a)
        assert y.shape == (BATCH, W)
    ts.sort()
    med = ts[len(ts) // 2]
    wps = BATCH / med
    gbs = BATCH * W / med / 1e9
    print(f"A\tW {W}\tbatch {BATCH}\tmedian {med*1e3:.2f} ms\t{wps:.0f} windows/s\t{gbs:.3f} GB/s\tmin {min(ts)*1e3:.2f}\tmax {max(ts)*1e3:.2f}", flush=True)
    out["rows"].append({"W": W, "median_s": med, "min_s": min(ts), "max_s": max(ts), "windows_per_s": wps, "gb_per_s": gbs})
json.dump(out, open(f"{outp}.json", "w"), indent=1)
