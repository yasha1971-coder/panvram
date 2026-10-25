#!/usr/bin/env python3
"""bench_b.py <out prefix> - GenVarLoader: T2T reference + the HPRC v1.1 minigraph-cactus variants of the same three
samples, haplotype sequences reconstructed on the fly. Speed only: the contents are reference-plus-variants, not the
assembly's own bases (see the report).
Dataset: chr1 of CHM13, 3 samples x 2 haplotypes, windows of W bases at the same batch size as A and C.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

import genvarloader as gvl

BATCH, RUNS, WS = 1024, 9, (8192, 32768)
outp = sys.argv[1]
D = Path(__import__("os").environ.get("D1_DIR", "."))
REF, VCF = D / "t2t.chr.fa.gz", D / "chr1.sorted.vcf.gz"
CHR1 = 248387328
rng = np.random.default_rng(20261003)

out = {"tool": f"genvarloader {gvl.__version__}", "device": "cpu", "batch": BATCH, "runs": RUNS,
       "scope": "chr1 of CHM13, 3 HPRC samples x 2 haplotypes, variants hprc-v1.1-mc-chm13 normalised", "rows": []}

for W in WS:
    ds_path = D / f"gvl_W{W}.gvl"
    n_reg = 4096
    starts = rng.integers(0, CHR1 - W - 1, n_reg)
    bed = pl.DataFrame({"chrom": ["chr1"] * n_reg, "chromStart": starts, "chromEnd": starts + W})
    t0 = time.perf_counter()
    if not ds_path.exists():
        gvl.write(path=ds_path, bed=bed, variants=VCF, samples=["HG00438", "HG00621", "HG00673"], max_jitter=0, overwrite=True)
    twrite = time.perf_counter() - t0

    t0 = time.perf_counter()
    ds = gvl.Dataset.open(ds_path, REF).with_len(W).with_seqs("haplotypes")
    topen = time.perf_counter() - t0
    NP = BATCH // 2            # one (region, sample) pair yields 2 haplotype sequences of W bases
    print(f"B\tW {W}\twrite {twrite:.1f} s\topen {topen:.3f} s\tregions {ds.n_regions}\tsamples {ds.n_samples}\tpairs per batch {NP}", flush=True)

    def pick():
        return rng.integers(0, ds.n_regions, NP), rng.integers(0, ds.n_samples, NP)

    for _ in range(3):
        r, sm = pick(); ds[r, sm]
    ts = []
    for _ in range(RUNS):
        r, sm = pick()
        a = time.perf_counter()
        x = ds[r, sm]
        ts.append(time.perf_counter() - a)
        assert x.shape == (NP, 2, W), x.shape
    ts.sort()
    med = ts[len(ts) // 2]
    print(f"B\tW {W}\tbatch {BATCH}\tmedian {med*1e3:.2f} ms\t{BATCH/med:.0f} windows/s\t{BATCH*W/med/1e9:.3f} GB/s\tmin {min(ts)*1e3:.2f}\tmax {max(ts)*1e3:.2f}", flush=True)
    out["rows"].append({"W": W, "write_s": twrite, "open_s": topen, "median_s": med, "min_s": min(ts), "max_s": max(ts),
                        "windows_per_s": BATCH / med, "gb_per_s": BATCH * W / med / 1e9})
json.dump(out, open(f"{outp}.json", "w"), indent=1)
