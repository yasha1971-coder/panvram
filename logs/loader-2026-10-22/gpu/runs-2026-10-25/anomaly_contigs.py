#!/usr/bin/env python3
"""Why is the window set that sample() draws at W = 8192 about 2.4x more expensive than the set it draws at
W = 32768, at the same W and the same number of decoded blocks?

sample() weights every contig by (length - W + 1) and drops contigs shorter than W, so at W = 8192 short contigs
become eligible. In a release-2 HPRC assembly those are the unplaced scaffolds, which align to T2T poorly, so their
blocks should carry more literals and more absolute copies - a longer stream for lane 0 to decode.

This script measures the same 1024 windows of 8192 bases drawn three ways: from every eligible contig (what sample()
does), from contigs of at least 32768 bases only (the set the 32 KiB sweep could see), and from contigs under 32768
only. If the long-contig set is the fast one, the sampling filter explains the v4 figures.
"""
import json
import subprocess
import sys
import time


def sh(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode:
        print(f"$ {cmd}\n{r.stdout[-1200:]}\n{r.stderr[-1200:]}", flush=True)
        raise SystemExit(f"failed: {cmd}")
    return r.stdout


print(sh("nvidia-smi --query-gpu=name,clocks.max.sm --format=csv,noheader").strip(), flush=True)
sh("apt-get -qq install -y libzstd-dev")
sh("git clone -q https://github.com/yasha1971-coder/panvram.git /content/pvrepo || true")
sh("cd /content/pvrepo && pip install -q -e .")

import importlib  # noqa: E402
import os  # noqa: E402

os.chdir("/tmp")
sys.path.insert(0, "/content/pvrepo")
importlib.invalidate_caches()

import torch  # noqa: E402

import panvram  # noqa: E402

assert "/content/pvrepo/panvram/" in panvram.__file__, panvram.__file__
assert panvram.with_cuda
print("panvram", panvram.__file__, "with_cuda", panvram.with_cuda, flush=True)

sys.path.insert(0, "/content/pvrepo/logs/loader-2026-10-22/scripts")
from colab_prepare import prepare_cohort, _sha256  # noqa: E402

REF_URL = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/"
           "GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz")
sh("mkdir -p /content/cohort")
sh(f'curl -fsSL -o /content/t2t.fna.gz "{REF_URL}"')
sh("gunzip -f /content/t2t.fna.gz && mv /content/t2t.fna /content/cohort/reference.fa")
assert _sha256("/content/cohort/reference.fa") == \
    "51ab7d3c2a375ecff53ee2a1be03e58d57cc1601c948a860c03a78d724f26e4c"
prepare_cohort(5, dest="/content/cohort", work="/content/dl")

BATCH, SEED, RUNS, W = 1024, 20261003, 9, 8192
c = panvram.Cohort.open("/content/cohort", device="cuda", dataset="q16k")
torch.cuda.synchronize()
Q = int(c._core.Q)

g = torch.Generator(device="cuda").manual_seed(SEED)
t0 = time.perf_counter()
while time.perf_counter() - t0 < 10:
    c.sample(BATCH, 32768, generator=g)
torch.cuda.synchronize()

# every contig of every assembly: (assembly, global start of the contig, length)
ctg = [(a, rec[3], int(rec[2])) for a in range(len(c.names)) for rec in c._core.contigs(a)]
print(f"contigs {len(ctg)}; >= 32768 bases: {sum(1 for x in ctg if x[2] >= 32768)}; "
      f"< 32768: {sum(1 for x in ctg if x[2] < 32768)}", flush=True)
tot = sum(x[2] for x in ctg)
print(f"bases in contigs >= 32768: {sum(x[2] for x in ctg if x[2] >= 32768)/tot*100:.2f} % of the cohort", flush=True)


def draw(pool, rng):
    """1024 starts drawn uniformly over the bases of `pool`, as sample() does, with start + W <= contig length."""
    w = torch.tensor([max(0, x[2] - W + 1) for x in pool], dtype=torch.float64)
    idx = torch.multinomial(w / w.sum(), BATCH, replacement=True)
    asm, start = [], []
    for i in idx.tolist():
        a, boff, ln = pool[i]
        asm.append(a)
        start.append(boff + int(rng.integers(0, ln - W + 1)))
    return (torch.tensor(asm, dtype=torch.int64, device="cuda"),
            torch.tensor(start, dtype=torch.int64, device="cuda"))


def timeit(asm, start):
    for _ in range(3):
        c.windows(asm, start, W)
    torch.cuda.synchronize()
    ts = []
    for _ in range(RUNS):
        torch.cuda.synchronize()
        a = time.perf_counter()
        c.windows(asm, start, W)
        torch.cuda.synchronize()
        ts.append(time.perf_counter() - a)
    ts.sort()
    return ts[len(ts) // 2]


import numpy as np  # noqa: E402

out = {"W": W, "batch": BATCH, "rows": []}
for name, pool in (("all eligible contigs (what sample() does)", [x for x in ctg if x[2] >= W]),
                   ("contigs >= 32768 only", [x for x in ctg if x[2] >= 32768]),
                   ("contigs < 32768 only", [x for x in ctg if W <= x[2] < 32768])):
    if not pool:
        print(f"{name}: empty pool", flush=True)
        continue
    rng = np.random.default_rng(SEED)
    asm, start = draw(pool, rng)
    nb = int(((start + W - 1) // Q - start // Q + 1).sum())
    med = timeit(asm, start)
    print(f"{name:42s}  {med*1e3:7.3f} ms  blocks {nb:6d}  {med*1e6/nb:6.3f} us per block  "
          f"{BATCH/med:9.0f} windows/s  pool {len(pool)} contigs", flush=True)
    out["rows"].append({"pool": name, "contigs": len(pool), "median_s": med, "blocks": nb,
                        "us_per_block": med * 1e6 / nb, "windows_per_s": BATCH / med})

print("\n" + json.dumps(out, indent=1), flush=True)
print("DONE", flush=True)
