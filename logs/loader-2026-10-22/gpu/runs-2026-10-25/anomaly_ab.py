#!/usr/bin/env python3
"""Does the 7.033 ms row reproduce?

In the first session, 1024 windows of 8192 bases drawn by `Cohort.sample` took 7.033 ms, while the same count of
windows derived from a 32768 draw took 2.905 ms and the same number of blocks was decoded in both. In the second
session an equivalent manual draw took 2.205 ms. Either `sample`'s draw really lands on expensive loci, or the
7.033 ms was a transient on the machine.

This interleaves the two draws in one session, three times each, and prints the clocks beside every measurement, so a
transient cannot hide behind an average.
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


sh("apt-get -qq install -y libzstd-dev")
sh("git clone -q https://github.com/yasha1971-coder/panvram.git /content/pvrepo || true")
sh("cd /content/pvrepo && pip install -q -e .")
import importlib  # noqa: E402
import os  # noqa: E402

os.chdir("/tmp")
sys.path.insert(0, "/content/pvrepo")
importlib.invalidate_caches()
import numpy as np  # noqa: E402
import torch  # noqa: E402

import panvram  # noqa: E402

assert "/content/pvrepo/panvram/" in panvram.__file__ and panvram.with_cuda
sys.path.insert(0, "/content/pvrepo/logs/loader-2026-10-22/scripts")
from colab_prepare import prepare_cohort, _sha256  # noqa: E402

REF = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/"
       "GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz")
sh("mkdir -p /content/cohort")
sh(f'curl -fsSL -o /content/t2t.fna.gz "{REF}"')
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

ctg = [(a, rec[3], int(rec[2])) for a in range(len(c.names)) for rec in c._core.contigs(a)]
pool = [x for x in ctg if x[2] >= W]
wt = torch.tensor([x[2] - W + 1 for x in pool], dtype=torch.float64)
wt /= wt.sum()


def clocks():
    return sh("nvidia-smi --query-gpu=clocks.sm,clocks_throttle_reasons.active --format=csv,noheader").strip()


def by_sample(seed):
    gg = torch.Generator(device="cuda").manual_seed(seed)
    _, co = c.sample(BATCH, W, generator=gg, return_coords=True)
    asm = co[:, 0].contiguous()
    boff = torch.tensor([c._core.contigs(int(a))[int(ci)][3] for a, ci in zip(co[:, 0], co[:, 1])],
                        dtype=torch.int64, device=co.device)
    return asm, (boff + co[:, 2]).contiguous()


def by_hand(seed):
    rng = np.random.default_rng(seed)
    idx = torch.multinomial(wt, BATCH, replacement=True).tolist()
    asm, st = [], []
    for i in idx:
        a, boff, ln = pool[i]
        asm.append(a)
        st.append(boff + int(rng.integers(0, ln - W + 1)))
    return (torch.tensor(asm, dtype=torch.int64, device="cuda"),
            torch.tensor(st, dtype=torch.int64, device="cuda"))


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
    return ts[len(ts) // 2], min(ts), max(ts)


out = {"W": W, "batch": BATCH, "rows": []}
for rep in range(3):
    for name, fn in (("Cohort.sample", by_sample), ("manual draw", by_hand)):
        asm, start = fn(SEED + rep)
        nb = int(((start + W - 1) // Q - start // Q + 1).sum())
        uniq_asm = int(torch.unique(asm).numel())
        spread = int((start.max() - start.min()).item())
        med, lo, hi = timeit(asm, start)
        print(f"rep {rep} {name:14s} {med*1e3:7.3f} ms  min {lo*1e3:6.3f} max {hi*1e3:6.3f}  blocks {nb:6d}  "
              f"{med*1e6/nb:6.3f} us/block  assemblies {uniq_asm}  start spread {spread/1e9:.2f} G  | {clocks()}",
              flush=True)
        out["rows"].append({"rep": rep, "draw": name, "median_s": med, "min_s": lo, "max_s": hi,
                            "blocks": nb, "us_per_block": med * 1e6 / nb, "assemblies": uniq_asm})

print("\n" + json.dumps(out, indent=1), flush=True)
print("DONE", flush=True)
