#!/usr/bin/env python3
"""The decisive test for the 8 KiB anomaly of logs/loader-2026-10-22/gpu/ANOMALY.md.

The two v4 rows sampled their own windows, so the 8 KiB batch and the 32 KiB batch touched different loci. This script
samples the windows ONCE at W = 32768 and then times the very same starts truncated to 8192, so the loci, the
assemblies and the locality are identical and only W differs. If microseconds per decoded block converge, locality
between the warps of one window explains the gap; if the 8 KiB number stays near its v4 value, it does not.

Run: colab run --gpu A100 anomaly_locality.py
"""
import json
import subprocess
import sys
import time


def sh(cmd, **kw):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)
    if r.returncode:
        print(f"$ {cmd}\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}", flush=True)
        raise SystemExit(f"failed: {cmd}")
    return r.stdout


print(sh("nvidia-smi --query-gpu=name,memory.total,driver_version,clocks.max.sm --format=csv,noheader").strip(),
      flush=True)
print("--- installing", flush=True)
sh("apt-get -qq install -y libzstd-dev")
sh("git clone -q https://github.com/yasha1971-coder/panvram.git /content/pvrepo || true")   # NOT a directory called
print(sh("cd /content/pvrepo && git log --oneline -1").strip(), flush=True)                 # panvram: that shadows it
sh("cd /content/pvrepo && pip install -q -e .")

import importlib
import os
os.chdir("/tmp")                        # leave /content so no sibling directory can shadow an import
sys.path.insert(0, "/content/pvrepo")   # an editable install writes a .pth, and a .pth is only read at interpreter
importlib.invalidate_caches()           # start: this kernel is already running, so put the root on the path by hand

import torch  # noqa: E402

import panvram  # noqa: E402

print("panvram", panvram.__file__, "| torch", torch.__version__, "| with_cuda",
      getattr(panvram, "with_cuda", None), flush=True)
assert "/content/pvrepo/panvram/" in panvram.__file__, f"shadowed or wrong import: {panvram.__file__}"
assert panvram.with_cuda, "the CUDA extension did not build"
print("--- import verified before any download", flush=True)

sys.path.insert(0, "/content/pvrepo/logs/loader-2026-10-22/scripts")
from colab_prepare import prepare_cohort, _sha256  # noqa: E402

print("--- reference", flush=True)
REF_URL = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/"
           "GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz")
sh("mkdir -p /content/cohort")
sh(f'curl -fsSL -o /content/t2t.fna.gz "{REF_URL}"')
sh("gunzip -f /content/t2t.fna.gz && mv /content/t2t.fna /content/cohort/reference.fa")
h = _sha256("/content/cohort/reference.fa")
print("reference sha256", h, flush=True)
assert h == "51ab7d3c2a375ecff53ee2a1be03e58d57cc1601c948a860c03a78d724f26e4c", "wrong reference"

print("--- archives", flush=True)
paths = prepare_cohort(5, dest="/content/cohort", work="/content/dl")
print(f"{len(paths)} archives", flush=True)

import torch  # noqa: E402

BATCH, SEED, RUNS = 1024, 20261003, 9
t0 = time.perf_counter()
c = panvram.Cohort.open("/content/cohort", device="cuda", dataset="q16k")
torch.cuda.synchronize()
print(f"OPEN {len(c.names)} assemblies {time.perf_counter()-t0:.3f} s, "
      f"{torch.cuda.memory_allocated()/1e9:.3f} GB", flush=True)

Q = int(c._core.Q)
print("ramping", flush=True)
g = torch.Generator(device="cuda").manual_seed(SEED)
t0 = time.perf_counter()
while time.perf_counter() - t0 < 10:
    c.sample(BATCH, 32768, generator=g)
torch.cuda.synchronize()
print(sh("nvidia-smi --query-gpu=clocks.sm,clocks_throttle_reasons.active --format=csv,noheader").strip(), flush=True)


def blocks_of(asm, start, W):
    return int(((start + W - 1) // Q - start // Q + 1).sum())


def timeit(asm, start, W):
    for _ in range(3):
        c.windows(asm, start, W)
    torch.cuda.synchronize()
    ts = []
    for _ in range(RUNS):
        torch.cuda.synchronize()
        a = time.perf_counter()
        x = c.windows(asm, start, W)
        torch.cuda.synchronize()
        ts.append(time.perf_counter() - a)
        assert x.shape == (BATCH, W)
    ts.sort()
    return ts[len(ts) // 2], min(ts), max(ts)


# the windows are drawn ONCE at 32768; the same (assembly, start) pairs are then used for both window sizes
g = torch.Generator(device="cuda").manual_seed(SEED)
_, coords = c.sample(BATCH, 32768, generator=g, return_coords=True)
asm = coords[:, 0].contiguous()
boff = torch.tensor([c._core.contigs(int(a))[int(ci)][3] for a, ci in zip(coords[:, 0], coords[:, 1])],
                    dtype=torch.int64, device=coords.device)
start = (boff + coords[:, 2]).contiguous()

out = {"gpu": torch.cuda.get_device_name(0), "batch": BATCH, "runs": RUNS, "Q": Q,
       "assemblies": len(c.names), "rows": []}
print("\n=== same windows, two sizes (the decisive comparison)", flush=True)
for W in (32768, 8192, 32768, 8192):
    med, lo, hi = timeit(asm, start, W)
    nb = blocks_of(asm, start, W)
    print(f"same-windows  W {W:6d}  {med*1e3:7.3f} ms  blocks {nb:6d}  {med*1e6/nb:6.3f} us per block  "
          f"{BATCH/med:9.0f} windows/s  {BATCH*W/med/1e9:6.3f} GB/s  min {lo*1e3:.3f} max {hi*1e3:.3f}", flush=True)
    out["rows"].append({"mode": "same_windows", "W": W, "median_s": med, "min_s": lo, "max_s": hi,
                        "blocks": nb, "us_per_block": med * 1e6 / nb})

print("\n=== each size draws its own windows (what v4 did)", flush=True)
for W in (32768, 8192):
    g2 = torch.Generator(device="cuda").manual_seed(SEED)
    _, cc = c.sample(BATCH, W, generator=g2, return_coords=True)
    a2 = cc[:, 0].contiguous()
    b2 = torch.tensor([c._core.contigs(int(a))[int(ci)][3] for a, ci in zip(cc[:, 0], cc[:, 1])],
                      dtype=torch.int64, device=cc.device)
    s2 = (b2 + cc[:, 2]).contiguous()
    med, lo, hi = timeit(a2, s2, W)
    nb = blocks_of(a2, s2, W)
    print(f"own-windows   W {W:6d}  {med*1e3:7.3f} ms  blocks {nb:6d}  {med*1e6/nb:6.3f} us per block  "
          f"{BATCH/med:9.0f} windows/s  {BATCH*W/med/1e9:6.3f} GB/s", flush=True)
    out["rows"].append({"mode": "own_windows", "W": W, "median_s": med, "min_s": lo, "max_s": hi,
                        "blocks": nb, "us_per_block": med * 1e6 / nb})

print("\n=== one warp per window: W = 16384, every window touches 1 or 2 blocks with bpw = 2", flush=True)
for W in (16384,):
    med, lo, hi = timeit(asm, start, W)
    nb = blocks_of(asm, start, W)
    print(f"same-windows  W {W:6d}  {med*1e3:7.3f} ms  blocks {nb:6d}  {med*1e6/nb:6.3f} us per block", flush=True)
    out["rows"].append({"mode": "same_windows", "W": W, "median_s": med, "blocks": nb,
                        "us_per_block": med * 1e6 / nb})

print("\n" + json.dumps(out, indent=1), flush=True)
open("/content/anomaly_result.json", "w").write(json.dumps(out, indent=1))
print("DONE", flush=True)
