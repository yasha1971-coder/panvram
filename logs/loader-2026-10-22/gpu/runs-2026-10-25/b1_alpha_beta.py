#!/usr/bin/env python3
"""T = alpha(W) + B * beta(W): where the time of a batch goes.

ONE fixed list of window starts, all valid for 32768 bases, is drawn once and reused for every batch size and both
window sizes, so nothing moves between rows except B and W. For each row three times are taken apart:
  host    - the python call with check=False: coordinate preparation and the kernel launch, no synchronisation
  kernel  - CUDA events around the same call
  status  - the extra cost of check=True, which reads the status word back and therefore synchronises
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

SEED, RUNS = 20261003, 7
BATCHES = (256, 1024, 4096, 16384)
WS = (8192, 32768)
NMAX = max(BATCHES)

c = panvram.Cohort.open("/content/cohort", device="cuda", dataset="q16k")
torch.cuda.synchronize()
Q = int(c._core.Q)
print(f"OPEN {len(c.names)} assemblies, {torch.cuda.memory_allocated()/1e9:.3f} GB, Q {Q}", flush=True)

g = torch.Generator(device="cuda").manual_seed(SEED)
t0 = time.perf_counter()
while time.perf_counter() - t0 < 5:
    c.sample(1024, 32768, generator=g)
torch.cuda.synchronize()
print(sh("nvidia-smi --query-gpu=clocks.sm,clocks_throttle_reasons.active --format=csv,noheader").strip(), flush=True)

# one fixed list of NMAX windows, every start valid for 32768 bases; every row below takes its first B entries
gg = torch.Generator(device="cuda").manual_seed(SEED)
_, co = c.sample(NMAX, 32768, generator=gg, return_coords=True)
ASM = co[:, 0].contiguous()
boff = torch.tensor([c._core.contigs(int(a))[int(ci)][3] for a, ci in zip(co[:, 0], co[:, 1])],
                    dtype=torch.int64, device=co.device)
START = (boff + co[:, 2]).contiguous()
print(f"fixed window list: {NMAX} starts, all valid for 32768", flush=True)

ev0, ev1 = torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)
out = {"seed": SEED, "runs": RUNS, "Q": Q, "assemblies": len(c.names), "rows": []}
print(f"\n{'W':>6} {'B':>6} {'blocks':>8} {'host us':>9} {'kernel us':>10} {'status us':>10} {'total us':>9} "
      f"{'us/window':>10} {'GB/s':>7}", flush=True)
for W in WS:
    for B in BATCHES:
        asm, start = ASM[:B].contiguous(), START[:B].contiguous()
        nb = int(((start + W - 1) // Q - start // Q + 1).sum())
        for _ in range(2):
            c.windows(asm, start, W)
        torch.cuda.synchronize()
        hs, ks, ts = [], [], []
        for _ in range(RUNS):
            torch.cuda.synchronize()
            a = time.perf_counter()
            ev0.record()
            c.windows(asm, start, W, check=False)      # no synchronisation: host-side cost only
            ev1.record()
            h = time.perf_counter() - a
            torch.cuda.synchronize()
            hs.append(h)
            ks.append(ev0.elapsed_time(ev1) / 1e3)
            torch.cuda.synchronize()
            a = time.perf_counter()
            c.windows(asm, start, W)                   # check=True: reads the status word back
            torch.cuda.synchronize()
            ts.append(time.perf_counter() - a)
        hs.sort(); ks.sort(); ts.sort()
        h, k, t = hs[RUNS // 2], ks[RUNS // 2], ts[RUNS // 2]
        print(f"{W:>6} {B:>6} {nb:>8} {h*1e6:>9.1f} {k*1e6:>10.1f} {max(0.0,t-k)*1e6:>10.1f} {t*1e6:>9.1f} "
              f"{t*1e6/B:>10.3f} {B*W/t/1e9:>7.3f}", flush=True)
        out["rows"].append({"W": W, "B": B, "blocks": nb, "host_s": h, "kernel_s": k, "total_s": t,
                            "status_s": max(0.0, t - k), "us_per_window": t * 1e6 / B,
                            "gb_per_s": B * W / t / 1e9})

print("\nfit T = alpha(W) + B * beta(W), least squares over the four batch sizes:", flush=True)
for W in WS:
    r = [x for x in out["rows"] if x["W"] == W]
    B = np.array([x["B"] for x in r], dtype=float)
    T = np.array([x["total_s"] for x in r]) * 1e6
    K = np.array([x["kernel_s"] for x in r]) * 1e6
    for label, Y in (("total", T), ("kernel only", K)):
        A = np.vstack([np.ones_like(B), B]).T
        al, be = np.linalg.lstsq(A, Y, rcond=None)[0]
        print(f"  W {W:5d} {label:12s}: alpha {al:8.1f} us   beta {be:7.4f} us per window   "
              f"(alpha is {100*al/(al+be*1024):.0f} % of a 1024 batch)", flush=True)
        out.setdefault("fit", []).append({"W": W, "series": label, "alpha_us": al, "beta_us": be})

print("\n" + json.dumps(out, indent=1), flush=True)
print("DONE", flush=True)
