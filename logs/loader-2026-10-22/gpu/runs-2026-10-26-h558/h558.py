#!/usr/bin/env python3
"""558 assemblies resident. The headline GPU number at the batch that actually saturates an A100: 16 384 windows, ten independent draws per
window size, reported as a distribution. Earlier runs timed one draw nine times, which measured the draw.
"""
import json, subprocess, sys, time
from pathlib import Path


def sh(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode:
        print(f"$ {cmd}\n{r.stdout[-1200:]}\n{r.stderr[-1200:]}", flush=True)
        raise SystemExit(f"failed: {cmd}")
    return r.stdout


T00 = time.perf_counter()
import threading
def _download():
    # the 16 parts, 4 at a time, while the extension builds; prepare_cohort() verifies every part's SHA-256
    REC = "https://zenodo.org/records/23232317/files"
    os.makedirs("/content/dl", exist_ok=True)
    t = time.perf_counter()
    r = subprocess.run(f"seq -w 0 15 | xargs -P 4 -I{{}} curl -fsSL --retry 5 -o /content/dl/cohort558_q16k.tar.part{{}} "
                       f"{REC}/cohort558_q16k.tar.part{{}}", shell=True, capture_output=True, text=True)
    DL["s"] = time.perf_counter() - t; DL["rc"] = r.returncode; DL["err"] = r.stderr[-500:]
DL = {}
import os
_th = threading.Thread(target=_download); _th.start()
sh("apt-get -qq install -y libzstd-dev")
sh("git clone -q https://github.com/yasha1971-coder/panvram.git /content/pvrepo || true")
sh("cd /content/pvrepo && pip install -q -e .")
import importlib, os  # noqa: E402

os.chdir("/tmp"); sys.path.insert(0, "/content/pvrepo"); importlib.invalidate_caches()
import numpy as np, torch, panvram  # noqa: E402

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
_th.join()
print(f"DOWNLOAD 16 parts, 4 streams: {DL['s']:.1f} s, rc {DL['rc']} {DL['err']}", flush=True)
print(f"setup so far {time.perf_counter()-T00:.0f} s", flush=True)
T_PREP = time.perf_counter()
prepare_cohort(558, dest="/content/cohort", work="/content/dl")
print(f"PREPARED 558 in {time.perf_counter()-T_PREP:.0f} s", flush=True)
sh("rm -rf /content/dl")

SEED, BATCH, RUNS, NDRAW = 20261003, 16384, 3, 10
t_open = time.perf_counter()
c = panvram.Cohort.open("/content/cohort", device="cuda", dataset="q16k")
torch.cuda.synchronize()
print(f"OPEN {len(c.names)} assemblies in {time.perf_counter()-t_open:.3f} s, resident {torch.cuda.memory_allocated()/1e9:.3f} GB", flush=True)
assert len(c.names) == 558
Q = int(c._core.Q)
g = torch.Generator(device="cuda").manual_seed(SEED)
t0 = time.perf_counter()
while time.perf_counter() - t0 < 5:
    c.sample(4096, 32768, generator=g)
torch.cuda.synchronize()
print(sh("nvidia-smi --query-gpu=clocks.sm,clocks_throttle_reasons.active --format=csv,noheader").strip(), flush=True)
print(f"batch {BATCH}, {NDRAW} draws per W, {RUNS} timings per draw", flush=True)

out = {"seed": SEED, "batch": BATCH, "draws": NDRAW, "runs_per_draw": RUNS, "Q": Q,
       "assemblies": len(c.names), "open_s": None, "gpu": torch.cuda.get_device_name(0), "dist": {}}
for W in (8192, 32768):
    per = []
    for d in range(NDRAW):
        gg = torch.Generator(device="cuda").manual_seed(SEED + d)
        _, co = c.sample(BATCH, W, generator=gg, return_coords=True)
        asm = co[:, 0].contiguous()
        boff = torch.tensor([c._core.contigs(int(a))[int(ci)][3] for a, ci in zip(co[:, 0], co[:, 1])],
                            dtype=torch.int64, device=co.device)
        start = (boff + co[:, 2]).contiguous()
        nb = int(((start + W - 1) // Q - start // Q + 1).sum())
        if d == 0:
            x = c.windows(asm[:256], start[:256], W)
            cpu = c._core.windows(asm[:256].cpu(), start[:256].cpu(), W, None, False, True, False, 0)
            ok = int((x.cpu() == cpu).all(dim=1).sum())
            print(f"CHECK W {W} draw 0: {ok}/256 windows identical to the CPU decoder", flush=True)
            assert ok == 256
        for _ in range(2):
            c.windows(asm, start, W)
        torch.cuda.synchronize()
        ts = []
        for _ in range(RUNS):
            torch.cuda.synchronize(); a = time.perf_counter()
            c.windows(asm, start, W)
            torch.cuda.synchronize(); ts.append(time.perf_counter() - a)
        ts.sort(); med = ts[len(ts) // 2]
        per.append({"draw": d, "seed": SEED + d, "median_s": med, "blocks": nb,
                    "windows_per_s": BATCH / med, "gb_per_s": BATCH * W / med / 1e9,
                    "us_per_block": med * 1e6 / nb})
        print(f"  W {W:6d} draw {d}  {med*1e3:8.3f} ms  {BATCH/med:9.0f} windows/s  "
              f"{BATCH*W/med/1e9:6.3f} GB/s  blocks {nb:7d}  {med*1e6/nb:.3f} us/block", flush=True)
    t = np.array([x["median_s"] for x in per]) * 1e3
    wps = np.array([x["windows_per_s"] for x in per])
    gbs = np.array([x["gb_per_s"] for x in per])
    print(f"  W {W:6d} over {NDRAW} draws: time median {np.median(t):.3f} ms  p90 {np.percentile(t,90):.3f}  "
          f"max {t.max():.3f}  min {t.min():.3f}  max/min {t.max()/t.min():.2f}x", flush=True)
    print(f"  W {W:6d}            windows/s median {np.median(wps):.0f}  p90 {np.percentile(wps,10):.0f} (slow end)  "
          f"worst {wps.min():.0f} | GB/s median {np.median(gbs):.3f}  worst {gbs.min():.3f}", flush=True)
    out["dist"][str(W)] = {"per_draw": per, "time_median_ms": float(np.median(t)),
                           "time_p90_ms": float(np.percentile(t, 90)), "time_max_ms": float(t.max()),
                           "time_min_ms": float(t.min()),
                           "windows_per_s_median": float(np.median(wps)), "windows_per_s_min": float(wps.min()),
                           "gb_per_s_median": float(np.median(gbs)), "gb_per_s_min": float(gbs.min())}
Path("/content/h558.json").write_text(json.dumps(out, indent=1))
print("\n" + json.dumps(out["dist"], indent=1)[:1200], flush=True)
print("DONE", flush=True)
