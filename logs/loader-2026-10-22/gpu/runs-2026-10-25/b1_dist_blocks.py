#!/usr/bin/env python3
"""Two things in one session.

PART 1 - the distribution over draws. Ten independent window sets per W, timed separately, so the spread between
draws is reported as a distribution (median, p90, max) instead of nine repeats of one draw, which is what the earlier
runs measured and why their 2 % "spread" meant nothing.

PART 2 - does a batch wait for its slowest block? The kernel is rebuilt with a per-block record: lane 0 writes the
cycles it spent decoding its block, the number of operations it emitted and how many bases came from literals, from
self copies and from the reference. Then the draw that takes 7.033 ms (seed 20261003, W = 8192) is run and its
slowest blocks are listed with their content.
"""
import json
import subprocess
import sys
import textwrap
import time
from pathlib import Path


def sh(cmd, cwd=None):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
    if r.returncode:
        print(f"$ {cmd}\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}", flush=True)
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

assert "/content/pvrepo/panvram/" in panvram.__file__ and panvram.with_cuda, panvram.__file__
print("panvram", panvram.__file__, "with_cuda", panvram.with_cuda, flush=True)
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

SEED, BATCH, RUNS, NDRAW = 20261003, 1024, 5, 10
c = panvram.Cohort.open("/content/cohort", device="cuda", dataset="q16k")
torch.cuda.synchronize()
Q = int(c._core.Q)
g = torch.Generator(device="cuda").manual_seed(SEED)
t0 = time.perf_counter()
while time.perf_counter() - t0 < 5:
    c.sample(1024, 32768, generator=g)
torch.cuda.synchronize()
print(sh("nvidia-smi --query-gpu=clocks.sm,clocks_throttle_reasons.active --format=csv,noheader").strip(), flush=True)


def coords_of(seed, W):
    gg = torch.Generator(device="cuda").manual_seed(seed)
    _, co = c.sample(BATCH, W, generator=gg, return_coords=True)
    asm = co[:, 0].contiguous()
    boff = torch.tensor([c._core.contigs(int(a))[int(ci)][3] for a, ci in zip(co[:, 0], co[:, 1])],
                        dtype=torch.int64, device=co.device)
    return asm, (boff + co[:, 2]).contiguous()


def timeit(asm, start, W, runs=RUNS):
    for _ in range(2):
        c.windows(asm, start, W)
    torch.cuda.synchronize()
    ts = []
    for _ in range(runs):
        torch.cuda.synchronize()
        a = time.perf_counter()
        c.windows(asm, start, W)
        torch.cuda.synchronize()
        ts.append(time.perf_counter() - a)
    ts.sort()
    return ts[len(ts) // 2]


print(f"\n=== PART 1: {NDRAW} independent draws per W, batch {BATCH}", flush=True)
out = {"seed": SEED, "batch": BATCH, "draws": NDRAW, "runs_per_draw": RUNS, "Q": Q, "dist": {}}
for W in (8192, 32768):
    per = []
    for d in range(NDRAW):
        asm, start = coords_of(SEED + d, W)
        nb = int(((start + W - 1) // Q - start // Q + 1).sum())
        med = timeit(asm, start, W)
        per.append({"draw": d, "seed": SEED + d, "median_s": med, "blocks": nb,
                    "us_per_block": med * 1e6 / nb, "windows_per_s": BATCH / med})
        print(f"  W {W:6d} draw {d} seed {SEED+d}  {med*1e3:7.3f} ms  blocks {nb:6d}  "
              f"{med*1e6/nb:6.3f} us/block", flush=True)
    t = np.array([x["median_s"] for x in per]) * 1e3
    b = np.array([x["blocks"] for x in per], dtype=float)
    print(f"  W {W:6d}  over {NDRAW} draws: median {np.median(t):.3f} ms  p90 {np.percentile(t,90):.3f}  "
          f"max {t.max():.3f}  min {t.min():.3f}  max/min {t.max()/t.min():.2f}x   "
          f"blocks {b.mean():.0f}+-{b.std():.0f}", flush=True)
    out["dist"][str(W)] = {"per_draw": per, "median_ms": float(np.median(t)), "p90_ms": float(np.percentile(t, 90)),
                           "max_ms": float(t.max()), "min_ms": float(t.min()),
                           "blocks_mean": float(b.mean()), "blocks_std": float(b.std())}

Path("/content/part1.json").write_text(json.dumps(out, indent=1))
print("PART1_DONE", flush=True)

# ---------------------------------------------------------------- PART 2
print("\n=== PART 2: per-block cycles, instrumented kernel", flush=True)
cu = Path("/content/pvrepo/csrc/pv_cuda.cu")
core = Path("/content/pvrepo/csrc/pv_core.cpp")
s, cc = cu.read_text(), core.read_text()
A0 = "template <bool PACKED, bool COMPACT>"
A1 = "        s_nk = r3_decode_stream(src, slen, D.T + a, ref, D.ref_n, blen, g0, sk);"
A2 = "    const uint32_t lane = threadIdx.x;"
for i, a in enumerate((A0, A1, A2)):
    assert s.count(a) == 1, f"anchor {i} appears {s.count(a)} times - upstream changed"
s = s.replace(A0, "#define PV_NBLK 65536\n__device__ unsigned long long pv_blk[PV_NBLK * 5];\n\n" + A0, 1)
s = s.replace(A1, """        const unsigned long long _t0 = clock64();
""" + A1 + """
        if (blockIdx.x < PV_NBLK) { const unsigned long long _d = clock64() - _t0;
            pv_blk[blockIdx.x * 5 + 0] = _d; pv_blk[blockIdx.x * 5 + 1] = (unsigned long long)s_nk;
            pv_blk[blockIdx.x * 5 + 2] = (unsigned long long)blen; pv_blk[blockIdx.x * 5 + 3] = (unsigned long long)slen;
            pv_blk[blockIdx.x * 5 + 4] = (unsigned long long)a; }""", 1)
cu.write_text(s)
hook = cc.index("PYBIND11_MODULE")
brace = cc.index("{", hook) + 1
cc = cc[:brace] + """
#ifdef PV_WITH_CUDA
    { extern std::vector<unsigned long long> pv_read_blocks(int64_t n);
      m.def("read_blocks", &pv_read_blocks); }
#endif
""" + cc[brace:]
core.write_text(cc)
cu.write_text(cu.read_text() + """
#include <vector>
std::vector<unsigned long long> pv_read_blocks(int64_t n) {
    if (n > PV_NBLK) n = PV_NBLK;
    std::vector<unsigned long long> h((size_t)n * 5, 0);
    cudaMemcpyFromSymbol(h.data(), pv_blk, h.size() * sizeof(unsigned long long));
    unsigned long long z = 0;
    for (size_t i = 0; i < h.size(); i++) cudaMemcpyToSymbol(pv_blk, &z, sizeof(z), i * sizeof(z));
    return h;
}
""")
print(sh("cd /content/pvrepo && git diff --stat -- csrc/").strip(), flush=True)
r = subprocess.run([sys.executable, "-m", "pip", "install", "-q", "-e", "."], cwd="/content/pvrepo",
                   capture_output=True, text=True)
if r.returncode:
    print((r.stdout + r.stderr)[-2000:], flush=True)
    print("PART 2 DID NOT BUILD - report part 1 only", flush=True)
    print("DONE", flush=True)
    raise SystemExit(0)

script = textwrap.dedent('''
    import json, sys, time
    sys.path.insert(0, "/content/pvrepo")
    import numpy as np, torch, panvram
    from panvram import _C
    assert "/content/pvrepo/panvram/" in panvram.__file__
    if not hasattr(_C, "read_blocks"):
        print("NO_READ_BLOCKS", file=sys.stderr); sys.exit(3)
    SEED, BATCH, W = 20261003, 1024, 8192
    c = panvram.Cohort.open("/content/cohort", device="cuda", dataset="q16k")
    Q = int(c._core.Q)
    gg = torch.Generator(device="cuda").manual_seed(SEED)
    _, co = c.sample(BATCH, W, generator=gg, return_coords=True)
    asm = co[:, 0].contiguous()
    boff = torch.tensor([c._core.contigs(int(a))[int(ci)][3] for a, ci in zip(co[:, 0], co[:, 1])],
                        dtype=torch.int64, device=co.device)
    start = (boff + co[:, 2]).contiguous()
    bpw = (W + Q - 1) // Q + 1
    for _ in range(2):
        c.windows(asm, start, W)
    torch.cuda.synchronize()
    _C.read_blocks(BATCH * bpw)                     # clear
    t = time.perf_counter(); c.windows(asm, start, W); torch.cuda.synchronize()
    wall = time.perf_counter() - t
    raw = np.array(_C.read_blocks(BATCH * bpw), dtype=np.uint64).reshape(-1, 5)
    cyc = raw[:, 0].astype(np.float64); nops = raw[:, 1].astype(np.float64); slen = raw[:, 3].astype(np.float64)
    live = cyc > 0
    print(f"batch wall {wall*1e3:.3f} ms; warps launched {len(cyc)}, decoding {int(live.sum())}")
    cl = cyc[live]
    print(f"cycles per block: median {np.median(cl):.0f}  p90 {np.percentile(cl,90):.0f}  "
          f"p99 {np.percentile(cl,99):.0f}  max {cl.max():.0f}  min {cl.min():.0f}  max/median {cl.max()/np.median(cl):.1f}x")
    print(f"at 1410 MHz that is median {np.median(cl)/1410:.1f} us, max {cl.max()/1410:.1f} us")
    print(f"sum of all block cycles / wall: {cl.sum()/1410/1e3/wall/1e3:.1f}x (concurrency actually achieved)")
    o = np.argsort(cl)[::-1][:10]
    sl = slen[live]; no = nops[live]
    print("slowest ten blocks: cycles, ops emitted, payload bytes of the block")
    for i in o:
        print(f"   {cl[i]:12.0f}  ops {no[i]:6.0f}  payload {sl[i]:7.0f} B")
    print(f"ops emitted: median {np.median(no):.0f} max {no.max():.0f}; payload bytes: median {np.median(sl):.0f} max {sl.max():.0f}")
    import numpy.linalg as la
    for name, v in (("ops", no), ("payload bytes", sl)):
        cc2 = float(np.corrcoef(v, cl)[0, 1])
        print(f"correlation of cycles with {name}: {cc2:+.3f}")
    json.dump({"wall_s": wall, "cycles": cl.tolist()[:4096], "nops": no.tolist()[:4096], "payload": sl.tolist()[:4096]},
              open("/content/part2.json", "w"))
    print("PART2_DONE")
''')
Path("/tmp/part2.py").write_text(script)
r = subprocess.run([sys.executable, "-P", "/tmp/part2.py"], cwd="/tmp", capture_output=True, text=True)
print(r.stdout[-4000:], flush=True)
if r.returncode:
    print(r.stderr[-2500:], flush=True)
print("DONE", flush=True)
