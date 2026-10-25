#!/usr/bin/env python3
"""#S2 - heavy-first block scheduling, no format change.

The kernel launches one warp per (window, block) slot in slot order. Here a launch-order permutation is added: slot
i of the grid runs logical block order[i]; order is argsort of the block's compressed stream length (from the block
table, off[b+1] - off[b]) descending, empty slots (b > b1) last. Computed on the GPU per call and timed inside the call.
Same build for both arms: order = none (before) and order = heavy-first (after).

Check: 256/256 windows of the heavy-first output against the CPU decoder, plus every window of every draw against the
before output (bytes equal).
Timing: the same 10 draws (seeds SEED+0..9, Cohort.sample), batch 1024 and 16384, W 8192.
Achieved parallelism: one instrumented run per draw and arm, %globaltimer at warp start and end;
sum of warp durations / (last end - first start), over warps that decoded a block.
"""
import json, subprocess, sys, textwrap
from pathlib import Path


def sh(cmd, cwd=None):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
    if r.returncode:
        print(f"$ {cmd}\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}", flush=True)
        raise SystemExit(f"failed: {cmd}")
    return r.stdout


sh("apt-get -qq install -y libzstd-dev")
sh("git clone -q https://github.com/yasha1971-coder/panvram.git /content/pvrepo || true")
print(sh("cd /content/pvrepo && git log --oneline -1").strip(), flush=True)

cu = Path("/content/pvrepo/csrc/pv_cuda.cu")
core = Path("/content/pvrepo/csrc/pv_core.cpp")
s, cc = cu.read_text(), core.read_text()
A0 = "template <bool PACKED, bool COMPACT>"
A1 = "    const uint64_t w = blockIdx.x / bpw; const uint32_t j = blockIdx.x % bpw;"
A2 = "    for (uint64_t x = x0 + lane; x < x1; x += 32) out[w * W + (x - s)] = blk[x - bst];\n}"
for i, a in enumerate((A0, A1, A2)):
    assert s.count(a) == 1, f"anchor {i} appears {s.count(a)} times - upstream changed"
s = s.replace(A0, """#define PV_NBLK 65536
__device__ const int32_t* pv_order = nullptr;      // launch slot -> logical (window, block) slot; null = identity
__device__ int pv_instr = 0;
__device__ unsigned long long pv_rec[PV_NBLK * 3];

""" + A0, 1)
s = s.replace(A1, """    const uint32_t _bid = pv_order ? (uint32_t)pv_order[blockIdx.x] : blockIdx.x;
    unsigned long long _g0 = 0; if (pv_instr && lane == 0) asm volatile("mov.u64 %0, %%globaltimer;" : "=l"(_g0));
    const uint64_t w = _bid / bpw; const uint32_t j = _bid % bpw;""", 1)
s = s.replace(A2, A2[:-1] + """    if (pv_instr) { __syncwarp();
        if (lane == 0 && blockIdx.x < PV_NBLK) { unsigned long long _g1; asm volatile("mov.u64 %0, %%globaltimer;" : "=l"(_g1));
            pv_rec[blockIdx.x * 3 + 0] = _g0; pv_rec[blockIdx.x * 3 + 1] = _g1; pv_rec[blockIdx.x * 3 + 2] = slen; } }
}""", 1)
s += """
#include <vector>
void pv_set_order(torch::Tensor t) {
    const int32_t* p = t.numel() ? t.data_ptr<int32_t>() : nullptr;
    C10_CUDA_CHECK(cudaMemcpyToSymbolAsync(pv_order, &p, sizeof(p), 0, cudaMemcpyHostToDevice, c10::cuda::getCurrentCUDAStream()));
    C10_CUDA_CHECK(cudaStreamSynchronize(c10::cuda::getCurrentCUDAStream()));
}
void pv_set_instr(int64_t on) { const int v = (int)on; C10_CUDA_CHECK(cudaMemcpyToSymbol(pv_instr, &v, sizeof(v))); }
std::vector<unsigned long long> pv_read_rec(int64_t n) {
    if (n > PV_NBLK) n = PV_NBLK;
    std::vector<unsigned long long> h((size_t)n * 3, 0);
    C10_CUDA_CHECK(cudaMemcpyFromSymbol(h.data(), pv_rec, h.size() * sizeof(unsigned long long)));
    void* a = nullptr; C10_CUDA_CHECK(cudaGetSymbolAddress(&a, pv_rec));
    C10_CUDA_CHECK(cudaMemset(a, 0, sizeof(unsigned long long) * PV_NBLK * 3));
    return h;
}
"""
cu.write_text(s)
hook = cc.index("PYBIND11_MODULE")
brace = cc.index("{", hook) + 1
cc = cc[:brace] + """
#ifdef PV_WITH_CUDA
    { extern void pv_set_order(torch::Tensor); extern void pv_set_instr(int64_t);
      extern std::vector<unsigned long long> pv_read_rec(int64_t);
      m.def("set_order", &pv_set_order); m.def("set_instr", &pv_set_instr); m.def("read_rec", &pv_read_rec); }
#endif
""" + cc[brace:]
core.write_text(cc)
print(sh("cd /content/pvrepo && git diff --stat -- csrc/").strip(), flush=True)
r = subprocess.run([sys.executable, "-m", "pip", "install", "-q", "-e", "."], cwd="/content/pvrepo",
                   capture_output=True, text=True)
if r.returncode:
    print((r.stdout + r.stderr)[-3000:], flush=True)
    raise SystemExit("BUILD FAILED")
print("built", flush=True)

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

script = textwrap.dedent('''
    import json, subprocess, sys, time
    sys.path.insert(0, "/content/pvrepo")
    import numpy as np, torch, panvram
    from panvram import _C
    assert "/content/pvrepo/panvram/" in panvram.__file__ and panvram.with_cuda
    assert hasattr(_C, "set_order") and hasattr(_C, "read_rec")
    print("panvram", panvram.__file__, flush=True)
    SEED, NDRAW, W = 20261003, 10, 8192
    c = panvram.Cohort.open("/content/cohort", device="cuda", dataset="q16k")
    torch.cuda.synchronize()
    Q = int(c._core.Q); d = c._dev; dev = d["off"].device
    assert not (c._core.packed_reference or c._core.compact_blocks)
    OFF = d["off"].long() & 0xFFFFFFFF; BB = d["blk_base"]
    bpw = (W + Q - 1) // Q + 1
    JJ = torch.arange(bpw, device=dev)
    _C.set_order(torch.empty(0, dtype=torch.int32, device=dev)); _C.set_instr(0)
    g = torch.Generator(device="cuda").manual_seed(SEED)
    t0 = time.perf_counter()
    while time.perf_counter() - t0 < 5:
        c.sample(4096, 32768, generator=g)
    torch.cuda.synchronize()
    print(subprocess.run("nvidia-smi --query-gpu=name,clocks.sm,clocks_throttle_reasons.active --format=csv,noheader",
                         shell=True, capture_output=True, text=True).stdout.strip(), flush=True)

    def heavy_order(asm, start):
        b0 = start // Q; b1 = (start + W - 1) // Q
        b = b0[:, None] + JJ
        valid = b <= b1[:, None]
        bb = BB[asm][:, None] + torch.minimum(b, b1[:, None])
        sl = (OFF[bb + 1] - OFF[bb]) % (1 << 32)
        sl = torch.where(valid, sl, torch.full_like(sl, -1)).flatten()
        return torch.argsort(sl, descending=True, stable=True).to(torch.int32)

    NONE = torch.empty(0, dtype=torch.int32, device=dev)
    keep = []
    def run(asm, start, heavy):
        if heavy:
            o = heavy_order(asm, start); keep.append(o); _C.set_order(o)
        else:
            _C.set_order(NONE)
        x = c.windows(asm, start, W)
        if heavy: keep.clear()
        return x

    def coords_of(seed, B):
        gg = torch.Generator(device="cuda").manual_seed(seed)
        _, co = c.sample(B, W, generator=gg, return_coords=True)
        asm = co[:, 0].contiguous()
        boff = torch.tensor([c._core.contigs(int(a))[int(ci)][3] for a, ci in zip(co[:, 0], co[:, 1])],
                            dtype=torch.int64, device=co.device)
        return asm, (boff + co[:, 2]).contiguous()

    def timeit(asm, start, heavy, runs):
        for _ in range(2): run(asm, start, heavy)
        torch.cuda.synchronize(); ts = []
        for _ in range(runs):
            torch.cuda.synchronize(); a = time.perf_counter()
            run(asm, start, heavy)
            torch.cuda.synchronize(); ts.append(time.perf_counter() - a)
        ts.sort(); return ts[len(ts) // 2]

    def parallelism(asm, start, heavy):
        n = asm.numel() * bpw
        _C.read_rec(n); _C.set_instr(1)
        run(asm, start, heavy); torch.cuda.synchronize()
        _C.set_instr(0)
        r = np.array(_C.read_rec(n), dtype=np.uint64).reshape(-1, 3)
        live = r[:, 1] > 0
        t0, t1 = r[live, 0].astype(np.float64), r[live, 1].astype(np.float64)
        span = t1.max() - t0.min(); dur = t1 - t0
        return {"warps_decoding": int(live.sum()), "span_us": span / 1e3, "sum_dur_us": dur.sum() / 1e3,
                "parallelism": dur.sum() / span, "max_warp_us": dur.max() / 1e3,
                "last_start_us": (t0.max() - t0.min()) / 1e3}

    # ---- check: 256/256 windows heavy-first vs the CPU decoder; every window heavy-first vs before
    asm, start = coords_of(SEED, 1024)
    xh = run(asm, start, True); xb = run(asm, start, False); torch.cuda.synchronize()
    cpu = c._core.windows(asm[:256].cpu(), start[:256].cpu(), W, None, False, True, False, 0)
    ok_cpu = int((xh[:256].cpu() == cpu).all(dim=1).sum())
    ok_all = int((xh == xb).all(dim=1).sum())
    print(f"CHECK heavy-first vs CPU decoder: {ok_cpu}/256 windows identical; heavy-first vs before: {ok_all}/1024", flush=True)
    assert ok_cpu == 256 and ok_all == 1024
    t = torch.cuda.Event(enable_timing=True); u = torch.cuda.Event(enable_timing=True)

    out = {"seed": SEED, "draws": NDRAW, "W": W, "Q": Q, "check_cpu": f"{ok_cpu}/256", "rows": {}}
    for B, RUNS in ((1024, 5), (16384, 3)):
        per = []
        for k in range(NDRAW):
            asm, start = coords_of(SEED + k, B)
            xb = run(asm, start, False); xh = run(asm, start, True); torch.cuda.synchronize()
            same = int((xh == xb).all(dim=1).sum())
            assert same == B, f"draw {k}: heavy-first differs from before in {B - same} windows"
            for _ in range(3): heavy_order(asm, start)
            torch.cuda.synchronize(); t.record()
            for _ in range(10): heavy_order(asm, start)
            u.record(); torch.cuda.synchronize(); ord_us = t.elapsed_time(u) / 10 * 1e3
            tb = timeit(asm, start, False, RUNS); th = timeit(asm, start, True, RUNS)
            tb2 = timeit(asm, start, False, RUNS)                       # before again: drift check
            pb = parallelism(asm, start, False); ph = parallelism(asm, start, True)
            per.append({"draw": k, "seed": SEED + k, "before_ms": tb * 1e3, "before_again_ms": tb2 * 1e3,
                        "heavy_ms": th * 1e3, "order_us": ord_us, "identical": f"{same}/{B}",
                        "par_before": pb, "par_heavy": ph})
            print(f"  B {B:5d} draw {k}  before {tb*1e3:7.3f} ms (again {tb2*1e3:7.3f})  heavy-first {th*1e3:7.3f} ms "
                  f"(order {ord_us:6.1f} us)  parallelism {pb['parallelism']:6.1f} -> {ph['parallelism']:6.1f}  "
                  f"span {pb['span_us']:8.1f} -> {ph['span_us']:8.1f} us  max warp {pb['max_warp_us']:7.1f} us  "
                  f"identical {same}/{B}", flush=True)
        def st(v):
            v = np.array(v); return {"median": float(np.median(v)), "p90": float(np.percentile(v, 90)),
                                     "max": float(v.max()), "min": float(v.min())}
        row = {"before_ms": st([x["before_ms"] for x in per]), "heavy_ms": st([x["heavy_ms"] for x in per]),
               "before_again_ms": st([x["before_again_ms"] for x in per]),
               "par_before": st([x["par_before"]["parallelism"] for x in per]),
               "par_heavy": st([x["par_heavy"]["parallelism"] for x in per]),
               "span_before_us": st([x["par_before"]["span_us"] for x in per]),
               "span_heavy_us": st([x["par_heavy"]["span_us"] for x in per]),
               "order_us": st([x["order_us"] for x in per]), "per_draw": per}
        out["rows"][str(B)] = row
        for key, lab in (("before_ms", "before     "), ("before_again_ms", "before again"), ("heavy_ms", "heavy-first")):
            s = row[key]; print(f"  B {B:5d} {lab} over {NDRAW} draws: median {s['median']:.3f} ms  p90 {s['p90']:.3f}  "
                                f"max {s['max']:.3f}  min {s['min']:.3f}", flush=True)
        for key, lab in (("par_before", "parallelism before"), ("par_heavy", "parallelism heavy-first")):
            s = row[key]; print(f"  B {B:5d} {lab}: median {s['median']:.1f}  p90 {s['p90']:.1f}  min {s['min']:.1f}  max {s['max']:.1f}", flush=True)
    json.dump(out, open("/content/s2_heavy.json", "w"), indent=1)
    print("S2_DONE", flush=True)
''')
Path("/tmp/s2run.py").write_text(script)
r = subprocess.run([sys.executable, "-P", "/tmp/s2run.py"], cwd="/tmp", capture_output=True, text=True)
print(r.stdout[-12000:], flush=True)
if r.returncode:
    print(r.stderr[-3000:], flush=True)
print("DONE", flush=True)
