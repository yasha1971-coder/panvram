t0 = time.time()
import os, sys, importlib
os.chdir("/tmp"); sys.path.insert(0, "/content/pvrepo"); importlib.invalidate_caches()
import numpy as np, torch, panvram
assert "/content/pvrepo/panvram/" in panvram.__file__ and panvram.with_cuda
c = panvram.Cohort.open("/content/cohort", device="cuda", dataset="q16k"); torch.cuda.synchronize()
print(f"OPEN {len(c.names)} assemblies, resident {torch.cuda.memory_allocated()/1e9:.3f} GB", flush=True)
assert len(c.names) == 558
stage("cohort_open", t0); t0 = time.time()
Q = int(c._core.Q); SEED, BATCH, RUNS = 20261003, 16384, 3
g = torch.Generator(device="cuda").manual_seed(SEED); tt = time.perf_counter()
while time.perf_counter() - tt < 5: c.sample(4096, 32768, generator=g)
torch.cuda.synchronize()
print(sh("nvidia-smi --query-gpu=clocks.sm,clocks_throttle_reasons.active --format=csv,noheader").strip(), flush=True)
REF0 = {8192: 12.255, 32768: 20.460}           # h558.log, draw 0 (seed 20261003), same batch
for W in (8192, 32768):
    gg = torch.Generator(device="cuda").manual_seed(SEED)
    _, co = c.sample(BATCH, W, generator=gg, return_coords=True)
    asm = co[:, 0].contiguous()
    boff = torch.tensor([c._core.contigs(int(a))[int(ci)][3] for a, ci in zip(co[:, 0], co[:, 1])], dtype=torch.int64, device=co.device)
    start = (boff + co[:, 2]).contiguous()
    x = c.windows(asm[:256], start[:256], W)
    cpu = c._core.windows(asm[:256].cpu(), start[:256].cpu(), W, None, False, True, False, 0)
    ok = int((x.cpu() == cpu).all(dim=1).sum()); print(f"CHECK W {W}: {ok}/256 windows identical to the CPU decoder", flush=True); assert ok == 256
    for _ in range(2): c.windows(asm, start, W)
    torch.cuda.synchronize(); ts = []
    for _ in range(RUNS):
        torch.cuda.synchronize(); a = time.perf_counter(); c.windows(asm, start, W); torch.cuda.synchronize(); ts.append(time.perf_counter() - a)
    ts.sort(); med = ts[len(ts)//2] * 1e3
    print(f"W {W} draw 0 seed {SEED}: runs {' '.join(f'{t*1e3:.3f}' for t in ts)} ms, median {med:.3f} ms, {BATCH/med*1e3:.0f} windows/s, "
          f"{BATCH*W/med/1e6:.3f} GB/s; h558.log draw 0 {REF0[W]:.3f} ms -> {100*(med/REF0[W]-1):+.1f} %", flush=True)
stage("measure", t0)
print("DONE", flush=True)
