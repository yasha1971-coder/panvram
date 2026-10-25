#!/usr/bin/env python3
"""m3_bench.py <assembly> - samples/s of the panvram loader (CPU cohort, 16 decoder threads) against torch DataLoader
over the .bin GPTDataset (what build_pretraining_data_loader returns for .bin), same sampler, same sample IDs,
micro-batch 8, 128 micro-batches after 8 of warm-up, median of 5 passes; silence gate (1-min load average < 0.5)
before every pass, the load at the start of each pass printed."""
import os, sys, time
import numpy as np, torch
sys.argv = [sys.argv[0], sys.argv[1] if len(sys.argv) > 1 else "HG00438.1", "8"]
src = open("/tmp/claude-1002/d2/m3_judge.py").read().split("FIELDS = (")[0]
g = {"__name__": "bench"}
exec(compile(src, "m3_judge.py", "exec"), g)
ds_bin, loader, MBS, SEQ = g["ds_bin"], g["loader"], g["MBS"], g["SEQ"]
NB, WARM = 128, 8
it = iter(loader.batch_sampler)
batches = [next(it) for _ in range(NB + WARM)]
def run_pan():
    for b in batches[:WARM]: loader.batch(b)
    t = time.perf_counter()
    for b in batches[WARM:]: loader.batch(b)
    return NB * MBS / (time.perf_counter() - t)
def run_bin(workers):
    dl = torch.utils.data.DataLoader(ds_bin, batch_sampler=batches, num_workers=workers,
                                     persistent_workers=False)
    t = None; k = 0
    for x in dl:
        k += 1
        if k == WARM: t = time.perf_counter()
    return NB * MBS / (time.perf_counter() - t)
def gate():
    t0 = time.time()
    while os.getloadavg()[0] >= 0.5:
        if time.time() - t0 > 3600:
            raise SystemExit("silence gate: load stayed >= 0.5 for 1 h")
        time.sleep(15)
    return os.getloadavg()[0]
for name, fn in (("panvram loader", run_pan), ("bin DataLoader workers 0", lambda: run_bin(0)),
                 ("bin DataLoader workers 2", lambda: run_bin(2)), ("bin DataLoader workers 8", lambda: run_bin(8))):
    raw, loads = [], []
    for _ in range(5):
        loads.append(gate()); raw.append(fn())
    r = sorted(raw)
    print(f"  SEQ {SEQ}  {name:26s} {r[2]:9.1f} samples/s  {r[2]*SEQ/1e6:7.2f} Mtoken/s  (min {r[0]:.1f} max {r[4]:.1f})  "
          f"raw {' '.join(f'{x:.1f}' for x in raw)}  load {' '.join(f'{x:.2f}' for x in loads)}", flush=True)
