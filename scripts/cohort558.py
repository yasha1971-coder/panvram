# SPDX-License-Identifier: MIT
"""cohort558.py <cohort dir> <manifest_v1.tsv> [dataset] - the whole cohort resident on the GPU, measured and checked
without source files. Prints C558 lines:
  open        Cohort.open time, resident bytes by component, torch allocated, nvidia-smi memory used before / after
  sample      sample(1024, 8192) over the cohort: median of 5 calls (after a warm-up), assemblies hit; == the CPU
              decoder's windows at the same coordinates
  fetch       1000 random coordinates (assembly uniform, contig by length, length 1 .. 2^20): GPU time per call
              (sync), == the CPU decoder's windows
  full        every assembly decoded in full on the GPU, copied to the host, the FASTA rebuilt from the contig table
              and hashed: XXH3 == fasta_xxh3 of the manifest (the source file's); time per assembly
  RESULT      PASS / FAIL"""
import csv
import random
import statistics
import subprocess
import sys
import time

import panvram
import torch


def smi():
    try:
        return int(subprocess.run(["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits", "-i", "0"],
                                  capture_output=True, text=True).stdout.split()[0])
    except Exception:
        return -1


def main():
    path, man = sys.argv[1], sys.argv[2]
    ds = sys.argv[3] if len(sys.argv) > 3 else "q4k"
    want = {r["name"]: r for r in csv.DictReader(open(man), delimiter="\t")}
    m0 = smi()
    t = time.time()
    c = panvram.Cohort.open(path, device="cuda", dataset=ds)
    torch.cuda.synchronize()
    t_open = time.time() - t
    m1 = smi()
    rb = c.resident_bytes()
    missing = sorted(set(want) - set(c.names))
    print(f"C558 open\t{len(c)} assemblies ({len(missing)} of the manifest missing)\t{ds} Q {c.block_size}\t{t_open:.1f} s\t"
          f"resident {rb['total'] / 1e9:.3f} GB (reference {rb['ref'] / 1e9:.3f}, payload {rb['P'] / 1e9:.3f}, block table "
          f"{(rb['off'] + rb['st']) / 1e9:.3f}, model tables {rb['tab'] / 1e9:.3f}, case runs {(rb['low_s'] + rb['low_l']) / 1e9:.3f})\t"
          f"torch allocated {torch.cuda.memory_allocated() / 1e9:.3f} GB\tnvidia-smi used {m0} -> {m1} MiB", flush=True)
    ok = not missing

    # sample
    g = torch.Generator(device="cuda").manual_seed(20261005)
    x, co = c.sample(1024, 8192, generator=g, return_coords=True)
    ts = []
    for _ in range(5):
        torch.cuda.synchronize()
        t = time.perf_counter()
        x, co = c.sample(1024, 8192, generator=g, return_coords=True)
        torch.cuda.synchronize()
        ts.append(time.perf_counter() - t)
    cl = co.cpu().tolist()
    asm = torch.tensor([r[0] for r in cl])
    start = torch.tensor([c.contigs(r[0])[r[1]][3] + r[2] for r in cl])
    y = c._core.windows(asm, start, 8192, None, False, True, True, 0)
    eq = bool(torch.equal(x.cpu(), y))
    ok &= eq
    t = statistics.median(ts)
    print(f"C558 sample\t1024 x 8192\tmedian {t * 1e3:.2f} ms per call ({1024 / t:,.0f} windows/s, {1024 * 8192 / t / 1e9:.2f} GB/s), "
          f"5 calls\t{len(set(r[0] for r in cl))} assemblies hit\t{'== CPU decoder' if eq else 'DIFFERS from the CPU decoder'}", flush=True)

    # fetch
    rng = random.Random(20261005)
    bad, ts = 0, []
    for _ in range(1000):
        a = rng.randrange(len(c))
        ctgs = c.contigs(a)
        p = rng.randrange(sum(r[2] for r in ctgs))
        ci = next(i for i, r in enumerate(ctgs) if (p := p - r[2]) < 0)
        n = ctgs[ci][2]
        L = min(n, 1 + int(2 ** (rng.random() * 20)))
        s = rng.randrange(n - L + 1)
        torch.cuda.synchronize()
        t = time.perf_counter()
        z = c.fetch(a, ci, s, L)
        torch.cuda.synchronize()
        ts.append(time.perf_counter() - t)
        w = c._core.windows(torch.tensor([a]), torch.tensor([ctgs[ci][3] + s]), L, None, False, True, True, 0)[0]
        bad += not torch.equal(z.cpu(), w)
    ok &= bad == 0
    print(f"C558 fetch\t1000 coordinates\tmedian {statistics.median(ts) * 1e3:.2f} ms, max {max(ts) * 1e3:.1f} ms per call (sync)\t"
          f"{1000 - bad} == CPU decoder, {bad} differ", flush=True)

    # full decode of every assembly on the GPU against the manifest's source XXH3
    buf = torch.empty(max(c.info(i)["bases"] for i in range(len(c))), dtype=torch.uint8, pin_memory=True)
    t0, nbad, times, tdec, nb = time.time(), 0, [], [], 0
    for i, name in enumerate(c.names):
        t = time.time()
        hb = c.decode(i, out=buf)
        td = time.time() - t
        h = c._core.fasta_xxh3(i, hb)
        dt = time.time() - t
        times.append(dt)
        tdec.append(td)
        nb += c.info(i)["bases"]
        good = name in want and h == want[name]["fasta_xxh3"]
        nbad += not good
        print(f"C558 full\t{name}\t{h}\t{'== manifest' if good else 'DIFFERS ' + want.get(name, {}).get('fasta_xxh3', '-')}\t{dt:.2f} s (GPU decode + copy {td:.2f})", flush=True)
    tt = time.time() - t0
    ok &= nbad == 0
    print(f"C558 full total\t{len(c)} assemblies, {nb / 1e12:.3f} T bases decoded on the GPU\t{len(c) - nbad} == manifest, {nbad} differ\t"
          f"{tt / 60:.1f} min (median per assembly {statistics.median(times):.2f} s; of it GPU decode + copy to host {statistics.median(tdec):.2f} s, sum {sum(tdec) / 60:.1f} min; {nb / sum(tdec) / 1e9:.2f} G bases/s)", flush=True)
    print(f"C558 RESULT\t{'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
