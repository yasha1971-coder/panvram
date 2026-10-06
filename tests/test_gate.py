"""The public gate (run before anything of panvram goes out):
  import panvram -> Cohort.open(cohort) -> sample(1024, 8192) on the GPU -> the same windows cut from the FASTA that the
  CPU decoder rebuilds (every block XXH3 and the FASTA XXH3 == the source file's, from the header) - byte for byte;
  the 1024 windows also == the CPU decoder's own windows; fetch at 1000 random coordinates == the FASTA.
Env: PANVRAM_COHORT (dir, required - else skipped), PANVRAM_DATASET (q4k), PANVRAM_REFERENCE, PANVRAM_DEVICE (cuda),
PANVRAM_SOURCES (dir with <assembly>.fa: the FASTA rebuilt is also compared with the source file byte for byte),
PANVRAM_SEED (20261004). Prints GATE lines for the log."""
import os
import random
import time

import panvram  # first, as a user would: it loads torch itself
import pytest
import torch
from fasta_cut import FastaIndex, transform

COHORT = os.environ.get("PANVRAM_COHORT")
pytestmark = pytest.mark.skipif(not COHORT, reason="PANVRAM_COHORT not set")


def test_gate():
    ds = os.environ.get("PANVRAM_DATASET", "q4k")
    dev = os.environ.get("PANVRAM_DEVICE", "cuda")
    seed = int(os.environ.get("PANVRAM_SEED", "20261004"))
    src = os.environ.get("PANVRAM_SOURCES")
    t0 = time.time()
    packed = os.environ.get("PANVRAM_PACKED_REF", "0") == "1"
    compact = os.environ.get("PANVRAM_COMPACT_BLOCKS", "0") == "1"
    c = panvram.Cohort.open(COHORT, device=dev, dataset=ds, reference=os.environ.get("PANVRAM_REFERENCE"), packed_reference=packed, compact_blocks=compact)
    t_open = time.time() - t0
    rb = c.resident_bytes()
    print(f"\nGATE form\tpacked_reference {int(packed)}\tcompact_blocks {int(compact)}")
    print(f"\nGATE open\t{len(c)} assemblies\tdataset {ds} (Q {c.block_size})\tdevice {c.device}\t{t_open:.1f} s\tresident {rb['total'] / 1e9:.3f} GB "
          f"(reference {sum(rb.get(k, 0) for k in ('ref', 'pw', 'pes', 'pel', 'peb')) / 1e9:.4f}, payload {rb['P'] / 1e9:.4f}, "
          f"block table {sum(rb.get(k, 0) for k in ('off', 'st', 'grp_off', 'grp_st', 'blen', 'bcode', 'exc_b', 'exc_st')) / 1e9:.4f}, model tables {rb['tab'] / 1e9:.4f})")
    g = torch.Generator(device=dev).manual_seed(seed)
    x, co = c.sample(1024, 8192, generator=g, return_coords=True)
    assert x.shape == (1024, 8192) and x.dtype == torch.uint8 and x.device.type == torch.device(dev).type
    xs, co = x.cpu().numpy(), co.cpu().tolist()
    # the CPU decoder's own windows for the same coordinates
    asm = torch.tensor([r[0] for r in co])
    start = torch.tensor([c.contigs(r[0])[r[1]][3] + r[2] for r in co])
    y = c._core.windows(asm, start, 8192, None, False, True, True, 0)       # CPU path, block XXH3 checked
    cpu_equal = bool(torch.equal(x.cpu(), y))
    # 1000 fetch coordinates, seeded
    rng = random.Random(seed)
    fetches = []
    for _ in range(1000):
        a = rng.randrange(len(c))
        ctgs = c.contigs(a)
        tot = sum(t[2] for t in ctgs)
        p = rng.randrange(tot)
        ci = next(i for i, t in enumerate(ctgs) if (p := p - t[2]) < 0)  # contig weighted by length
        n = ctgs[ci][2]
        L = min(n, 1 + int(2 ** (rng.random() * 20)))
        s = rng.randrange(n - L + 1)
        fetches.append((a, ci, s, L))
    win_bad = fetch_bad = 0
    touched = sorted({r[0] for r in co} | {f[0] for f in fetches})
    t1 = time.time()
    for a in touched:
        fa = c.fasta(a)                                                     # CPU full decode, all hashes checked
        if src:
            with open(os.path.join(src, c.names[a] + ".fa"), "rb") as f:
                assert f.read() == fa, f"{c.names[a]}: rebuilt FASTA differs from the source file"
        idx = FastaIndex(fa)
        names = [t[0] for t in c.contigs(a)]
        for k, r in enumerate(co):
            if r[0] == a and xs[k].tobytes() != idx.cut(names[r[1]], r[2], 8192):
                win_bad += 1
        for (fa_, ci, s, L) in fetches:
            if fa_ == a and c.fetch(a, ci, s, L).cpu().numpy().tobytes() != idx.cut(names[ci], s, L):
                fetch_bad += 1
        del fa, idx
    print(f"GATE sample\t1024 x 8192 on {c.device}\t{1024 - win_bad} == FASTA\t{win_bad} differ\tCPU decoder windows {'==' if cpu_equal else 'DIFFER'}")
    print(f"GATE fetch\t1000 coordinates (length 1 .. 2^20, contig by length)\t{1000 - fetch_bad} == FASTA\t{fetch_bad} differ")
    print(f"GATE fasta\t{len(touched)} assemblies rebuilt by the CPU decoder, XXH3 == source{' and == source files' if src else ''}\t{time.time() - t1:.0f} s")
    ok = win_bad == 0 and fetch_bad == 0 and cpu_equal
    print(f"GATE RESULT\t{'PASS' if ok else 'FAIL'}")
    assert ok
