# SPDX-License-Identifier: MIT
"""cohort_region.py <cohort dir> <requests.tsv> <truth.tsv> [dataset] - PROTOCOL_COHORT_REGION (GPU row): the task =
every request of requests.tsv (columns id, sample = manifest name, contig, start 0-based, length) answered by
Cohort.fetch on the GPU; 3 warm-ups + 9 timed runs of the whole task (one device sync per run); after every run (outside
timing) SHA-256 of every answer == truth.tsv (id, SHA-256). Resident form from PANVRAM_PACKED_REF / PANVRAM_COMPACT_BLOCKS.
Prints CR lines: open, every run, RESULT."""
import hashlib
import os
import statistics
import sys
import time

import panvram
import torch


def main():
    path, reqf, truthf = sys.argv[1:4]
    ds = sys.argv[4] if len(sys.argv) > 4 else "q4k"
    packed = os.environ.get("PANVRAM_PACKED_REF", "0") == "1"
    compact = os.environ.get("PANVRAM_COMPACT_BLOCKS", "0") == "1"
    reqs = [l.rstrip("\n").split("\t") for l in open(reqf) if l.strip() and not l.startswith("id\t")]
    truth = dict(l.rstrip("\n").split("\t")[:2] for l in open(truthf) if l.strip() and not l.startswith("id\t"))
    t = time.time()
    c = panvram.Cohort.open(path, device="cuda", dataset=ds, packed_reference=packed, compact_blocks=compact)
    torch.cuda.synchronize()
    rb = c.resident_bytes()
    print(f"CR open\t{len(c)} assemblies\t{ds}\tpacked_reference {int(packed)} compact_blocks {int(compact)}\t{time.time() - t:.1f} s\t"
          f"resident {rb['total'] / 1e9:.3f} GB\t{len(reqs)} requests, {sum(int(r[4]) for r in reqs)} bases", flush=True)
    times, bad_runs = [], 0
    for run in range(12):
        torch.cuda.synchronize()
        t0 = time.perf_counter()
        outs = [c.fetch(s, ctg, int(st), int(n)) for (_, s, ctg, st, n) in reqs]
        torch.cuda.synchronize()
        dt = time.perf_counter() - t0
        ok = sum(hashlib.sha256(o.cpu().numpy().tobytes()).hexdigest() == truth[r[0]] for o, r in zip(outs, reqs))
        del outs
        bad_runs += ok != len(reqs)
        if run >= 3:
            times.append(dt)
        print(f"CR run\t{run}\t{'warmup' if run < 3 else 'timed'}\t{dt:.4f} s\tsha256 {ok}/{len(reqs)}", flush=True)
    m = statistics.median(times)
    print(f"CR RESULT\t{'PASS' if bad_runs == 0 else 'FAIL'}\tmedian {m:.4f} s per task\tmin {min(times):.4f}\tmax {max(times):.4f}\t"
          f"raw {','.join(f'{x:.4f}' for x in times)}\t{len(reqs) / m:.1f} answers/s")
    sys.exit(0 if bad_runs == 0 else 1)


if __name__ == "__main__":
    main()
