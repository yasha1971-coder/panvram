"""bench.py <cohort dir> [dataset] [device] [quick] - sample() throughput on a resident cohort: windows/s and output GB/s per
(W, n, tokens), each the median of 5 timed runs of 10 back-to-back calls after a warm-up (one device sync per run, the
decode status checked after); fetch latency (one call, sync) by length. Prints BENCH lines."""
import statistics
import sys
import time

import torch

import panvram


def sync(dev):
    if dev.type == "cuda":
        torch.cuda.synchronize(dev)


def main():
    path = sys.argv[1]
    ds = sys.argv[2] if len(sys.argv) > 2 else "q4k"
    dev = torch.device(sys.argv[3] if len(sys.argv) > 3 else "cuda")
    quick = len(sys.argv) > 4
    t0 = time.time()
    c = panvram.Cohort.open(path, device=dev, dataset=ds)
    sync(c.device)
    rb = c.resident_bytes()
    name = torch.cuda.get_device_name(c.device) if c.device.type == "cuda" else "cpu"
    print(f"BENCH open\t{name}\t{len(c)} assemblies\t{ds} Q {c.block_size}\t{time.time() - t0:.1f} s\tresident {rb['total'] / 1e9:.3f} GB")
    g = torch.Generator(device=c.device).manual_seed(1)
    for W in ((1024, 8192) if quick else (1024, 8192, 65536)):
        for n in ((64, 256) if quick else (1024, 16384, 65536)):
            if n * W > (1 << 30):
                continue
            for tok in (False, True):
                c.sample(n, W, generator=g, tokens=tok)
                runs = []
                for _ in range(5):
                    sync(c.device)
                    t = time.perf_counter()
                    for _ in range(10):
                        c.sample(n, W, generator=g, tokens=tok, check=False)
                    sync(c.device)
                    runs.append(time.perf_counter() - t)
                t = statistics.median(runs) / 10
                c.sample(n, W, generator=g, tokens=tok)                     # status checked
                print(f"BENCH sample\tW {W}\tn {n}\ttokens {int(tok)}\t{n / t:,.0f} windows/s\t{n * W / t / 1e9:.2f} GB/s out\t{t * 1e3:.2f} ms per call")
    for L in ((100, 10_000) if quick else (100, 10_000, 1_000_000, 10_000_000)):
        ts = []
        for k in range(5):
            name0, _, n0, _, _ = c.contigs(k % len(c))[0]
            sync(c.device)
            t = time.perf_counter()
            y = c.fetch(k % len(c), name0, 1000 * k, min(L, n0 - 1000 * k))
            sync(c.device)
            ts.append(time.perf_counter() - t)
        print(f"BENCH fetch\tlength {L}\tmedian {statistics.median(ts) * 1e3:.2f} ms (one call, sync)")
    print("BENCH done")


if __name__ == "__main__":
    main()
