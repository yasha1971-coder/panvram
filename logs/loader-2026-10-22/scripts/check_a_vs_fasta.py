#!/usr/bin/env python3
"""check_a_vs_fasta.py - every window of the A request sets decoded by panvram must equal the same bases read from the
source FASTA (bgzip of the files hprc50 was built from), byte for byte, case kept."""
import json

import numpy as np
import torch
from pyfaidx import Fasta

import panvram

import os
COHORT = os.environ.get("COHORT", "cohort5")
BGZ = os.environ.get("BGZ", "bgz")

PAN = [l.strip() for l in open("requests/names5_panvram.txt") if l.strip()]
FAN = [l.strip() for l in open("requests/names5_fa.txt") if l.strip()]
c = panvram.Cohort.open(COHORT, device="cpu", dataset="q16k", threads=16)
assert list(c.names) == PAN, (list(c.names), PAN)
fa = [Fasta(f"{BGZ}/{n}.fa.gz", as_raw=True, sequence_always_upper=False, rebuild=False) for n in FAN]
ctg = [[l.split("\t")[0] for l in open(f"{BGZ}/{n}.fa.gz.fai")] for n in FAN]

for W in (8192, 32768):
    co = json.load(open(f"requests/A5.coords.W{W}.json"))
    asm = torch.tensor([r[0] for r in co], dtype=torch.int64)
    boff = torch.tensor([c._core.contigs(r[0])[r[1]][3] for r in co], dtype=torch.int64)
    start = boff + torch.tensor([r[2] for r in co], dtype=torch.int64)
    x = c.windows(asm, start, W).numpy()
    bad = 0
    for i, (a, cc, s) in enumerate(co):
        ref = np.frombuffer(str(fa[a][ctg[a][cc]][s:s + W]).encode(), dtype=np.uint8)
        if not np.array_equal(x[i], ref):
            bad += 1
            if bad == 1:
                d = int(np.flatnonzero(x[i] != ref)[0])
                print(f"  first difference: window {i} {PAN[a]} {ctg[a][cc]} start {s} offset {d}: panvram {x[i][d]} fasta {ref[d]}")
    print(f"CHECK\tW {W}\twindows {len(co)}\tidentical {len(co)-bad}\tdiffer {bad}", flush=True)
