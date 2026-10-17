#!/usr/bin/env python3
"""window_requests.py <hw-apex checkout> <contigs.tsv> <names.txt> <out.tsv>

Window requests from the frozen generator of the hw-apex window-law job (protocol v1.1):
tools/axis3_window_engine.sample_coordinates of the given hw-apex-bench checkout, imported unchanged; 10 000 requests
per W, seed 20261003, W = 1 (c0 probe) and the verdict windows 1024 / 8192 / 65536; canonical domain = upper-case bases
of the contig record, half-open coordinates. contigs.tsv: "<assembly> <contig name> <length>" per line, archive record
order within an assembly; names.txt fixes the assembly order. Output TSV: W id assembly contig start0 end0."""
import sys

hwapex, contigs_path, names_path, out = sys.argv[1:5]
sys.path.insert(0, hwapex)
from tools.axis3_window_engine import sample_coordinates  # noqa: E402

order = [l.strip() for l in open(names_path) if l.strip()]
ctg = {}
for l in open(contigs_path):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 3:
        ctg.setdefault(f[0], []).append((f[1], int(f[2])))
contigs = [(a, n, L) for a in order for (n, L) in ctg[a]]
with open(out, "w") as o:
    o.write("W\tid\tassembly\tcontig\tstart0\tend0\n")
    for w in (1, 1024, 8192, 65536):
        for (i, a, c, s, e) in sample_coordinates(contigs, w, count=10000, seed=20261003):
            o.write(f"{w}\t{i}\t{a}\t{c}\t{s}\t{e}\n")
print(f"{len(order)} assemblies, {len(contigs)} contigs, 4 x 10000 requests -> {out}")
