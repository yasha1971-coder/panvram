# SPDX-License-Identifier: MIT
"""samtools_truth.py <work dir> <cohort dir> <sources.tsv> [n per assembly] [only assembly]
Independent truth for panvram fetch: samtools faidx on the SOURCE FASTA of every assembly, compared byte for byte with
panvram fetch (CPU path) for q4k and q16k.

Coordinates - hw-apex protocol v1.1 (review/axis3/PROTOCOL_AXIS3.md, SHA-256 f0bb442a...), "Canonical sequence bytes and
coordinates": a request is (assembly_id, contig_id, start0, end0), zero-based half-open [start0, end0), valid iff
0 <= start0 < end0 <= contig_length; the answer is the uppercase sequence bytes (no header, no line separators).
Translation at each library boundary:
  samtools faidx region string (one-based inclusive)  contig:{start0 + 1}-{end0}
  panvram Cohort.fetch (zero-based start, length)      fetch(assembly, contig, start0, end0 - start0)
Both answers are uppercased; SHA-256 of the samtools answer is the truth.

Requests: contig lengths from the .fai that samtools builds on the source (not from the archives);
random.Random(20261009 + i) for the i-th assembly in sources.tsv order; contig weighted by length; length
L = min(contig_length, 1 + floor(2 ** (u * 20))) (1 .. 1 048 576); start0 uniform over every valid start.
sources.tsv: assembly_id, source FASTA path (plain or BGZF), FASTA XXH3 expected (MANIFEST.tsv) - checked by the caller.
Output in <work dir>: requests.tsv, truth.tsv (id, sha256 of samtools answer), compare.tsv, and a RESULT line."""
import hashlib
import os
import random
import subprocess
import sys

import panvram

SAMTOOLS = os.environ.get("SAMTOOLS", "samtools")
W, COH, SRC = sys.argv[1], sys.argv[2], sys.argv[3]
N = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
ONLY = sys.argv[5] if len(sys.argv) > 5 else None
os.makedirs(W, exist_ok=True)
sha = lambda b: hashlib.sha256(b).hexdigest()


def fai(path):
    if not os.path.exists(path + ".fai"):
        subprocess.run([SAMTOOLS, "faidx", path], check=True)
    return [(l.split("\t")[0], int(l.split("\t")[1])) for l in open(path + ".fai")]


def requests(i, asm, path):
    ctgs = fai(path); rng = random.Random(20261009 + i); tot = sum(n for _, n in ctgs); out = []
    for r in range(N):
        p = rng.randrange(tot); k = 0
        while p >= ctgs[k][1]: p -= ctgs[k][1]; k += 1
        name, n = ctgs[k]
        L = min(n, 1 + int(2 ** (rng.random() * 20)))
        s0 = rng.randrange(n - L + 1); e0 = s0 + L
        assert 0 <= s0 < e0 <= n
        out.append((f"{asm}#{r}", asm, name, s0, e0))
    return out


def samtools(path, reqs):
    rf = os.path.join(W, "regions.tmp")
    with open(rf, "w") as f:
        for _, _, c, s0, e0 in reqs: f.write(f"{c}:{s0 + 1}-{e0}\n")
    out = subprocess.run([SAMTOOLS, "faidx", "-n", "0", "-r", rf, path], check=True, capture_output=True).stdout
    recs, cur = [], None
    for line in out.split(b"\n"):
        if line.startswith(b">"): cur = []; recs.append(cur)
        elif cur is not None: cur.append(line)
    assert len(recs) == len(reqs), (len(recs), len(reqs))
    return [b"".join(r).upper() for r in recs]


def main():
    src = [l.rstrip("\n").split("\t") for l in open(SRC) if l.strip() and not l.startswith("#")]
    allreq, truth = [], {}
    for i, (asm, path, _x3) in enumerate(src):
        if ONLY and asm != ONLY: continue
        reqs = requests(i, asm, path); ans = samtools(path, reqs)
        for q, a in zip(reqs, ans):
            assert len(a) == q[4] - q[3], q
            truth[q[0]] = sha(a)
        allreq += reqs
    with open(os.path.join(W, "requests.tsv"), "w") as f:
        f.write("id\tassembly_id\tcontig_id\tstart0\tend0\n")
        for q in allreq: f.write("%s\t%s\t%s\t%d\t%d\n" % q)
    with open(os.path.join(W, "truth.tsv"), "w") as f:
        f.write("id\tsha256\n")
        for q in allreq: f.write(f"{q[0]}\t{truth[q[0]]}\n")
    rsha = sha(open(os.path.join(W, "requests.tsv"), "rb").read())
    res = {}
    with open(os.path.join(W, "compare.tsv"), "w") as f:
        f.write("dataset\tid\tequal\n")
        for ds in ("q4k", "q16k"):
            c = panvram.Cohort.open(COH, device="cpu", dataset=ds, assemblies=sorted({q[1] for q in allreq}))
            ok = 0
            for q in allreq:
                b = bytes(c.fetch(q[1], q[2], q[3], q[4] - q[3]).cpu().numpy().tobytes()).upper()
                eq = sha(b) == truth[q[0]]; ok += eq
                f.write(f"{ds}\t{q[0]}\t{int(eq)}\n")
            res[ds] = ok
            del c
    n = len(allreq)
    print(f"requests {n} ({len({q[1] for q in allreq})} assemblies) requests.tsv sha256 {rsha}")
    print(f"truth.tsv sha256 {sha(open(os.path.join(W, 'truth.tsv'), 'rb').read())}")
    for ds in res: print(f"{ds}: panvram fetch == samtools {res[ds]}/{n}")
    print("RESULT", "PASS" if all(v == n for v in res.values()) else "FAIL")


if __name__ == "__main__":
    main()
