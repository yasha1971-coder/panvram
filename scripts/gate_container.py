# SPDX-License-Identifier: MIT
"""gate_container.py <cohort dir> <source_sha256.tsv> - the 8-assembly gate run inside the CPU container (--network none):
for q4k and q16k, every assembly: check_full (full decode, FASTA XXH3 rebuilt == header == MANIFEST) and the CPU full
decode (all block XXH3 checked) == the source file byte for byte (SHA-256 of the uncompressed source, computed outside).
sources tsv: assembly, SHA-256 of the uncompressed source FASTA, MANIFEST fasta_xxh3."""
import hashlib, sys, time
import panvram

coh, tsv = sys.argv[1], sys.argv[2]
src = [l.rstrip("\n").split("\t") for l in open(tsv) if l.strip()]
names = [s[0] for s in src]
ok = n = 0; t0 = time.time()
for ds in ("q4k", "q16k"):
    c = panvram.Cohort.open(coh, device="cpu", dataset=ds, assemblies=names)
    for a, sha, x3 in src:
        t = time.time(); eq, h = c.check_full(a); fa = c.fasta(a)
        good = eq and h == x3 and hashlib.sha256(fa).hexdigest() == sha
        ok += good; n += 1
        print(f"GATE {ds} {a}\tcheck_full {'==' if eq else '!='} header ({h}, MANIFEST {x3})\tFASTA {len(fa)} B sha256 {'== source' if hashlib.sha256(fa).hexdigest() == sha else 'DIFFERS'}\t{time.time() - t:.1f} s", flush=True)
        del fa
    del c
print(f"GATE RESULT {ok}/{n} {'PASS' if ok == n else 'FAIL'}\t{time.time() - t0:.0f} s")
sys.exit(0 if ok == n else 1)
