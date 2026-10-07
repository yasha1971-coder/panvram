# SPDX-License-Identifier: MIT
"""check_full_all.py <cohort dir> <manifest_v1.tsv> <reference> - every archive of the cohort (q4k and q16k), CPU path:
Cohort.check_full (full decode, FASTA rebuilt from the contig table, XXH3) == the header's == fasta_xxh3 of the manifest.
Opens one archive at a time (reference loaded once per dataset). Prints CF lines and RESULT."""
import csv, os, sys, time
import panvram

coh, man, ref = sys.argv[1:4]
want = {r["name"]: r["fasta_xxh3"] for r in csv.DictReader(open(man), delimiter="\t")}
ok = n = 0; t0 = time.time()
for ds in ("q4k", "q16k"):
    names = sorted(want)
    c = panvram.Cohort.open(coh, device="cpu", dataset=ds, reference=ref, assemblies=names)
    for a in names:
        t = time.time()
        try: eq, h = c.check_full(a); good = eq and h == want[a]
        except Exception as e: good, h = False, "ERROR " + str(e)
        ok += good; n += 1
        print(f"CF\t{ds}\t{a}\t{'ok' if good else 'FAIL'}\t{h}\t{want[a]}\t{time.time() - t:.1f} s", flush=True)
    del c
print(f"RESULT\t{ok}/{n}\t{'PASS' if ok == n == 2 * len(want) else 'FAIL'}\t{time.time() - t0:.0f} s", flush=True)
