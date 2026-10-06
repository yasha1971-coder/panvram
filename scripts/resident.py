# SPDX-License-Identifier: MIT
"""resident.py <cohort dir> <dataset> - bytes held by Cohort.open on the CPU path for the four resident forms
(packed_reference x compact_blocks), by component; prints RES lines. Opening checks every archive (FORMAT.md 5)."""
import sys
import time

import panvram

path, ds = sys.argv[1], sys.argv[2]
for packed, compact in ((False, False), (True, False), (False, True), (True, True)):
    t = time.time()
    c = panvram.Cohort.open(path, device="cpu", dataset=ds, packed_reference=packed, compact_blocks=compact)
    rb = c.resident_bytes()
    ref = sum(rb.get(k, 0) for k in ("ref", "pw", "pes", "pel", "peb"))
    tab = sum(rb.get(k, 0) for k in ("off", "st", "grp_off", "grp_st", "blen", "bcode", "exc_b", "exc_st"))
    print(f"RES\t{ds}\t{len(c)} assemblies\tpacked_reference {int(packed)}\tcompact_blocks {int(compact)}\ttotal {rb['total']}\treference {ref}\t"
          f"payload {rb['P']}\tblock table {tab}\tmodel tables {rb['tab']}\tcase runs {rb.get('low_s', 0) + rb.get('low_l', 0)}\topen {time.time() - t:.1f} s", flush=True)
    del c
