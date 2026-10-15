# Order effect of 64 KiB windows at N = 558 (v1 q16k) - perf counters (2026-10-09)

Copy of aceapex refrel-v2 `research/refrel_v2/results/mech558/v1_*` (the v1 rows). The same 10 000 windows of 65 536
bases in request order (`random`) and sorted by assembly/contig/start (`sorted`), one thread (taskset -c 3), every answer
SHA-256 == truth (`*.log`). `perf stat -x,` user-space counters (`*.csv`: value, unit, event, running ns, % running);
`n1` = open + 1 timed pass, `n11` = open + 11 passes; per-window values are (n11 - n1) / (10 x 10 000).
Clock during the counted region: cycles(n11 - n1) / enabled time(n11 - n1) = 5.23-5.28 GHz over the six variants (g1-g3 x random/sorted)
(enabled = running / 0.8333).

## Per window (derived from the CSVs by the formula above)

| group | events | cycles random | cycles sorted | extra cycles | clock GHz random / sorted |
|---|---|---:|---:|---:|---|
| g1 | cycles, instructions, demand fills by source | 453358 | 341536 | 111822 | 5.28 / 5.24 |
| g2 | cycles, dTLB misses, table walks, page faults | 449992 | 344994 | 104998 | 5.26 / 5.25 |
| g3 | cycles, L1d, L2 misses, L2 prefetch misses, any fills from DRAM | 449541 | 346763 | 102779 | 5.26 / 5.23 |

Group g1, demand fills from local DRAM (`ls_dmnd_fills_from_sys.mem_io_local`) per window: random 2775, sorted 2485, extra 289.
Extra cycles per extra DRAM demand fill (g1): 111822 / 289 = 387 cycles = 73.5 ns at 5.26 GHz.
Group g3, any fills from local DRAM per window: random 3129, sorted 2896. Group g2, table walks per window (`ls_tablewalker.dside`): random 131, sorted 91.
