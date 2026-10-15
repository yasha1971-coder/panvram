# Order effect of 64 KiB windows at N = 558 (v1 q16k) - perf counters (2026-10-09)

Copy of aceapex refrel-v2 `research/refrel_v2/results/mech558/v1_*` (the v1 rows). The same 10 000 windows of 65 536
bases in request order (`random`) and sorted by assembly/contig/start (`sorted`), one thread (taskset -c 3), every answer
SHA-256 == truth (`*.log`). `perf stat -x,` user-space counters (`*.csv`: value, unit, event, running ns, % running);
`n1` = open + 1 timed pass, `n11` = open + 11 passes; per-window values are (n11 - n1) / (10 x 10 000).
Clock during the counted region: cycles(n11 - n1) / enabled time(n11 - n1) = 5.23-5.28 GHz over the six variants (g1-g3 x random/sorted)
(enabled = running / 0.8333).
