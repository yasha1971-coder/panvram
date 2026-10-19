# Time and memory, 558 assemblies

Host in `host.txt`: AMD EPYC 4344P, 8 cores and 16 threads, Ubuntu 22.04 with kernel 5.15.0-163, g++ 11.4.0.
All runs used 16 threads behind a silence gate (`silence_n558.txt`: one-minute load average below 0.5 and no process
above 5 % CPU before the run starts). CPU seconds are summed over the 16 threads; wall is the elapsed time of the
pass. Raw tool lines in `raw_lines.txt`.

## Exact search, 2 993 arbitrary 23-mers

| stage | full decode of every window | transfer plus residual |
|---|---:|---:|
| token decode, CPU s | 320 | 515 |
| copy execution, CPU s | 321 | 464 |
| window scan, CPU s | 6 588 | 658 |
| transfer and sort, CPU s | — | 514 |
| reference scan, wall | — | 0.9 s |
| **wall** | **476 s** | **156 s** |
| peak RSS | 5.9 GB | 6.9 GB |

The reference scan happens once for the whole cohort and yields 6 824 052 occurrences, 0.11 GB.

## CRISPR-style search, 120 guides, PAM NGG, up to 6 substitutions

| stage | full decode | transfer plus residual |
|---|---:|---:|
| window scan, CPU s | 37 764 | 1 660 |
| **wall** | **2 457 s** | **269 s** |
| peak RSS | 8.4 GB | 9.6 GB |

## The same decomposition with the per-guide origin split, 12 guides

145 014 299 sites, wall 106 s, of which the reference scan is 1.0 s, peak RSS 6.0 GB. Collecting the unique residual
sequences took a further 71 s, and checking all 118 078 of them against the reference 22 s.

## Why the ratios differ

Token decode and copy execution are needed by both paths and do not shrink. What shrinks is the window scan: 4.5 % of
the windows are residual, the rest are answered by the single reference scan. So the cheaper the per-window test, the
smaller the end-to-end gain: for one hash lookup per window the whole pass is 3.1 times faster, while for 120 guides
compared against each window it is 9.1 times faster and the matching step alone 22.7 times.

The variant that avoids executing copies and resolves residual bytes on demand (`parse2` in the tool) is byte-exact
but slower, 1 472 s wall at this scale, because self copies form chains up to the full block length of 16 384 and each
resolved byte then costs many hops. Numbers for it are in `raw_lines.txt`; it is kept as a measured negative result.
