# Window-cost verdict of refrel3 v1 (q16k) at N = 558 (2026-10-17)

Model B of the hw-apex window-law job (protocol v1.1) applied to the 558 v1 archives of the manifest, measured by
`tools/window_cost/window_cost` of this repository: a one-file harness built on the public reader only
(`csrc/pv_v1.h`: `open_v1` with every check of FORMAT.md section 5, `decode_block`). Run by the author on ace-core
(AMD EPYC 4344P, Ubuntu 22.04.5, g++ 11.4.0), one thread pinned to CPU 3, behind a silence gate (load < 0.5, no
process above 5 % CPU; `silence_n558_v1.txt`).

Rules (unchanged from the N = 4 verdict): 10 000 requests per W from the frozen hw-apex generator (seed 20261003,
`req558.tsv.gz`, produced by `tools/window_cost/window_requests.py`); D_Q = canonical bases / median of 9 full
sequential decodes after 3 warm-ups, output buffer allocated outside the timer, XXH3 of every decoded assembly
checked after its interval; per request one timed fetch (contig lookup by name, the touched blocks decoded, copy),
SHA-256 of the answer checked after the timer against `truth558.tsv.gz`; p50 nearest rank; c0 = p50(W = 1) - Q / D_Q;
predicted t(W) = c0 + (W + Q - 1) / D_Q; error = 100 (p50 - pred) / pred; PASS requires c0 >= 0 and |error| <= 20 %
at W = 1024, 8192 and 65536. Block hashes were not checked inside the timer (`--verify-blocks` off); the truth table
was produced with every block hash and the header XXH3 checked.

## Result (`n558_v1.log`, 2026-10-17 13:20:37 - 15:03:22 UTC, 1 h 42 min, peak RSS 15.2 GB)

| quantity | value |
|---|---|
| assemblies / canonical bases | 558 / 1 680 437 743 110 |
| Q | 16 384 |
| D_Q | 3.543753 Gbases/s (median 474 197 192 697 ns per full decode; min 474 108 461 248, max 474 243 982 581) |
| p50 at W = 1 (calibration) | 10.559 us |
| c0 | 5.936 us |
| W = 1 024 | p50 10.199 us, predicted 10.848 us, error -5.98 %, PASS |
| W = 8 192 | p50 13.245 us, predicted 12.870 us, error +2.91 %, PASS |
| W = 65 536 | p50 43.282 us, predicted 29.052 us, error +48.98 %, FAIL |
| verdict | FAIL (64 KiB) |
| requests | 10 000 per W, every answer SHA-256 == truth, bad 0 |

Every timed request is in `raw_n558_v1.tsv.gz` (W, request id, nanoseconds; 40 000 rows).

## Inputs and binary (`n558_v1_inputs.sha256`)

| file | SHA-256 |
|---|---|
| req558.tsv (requests; here gzip-compressed) | 4f89a4d4f70c5b95eb3afe3a6140e9abec306ab9b7bd4c95c520224c0de2ce43 |
| truth558.tsv (truth table; here gzip-compressed) | 0193fa4c0a1ced3ef511d5d665f3be86e43d324925da9962518a6b178ed5a3eb |
| names558.txt (assembly order) | 22939508bc9cdaa1b80d03daa3ddafc8ae01a8adbcbeab71b5cb528dd483edc5 |
| window_cost (binary, `tools/window_cost/build.sh`) | caca51c322ee5a913bd6abd0b3ebf0cd9816d4c2f6625d42ca108012c328ce36 |

Archives: the q16k set of `MANIFEST.tsv` (`SHA256SUMS.q16k`), files `<name>.q16k.rr3`; reference: the T2T FASTA of
DATASET.md (3 156 259 565 bytes, SHA-256 51ab7d3c2a375ecff53ee2a1be03e58d57cc1601c948a860c03a78d724f26e4c), checked
by `open_v1` against every archive header. Build: `bash tools/window_cost/build.sh` (g++ -std=c++17 -O3 -march=native
-funroll-loops, libzstd); `truth_check.txt`: the truth table of this harness is byte-identical to the one produced
earlier by the author's private harness.

## Method check (`method-check-n50/`)

The same protocol on the first 50 assemblies, same requests and truth, same CPU, the private harness that produced the
2026-10-09 result (`old_n50.log`) and `window_cost` (`new_n50.log`), one after the other behind the gate:

| harness | D_Q, Gbases/s | c0, us | 1 KiB | 8 KiB | 64 KiB | verdict | wall |
|---|---:|---:|---|---|---|---|---|
| private (2026-10-09 harness) | 2.397214 | 1.611 | -8.65 % | -13.52 % | +3.14 % | PASS | 13 min 34 s |
| window_cost (this repository) | 2.397611 | 1.613 | -6.39 % | -10.56 % | +5.71 % | PASS | 13 min 31 s |

D_Q and c0 agree within 0.1 %; the window errors have the same sign at every W (the spread of 2-3 points between
two runs is the run-to-run variation of a p50 over 10 000 requests). The truth tables of both harnesses are identical
(`method-check-n50/truth_check.txt`).

For comparison, the private harness on the same 558 archives on 2026-10-09 gave D_Q 3.400 Gbases/s, c0 6.442 us and
+47.88 % at 64 KiB (FAIL); that run is superseded by this one and is not published.
