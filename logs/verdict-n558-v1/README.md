# Window-law verdict of refrel3 v1 (q16k) at N = 558 (2026-10-09)

Model B of the hw-apex window-law job (protocol v1.1), run by the author on ace-core (AMD EPYC 4344P, one thread,
`taskset -c 3`) with the in-process C++ harness `rv2 law` of the refrel v2 prototype, not the Python job used for the
N = 4 verdict. Rules are the same: 10 000 requests per W from the frozen hw-apex generator (seed 20261003);
D_Q = canonical bases / median of 9 full sequential decodes after 3 warm-ups (XXH3 of every decoded assembly checked);
c0 = p50(W = 1) - Q / D_Q; predicted t(W) = c0 + (W + Q - 1) / D_Q; error = 100 (p50 - pred) / pred, p50 nearest rank;
PASS if c0 >= 0 and |error| <= 20 % at W = 1024, 8192, 65536. Every request SHA-256 == truth (bad 0). Block-hash
checks of v1 were off in this run (the rv2 comparison condition).

Result (`n558_v1.log`): 558 assemblies, 1 680 437 743 110 canonical bases, Q = 16 384, D_Q 3.400120 Gbases/s,
c0 6.442 us; p50 10.710 us / -7.37 % at 1 KiB (PASS), 13.716 us / +0.34 % at 8 KiB (PASS), 45.155 us / +47.88 % at
64 KiB (predicted 30.535 us, FAIL). Verdict FAIL. Peak RSS 16.6 GB, 1 h 48 min wall.

Files:
- `n558_v1.log` - the run log (D_Q passes, per-W rows, `/usr/bin/time -v`); home paths replaced by `~` in the recorded
  command line (original SHA-256 7d50969da7773b757c87214efe284e24292309a81010a882b9d12892235316e8).
- `raw_n558_v1.tsv.gz` - every timed request: W, request id, nanoseconds (40 000 rows).
- `silence_n558_v1.txt` - host state before the run (load, free disk and memory, top processes).
- `n558_inputs.sha256` - SHA-256 of the request file, truth table, assembly list and the `rv2` binary used.

Harness source at commit 8baee1f of the private research branch (not published with this log): `rv2.cpp`
SHA-256 eaafe7a059a3eb84eb58866619c26357439fc2f86f2e3f73423e2fd9197e46a1, `law_requests.py`
5b34f6dca59c2d52a6eac094dde6a5c2bcc2cd46f5e5233e5b5b0f0a48213f2b, `run_law.sh`
302a38bc114d86856ba2e2a3639703346ed00c51c61320b995825f8d5f080106. The request and truth files (hashes above) are
derived from the published v1 archives and the T2T reference and are not copied here.
