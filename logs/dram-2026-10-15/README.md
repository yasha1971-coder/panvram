# DRAM load latency on ace-core (2026-10-15)

Pointer-chasing benchmark `dram_lat.c` (gcc -O2 -std=c11; random cyclic permutation of 64-B nodes; seed 20261012;
3 warm-ups + 9 measurements per working set, median reported) run by the author on ace-core (AMD EPYC 4344P, CPU 0;
`run/lscpu.txt`, `run/uname.txt`, silence state `run/gate.txt`, processes before/after). Results:
`dram_results_2026-10-15.tsv` (same table as `run/results.tsv`): working sets 32 KiB .. 4 GiB; 103.3 ns per load at
256 MiB rising to 115.8 ns at 4 GiB (4 KiB pages). `run/STATUS` = PASS.

- `SHA256SUMS` is the package's own list. Its entry for `run.sh` does not match the `run.sh` shipped next to it
  (eac193e9... listed, b805687d... present); the run's `run/command.txt` and `run/source_sha256.txt` record what ran.
- The binary `run/dram_lat` is not copied (its SHA-256 is in `PROVENANCE.tsv`); local home paths replaced by `~` in
  2 files (`PROVENANCE.tsv` lists original and copy SHA-256).
