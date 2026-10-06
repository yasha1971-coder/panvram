# Benchmarks - how every number in README was measured

Only runs with committed logs (or a named log outside this repository) are listed. CPU and GPU numbers come from
different hardware and are never put in one column.

## CPU: refrel3 against AGC 3.2.4, HPRC N = 50 (in process)

- Where: aceapex branch `refrel`, `research/agc_vs_refrel3/` (commit f526bc7, local until published): `README.md`,
  `SUMMARY_runs.md`, `logs/`, tools (`truth.cpp`, `agcbench.cpp`, `rr3bench.cpp`, `agc_cli_batch.py`,
  `make_requests.py`, `run_all.sh`, `summarize.py`), `binaries.sha256`.
- Hardware: ace-core, AMD EPYC 4344P (8 cores / 16 threads), 125 GB RAM, Ubuntu 22.04; silence gate before the runs
  (load < 0.5; `logs/silence_gate*.txt`).
- Data: the first 50 rows of MANIFEST.tsv; truth from the 50 source `.fa.gz` (SHA-256 == the HPRC index).
- refrel3: aceapex research/refrel @ 5b6d5ce (`refrel3v1.cpp`), g++ 11.4 -O3 -march=native -funroll-loops; the v1
  archives of MANIFEST.tsv.
- AGC: 3.2.4; in process via the C API of `libagc.a` built from source tag v3.2.4 (commit e67e3fc; g++ -O3
  -march=native -std=c++20), one handle per thread, `prefetching = 1`; CLI rows with the release binary (build
  20260722.1), all requests in one `agc getctg` process. Archives: `t2t_50.agc` (one `create`, T2T + 50) and
  `noref_50.agc`.
- Protocol: 10 000 windows of 8 192 bases (seed 20261003), 1 000 regions (seed 20261005), whole samples (first 4; all
  50 with one run); 3 warm-ups + 9 timed runs; threads 1 (pinned) and 16; every answer SHA-256 == truth.

## CPU: one-thread full decode D_Q (window law), HPRC N = 4

- Where: hw-apex-bench PR #63, `review/axis3/results/ace-core-2026-10-05-407088c-dq/` (commit 1351dcf): `RUN.md`,
  `SUMMARY.md`, logs, `dq_refrel3.cpp`, `dq_others.c`.
- Same host and gate procedure; refrel3 q4k / q16k, BGZF (htslib 1.24 + libdeflate 1.19, level 6), zstd seekable
  (1.5.7, level 3, 16 KiB frames), lz4 frame (1.10.0, 4 MiB independent blocks); OpenZL v0.3.0 failed at build
  (input limit ~1.5 GiB).

## CPU path of panvram: gate and resident bytes

- `logs/tests-cpu-2026-10-04.log`, `logs/gate-cpu-2026-10-04.log`, `logs/gate-cpu-compact-2026-10-05.log`,
  `logs/resident-558-2026-10-06.log`; ace-core, torch 2.14.1+cpu (requirements.lock); scripts `tests/test_gate.py`,
  `scripts/resident.py`.

## GPU

- Public gate on 4 HPRC assemblies, q4k and q16k, archives == MANIFEST.tsv: `run_colab_panvram.sh` on Colab; log
  `MyDrive/aceapex_logs/panvram_<date>.txt` (not in this repository): PASS. No throughput number of the resident cohort
  is published yet.
- refrel3 research kernel (the queue kernel's origin), 4 assemblies, RTX PRO 6000 Blackwell: aceapex
  `research/refrel/logs/colab-blackwell-2026-10-03-excerpt.txt`.
- The compact resident forms (`packed_reference`, `compact_blocks`) are compiled for sm_80 but have not run on a GPU.
