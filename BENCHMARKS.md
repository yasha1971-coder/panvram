# BENCHMARKS - panvram (one structure per axis; numbers only from the named evidence)

Each axis: question - protocol - hardware - versions and flags - table - evidence; "Where we lose" separately. CPU and GPU
numbers come from different hardware and are never put in one column. Measurements of the format and of AGC are
kept in aceapex `BENCHMARKS.md` (branch `refrel`); this file lists the panvram paths and points there.

## 1. Correctness gates

| gate | scope | result | evidence |
|---|---|---|---|
| synthetic tests (Q 1024 / 4096 / 16384, all four resident forms) | CPU path | 17 passed, 1 skipped (no CUDA) | `tests/test_synth.py`; `logs/tests-cpu-2026-10-04.log` (default form) |
| public gate, 8 HPRC + T2T, q4k and q16k | CPU path, default form | 1024 / 1024 windows, 1000 / 1000 fetches == FASTA | `logs/gate-cpu-2026-10-04.log` |
| same, all four resident forms | CPU path | PASS in all 8 runs | `logs/gate-cpu-compact-2026-10-05.log` |
| public gate, 4 HPRC + T2T, q4k and q16k, archives == MANIFEST.tsv | GPU (Colab) | PASS | Colab log `MyDrive/aceapex_logs/panvram_<date>.txt` (not in this repository) |

- Hardware (CPU rows): ace-core, AMD EPYC 4344P 8C/16T, 125 GB, Ubuntu 22.04; torch 2.14.1+cpu (requirements.lock).

## 2. Resident bytes, 558 HPRC + T2T (CPU path pools = what goes to the device)

| form | q4k | q16k |
|---|---:|---:|
| default (byte reference, 12 B per block) | 16.64 GB | 11.56 GB |
| packed_reference | 14.30 GB | 9.23 GB |
| compact_blocks | 13.97 GB | 10.90 GB |
| both | 11.63 GB | 8.56 GB |

Evidence: `logs/resident-558-2026-10-06.log` (`scripts/resident.py`). GPU VRAM of the forms: not measured yet (C1 c).

## 3. CPU random access of the format (refrel3 v1) against AGC

See aceapex `BENCHMARKS.md` sections 2-3 (N = 50 measured; N = 558 under PROTOCOL_AGC558).

## 4. GPU throughput of the resident cohort

Not measured yet (`run_colab_cohort558.sh`, C1 b/c prepared, not run).

## Where we lose

- **Size against AGC** (aceapex BENCHMARKS section 1): 13.04 - 16.22 MB per assembly without block hashes against 5.40 MB
  for AGC with T2T at N = 558 (2.4 - 3.0 x); at N = 50 the gap was 1.5 - 1.9 x (8.65 against 13.14 - 16.33 MB): with a
  larger cohort AGC wins by more (two measurements; the cause is not stated). Plus the decoded T2T resident (3.117 G
  bases; 0.779 GB with `packed_reference`).
- **GPU:** CUDA only, sm_80 or newer; compact forms not yet run on a GPU.
