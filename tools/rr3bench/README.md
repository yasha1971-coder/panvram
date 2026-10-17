# rr3bench - the refrel3 v1 in-process benchmark of M6 and mech558

`rr3bench` produced the refrel3 rows of M6 (`logs/m6-2026-10-07/`, e.g. 145 222 requests/s q4k and 73 690 q16k at
N = 558, one thread, 4 KiB windows) and the v1 rows of the 64 KiB order-effect counters (`logs/mech558-2026-10-09/`).
It reads only refrel3 v1 archives: the reference and every archive are opened with all checks of FORMAT.md section 5
(`open_v1`) before timing; a request decodes the blocks it touches (`block_v1`, block XXH3 not checked) and copies its
bases; every answer is compared by SHA-256 with the truth table after the timer.

Sources (unchanged copies):

| file | origin |
|---|---|
| `rr3bench.cpp`, `bench_common.h` | aceapex `research/agc_vs_refrel3` @ 2cc43cc |
| `refrel3v1.cpp` (v1 container), `refrel3.cpp`, `refrel3.h`, `refrel_format.h`, `refrel.cpp`, `refrel2.cpp`, `refrel_v2.h`, `refrel_io.h`, `rr_sha256.h` | aceapex `research/refrel` @ 5b6d5ce |
| `src/` of the aceapex library (`aceapex_api.cpp`, linked; `xxhash.h`) | public aceapex repository @ e902934 (same `src` tree as 5b6d5ce) |

`refrel2.cpp` / `refrel_v2.h` are the October 2-4 token-stream experiments that `refrel3.cpp` includes for shared
code (the "v2" in their names is the refrel token-stream variant, not a separate archive format); the decode path used
here is `open_v1` / `block_v1` of `refrel3v1.cpp`.

## Build

    git clone https://github.com/yasha1971-coder/aceapex && git -C aceapex checkout e902934
    bash tools/rr3bench/build.sh aceapex [out dir]

`build.sh` derives `refrel3v1_nomain.cpp` (the `main` of `refrel3v1.cpp` renamed by `sed`) and compiles with
`g++ -std=c++17 -O3 -march=native -funroll-loops ... -lzstd -lpthread libcrypto.so.3`.

## Binary identity and method check (ace-core, 2026-10-17)

| binary | SHA-256 |
|---|---|
| used for M6 and mech558 (recorded in `logs/mech558-2026-10-09/inputs_counters.sha256`) | f2f84fe9f396d3a67ede44d81eafb4959330baa7069584953c47d194a6e21673 |
| rebuilt from this directory with `build.sh` (g++ 11.4.0, Ubuntu 22.04.5, AMD EPYC 4344P) | f2f84fe9f396d3a67ede44d81eafb4959330baa7069584953c47d194a6e21673 |

The rebuilt binary is byte-identical to the one that produced the published logs.

Short re-run of the M6 N = 50 one-thread window rows with the M6 inputs (`windows_4096.tsv` abdc0460...,
`truth.tsv` 938f0359..., `names.tsv` 9df94ca1...; 3 warm-ups + 9 timed runs, `taskset -c 3`, silence gate), logs in
`method-check-n50/`:

| row | published median s (min - max) | re-run median s (min - max) | difference | answers |
|---|---|---|---:|---|
| rr3_q4k_win_t1 | 0.0498 (0.0486 - 0.0519) | 0.0517 (0.0515 - 0.0519) | +3.8 % | SHA-256 10 000 / 10 000 in all 12 runs |
| rr3_q16k_win_t1 | 0.1070 (0.1062 - 0.1112) | 0.1050 (0.1043 - 0.1057) | -1.9 % | SHA-256 10 000 / 10 000 in all 12 runs |

Both within 5 % of the published medians (published: `logs/m6-2026-10-07/raw/RESULTS_M6.md`, N = 50 rows).
