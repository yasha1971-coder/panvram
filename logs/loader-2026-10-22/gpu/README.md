# GPU path of A: measured distributions (2026-10-25)

NVIDIA A100-SXM4-40GB on Colab, 1 410 MHz with no throttle reason in every run of `runs-2026-10-25/`. Those runs use a
cohort of 5 assemblies (`colab_prepare.prepare_cohort(5)`, q16k, Q = 16 384) and report **distributions over 10
independent draws of windows** (`Cohort.sample`, seeds 20261003 .. 20261012); each draw is timed several times and its
median taken. The CPU rows of the parent README were measured on ace-core and are not changed by anything here.

## Headline: batch 16 384 (`runs-2026-10-25/b2head.log`, 3 timings per draw)

| W | time median | p90 | max | min | windows/s median | windows/s worst | GB/s median | GB/s worst |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 192 | 11.938 ms | 13.681 | 14.766 | 10.541 | 1 379 788 | 1 109 551 | 11.303 | 9.089 |
| 32 768 | 17.427 ms | 18.781 | 19.332 | 15.644 | 940 250 | 847 528 | 30.810 | 27.772 |

## Batch 1 024 (`runs-2026-10-25/b1dist.log`, 5 timings per draw)

| W | time median | p90 | max | min | max/min |
|---:|---:|---:|---:|---:|---:|
| 8 192 | 3.125 ms | 4.175 | 6.899 | 1.270 | 5.43x |
| 32 768 | 3.216 ms | 3.642 | 4.249 | 2.202 | 1.93x |

A batch of 1 024 does not fill the card and waits for its slowest block; why, with the per-block record and the
refuted hypotheses: `ANOMALY.md`.

## Heavy-first block order (#S2, `runs-2026-10-25/s2heavy_b1024.log`, `s2h16k.log`)

No format change. A launch-order permutation is passed to the kernel: grid slot i runs logical (window, block) slot
order[i], order = argsort of the block's compressed stream length `off[b+1] - off[b]` from the block table,
descending, empty slots last; computed on the GPU per call. Same build for both arms (order none / heavy-first),
W = 8 192, the same 10 draws. Check: heavy-first against the CPU decoder 256/256 windows identical (at batch 1 024 and
at 16 384), and heavy-first against before every window of every draw identical (1 024/1 024, 16 384/16 384).
Achieved parallelism = sum of warp durations / (last warp end - first warp start), `%globaltimer`, one instrumented
run per draw and arm.

| batch | arm | median | p90 | max | min | parallelism median (min .. max) |
|---:|---|---:|---:|---:|---:|---:|
| 1 024 | before | 3.161 ms | 4.235 | 7.069 | 1.305 | 113.9 (49.6 .. 268.3) |
| 1 024 | heavy-first, order in the call | 3.662 ms | 4.504 | 7.559 | 1.686 | 111.1 (48.7 .. 270.0) |
| 16 384 | before | 10.973 ms | 11.839 | 11.918 | 9.741 | 561.0 (520.5 .. 634.4) |
| 16 384 | heavy-first, order in the call | 10.212 ms | 10.634 | 11.174 | 8.809 | 612.2 (559.6 .. 721.3) |
| 16 384 | heavy-first, order precomputed | 9.720 ms | 10.138 | 10.646 | 8.298 | - |

Computing the order costs 343 - 357 us at batch 1 024 and 410 - 455 us at 16 384 (409.9 - 454.7). "Before" timed again after the
heavy-first arm: medians 3.138 / 11.041 ms, so there is no drift between arms.

- **Batch 1 024: no gain.** The slowest warp alone is 94 - 100 % of the batch span in every draw (e.g. 6 785 of
  6 964 us), so reordering cannot shorten the batch and the order computation is pure cost (+0.50 ms median).
- **Batch 16 384: -6.9 % median, -10.2 % p90, -6.2 % max with the order in the call; -11.4 % median with it
  precomputed.** Parallelism rises 561 -> 612 (median). The floor is again the slowest block (longest warp 6.4 -
  10.7 ms per draw).

The tail is set by serial rANS on lane 0 (P1 in `ANOMALY.md`); scheduling only moves it.

## Earlier runs: one draw each

`v3_run.log` / `v3_result.json` (`panvram_loader_A100_v3.ipynb`) and `v4_run.log` / `v4_result.json`
(`panvram_loader_A100_v4.ipynb`), 558 assemblies resident, batch 1 024. **One draw** of windows per W, timed 9
times; the spread they report is the clock's, not the workload's.

| run | W | ms per batch (one draw) | windows/s | GB/s |
|---|---:|---:|---:|---:|
| v3 | 8 192 | 3.872 | 264 475 | 2.167 |
| v3 | 32 768 | 2.575 | 397 729 | 13.033 |
| v4 | 8 192 | 3.689 | - | - |
| v4 | 32 768 | 2.534 | - | - |

Two figures from v3 are not draw-dependent and stand: the 558-assembly q16k cohort is resident in **11.567 GB** of
VRAM (`nvidia-smi` 11 466 MiB), and opening it took **24.499 s**. The CPU decoder check of v3:
`CHECK W 8192 windows 256 identical 256 differ 0`.
