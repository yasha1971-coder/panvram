# The 8 KiB row on the GPU: under-occupancy and the spread between draws (2026-10-25)

A100-SXM4-40GB on Colab, 1 410 MHz with no throttle reason in every run below, cohort of 5 assemblies (`q16k`,
Q = 16 384). Runs and their logs are in `runs-2026-10-25/` (`v4_run.log` in this directory); every number below carries its log.

## What it is

The v3/v4 rows (558 assemblies resident; 8 KiB 3.689 ms against 32 KiB 2.534 ms per batch of 1 024, "2.9x more per
decoded block") were
**one draw of windows timed nine times**. Their 2 % spread was the spread of the clock, not of the workload. Two
effects explain them, and both are measured:

**1. Under-occupancy at batch 1 024.** `T = alpha + B * beta` fitted over B = 256 .. 16 384 on one fixed list of
windows (`b1.log`):

| W | alpha, kernel only | beta per window | alpha in a batch of 1 024 | GB/s at 1 024 | GB/s at 16 384 |
|---:|---:|---:|---:|---:|---:|
| 8 192 | 2 013 us | 0.466 us | 81 % | 2.945 | 14.393 |
| 32 768 | 1 896 us | 0.878 us | 68 % | 10.424 | 33.431 |

Host time per call is 71 - 112 us and reading the status word 0 - 99 us, so alpha is inside the kernel: it is the
latency of one decoding warp (median 170.8 us of rANS on lane 0, see 2.), not a launch cost. A batch of 1 024 at
8 KiB has 1 534 decoding warps against 972 resident slots (17 152 B of shared memory per block, 9 per SM, 108 SMs),
so the card is mostly waiting for the last warps.

**2. A batch waits for its slowest block, and the slowest block depends on the draw.** Ten independent draws per W,
batch 1 024, five timings each (`b1dist.log`):

| W | median | p90 | max | min | max/min | decoded blocks |
|---:|---:|---:|---:|---:|---:|---:|
| 8 192 | 3.125 ms | 4.175 | 6.899 | 1.270 | 5.43x | 1 544 +- 10 |
| 32 768 | 3.216 ms | 3.642 | 4.249 | 2.202 | 1.93x | 3 072 |

Over draws the two window sizes cost the same median time per batch. Per-block instrumentation of the slowest draw
(seed 20261003, W 8 192, clock64 around `r3_decode_stream` on lane 0):

- batch wall 7.740 ms; 2 048 warps launched, 1 534 decoding;
- decode cycles per block: median 240 833 (170.8 us), p90 521 529, p99 1 230 577, **max 9 696 905 (6.88 ms)**,
  max/median 40.3x;
- **one block is 6.88 ms of the 7.74 ms batch**; achieved concurrency (sum of block cycles / wall) 40.2x of 972 slots;
- slowest blocks carry 3 919 / 2 258 / 977 B of payload against a median of 38 B; correlation of cycles with payload
  bytes **+0.966**, with emitted operations +0.676. Payload is literal-heavy: rANS on one lane is the cost.

At batch 16 384 (the headline, `b2head.log`, 10 draws, 3 timings each) the tail is diluted but still visible:

| W | median | p90 | max | min | max/min | windows/s median | GB/s median / worst |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 192 | 11.938 ms | 13.681 | 14.766 | 10.541 | 1.40x | 1 379 788 | 11.303 / 9.089 |
| 32 768 | 17.427 ms | 18.781 | 19.332 | 15.644 | 1.24x | 940 250 | 30.810 / 27.772 |

## Refuted hypotheses

| hypothesis | test | result | log |
|---|---|---|---|
| cold card / clocks | 10 s ramp, sweep order shuffled, clocks read per row | 1 410 MHz in every row, no throttle reason; repeats spread 2.0 % (8 KiB) / 0.8 % (32 KiB); difference survived | `v4_run.log` |
| locality between the warps of one window | same 1 024 starts drawn once at 32 KiB, timed at 32 KiB and truncated to 8 KiB | same windows: 8 KiB **2.905 / 2.886 ms**, 32 KiB 3.198 / 3.195 ms - 8 KiB is *faster*; the 7.033 ms came only with 8 KiB's own draw | `anom3.log` |
| short contigs (unplaced scaffolds) entering the 8 KiB draw | draw from contigs >= 32 768 only vs < 32 768 only | 428 of 430 contigs are >= 32 768 and hold 100.00 % of the bases; long-contig pool 1.964 ms vs all 2.205 ms; short pool 0.383 ms (fastest) | `anomctg.log` |

The seed-20261003 `Cohort.sample` draw reproduces at 7.033 ms (min 6.852, max 7.216) in a separate session, while
other draws of the same size take 2.266 - 3.911 ms (`anomab.log`): the draw, not the card, sets the time.

"Microseconds per decoded block" (batch time / blocks) was the wrong unit for a batch dominated by alpha and by one
slow block; it is not used any more.

## Open direction (P1, not started)

Lane 0 decodes a block's rANS stream serially: median 170.8 us per block, alpha about 2 ms per batch, and the
slowest block 6.88 ms of a 7.74 ms batch (`b1dist.log`, `b1.log`). A multi-way interleaved rANS (several states per block decoded by
several lanes) would cut both alpha and the tail; it changes the stream layout, so it is a format decision and is
recorded here only.
