# GPU path of A: what is measured and what is still open (2026-10-25)

The CPU rows of this directory were measured by the author on ace-core. The GPU numbers here come from one Colab run
of `panvram_loader_A100_v3.ipynb` on 2026-10-25, **NVIDIA A100-SXM4-40GB**, with the 558-assembly q16k cohort
resident on the card. Raw output as received: `v3_run.log`; the same with the derived block figures: `v3_result.json`.

| | value |
|---|---|
| cohort opened | 558 assemblies, **24.499 s** |
| resident on the card | **11.567 GB** allocated and peak; `nvidia-smi` 11 466 MiB |
| batch / timed runs | 1 024 windows / 9, median reported |
| correctness | `CHECK W 8192 windows 256 identical 256 differ 0` - every window equals the CPU decoder's answer |

| W | ms per batch | min / max, ms | windows/s | GB/s delivered | blocks decoded per batch | us per block decode |
|---:|---:|---|---:|---:|---:|---:|
| 8 192 | 3.872 | 3.693 / 4.007 | 264 475 | 2.167 | 1 536 | **2.521** |
| 32 768 | 2.575 | 2.501 / 2.766 | 397 729 | 13.033 | 3 072 | **0.838** |

**The 8 KiB row is an open anomaly; it is not a window-size result and nothing here supersedes the CPU table of the
parent README. A v4 run is awaited.** Why:

The last column is the problem. Each warp decodes exactly one 16 KiB block (`csrc/pv_cuda.cu`: one warp per block,
lane 0 runs `r3_decode_stream`, lanes 1..31 execute the copies), so microseconds per decoded block must not depend on
W. Here it differs by **3.0x**.

The block counts are not a guess. `csrc/pv_cuda.cu:128` launches `bpw = (W + Q - 1) / Q + 1` warps per window, and a
window of W bases starting at an arbitrary offset touches `1 + (W-1)/Q` blocks on average: 1.50 blocks at W = 8 192,
3.00 at W = 32 768 with Q = 16 384. So the 8 KiB batch decodes **fewer** blocks (1 536 against 3 072) and must
therefore be the **cheaper** of the two. With equal cost per block, the 32 KiB batch should have taken
3.872 x 3 072 / 1 536 = 7.744 ms; it took 2.575 ms, a factor of 3.0 the other way.

A second consequence of the same formula, which is a real property and not an artefact: at q16k a 8 KiB window
decodes 24 575 bases to deliver 8 192, a 3.00x overhead, while a 32 KiB window decodes 49 151 to deliver 32 768,
1.50x. Per delivered base the small window costs twice the decode. That is why the GB/s column rises with W on the
CPU rows as well (3.37 -> 5.92 GB/s at N = 558).

## The hypothesis for the 3.0x

W = 8 192 was measured **first**, immediately after `Cohort.open` pushed ~14 GB to the card - a memcpy-bound phase
during which the SM clocks sit low - and after only three warm-up batches of 8 MB each, far too little work to bring
an A100 from its idle clock to boost. By the time W = 32 768 was measured the card had been busy for seconds. An idle
to boost ratio on A100 is of the order of 1 410 / 300 MHz, which is enough to cover a factor of 3 on its own. Nothing
in the kernel, the block arithmetic or the output allocation predicts a 3.0x advantage for the larger window; the
order of the sweep does.

This is a hypothesis, not a measurement: the v3 run recorded no clock readings, so **the cause is NOT MEASURED**.

## What v4 changes, so the next run settles it

`panvram_loader_A100_v4.ipynb`:

- a 10-second ramp at the largest W before anything is timed, with `nvidia-smi` clocks printed before and after;
- the W order **shuffled**, the whole sweep repeated 3 times, every repeat reported, and the spread over repeats
  printed as a percentage - an order effect cannot hide in a median any more;
- blocks decoded per batch computed from the drawn coordinates, so **microseconds per decoded block** is printed
  directly; that number must be flat in W, and if it is not, the cause is not the window size;
- the cost of allocating the output tensor measured on its own and printed beside each batch time, so it can be
  subtracted;
- clocks and throttle reasons recorded with every row, and the whole result written to `gpu_result.json`.

Until that run exists, the two rows above are the only GPU numbers in this directory, the 8 KiB one is flagged as an
anomaly, and the window-size comparison on the GPU stays open.

Two figures from this run are measurements in their own right and are not in doubt: the 558-assembly q16k cohort is
resident in **11.567 GB** of VRAM, and opening it took **24.499 s**. The 2-bit and compact forms measured elsewhere in
this repository are not used here.
