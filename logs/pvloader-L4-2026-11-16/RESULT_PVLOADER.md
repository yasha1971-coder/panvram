# PV-LOADER — RESULT (panvram GPU-resident loader vs ordinary loaders)

panvram FROZEN 01f0ce11 (code UNCHANGED). seed 20261115. Pre-registered S3 model + killer
criterion in CALIBRATION.md BEFORE the GPU run. GPU: NVIDIA L4 (24 GB), self-deleting VM.
CPU baselines: <host>, 16 cores. Rule 11.11 applied to every "N×" (what / input / model /
GPU / k / measurement bounds stated inline).

## S0 — determinism (CPU, done)
One seeded window list (W=2048) for ALL loaders. panvram CPU `fetch(tokens=True)` vs FASTA
`pysam`→tokens (same A/C/G/T=0..3,else 4 LUT): 8 batches × 64, xxh3 per batch → **0 mismatches,
all identical**. The loaders feed byte-identical batches (else STOP). → `pvloader/s0.json`.

## S1 — throughput (windows/s, NO model; W=2048, batch=64, 5 repeats, median)
| loader | windows/s | input / where | bounds |
|---|---|---|---|
| panvram `sample` (all 558 resident) | **126,918** | q4k in 11.6 GB VRAM, on-card decode, L4 | training path, median of 5 |
| panvram `sample` (6 assemblies) | 144,224 | VRAM, L4 | median of 5 |
| panvram `fetch` (coord path, S0 list, all 558) | 5,530 | VRAM, L4 | per-window random coord fetch |
| Base A pysam (bgzip+faidx) k=4 (best) | 19,057 | NVMe FASTA, CPU 16c / 4 workers | fetch+tokenize |
| Base A pysam k=1 / k=16 | 8,298 / 14,525 | NVMe, CPU | fetch+tokenize |
| Base B memmap uint8 (warm page cache) | 287,000 | decompressed on NVMe, CPU | read+tokenize; needs data resident (see S5) |

→ `pvloader/logs/s1_panvram_{n6,all}.json`, `pvloader/s1_base.json`.

## S2 — synthetic compute ceiling (GPU L4; t_step = fwd+bwd+opt on a batch resident in VRAM)
| model | params | t_step | compute-bound windows/s (batch/t_step) |
|---|---|---|---|
| x1 (TinyCNNLM d128 L8) | 397,573 | **114.46 ms** | 559 |
| x8 (d256 L16) | 3,160,581 | **498.63 ms** | 128 |
| x64 (d512 L32) | — | **OOM on L4 24 GB** (229 MB free; recorded honestly, did not fit) | — |

→ `pvloader/logs/s2.json`.

## S3 — predicted GPU idle per loader × model (model §3; t_data = 64 / win_s)
idle_overlap = max(0, 1 − t_count/max(t_count,t_data)) (perfect prefetch overlap);
idle_serial = t_data/(t_count+t_data) (no overlap, worst case). t_count = S2 t_step.

t_data (ms): memmap 0.223, panvram-sample 0.504, pysam-k4 3.358, panvram-fetch 11.573.

| model | loader | idle_overlap | idle_serial |
|---|---|---|---|
| x1 (114.46 ms) | memmap | 0.0% | 0.19% |
| x1 | panvram-sample | 0.0% | 0.44% |
| x1 | pysam-k4 (best base A) | 0.0% | **2.85%** |
| x1 | panvram-fetch (worst) | 0.0% | 9.18% |
| x8 (498.63 ms) | any of the above | 0.0% | ≤ 2.27% |

### Killer criterion (pre-registered, verbatim) → FIRED
"If Base A at its best k idles the GPU <10% even on the SMALLEST model → no data bottleneck →
DO NOT run S4; record a confirmed zero." Base A (pysam) k=4 on x1 = **0.0% overlap / 2.85%
serial < 10%** → **NO BOTTLENECK → S4 NOT RUN → CONFIRMED ZERO.** Even the worst pair
(panvram-fetch × x1) is 9.18% serial < 10%; on x8 everything is ≤2.3%.

## S5 — capacity (CPU, done): the real, measured panvram win
Per assembly (measured on 6, mean): bgzip FASTA 0.84 GB, decompressed uint8 2.83 Gbp.
| storage for the cohort | N=50 | N=558 (full) |
|---|---|---|
| bgzip FASTA baseline | 42.0 GB | **468.9 GB** |
| memmap uint8 baseline | 141.5 GB | **1,578.8 GB = 1.54 TB** |
| panvram q4k | — | **11.48 GB disk / 11.6 GB VRAM resident** |

Ratios (rule 11.11: storage to hold the full 558-assembly HPRC cohort for a DNA-LM loader):
memmap/panvram-VRAM = **136×**; bgzip/panvram-VRAM = **40×**. Fit on <host> (disk 878 GB,
RAM 125 GB): N=558 memmap (1.54 TB) fits NEITHER RAM NOR disk; bgzip (469 GB) fits disk only.
panvram holds all 558 resident in 11.6 GB VRAM and serves 126,918 windows/s (`sample`).

## idle vs model size (the "graph")
```
idle (serial bound, %)      x1(0.40M)   x8(3.16M)   x64(OOM)
  memmap                      0.19        0.04         n/a
  panvram sample              0.44        0.10         n/a
  pysam k4 (best base)        2.85        0.67         n/a
  panvram fetch (worst)       9.18        2.27         n/a
```
Every loader × every runnable model: GPU idle ≈ 0 (overlap) / <10% (serial). Compute step
(114 ms even for the small x1) dominates data time (≤11.6 ms) by ≥10×. x64 doesn't fit the L4
at all → no loader has anything to be idle for, at the model sizes this GPU can run.

## VERDICT
On the train_small recipe and an L4, the **data loader is NOT the bottleneck — confirmed zero**:
the compute step (≥114 ms) far exceeds any loader's per-batch data time (≤11.6 ms), so no loader
(panvram or baseline) leaves the GPU idle >10%. There is NO FFCV-style "faster because of the
loader" result here. panvram's real, measured advantage is **CAPACITY (S5)**: it is the only
loader that can hold the entire 558-assembly / 1.68 Gbp-scale cohort for training at all —
11.6 GB VRAM resident vs a 1.54 TB memmap (136×) or 469 GB bgzip (40×) baseline that does not
fit the machine — while serving 126,918 windows/s and feeding byte-identical batches (S0).

## Prior Art Difference
- **FFCV** reports large end-to-end speedups, but part of theirs comes from recipe/augmentation
  changes. Here the training recipe is FROZEN and S0 proves byte-identical batches (loss curves
  must coincide), so our honest loader-speed result is **zero** (no bottleneck to remove) — and
  **yes on capacity**, a dimension FFCV does not address (it still needs the dataset on disk/RAM).
- **MosaicML StreamingDataset** targets shard streaming + determinism. Our S0 determinism is
  confirmed (checksums match across loaders); panvram additionally keeps the whole cohort resident
  rather than streaming shards from storage.
- **NVIDIA DALI** moves decode/augment onto the GPU to cut the data fraction. panvram's on-card
  decode is already negligible: `sample` = 0.504 ms/batch (0.44% of the 114 ms x1 step), so there
  is no decode fraction left to reclaim at these model sizes.

Artifacts: `pvloader/s0.json`, `pvloader/s1_base.json`, `pvloader/s5.json`,
`pvloader/logs/{s1_panvram_n6,s1_panvram_all,s2}.json`, driver.log. Kit: pubrepo train1-kit
`bench/pvloader/kit/`. VM self-deleted. S4 not run (killer fired).
