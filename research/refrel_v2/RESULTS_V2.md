# refrel v2 prototype (cohort anchors): window law at N = 558

**Status: research prototype. The v2 container ("bt3", version 3) is NOT frozen and is unstable: archives written now
may not be readable by a later v2.** refrel3 v1 (FORMAT.md, the 558 archives of the Zenodo dataset, `csrc/`) is not
changed by anything here. Published: this text, the logs in
[`logs/refrel-v2-2026-10-21/`](../../logs/refrel-v2-2026-10-21/) and the decoder (`rv2dec.cpp`, `rv2.h`,
`build.sh`). **The encoder is not published**, so the compression cannot be repeated from this repository; decoding
and the window law can.

Peak RSS is the `Maximum resident set size` of `/usr/bin/time -v` in kbytes, divided by 10^6 for GB, as in the v1
verdict. Every number below comes from a file of `logs/refrel-v2-2026-10-21/` (named as `L/<path>`) or from an earlier log of
this repository named in full.

## Result: the window law holds at N = 558 for v2, and fails for v1 at 64 KiB

Model B of the hw-apex window-law job (protocol v1.1), the same rules as the v1 verdict at N = 558
(`logs/verdict-n558-v1/README.md`): the same 10 000 requests per W (seed 20261003) and the same truth table (SHA-256
`4f89a4d4...` and `0193fa4c...`, identical files), D_Q = canonical bases / median of 9 full sequential decodes after 3
warm-ups with the XXH3 of every decoded assembly checked, one timed fetch per request with its SHA-256 checked after
the timer, c0 = p50(W = 1) - Q / D_Q, predicted t(W) = c0 + (W + Q - 1) / D_Q, PASS iff c0 >= 0 and |error| <= 20 % at
1 / 8 / 64 KiB. One thread pinned to CPU 3, silence gate before the run; ace-core (AMD EPYC 4344P, Ubuntu 22.04.5,
g++ 11.4.0). 1 680 437 743 110 canonical bases, Q = 16 384.

| N = 558 | D_Q, Gbases/s | c0, us | 1 KiB p50 / error | 8 KiB p50 / error | 64 KiB p50 / error | verdict | peak RSS | log |
|---|---:|---:|---|---|---|---|---:|---|
| refrel3 v1 q16k | 3.544 | 5.936 | 10.199 us / -5.98 % | 13.245 us / +2.91 % | **43.282 us / +48.98 %** | **FAIL (64 KiB)** | 15.2 GB | `logs/verdict-n558-v1/n558_v1.log` |
| **refrel v2, K = 8 anchors, version 3, decoder of this directory** (2026-10-21) | **3.400** | **3.266** | 7.755 us / -7.52 % | 9.968 us / -5.02 % | 32.772 us / **+19.77 %** | **PASS** | **36.6 GB** | `L/law-2026-10-21-rv2dec/n558v3_rv2dec_law.log` |
| same archives, earlier prototype run (2026-10-10, prototype binary with the encoder) | 3.433 | 3.542 | 8.245 us / -4.27 % | 10.399 us / -2.82 % | 32.701 us / +19.32 % | PASS | 36.6 GB | `L/law-2026-10-10-v3/n558v3_rv2.log` |

Every request of every W: SHA-256 == truth, 10 000 / 10 000, bad 0; every full decode XXH3 == truth, 12 x 558 per run
(the `DQ` lines of each log).

What changes from v1 to v2 is c0, the fixed cost of a request: 5.936 -> 3.266 us. D_Q, the sequential decode rate, is
4 % lower in v2 (3.400 against 3.544). The 64 KiB error falls from +48.98 % to +19.77 %: **the margin to the 20 %
threshold at 64 KiB is 0.23 percentage points.**

**The two v2 runs differ** (same archives, requests, truth and rules; binaries `6e1f91b7...` and `7b4152f1...`, in the
`*_inputs.sha256` of each run):

| | 2026-10-21 (decoder of this directory) | 2026-10-10 (prototype binary) | difference |
|---|---:|---:|---:|
| D_Q, Gbases/s | 3.399679 | 3.432540 | -0.96 % |
| c0, us | 3.266 | 3.542 | -0.276 |
| error 1 KiB | -7.52 % | -4.27 % | -3.25 points |
| error 8 KiB | -5.02 % | -2.82 % | -2.20 points |
| error 64 KiB | +19.77 % | +19.32 % | +0.45 points |
| min / max of the 9 timed full decodes, ns | 494 274 856 641 / 494 461 919 521 | 489 382 795 106 / 489 626 698 958 | |

The D_Q difference (0.96 %) is larger than the spread inside either run (0.04 %) and larger than the difference between
the 2026-10-09 version-2 run and the 2026-10-10 version-3 run (0.3 %, below). Its cause (the day, or the code layout of
the binary without the encoder) is NOT MEASURED. The verdict is the same in both runs.

**Costs, stated next to the result:**

- **Memory.** Peak RSS 36.6 GB for v2 against 15.2 GB for v1 in the same protocol, 2.4 times. The cause is the target
  space every v2 archive refers to: T2T (3.117 G bases, README "Limitations") plus the 8 anchor assemblies decoded (24.29 G bases in total,
  the `anchor ...` lines of the logs), about 27.4 GB at 1 byte per base, plus the lift tables. Decoding even one v2
  assembly needs all of it; the 5-assembly subset of `L/subset5/` needs the same memory.
- **The 64 KiB margin is nearly spent.** +19.77 % against the 20 % threshold: 0.23 points left. The mechanism measured for v1 (random
  order over a resident set larger than the caches turns L2 / L3 reuse into DRAM fills; `logs/mech558-2026-10-09/`)
  was not measured for v2 in this package; a larger cohort or a slower memory system can push v2 over the threshold.
- **Version 3 differs from version 2 only in how the block table is read at open.** Payload, grammar and window decode
  are the same code. The version-2 run (`L/law-2026-10-09-v2/n558_rv2.log`, the same 550 payloads with the version-2
  block table: D_Q 3.423, c0 3.159 us, -5.44 / -2.11 / +18.54 %, PASS) and the version-3 run differ by less than the
  difference between two v2 runs above (D_Q 0.3 %, window medians 0.3 - 0.6 us).
- **The decoder here is the prototype decoder with the encoder cut out.** Its run is the headline row; the earlier run of
  the prototype binary on the same archives is kept as a separate row with the difference shown above.

## Size: comparable to AGC 5.40, not a claim of being smaller

| N = 558, T2T not counted | members' bytes | MB per assembly | log |
|---|---:|---:|---|
| AGC 3.2.4, one create with T2T, T2T alone subtracted | 3 014 312 059 | 5.40 | `logs/sizes-558-2026-10-06.txt` |
| **refrel v2, K = 8, version 3** | anchors 116 291 305 + 550 v2 archives 2 856 573 432 = 2 972 864 737 | **5.3277** | `L/roundtrip-n558/n558_k8_m3_sizes.txt` |
| refrel3 v1 q16k without block hashes | 7 277 458 363 | 13.04 | `logs/sizes-558-2026-10-06.txt` |

The v2 size is comparable to AGC's on this cohort. The difference, 1.3 %, is smaller than what a change of cohort or of
the anchor choice can move, and it is not offered as a result. Peak RSS of the v2 decoder (36.6 GB) belongs next to this
size. The 8 anchors are counted as ordinary members at their v1 size without block hashes; they are v1 archives of the
Zenodo dataset (`L/archives/SHA256SUMS.anchors`).

Round trip of the 558 (`L/roundtrip-n558/n558_k8_m3_judge.tsv`): every assembly decoded to FASTA (anchors from their
v1 archives, the 550 others from v2); FASTA bytes and FASTA XXH3 == the source (MANIFEST.tsv) and SHA-256 of the
upper-case bases == the truth of the AGC comparison: 558 / 558 PASS; 5 768 s, peak RSS 40.2 GB.

## Encoder (not published)

A first version: one greedy pass, the v1 parse of refrel3 run over the target space (T2T + the 8 anchors) with one new
copy kind, LIFT (a copy into an anchor at the position that the anchor's own reference copies map from the current
diagonal), no tuning of its costs. Anchors: chosen greedily by the bytes of shared breaks against T2T
(`L/archives/anchors.txt`). None of this can be checked from this repository.

## Not done / open

Carried over in full from the prototype's working notes (section 6 there), with the later state of each item:

- f = 5 anchor lists were selected but not encoded; anchor choice by the measured v2 gain (not the H4 proxy) not tried.
  *(still open)*
- Anchors against anchors (depth 2), anchor block hashes, a GPU decoder, compact (2-bit) target space: not built.
  *(still open; a GPU decoder is a separate prototype, not part of this package)*
- Encoder is the first simple version (one greedy pass, LIFT window +- 8, no tuning of costs). *(still open)*
- Section 3 used one window size (8192) with a derived c0; the model-B window law with a measured c0 at 1 / 8 / 64 KiB
  is in section 8 (K = 4 only). *(superseded: the window law with a measured c0 at N = 558, K = 8, is the result above)*
- N = 558: only the H4 estimate (section 9); no 558 build. *(superseded: built, 558 / 558 round trip, above)*

Further open points of the prototype:

- K = 16 not built: its encoder peak was predicted above the 87.5 GB memory gate of ace-core.
- The format binds an anchor by name, base count and XXH3-64 of its bases; the lift tables come from the anchor's
  encoding, so an anchor re-encoded to the same bases by another build breaks every dependent archive while that check
  passes (only the FASTA XXH3 catches it). A frozen v2 must bind anchors by the SHA-256 of their archives or define the
  lift tables from bases only.
- A damaged anchor damages every assembly that depends on it; anchor block hashes are not implemented.
- The corruption check of this package (`L/negative_controls/`) breaks the decoder, not the archives: archive
  corruption needs the encoder, which is not published. Five deliberate one-base faults (lift, absolute copy into an
  anchor, absolute copy into the reference, self copy, block-table prediction) are each caught by the judge; a sixth,
  `break_lower`, was dropped because the five assemblies of the subset carry no lower-case runs, so it had nothing to
  shift (`L/negative_controls/`).

## Reproduce the decode and the window law

```
bash research/refrel_v2/build.sh <aceapex checkout at e902934> /tmp/rv2     # prints the SHA-256, see build.sh
```
Inputs: T2T (DATASET.md), the 8 anchors and the v1 archives as q16k files of the Zenodo dataset (DATASET.md, "Rebuilding
a tar from its parts"; check `L/archives/SHA256SUMS.anchors`), the v2 archives (5 of them as a release asset of this
repository, `L/subset5/`; all 550 listed in `L/archives/SHA256SUMS.rv2`). Commands and expected outputs:
`logs/refrel-v2-2026-10-21/README.md`.
