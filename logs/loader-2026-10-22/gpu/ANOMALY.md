# The 8 KiB row on the GPU: what has been excluded, and what is left (2026-10-25)

Two Colab runs on an A100-SXM4-40GB, 558 assemblies resident in 11.567 GB, batches of 1 024 windows:

| run | W | ms per batch | blocks decoded | us per decoded block | clocks | spread over repeats |
|---|---:|---:|---:|---:|---|---:|
| v3 | 8 192 | 3.872 | 1 536 (predicted) | 2.521 | not recorded | - |
| v3 | 32 768 | 2.575 | 3 072 (predicted) | 0.838 | not recorded | - |
| **v4** | 8 192 | 3.689 | **1 544 (counted)** | **2.389** | 1 410 MHz | 2.0 % |
| **v4** | 32 768 | 2.534 | **3 072 (counted)** | **0.826** | 1 410 MHz | 0.8 % |

## Excluded: the state of the card

v4 ramps the card for ten seconds before timing, shuffles the window-size order and repeats the whole sweep three
times. Clocks read **1 410 MHz in every row**, before and after, with no throttle reason set; the spread over repeats
is 2.0 % at 8 KiB and 0.8 % at 32 KiB, and the order of the sweep makes no difference. **The cold-card hypothesis of
the v3 run is refuted.** The cost of allocating the output tensor is 14 to 21 us, 0.4 to 0.6 % of a batch, so that is
excluded too.

The difference survives: 2.389 us against 0.826 us per decoded block, a factor of 2.9.

## What the warps actually do

`csrc/pv_cuda.cu:128` launches `bpw = (W + Q - 1) / Q + 1` warps per window and `:46` maps `blockIdx.x` to
`w = blockIdx.x / bpw`, `j = blockIdx.x % bpw`, so **bpw consecutive warps serve one window** and decode its
consecutive blocks. `:50` is `b = b0 + j; if (b > b1) return;` - a warp whose index is past the last block the window
touches returns there, before the shared-memory head/tail are initialised and before any decode.

| W | warps launched | warps that decode | warps that return at `:50` | waves of the 972-block residency |
|---:|---:|---:|---:|---:|
| 8 192 | 2 048 | 1 544 | **504** | 2.11 |
| 32 768 | 3 072 | 3 072 | 0 | 3.16 |

Residency: shared memory per block is `Q + QN * sizeof(QOp)` = 16 384 + 64 x 12 = **17 152 B** (`:129`, `:18-19`), so
an A100 SM with the 164 KiB opt-in carveout holds **9 blocks**, 972 on 108 SMs - `__launch_bounds__(32, 16)` asks for
16 and cannot get them. Tail efficiency is therefore 2.11/3 = 0.70 at 8 KiB against 3.16/4 = 0.79 at 32 KiB.

**So the 504 warps are nearly free, and they do not explain the gap.** They return before touching the queue, they do
not spin to `QSPIN` (that path is `:78`, reachable only by a consumer lane inside the loop), and they do not wait for
another warp's block: there is no cross-warp dependency in this kernel at all - each warp decodes its own block into
its own shared buffer. Nor is a block decoded twice by two windows except when two sampled windows happen to touch
the same block, which over 558 assemblies is rare. The worse tail accounts for about 13 %, not 190 %.

## What is left, and the one measurement that settles it

The work per decoding warp is **the same** whatever W is: a full 16 KiB block is decoded and its copies executed, and
only the window's slice is written out (`:84-85`). The slice is *larger* at 32 KiB (10 923 bases on average against
5 461), so if anything the 32 KiB warp should be the slower one. It is not. Therefore the difference cannot be in the
per-block work, and the surviving candidate is **locality between the warps of one window**:

- `w = blockIdx.x / bpw` puts the warps of one window at adjacent `blockIdx`, so they are co-resident and read
  neighbouring parts of the reference and of the same archive. At bpw = 3 each window has three such warps covering a
  contiguous 48 KiB span of its assembly; at bpw = 2 only two, and 504 windows have a single warp with no neighbour
  at all. The reference is 3.117 G bases in HBM and a block's copies are scattered through it, so L2 reuse between
  the warps of one window is the one quantity that grows with bpw.

This is a hypothesis. **No change to the kernel should be made on it until it is measured**, and two cheap
measurements would settle it, neither of them touching the kernel:

1. Sample the windows **once** at W = 32 768 and time the same starts truncated to 8 192. Same loci, same locality,
   different W. If microseconds per decoded block converge, locality explains the gap; if the 8 KiB number stays at
   2.4 us, it does not.
2. Nsight Compute on the two launches: `lts__t_sector_hit_rate`, `dram__bytes_read.sum` and
   `sm__warps_launched.sum` per kernel. A locality effect shows as a lower L2 hit rate and more DRAM traffic per
   decoded block at 8 KiB.

Until then the cause is **NOT MEASURED**, and the 8 KiB GPU row is not a statement about window size.
