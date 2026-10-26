# Headline at 558 assemblies, batch 16 384 (Colab, 2026-10-26)

Raw output: `h558.log`. Script: `h558.py`. Commands, from the directory holding the scripts:

```bash
colab run --timeout 600 -s h558cpu speed.py                 # CPU VM: download speed of Zenodo from Colab
colab run --gpu A100 --timeout 1500 -s h558 h558.py         # GPU VM, released when the script ends
```

## Card

`--gpu A100` was requested. The log records the clocks: `nvidia-smi --query-gpu=clocks.sm,clocks_throttle_reasons.active`
reports **1 410 MHz, throttle reasons 0x0** after a 5 s ramp. The card's name and memory size were **not printed** by
this script (they went only to `/content/h558.json` on the released VM), so 40 GB or 80 GB is NOT MEASURED here.

## Archive

The 558 archives of record 10.5281/zenodo.23232317, as `cohort558_q16k.tar` in 16 parts. `PARTS.sha256` and
`cohort558_q16k.tar.sha256` are the record's own files, fetched for this directory. In the run
(`colab_prepare.prepare_cohort(558)`): every part's SHA-256 matched `PARTS.sha256` (`part ok` x 16), the joined tar of
8 098 467 840 B matched `3b8c497e9193984724128e0fd18df02d131a05c341df41df44fe3ed89a8114d5`, and the 558 archives
inside matched the tar's `SHA256SUMS` (`VERIFIED 558 archives`). Reference: GCA_009914755.4, SHA-256 of the
decompressed FASTA `51ab7d3c2a375ecff53ee2a1be03e58d57cc1601c948a860c03a78d724f26e4c`, asserted in the script.

## Timeline (from the log)

| step | time |
|---|---:|
| Zenodo, part00 (524 288 000 B) in 4 range streams, CPU VM | 7.5 s, 69.9 MB/s (`speed.py`; the line it printed: `part00 524288000 B in 7.5 s = 69.9 MB/s, 4 range streams, sha256 ok True`) |
| 16 parts, 4 at a time, during the extension build | 109.1 s |
| setup until the parts were in | 197 s |
| verify, join, extract, verify 558 | 114 s |
| `Cohort.open`, 558 assemblies | 24.711 s, 11.567 GB resident |

## Check

`CHECK W 8192 draw 0: 256/256 windows identical to the CPU decoder`, the same at W 32 768.

## Result: batch 16 384, 10 independent draws (seeds 20261003 .. 20261012), 3 timings per draw, median per draw

| W | time median / p90 / max | windows/s median / p90 (slow end) / worst | GB/s median / p90 (slow end) / worst |
|---:|---|---|---|
| 8 192 | 13.357 / 15.916 / 16.838 ms | 1 226 694 / 1 029 738 / 973 039 | 10.049 / 8.435 / 7.971 |
| 32 768 | 20.582 / 22.494 / 22.577 ms | 796 046 / 728 381 / 725 682 | 26.085 / 23.867 / 23.779 |

The p90 of the GB/s column is computed from the ten per-draw rows of the log (the script prints its median and worst).
Compute units used, CPU probe and GPU run together: 0.44 (balance 18.01 -> 17.57).
