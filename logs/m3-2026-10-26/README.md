# M3: panvram as the data path of Megatron / Gengram (runs of 2026-10-25)

Code: `tools/megatron/` (patch against BGI-HangzhouAI/Gengram `c7053fa`, applied there as local branch commit
`4ee6879`). The scripts below carry the paths of the machine they ran on (`/tmp/claude-1002/...`); they are the
exact files that produced the logs.

## judge/ - fields against .bin/.idx, CPU, ace-core

`m3_judge.py <assembly> <samples>` builds two train datasets with `BlendedMegatronDatasetBuilder(GPTDataset, ...)`,
one config, seed 20261003, `seq_length` 8192, resets on: blend `<assembly>_text_document` (`.bin/.idx` from
`tools/preprocess_data.py --append-eod`, NullTokenizer vocab 12) and blend `panvram:<cohort>` (the patched
`build_low_level_dataset`). The loader is what the patched `build_pretraining_data_loader` returns. Compared:
tokens, labels, loss_mask, position_ids of every row against `GPTDataset.__getitem__` of the `.bin` dataset.

| assembly | documents | sampler order | fast / slow rows | boundary or end-of-document samples | all four fields |
|---|---:|---:|---|---:|---|
| y1_HG00438.1 | 276 | 4 096 | 4 092 / 4 | 551 | identical |
| y1_HG00438.2 | 258 | 4 096 | 4 095 / 1 | 515 | identical |
| y1_HG00621.1 | 287 | 4 096 | 4 093 / 3 | 573 | identical |
| y1_HG00621.2 | 253 | 4 096 | 4 092 / 4 | 505 | identical |
| y1_HG00673.1 | 420 | 4 096 | 4 091 / 5 | 839 | identical |
| y1_HG00438.2, resets off | 258 | 2 048 | 2 047 / 1 | 515 | identical |

Sequence lengths and the document, sample and shuffle index are asserted identical before any field is compared.
One defect was found on the way and fixed before these runs: GPTDataset caches its indices under
`<path_prefix>/cache` keyed by the data path string, which for `panvram:<dir>` does not name the assemblies, so a
second selection reused the first one's indices; `path_prefix` is now one directory per selection.

## cpu_gated/ - loader rate, CPU, ace-core

AMD EPYC 4344P (8 cores, 16 threads), 125 GB, Linux 5.15, torch 2.14.1+cpu. y1_HG00438.1, micro-batch 8, the same
136 micro-batches for every path (8 warm-up, 128 timed), median of 5 passes, one-minute load average below 0.5
before every pass (0.16 - 0.46, printed per pass in the log). panvram: CPU cohort, 16 decoder threads.

| path | 8 KiB samples/s | 32 KiB samples/s |
|---|---:|---:|
| panvram loader | **13 947** | **8 030** |
| `.bin` DataLoader, 0 workers | 12 159 | 4 859 |
| `.bin` DataLoader, 2 workers (Megatron default) | 4 581 | 2 084 |
| `.bin` DataLoader, 8 workers | 4 996 | 3 339 |

The `.bin` file was in the page cache. The passes are short (0.1 - 0.5 s), so the gate is read before a pass, not
during it.

## phaseB/ - 200 training steps on a GPU, Colab L4

`m3b.py` (and the same code as `panvram_gengram_M3_phaseB.ipynb`): fresh Colab VM, NVIDIA L4 23 034 MiB, the fork at
`c7053fa` plus the patch, one assembly (r2_HG00097_hap1, 75 contigs: alphabetically the first of the five archives `colab_prepare.prepare_cohort(5)` takes from the record). The
`.bin/.idx` is written on the VM from the same archive through Megatron's `IndexedDatasetBuilder` (one document per
contig, end-of-document id appended; sequence lengths asserted equal to the panvram dataset). Model: 4 layers, hidden
256, 4 heads, Gengram on layer 2 (n-grams 2..6, 64 per n-gram, short conv 4, window 12), RMSNorm, RoPE, sequence
2 048, micro and global batch 16, 200 steps, bf16, seed 1234, `--transformer-impl local`, `--reset-position-ids
--eod-mask-loss`. GPU utilisation from `nvidia-smi` every 200 ms.

| data path | ms per step, median after step 20 | tokens/s | GPU util median / mean | wall |
|---|---:|---:|---|---:|
| `.bin/.idx`, DataLoader 2 workers | 383.0 | 85 556 | 92 % / 63.1 % | 119 s |
| `panvram:` | 384.2 | 85 278 | 91 % / 60.6 % | 123 s |

Disk: `.bin` 6 068 262 842 B against the archive 13 279 654 B. lm loss every 10 steps is identical in the first five
logged values and then differs in the fifth decimal (1.344938 / 1.344950 at step 60); the data are the same by the
CPU judge, the difference is consistent with non-deterministic GPU reductions, which these runs do not switch off.
The mean utilisation includes start-up and shutdown.

Runs 1 - 4 failed on the fork's build without Transformer Engine / Apex and are kept as logs
(`m3b_run{1..4}_failed.log`: pybind11 missing; `apply_rope_fusion` needs TE; `helpers_cpp` built for another
interpreter; `scaled_masked_softmax_cuda` missing). `cpu_model_test.py` is the CPU check of model construction and one
forward and backward pass used between runs. Five runs cost 0.45 Colab compute units.

**Boundary:** at this model size both paths are bound by GPU compute; equal speed says the panvram path adds nothing,
not how it behaves when data loading is the bottleneck.
