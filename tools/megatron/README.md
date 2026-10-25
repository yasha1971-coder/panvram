# panvram as a Megatron data path (Gengram fork c7053fa)

Train a Megatron-LM / Gengram model straight from panvram archives: `--data-path panvram:<cohort directory>` instead
of a `.bin/.idx` prefix. No preprocessing step, no token file on disk; the cohort is resident on the GPU and one decode
call fills a micro-batch.

| file | what it is |
|---|---|
| `panvram_dataset.py` | `PanvramIndexedDataset`: Megatron's low-level dataset contract (`sequence_lengths`, `get(idx, offset, length)`, `document_indices`) over a panvram cohort. One document = one contig, archive record order; Gengram ids A=8 C=5 G=6 T=7 other=9; a virtual end-of-document id 12 (NullTokenizer `--vocab-size 12`) ends every document, as `tools/preprocess_data.py --append-eod` writes it. |
| `panvram_loader.py` | `PanvramGPTLoader`: replaces torch's `DataLoader` for that dataset. Megatron's own document / sample / shuffle index is kept; a sample inside one document goes through one `Cohort.windows` call per micro-batch on the device, every other sample (contig boundary, end-of-document id, padding, `create_attention_mask`) through `GPTDataset.__getitem__` itself. |
| `gengram-c7053fa-panvram.patch` | `git format-patch` against [BGI-HangzhouAI/Gengram](https://github.com/BGI-HangzhouAI/Gengram) `c7053faf087acf05bf3edb8ebe2a22c84620ed9d`: adds the two files above under `megatron/training/datasets/`, +5 lines in `GPTDataset.build_low_level_dataset` (`panvram:` prefix) and +8 lines in `build_pretraining_data_loader` (return the panvram iterator). `pretrain_gpt.py` is not touched. |

## Licence

The fork is NVIDIA Megatron-LM under its **BSD 3-Clause** licence (`LICENSE` of the fork, "Copyright (c) 2019-2025,
NVIDIA CORPORATION"), with Gengram's own module `megatron/core/transformer/gengram.py` under **Apache-2.0**. The patch's
changes to `gpt_dataset.py` and `data_samplers.py` are modifications of BSD-3-Clause files and stay under that licence.
`panvram_dataset.py` and `panvram_loader.py` are new files, MIT like the rest of this repository. Nothing of the fork is
copied here except the context lines of the patch.

## Reproduce

```bash
git clone https://github.com/BGI-HangzhouAI/Gengram.git && cd Gengram
git checkout c7053faf087acf05bf3edb8ebe2a22c84620ed9d
git apply /path/to/panvram/tools/megatron/gengram-c7053fa-panvram.patch
pip install -e /path/to/panvram regex einops pybind11
# Megatron's index helpers, built for the interpreter that will train
cd megatron/core/datasets && g++ -O3 -shared -std=c++11 -fPIC $(python -m pybind11 --includes) helpers.cpp \
    -o helpers_cpp$(python -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))') && cd -

# a cohort directory: reference.fa + *.q16k.rr3 (see DATASET.md); optional PANVRAM_ASSEMBLIES=a,b  PANVRAM_DEVICE=cuda|cpu
PANVRAM_ASSEMBLIES=r2_HG00097_hap1_hprc_r2_v1.0.1 \
torchrun --nproc_per_node 1 pretrain_gpt.py --data-path panvram:/content/cohort \
    --tokenizer-type NullTokenizer --vocab-size 12 --split 100,0,0 \
    --gengram-enabled --gengram-token-ids 8 5 6 7 9 ...   # the rest of the model and training arguments
```

The exact run measured below, with every argument, is `logs/m3-2026-10-26/phaseB/m3b.py` (also as a Colab notebook,
`panvram_gengram_M3_phaseB.ipynb`). Without Transformer Engine and Apex the fork additionally needs
`--transformer-impl local --no-rope-fusion --no-persist-layer-norm --no-masked-softmax-fusion
--no-gradient-accumulation-fusion`; with `--transformer-impl local`, `--no-create-attention-mask-in-dataloader` keeps
the data loader from building a full attention mask per sample (both paths were run that way).

## Prior Art Difference

Megatron trains from a `.bin/.idx` pair written by `tools/preprocess_data.py`: every token stored, 2 bytes each
(`uint16` for a vocabulary under 65 500). panvram serves the same tokens from the compressed archive.

| | `.bin/.idx` | panvram | source |
|---|---:|---:|---|
| disk, one assembly (r2_HG00097_hap1, 75 contigs) | **6 068 262 842 B** (6.07 GB) | **13 279 654 B** (13.3 MB) archive | `logs/m3-2026-10-26/phaseB/m3b_run5_ok.log` |
| ratio | | **457x** smaller | |
| preprocessing | a token file must be written first | none | |
| fields of `GPTDataset.__getitem__` on identical sample IDs, five assemblies | reference | **identical**, 4 fields, 23 463 samples (+ 2 563 with the masks not reset) | `logs/m3-2026-10-26/judge/` |
| training speed, L4, 200 steps (median ms per step after 20) | 383.0 ms, 85 556 tokens/s, GPU util median 92 % | 384.2 ms, 85 278 tokens/s, 91 % | `logs/m3-2026-10-26/phaseB/` |

**Boundary of that last row.** The model was small (4 layers, hidden 256, Gengram on layer 2, sequence 2 048, global
batch 16, bf16, one L4) and both paths keep the GPU at ~92 %: the step is bound by compute, not by data. Equal speed
here means the panvram path costs nothing at this size; it does **not** show what happens when the data path is the
bottleneck (larger batches per GPU, faster GPUs, many GPUs per host reading one disk). The CPU-side loader rate,
measured separately (`logs/m3-2026-10-26/cpu_gated/`), is 13 947 samples/s at 8 KiB against 12 159 for the `.bin`
DataLoader without workers and 4 581 with Megatron's default 2 workers.
