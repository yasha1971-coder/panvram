#!/usr/bin/env python3
"""m3_judge.py <assembly> [samples] - #M3 CPU judge through the fork's own entry points.

Two train datasets from BlendedMegatronDatasetBuilder(GPTDataset, ...) with one config and one seed: blend
`/tmp/claude-1002/d2/bin/<assembly>_text_document` (.bin/.idx of tools/preprocess_data.py) and blend
`panvram:/tmp/claude-1002/d2/cohort5` (PANVRAM_ASSEMBLIES=y1_<assembly>), i.e. build_low_level_dataset is the
patched one. The panvram loader is what build_pretraining_data_loader returns. Compared on identical sample IDs:
tokens, labels, loss_mask, position_ids of every row of every micro-batch against GPTDataset.__getitem__ of the
.bin dataset, (1) for the first <samples> samples of the sampler, (2) for every sample that crosses a document
boundary or reaches an end-of-document id.
"""
import argparse
import os
import sys
import time

import numpy as np
import torch

NAME = sys.argv[1] if len(sys.argv) > 1 else "HG00438.1"
NS = int(sys.argv[2]) if len(sys.argv) > 2 else 2048
SEQ = int(os.environ.get("SEQ_LENGTH", "8192"))
MBS = 8
RESETS = os.environ.get("RESETS", "1") == "1"
SEED = 20261003
os.environ["PANVRAM_ASSEMBLIES"] = f"y1_{NAME}"
os.environ["PANVRAM_DEVICE"] = "cpu"
os.environ.setdefault("PANVRAM_THREADS", "16")

sys.path.insert(0, "/tmp/claude-1002/d2/gg")
from megatron.core.datasets.blended_megatron_dataset_builder import BlendedMegatronDatasetBuilder  # noqa: E402
from megatron.core.datasets.gpt_dataset import GPTDataset, GPTDatasetConfig  # noqa: E402
from megatron.training import global_vars  # noqa: E402
from megatron.training.datasets.data_samplers import build_pretraining_data_loader  # noqa: E402
from megatron.training.datasets.panvram_loader import PanvramGPTLoader  # noqa: E402
from megatron.training.tokenizer.tokenizer import _NullTokenizer  # noqa: E402


def config(prefix):
    return GPTDatasetConfig(
        random_seed=SEED, sequence_length=SEQ, blend=([prefix], None), split="100,0,0",
        path_to_cache=None, tokenizer=_NullTokenizer(12), mmap_bin_files=True,
        reset_position_ids=RESETS, reset_attention_mask=RESETS, eod_mask_loss=RESETS,
        create_attention_mask=False,
    )


def build(prefix, n):
    t = time.perf_counter()
    tr, _, _ = BlendedMegatronDatasetBuilder(GPTDataset, [n, 0, 0], lambda: True, config(prefix)).build()
    return tr, time.perf_counter() - t


from megatron.core import parallel_state  # noqa: E402

parallel_state._MPU_DATA_PARALLEL_WORLD_SIZE = 1     # single process, no torch.distributed: data parallel 1, rank 0
parallel_state._MPU_DATA_PARALLEL_RANK = 0
global_vars.set_args(argparse.Namespace(
    full_validation=False, dataloader_type="single", hybrid_context_parallel=False, micro_batch_size=MBS,
    num_workers=0, exit_signal_handler=False, exit_signal=None))

N_TRAIN = 400_000                                   # >= every sample of one assembly at 8 KiB (~370 k)
ds_bin, t_bin = build(f"/tmp/claude-1002/d2/bin/{NAME}_text_document", N_TRAIN)
ds_pan, t_pan = build("panvram:/tmp/claude-1002/d2/cohort5", N_TRAIN)
print(f"{NAME}: build  bin {t_bin:.2f} s  panvram {t_pan:.2f} s;  documents {len(ds_bin.dataset.sequence_lengths)} / "
      f"{len(ds_pan.dataset.sequence_lengths)}", flush=True)
assert np.array_equal(ds_bin.dataset.sequence_lengths.astype(np.int64), ds_pan.dataset.sequence_lengths.astype(np.int64))
for k in ("document_index", "sample_index", "shuffle_index"):
    assert np.array_equal(getattr(ds_bin, k), getattr(ds_pan, k)), k
print(f"  sequence_lengths, document/sample/shuffle index identical; samples {len(ds_pan)}", flush=True)

loader = build_pretraining_data_loader(ds_pan, 0)
assert isinstance(loader, PanvramGPTLoader), type(loader)
FIELDS = ("tokens", "labels", "loss_mask", "position_ids")


def judge(batches):
    bad = {f: 0 for f in FIELDS}
    n = 0
    for ids in batches:
        out = loader.batch(ids)
        for r, i in enumerate(ids):
            ref = ds_bin[i]
            for f in FIELDS:
                if not torch.equal(out[f][r].cpu(), ref[f]):
                    bad[f] += 1
            n += 1
    return n, bad


# (1) the first NS samples in the order the sampler gives them
it = iter(loader.batch_sampler)
first = [next(it) for _ in range(NS // MBS)]
f0, s0 = loader.stats["fast"], loader.stats["slow"]
n1, bad1 = judge(first)
print(f"  sampler order: {n1} samples, fast {loader.stats['fast']-f0} slow {loader.stats['slow']-s0}; " +
      "  ".join(f"{f} {n1-bad1[f]}/{n1}" for f in FIELDS), flush=True)

# (2) every sample that crosses a document boundary or reaches the end-of-document id of its document
low = ds_pan.dataset
inv = np.empty_like(ds_pan.shuffle_index)
inv[ds_pan.shuffle_index] = np.arange(ds_pan.shuffle_index.shape[0], dtype=ds_pan.shuffle_index.dtype)
L = SEQ + 1
tgt = []
si, di = ds_pan.sample_index, ds_pan.document_index
for j in range(si.shape[0] - 1):
    beg, bo = si[j]
    end, _ = si[j + 1]
    i = int(inv[j])
    if i >= len(ds_pan):
        continue
    if beg != end or int(bo) + L > int(low.doc_bases[int(di[beg])]):
        tgt.append(i)
f0, s0 = loader.stats["fast"], loader.stats["slow"]
n2, bad2 = judge([tgt[k:k + MBS] for k in range(0, len(tgt), MBS)])
print(f"  boundary / EOD samples: {n2}, fast {loader.stats['fast']-f0} slow {loader.stats['slow']-s0}; " +
      "  ".join(f"{f} {n2-bad2[f]}/{n2}" for f in FIELDS), flush=True)
ok = not any(bad1.values()) and not any(bad2.values())
print(f"  {NAME}: {'ALL FIELDS IDENTICAL' if ok else 'DIFFER'}", flush=True)
sys.exit(0 if ok else 4)
