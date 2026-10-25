# Copyright (c) 2026, Yakiv Shavidze. MIT licence.
"""PanvramGPTLoader: the batch iterator for a GPTDataset over a PanvramIndexedDataset.

It replaces torch's DataLoader for that one case (data_samplers.build_pretraining_data_loader): DataLoader would call
GPTDataset.__getitem__ once per sample in worker processes with pin_memory, which cannot hold device tensors. Here the
index arithmetic stays Megatron's (document / sample / shuffle index of the GPTDataset, the same arrays as for
.bin/.idx), and only the fetch and the masks change:

  fast path  a sample inside one document that does not reach its end-of-document id -> one (assembly, start) pair;
             every such sample of the micro-batch is decoded by one `Cohort.windows` call where the cohort is
             resident, tokens 0..4 mapped to the Gengram ids on the device;
  slow path  every other sample (crossing a document boundary, carrying an end-of-document id, padded, idx None, or
             create_attention_mask) -> GPTDataset.__getitem__ itself, copied into its row.

The masks of the fast rows repeat _get_ltor_masks_and_position_ids and the padding rules of GPTDataset.__getitem__,
vectorised over the batch. Output: the dict the DataLoader's default collate would give, on the cohort's device.
"""
import numpy as np
import torch

from megatron.training.datasets.panvram_dataset import PanvramIndexedDataset


def maybe_panvram_loader(dataset, batch_sampler):
    """The panvram iterator if `dataset` is a GPTDataset over a PanvramIndexedDataset, else None."""
    low = getattr(dataset, "dataset", None)
    if not isinstance(low, PanvramIndexedDataset):
        return None
    return PanvramGPTLoader(dataset, batch_sampler)


class PanvramGPTLoader:
    def __init__(self, dataset, batch_sampler):
        self.ds = dataset
        self.batch_sampler = batch_sampler
        self.low = low = dataset.dataset
        cfg = dataset.config
        self.cfg = cfg
        self.dev = low.device
        self.extra = int(cfg.add_extra_token_to_sequence)
        self.S = int(cfg.sequence_length)
        self.L = self.S + self.extra
        self.eod = int(cfg.tokenizer.eod)
        self.pad = int(dataset._pad_token_id)
        self.lut5 = torch.tensor(low.token_ids, dtype=torch.int64, device=self.dev)   # panvram token 0..4 -> id
        self.ar = torch.arange(self.S, dtype=torch.int64, device=self.dev)
        self.stats = {"fast": 0, "slow": 0}

    def __len__(self):
        return len(self.batch_sampler)

    def __iter__(self):
        for idxs in self.batch_sampler:
            yield self.batch(idxs)

    # ---- one micro-batch ----
    def _plan(self, idxs):
        """Split the batch into fast rows (assembly, start) and slow rows."""
        ds, low = self.ds, self.low
        fast_r, fast_a, fast_s, slow_r = [], [], [], []
        if self.cfg.create_attention_mask:
            return fast_r, fast_a, fast_s, list(range(len(idxs)))
        for r, idx in enumerate(idxs):
            if idx is None:
                slow_r.append(r)
                continue
            j = ds.shuffle_index[idx]
            beg, bo = ds.sample_index[j]
            end, eo = ds.sample_index[j + 1]
            if beg == end and int(eo) - int(bo) + self.extra == self.L:
                doc = int(ds.document_index[beg])
                if int(bo) + self.L <= int(low.doc_bases[doc]):          # stays before the end-of-document id
                    fast_r.append(r)
                    fast_a.append(int(low.doc_asm[doc]))
                    fast_s.append(int(low.doc_start[doc]) + int(bo))
                    continue
            slow_r.append(r)
        return fast_r, fast_a, fast_s, slow_r

    def batch(self, idxs):
        idxs = list(idxs)
        n, S, dev = len(idxs), self.S, self.dev
        fast_r, fast_a, fast_s, slow_r = self._plan(idxs)
        self.stats["fast"] += len(fast_r)
        self.stats["slow"] += len(slow_r)

        tokens = torch.empty((n, S), dtype=torch.int64, device=dev)
        labels = torch.empty((n, S), dtype=torch.int64, device=dev)
        loss_mask = torch.empty((n, S), dtype=torch.float32, device=dev)
        position_ids = torch.empty((n, S), dtype=torch.int64, device=dev)
        attention_mask = None

        if fast_r:
            asm = torch.tensor(fast_a, dtype=torch.int64, device=dev)
            st = torch.tensor(fast_s, dtype=torch.int64, device=dev)
            raw = self.low.cohort.windows(asm, st, self.L, tokens=True)          # uint8 0..4, [k, L], on dev
            text = self.lut5[raw.long()]
            if self.extra:
                tk, lb = text[:, :-1], text[:, 1:]
            else:
                tk = text
                lb = torch.roll(text, shifts=-1, dims=1)
                lb[:, -1] = self.pad
            lm = torch.ones(tk.shape, dtype=torch.float32, device=dev)
            if self.cfg.eod_mask_loss:
                lm[tk == self.eod] = 0.0
            pos = self.ar.expand(tk.shape[0], S)
            if self.cfg.reset_position_ids:
                # position restarts at 0 right after every end-of-document id (as the loop in
                # _get_ltor_masks_and_position_ids does)
                mark = torch.zeros(tk.shape, dtype=torch.bool, device=dev)
                mark[:, 1:] = tk[:, :-1] == self.eod
                seg = torch.where(mark, self.ar.expand_as(mark), torch.zeros_like(pos)).cummax(dim=1).values
                pos = pos - seg
            lm[lb == self.pad] = 0.0
            tk = tk.masked_fill(tk == self.pad, 0)
            lb = lb.masked_fill(lb == self.pad, 0)
            rows = torch.tensor(fast_r, dtype=torch.int64, device=dev)
            tokens[rows] = tk
            labels[rows] = lb
            loss_mask[rows] = lm
            position_ids[rows] = pos

        if slow_r:
            masks = []
            for r in slow_r:
                d = self.ds[idxs[r]]
                tokens[r] = d["tokens"].to(dev)
                labels[r] = d["labels"].to(dev)
                loss_mask[r] = d["loss_mask"].to(dev)
                position_ids[r] = d["position_ids"].to(dev)
                if "attention_mask" in d:
                    masks.append(d["attention_mask"])
            if masks:
                assert len(masks) == n
                attention_mask = torch.stack(masks).to(dev)

        out = {"tokens": tokens, "labels": labels, "loss_mask": loss_mask, "position_ids": position_ids}
        if attention_mask is not None:
            out["attention_mask"] = attention_mask
        return out
