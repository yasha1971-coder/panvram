# Copyright (c) 2026, Yakiv Shavidze. MIT licence.
"""PanvramIndexedDataset: Megatron's low-level dataset contract served straight from panvram (refrel3 v1) archives,
with no .bin/.idx on disk.

Contract (what megatron/core/datasets/gpt_dataset.py needs of a low-level dataset):
  * `sequence_lengths` -> numpy int32 array, one entry per document;
  * `get(idx, offset=0, length=None)` -> numpy array of token ids of document `idx`;
  * `__len__`, `document_indices`.

One document = one contig, in the archive's record order, assemblies in the order given. Tokens are the ids of
`--gengram-token-ids` (A=8, C=5, G=6, T=7, anything else=9; megatron/core/transformer/transformer_config.py), and a
virtual end-of-document id (12 = `--vocab-size` of the NullTokenizer) ends every document, the token
`tools/preprocess_data.py --append-eod` writes. The end-of-document id is produced on read; nothing is stored for it.

Data path syntax (`--data-path`): `panvram:<cohort directory>`. Environment, all optional:
  PANVRAM_ASSEMBLIES  comma-separated assembly names (default: every archive in the directory)
  PANVRAM_DEVICE      cuda | cpu (default: cuda when available) - where the cohort is resident
  PANVRAM_DATASET     archive kind (default q16k)
  PANVRAM_THREADS     CPU decoder threads (default 0 = all)
"""
import hashlib
import json
import os
from typing import List, Optional

import numpy as np
import torch

PREFIX = "panvram:"
EOD_ID = 12
TOKEN_IDS = (8, 5, 6, 7, 9)                                  # A, C, G, T, other


def is_panvram_path(path: Optional[str]) -> bool:
    return isinstance(path, str) and path.startswith(PREFIX)


def _lut(token_ids=TOKEN_IDS):
    """byte -> token id for the bytes panvram returns (case kept, so both cases are mapped)."""
    a, c, g, t, n = token_ids
    lut = np.full(256, n, dtype=np.uint16)
    for ch, v in (('A', a), ('C', c), ('G', g), ('T', t)):
        lut[ord(ch)] = v
        lut[ord(ch.lower())] = v
    return lut


class PanvramIndexedDataset:
    def __init__(self, path: str, assemblies: Optional[List[str]] = None, dataset: str = "q16k",
                 device: Optional[str] = None, threads: int = 0, token_ids=TOKEN_IDS, eod_id: int = EOD_ID,
                 append_eod: bool = True):
        import panvram
        if device is None:
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.cohort = panvram.Cohort.open(path, device=device, dataset=dataset, threads=threads,
                                          assemblies=assemblies)
        # GPTDataset caches its indices under <path_prefix>/cache, keyed by the data path string only; the data path
        # does not name the assemblies, so the prefix must: one directory per (assemblies, kind, ids, eod) selection
        sel = json.dumps([list(self.cohort.names), dataset, list(token_ids), eod_id, append_eod])
        self.path_prefix = os.path.join(path, "panvram-" + hashlib.md5(sel.encode()).hexdigest()[:16])
        self.device = torch.device(device)
        self.token_ids = tuple(token_ids)
        self.lut = _lut(token_ids)
        self.eod_id = eod_id
        self.append_eod = append_eod
        self.dtype = np.uint16
        asm, ln, boff, name = [], [], [], []
        for a in range(len(self.cohort.names)):
            for rec in self.cohort._core.contigs(a):
                cname, _, clen, off, _ = rec
                asm.append(a)
                ln.append(int(clen))
                boff.append(int(off))
                name.append(cname)
        self.doc_asm = np.array(asm, dtype=np.int64)              # assembly of every document
        self.doc_start = np.array(boff, dtype=np.int64)           # offset of its first base in the assembly's stream
        self.doc_bases = np.array(ln, dtype=np.int64)             # bases, without the end-of-document id
        self._name = name
        self._lengths = (self.doc_bases + (1 if append_eod else 0)).astype(np.int32)
        self._doc_idx = np.arange(len(ln) + 1, dtype=np.int64)

    @classmethod
    def from_data_path(cls, data_path: str):
        assert is_panvram_path(data_path), data_path
        env = os.environ
        names = env.get("PANVRAM_ASSEMBLIES")
        return cls(data_path[len(PREFIX):],
                   assemblies=[x for x in names.split(",") if x] if names else None,
                   dataset=env.get("PANVRAM_DATASET", "q16k"),
                   device=env.get("PANVRAM_DEVICE") or None,
                   threads=int(env.get("PANVRAM_THREADS", "0")))

    # ---- the contract ----
    @property
    def sequence_lengths(self) -> np.ndarray:
        return self._lengths

    @property
    def document_indices(self) -> np.ndarray:
        return self._doc_idx

    def __len__(self) -> int:
        return len(self._lengths)

    def get(self, idx: int, offset: int = 0, length: Optional[int] = None) -> np.ndarray:
        idx, offset = int(idx), int(offset)                      # Megatron passes numpy ints here
        total = int(self._lengths[idx])
        length = total - offset if length is None else int(length)
        if offset < 0 or length < 0 or offset + length > total:
            raise IndexError(f"document {idx}: offset {offset} length {length} outside {total}")
        take = max(0, min(length, int(self.doc_bases[idx]) - offset))   # the part of the request that is bases
        out = np.empty(length, dtype=self.dtype)
        if take:
            asm = torch.tensor([int(self.doc_asm[idx])], dtype=torch.int64)
            st = torch.tensor([int(self.doc_start[idx]) + offset], dtype=torch.int64)
            raw = self.cohort.windows(asm, st, take)[0]
            out[:take] = self.lut[raw.cpu().numpy()]
        if take < length:                                        # the tail is the virtual end-of-document id
            out[take:] = self.eod_id
        return out

    def __getitem__(self, idx):
        return self.get(idx)

    def document_name(self, idx: int) -> str:
        return f"{self.cohort.names[int(self.doc_asm[idx])]}:{self._name[idx]}"
