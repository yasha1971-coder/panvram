"""panvram.Cohort - a pangenome cohort resident in memory (GPU or CPU) as refrel3 v1 archives against one decoded
reference; windows and coordinate slices decoded on demand by the CUDA queue kernel or the CPU decoder."""
import glob
import os
import re

import torch

from . import _C

DATASETS = ("q4k", "q16k")


def _arch_ok(device):
    cap = torch.cuda.get_device_capability(device)
    if cap < (8, 0):
        raise RuntimeError(f"panvram needs sm_80 or newer, {torch.cuda.get_device_name(device)} is sm_{cap[0]}{cap[1]}")


class Cohort:
    """Open with Cohort.open(dir, device=...). Coordinates are 0-based; a contig is its FASTA header up to the first
    space or tab, or its index in the assembly; an assembly is its name (archive file name without .<dataset>.rr3) or
    its index."""

    def __init__(self, core, names, device):
        self._core, self.names, self.device = core, list(names), torch.device(device)
        self._index = {n: i for i, n in enumerate(self.names)}
        pools = core.pools()
        self._ctg = {k: pools[k] for k in ("ctg_asm", "ctg_local", "ctg_boff", "ctg_len")}
        self._ctg_names = [{c[0]: k for k, c in enumerate(core.contigs(i))} for i in range(core.n)]
        if self.device.type == "cuda":
            if not _C.with_cuda:
                raise RuntimeError("panvram was built without CUDA (CPU path only): reinstall where nvcc and a CUDA torch are present")
            _arch_ok(self.device)
            keep = ("P", "pay_base", "off", "st", "tab", "blk_base", "nbases", "low_s", "low_l", "low_base", "low_cnt", "ref")
            self._dev = {k: pools[k].to(self.device) for k in keep}
            self._ctg = {k: v.to(self.device) for k, v in self._ctg.items()}
        else:
            self._dev = None
        self._valid = {}

    @classmethod
    def open(cls, path, device="cuda", dataset="q4k", reference=None, threads=0, assemblies=None):
        """Open the cohort directory `path`: every `*.<dataset>.rr3` in it (sorted by name; `assemblies` = a list of
        names to take only those) and the reference they were encoded against - `reference` if given, else
        `path/reference.fa`, else the file named in the archives' header in `path`. Every archive is checked before
        any decode (FORMAT.md section 5): header XXH3, reference SHA-256 and size, sections, tables. One block size
        per cohort. device: "cuda", "cuda:N" or "cpu"."""
        if dataset is not None and not re.fullmatch(r"q\d+k", dataset):
            raise ValueError(f"dataset: {DATASETS} (or another q<N>k suffix), or None for every *.rr3")
        pat = f"*.{dataset}.rr3" if dataset else "*.rr3"
        files = sorted(glob.glob(os.path.join(path, pat)))
        suffix = f".{dataset}.rr3" if dataset else ".rr3"
        names = [os.path.basename(f)[: -len(suffix)] for f in files]
        if assemblies is not None:
            want = list(assemblies)
            have = dict(zip(names, files))
            missing = [a for a in want if a not in have]
            if missing:
                raise FileNotFoundError(f"not in {path}: {missing[:5]}")
            names, files = want, [have[a] for a in want]
        if not files:
            raise FileNotFoundError(f"no {pat} in {path}")
        if reference is None:
            cand = os.path.join(path, "reference.fa")
            if not os.path.exists(cand):
                with open(files[0], "rb") as f:
                    h = f.read(138 + 65535)
                nl = h[136] | h[137] << 8
                cand = os.path.join(path, h[138:138 + nl].decode())
            reference = cand
        if not os.path.exists(reference):
            raise FileNotFoundError(f"reference {reference} not found (pass reference=)")
        device = torch.device(device)
        if device.type == "cuda" and device.index is None:
            device = torch.device("cuda", torch.cuda.current_device())
        core = _C.Core(os.fspath(reference), [os.fspath(f) for f in files], int(threads))
        return cls(core, names, device)

    # ------------------------------------------------------------------ description
    def __len__(self):
        return len(self.names)

    @property
    def block_size(self):
        return self._core.Q

    @property
    def reference_sha256(self):
        return self._core.ref_sha256

    def info(self, assembly):
        return self._core.info(self._asm(assembly))

    def contigs(self, assembly):
        """[(name, header, length, offset in the base stream, line width)] of an assembly."""
        return self._core.contigs(self._asm(assembly))

    def resident_bytes(self):
        """Bytes of the cohort on its device (reference included), by component."""
        src = self._dev if self._dev is not None else self._core.pools()
        d = {k: v.numel() * v.element_size() for k, v in src.items() if k in ("P", "off", "st", "tab", "low_s", "low_l", "ref")}
        d["total"] = sum(d.values())
        return d

    # ------------------------------------------------------------------ decode
    def _asm(self, a):
        if isinstance(a, str):
            if a not in self._index:
                raise KeyError(f"no assembly {a}")
            return self._index[a]
        a = int(a)
        if not 0 <= a < len(self.names):
            raise IndexError(f"assembly index {a}")
        return a

    def _contig(self, a, c):
        if isinstance(c, str):
            if c not in self._ctg_names[a]:
                raise KeyError(f"no contig {c} in {self.names[a]}")
            return self._ctg_names[a][c]
        c = int(c)
        if not 0 <= c < len(self._ctg_names[a]):
            raise IndexError(f"contig index {c}")
        return c

    def windows(self, asm, start, W, reverse_complement=None, tokens=False, check=True):
        """Bases [start, start + W) of assembly asm[i] (start = offset in the assembly's base stream: contigs
        concatenated in file order) for every i -> uint8 [n, W] on the cohort's device. reverse_complement: None or a
        bool tensor [n]. tokens: A/a 0, C/c 1, G/g 2, T/t 3, anything else 4; else the FASTA's bytes, case kept.
        check: wait for the kernel and raise if any block failed to decode (one device sync)."""
        W = int(W)
        if self._dev is None:
            return self._core.windows(asm, start, W, reverse_complement, tokens, True, False, 0)
        d = self._dev
        out, status = _C.windows_cuda(d["P"], d["pay_base"], d["off"], d["st"], d["tab"], d["blk_base"], d["nbases"], d["low_s"], d["low_l"],
                                      d["low_base"], d["low_cnt"], d["ref"], self._core.Q, asm, start, W, reverse_complement, tokens, True)
        if check:
            s = int(status.item())
            if s:
                raise RuntimeError(f"GPU decode failed (status {s}: 1 corrupt block, 4 bounded wait expired, 8 window outside an assembly)")
        return out

    def sample(self, n, W, generator=None, reverse_complement=False, tokens=False, return_coords=False, check=True):
        """n windows of W bases, uniform over every (assembly, contig, start) with start + W <= contig length (a
        window never crosses contigs) -> uint8 [n, W] on the cohort's device. generator: a torch.Generator (its device
        draws the numbers). reverse_complement: False, True, or p in (0, 1) - each window reverse-complemented with
        probability p. return_coords: also an int64 [n, 4] tensor (assembly, contig index, start in the contig,
        reverse-complemented 0/1)."""
        n, W = int(n), int(W)
        if W <= 0 or n < 0:
            raise ValueError("W must be positive, n non-negative")
        if W not in self._valid:
            v = (self._ctg["ctg_len"] - W + 1).clamp_min(0)
            self._valid[W] = (torch.cumsum(v, 0), v)
        cum, v = self._valid[W]
        total = int(cum[-1].item()) if cum.numel() else 0
        if total == 0:
            raise ValueError(f"no contig of {W} bases or more")
        gdev = generator.device if generator is not None else self.device
        u = torch.randint(0, total, (n,), generator=generator, device=gdev, dtype=torch.int64).to(self.device)
        k = torch.searchsorted(cum, u, right=True)
        pos = u - (cum[k] - v[k])
        asm = self._ctg["ctg_asm"][k]
        start = self._ctg["ctg_boff"][k] + pos
        if reverse_complement is True:
            rc = torch.ones(n, dtype=torch.bool, device=self.device)
        elif reverse_complement is False or reverse_complement is None or reverse_complement == 0:
            rc = None
        else:
            p = float(reverse_complement)
            if not 0 < p <= 1:
                raise ValueError("reverse_complement: False, True or a probability")
            rc = (torch.rand(n, generator=generator, device=gdev) < p).to(self.device)
        x = self.windows(asm, start, W, rc, tokens, check)
        if not return_coords:
            return x
        strand = rc.to(torch.int64) if rc is not None else torch.zeros(n, dtype=torch.int64, device=self.device)
        return x, torch.stack([asm, self._ctg["ctg_local"][k], pos, strand], 1)

    def fetch(self, assembly, contig, start, length, reverse_complement=False, tokens=False):
        """Bases [start, start + length) of a contig (0-based) -> uint8 [length] on the cohort's device."""
        a = self._asm(assembly)
        c = self._contig(a, contig)
        name, _, clen, boff, _ = self._core.contigs(a)[c]
        start, length = int(start), int(length)
        if start < 0 or length < 0 or start + length > clen:
            raise IndexError(f"{self.names[a]} {name}:{start}+{length} outside the contig ({clen} bases)")
        asm = torch.tensor([a], dtype=torch.int64)
        st = torch.tensor([boff + start], dtype=torch.int64)
        rc = torch.tensor([bool(reverse_complement)]) if reverse_complement else None
        return self.windows(asm, st, length, rc, tokens)[0]

    def fasta(self, assembly, verify=True):
        """Full decode of an assembly on the CPU -> the FASTA file's bytes. verify: every block XXH3 (if present) and
        the XXH3 of the rebuilt FASTA against the header's (the source file's)."""
        return self._core.fasta(self._asm(assembly), verify)
