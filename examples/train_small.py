# SPDX-License-Identifier: MIT
"""Train a small next-base DNA model straight from panvram windows, and compare it with the same model trained on the
reference alone.

Both models see the same number of training tokens. Evaluation is on windows from assemblies held out of training
(split by assembly, never by window), so a gain from the cohort is a gain on unseen people, not memorisation.

    python examples/train_small.py --cohort cohort/ --device cuda --holdout 50 --steps 2000
    python examples/train_small.py --cohort tests_cohort/ --reference ref.fa --device cpu --holdout 1 --steps 50 --W 512

Output: one JSON line per model with final held-out loss (bits per base) and training throughput.
"""
import argparse
import json
import random
import time

import torch
import torch.nn as nn
import torch.nn.functional as F

import panvram

V = 5  # A C G T N (panvram tokens=True)


class TinyCNNLM(nn.Module):
    """Causal dilated CNN: small enough to train in minutes, big enough to show a data effect."""

    def __init__(self, d=128, layers=8):
        super().__init__()
        self.emb = nn.Embedding(V, d)
        self.convs = nn.ModuleList(nn.Conv1d(d, d, 3, dilation=2 ** (i % 6)) for i in range(layers))
        self.norms = nn.ModuleList(nn.LayerNorm(d) for _ in range(layers))
        self.out = nn.Linear(d, V)

    def forward(self, x):  # x int64 [B, W] -> logits [B, W, V] predicting x[:, t+1] at t
        h = self.emb(x)
        for conv, norm in zip(self.convs, self.norms):
            pad = conv.dilation[0] * 2
            y = conv(F.pad(h.transpose(1, 2), (pad, 0))).transpose(1, 2)
            h = h + F.gelu(norm(y))
        return self.out(h)


def bits_per_base(model, x):
    logits = model(x[:, :-1])
    return F.cross_entropy(logits.reshape(-1, V), x[:, 1:].reshape(-1)).item() / torch.log(torch.tensor(2.0)).item()


def read_fasta_tokens(path):
    lut = torch.full((256,), 4, dtype=torch.uint8)
    for i, ch in enumerate(b"ACGT"):
        lut[ch] = i
        lut[ch + 32] = i
    seqs, cur = [], []
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if cur:
                    seqs.append(b"".join(cur))
                cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append(b"".join(cur))
    return [lut[torch.frombuffer(bytearray(s), dtype=torch.uint8).long()] for s in seqs if s]


def reference_sampler(contigs, n, W, device, gen):
    ok = [c for c in contigs if len(c) > W]
    lens = torch.tensor([len(c) - W for c in ok], dtype=torch.float64)
    pick = torch.multinomial(lens / lens.sum(), n, replacement=True, generator=gen)
    out = torch.empty(n, W, dtype=torch.long)
    for k, ci in enumerate(pick.tolist()):
        s = int(torch.randint(0, len(ok[ci]) - W, (1,), generator=gen))
        out[k] = ok[ci][s:s + W]
    return out.to(device)


def train(name, next_batch, evalset, args, device):
    torch.manual_seed(args.seed)
    model = TinyCNNLM(args.d, args.layers).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr)
    t0, toks = time.time(), 0
    for step in range(args.steps):
        x = next_batch()
        loss = F.cross_entropy(model(x[:, :-1]).reshape(-1, V), x[:, 1:].reshape(-1))
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        toks += x.numel()
    if device.type == "cuda":
        torch.cuda.synchronize()
    dt = time.time() - t0
    model.eval()
    with torch.no_grad():
        held = sum(bits_per_base(model, e) for e in evalset) / len(evalset)
    return {"model": name, "steps": args.steps, "train_tokens": toks, "train_s": round(dt, 2),
            "tokens_per_s": round(toks / dt), "heldout_bits_per_base": round(held, 4)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cohort", required=True)
    p.add_argument("--reference", default=None, help="reference FASTA (default: the one panvram finds in --cohort)")
    p.add_argument("--dataset", default="q4k")
    p.add_argument("--device", default="cuda")
    p.add_argument("--holdout", type=int, default=50, help="assemblies held out for evaluation")
    p.add_argument("--steps", type=int, default=2000)
    p.add_argument("--batch", type=int, default=64)
    p.add_argument("--W", type=int, default=2048)
    p.add_argument("--eval_batches", type=int, default=8)
    p.add_argument("--d", type=int, default=128)
    p.add_argument("--layers", type=int, default=8)
    p.add_argument("--lr", type=float, default=2e-3)
    p.add_argument("--seed", type=int, default=20261115)
    args = p.parse_args()
    device = torch.device(args.device)

    names = panvram.Cohort.open(args.cohort, device="cpu", dataset=args.dataset, reference=args.reference).names
    rnd = random.Random(args.seed)
    held = sorted(rnd.sample(names, args.holdout))
    train_names = [n for n in names if n not in held]
    tr = panvram.Cohort.open(args.cohort, device=args.device, dataset=args.dataset, reference=args.reference, assemblies=train_names)
    ev = panvram.Cohort.open(args.cohort, device=args.device, dataset=args.dataset, reference=args.reference, assemblies=held)

    g_ev = torch.Generator(device=args.device).manual_seed(args.seed + 1)
    evalset = [ev.sample(args.batch, args.W, generator=g_ev, tokens=True).long() for _ in range(args.eval_batches)]

    g_tr = torch.Generator(device=args.device).manual_seed(args.seed)
    cohort_batch = lambda: tr.sample(args.batch, args.W, generator=g_tr, tokens=True, reverse_complement=0.5).long()

    ref =read_fasta_tokens(args.reference) if args.reference else None
    results = [train("cohort", cohort_batch, evalset, args, device)]
    if ref is not None:
        g_ref = torch.Generator().manual_seed(args.seed)
        ref_batch = lambda: reference_sampler(ref, args.batch, args.W, device, g_ref)
        results.append(train("reference_only", ref_batch, evalset, args, device))
    meta = {"train_assemblies": len(train_names), "heldout_assemblies": held, "W": args.W, "batch": args.batch,
            "device": str(device), "panvram": getattr(panvram, "__version__", "?")}
    for r in results:
        print(json.dumps({**r, **meta}))


if __name__ == "__main__":
    main()
