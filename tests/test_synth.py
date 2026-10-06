"""Fast tests on the synthetic fixtures (tests/data, encoded from tests/synth.py by aceapex research/refrel/refrel3v1
@ 5b6d5ce): every path (CPU, and CUDA when present) against windows cut from the FASTA; refusals."""
import os
import random
import shutil

import pytest
import torch

import panvram
from fasta_cut import FastaIndex, transform
import synth

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
DEVICES = ["cpu"] + (["cuda"] if panvram.with_cuda and torch.cuda.is_available() else [])


@pytest.fixture(scope="session")
def gen(tmp_path_factory):
    d = tmp_path_factory.mktemp("synth")
    files = synth.write(str(d))
    out = {}
    for ds in ("q1k", "q4k", "q16k"):
        dd = d / ds
        dd.mkdir()
        for nm in ("synthA", "synthB"):
            shutil.copy(os.path.join(DATA, f"{nm}.{ds}.rr3"), dd / f"{nm}.{ds}.rr3")
        out[ds] = str(dd)
    fa = {nm: open(files[nm + ".fa"], "rb").read() for nm in ("synthA", "synthB")}
    return {"dirs": out, "ref": files["synth_ref.fa"], "fa": fa, "idx": {k: FastaIndex(v) for k, v in fa.items()}}


FORMS = [(False, False), (True, False), (False, True), (True, True)]   # (packed_reference, compact_blocks)


@pytest.mark.parametrize("form", FORMS, ids=["bytes-full", "packed-full", "bytes-compact", "packed-compact"])
@pytest.mark.parametrize("device", DEVICES)
@pytest.mark.parametrize("ds", ["q1k", "q4k", "q16k"])
def test_windows_fetch_fasta(gen, ds, device, form):
    c = panvram.Cohort.open(gen["dirs"][ds], device=device, dataset=ds, reference=gen["ref"], packed_reference=form[0], compact_blocks=form[1])
    assert c.names == ["synthA", "synthB"] and c.block_size == {"q1k": 1024, "q4k": 4096, "q16k": 16384}[ds]
    for i, nm in enumerate(c.names):
        assert c.fasta(i) == gen["fa"][nm]
    g = torch.Generator(device=device).manual_seed(7)
    for W in (1, 13, 1024, 4096, 5000, 16385):
        x, co = c.sample(300, W, generator=g, return_coords=True)
        assert x.shape == (300, W) and x.dtype == torch.uint8 and x.device.type == device
        xs, co = x.cpu().numpy(), co.cpu().tolist()
        for k, (a, ci, s, strand) in enumerate(co):
            name = c.contigs(a)[ci][0]
            assert strand == 0
            assert xs[k].tobytes() == gen["idx"][c.names[a]].cut(name, s, W), (W, a, name, s)
    for W in (100, 8192):
        x, co = c.sample(200, W, generator=g, reverse_complement=0.5, tokens=True, return_coords=True)
        xs, co = x.cpu().numpy(), co.cpu().tolist()
        assert 0 < sum(r[3] for r in co) < 200
        for k, (a, ci, s, strand) in enumerate(co):
            want = transform(gen["idx"][c.names[a]].cut(c.contigs(a)[ci][0], s, W), rc=bool(strand), tokens=True)
            assert xs[k].tobytes() == want
    rng = random.Random(11)
    for _ in range(200):
        a = rng.randrange(len(c))
        ctgs = c.contigs(a)
        ci = rng.randrange(len(ctgs))
        name, _, n, _, _ = ctgs[ci]
        L = min(n, int(2 ** (rng.random() * 15)))
        s = rng.randrange(n - L + 1)
        rc, tok = rng.random() < 0.3, rng.random() < 0.3
        y = c.fetch(c.names[a] if rng.random() < 0.5 else a, name if rng.random() < 0.5 else ci, s, L, reverse_complement=rc, tokens=tok)
        assert y.shape == (L,) and y.cpu().numpy().tobytes() == transform(gen["idx"][c.names[a]].cut(name, s, L), rc, tok)
    # whole contigs
    for a in range(len(c)):
        for name, _, n, _, _ in c.contigs(a):
            assert c.fetch(a, name, 0, n).cpu().numpy().tobytes() == gen["idx"][c.names[a]].cut(name, 0, n)


@pytest.mark.skipif(len(DEVICES) < 2, reason="no CUDA")
def test_cuda_equals_cpu(gen):
    cg = panvram.Cohort.open(gen["dirs"]["q4k"], device="cuda", reference=gen["ref"])
    cc = panvram.Cohort.open(gen["dirs"]["q4k"], device="cpu", reference=gen["ref"])
    g = torch.Generator(device="cuda").manual_seed(3)
    x, co = cg.sample(4096, 777, generator=g, reverse_complement=0.5, return_coords=True)
    asm = co[:, 0].cpu()
    start = torch.tensor([cc.contigs(int(a))[int(ci)][3] + int(s) for a, ci, s, _ in co.cpu().tolist()])
    y = cc.windows(asm, start, 777, co[:, 3].cpu().bool())
    assert torch.equal(x.cpu(), y)


def test_refusals(gen, tmp_path):
    good = os.path.join(gen["dirs"]["q4k"], "synthA.q4k.rr3")
    raw = open(good, "rb").read()
    # wrong reference: one base changed
    ref = bytearray(open(gen["ref"], "rb").read())
    ref[ref.index(b"\n") + 5] ^= 0x01 if ref[ref.index(b"\n") + 5] != ord("A") else 0x02
    (tmp_path / "bad_ref.fa").write_bytes(bytes(ref))
    d = tmp_path / "one"
    d.mkdir()
    shutil.copy(good, d / "synthA.q4k.rr3")
    with pytest.raises(panvram.FormatError, match="wrong reference"):
        panvram.Cohort.open(str(d), device="cpu", reference=str(tmp_path / "bad_ref.fa"))
    # header byte, header XXH3 field, truncation, unknown version
    for k, mut in enumerate((lambda b: b[:40] + bytes([b[40] ^ 1]) + b[41:], lambda b: b[:130] + bytes([b[130] ^ 1]) + b[131:],
                             lambda b: b[:-1], lambda b: b[:8] + b"\x02" + b[9:])):
        dd = tmp_path / f"bad{k}"
        dd.mkdir()
        (dd / "synthA.q4k.rr3").write_bytes(mut(raw))
        with pytest.raises(panvram.FormatError, match="refused"):
            panvram.Cohort.open(str(dd), device="cpu", reference=gen["ref"])
    # one block size per cohort
    mix = tmp_path / "mix"
    mix.mkdir()
    shutil.copy(good, mix / "a.rr3")
    shutil.copy(os.path.join(gen["dirs"]["q1k"], "synthB.q1k.rr3"), mix / "b.rr3")
    with pytest.raises(panvram.FormatError, match="block size"):
        panvram.Cohort.open(str(mix), device="cpu", dataset=None, reference=gen["ref"])
    # coordinates outside a contig
    c = panvram.Cohort.open(gen["dirs"]["q4k"], device="cpu", reference=gen["ref"])
    name, _, n, _, _ = c.contigs(0)[0]
    with pytest.raises(IndexError):
        c.fetch(0, name, n - 5, 6)
    with pytest.raises(ValueError):
        c.sample(4, 10 ** 9)


def test_corrupt_payload_caught_with_hashes(gen, tmp_path):
    """A payload byte flipped: the CPU full decode with verification never returns wrong bytes."""
    raw = bytearray(open(os.path.join(gen["dirs"]["q4k"], "synthA.q4k.rr3"), "rb").read())
    caught = 0
    for k in range(40):
        b = bytearray(raw)
        p = len(b) - 1 - (k * 97) % 3000
        b[p] ^= 1 << (k % 8)
        d = tmp_path / f"c{k}"
        d.mkdir()
        (d / "synthA.q4k.rr3").write_bytes(bytes(b))
        c = panvram.Cohort.open(str(d), device="cpu", reference=gen["ref"])
        try:
            fa = c.fasta(0)
        except panvram.FormatError:
            caught += 1
            continue
        assert fa == gen["fa"]["synthA"]
    assert caught > 0


@pytest.mark.parametrize("form", FORMS[1:], ids=["packed-full", "bytes-compact", "packed-compact"])
def test_compact_forms_resident_bytes(gen, form):
    """The compact forms hold fewer bytes and the same windows as the default."""
    base = panvram.Cohort.open(gen["dirs"]["q1k"], device="cpu", reference=gen["ref"], dataset="q1k")
    c = panvram.Cohort.open(gen["dirs"]["q1k"], device="cpu", reference=gen["ref"], dataset="q1k", packed_reference=form[0], compact_blocks=form[1])
    rb, rc = base.resident_bytes(), c.resident_bytes()
    assert rc["total"] < rb["total"]
    g = torch.Generator().manual_seed(5)
    x, co = base.sample(500, 3000, generator=g, return_coords=True)
    asm = co[:, 0]
    start = torch.tensor([base.contigs(int(a))[int(ci)][3] + int(s) for a, ci, s, _ in co.tolist()])
    assert torch.equal(c.windows(asm, start, 3000), x)
